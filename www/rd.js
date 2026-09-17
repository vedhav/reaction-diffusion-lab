/*
 * Reaction-Diffusion Lab -- client-side Gray-Scott solver.
 *
 * The whole simulation lives here so that interaction stays at 60fps: the
 * server is only told *about* the state, it never drives it. Communication
 * with Shiny is deliberately narrow:
 *
 *   JS -> R   input$sim_state  (debounced parameter snapshot)
 *             input$activity   (mean concentration of B, ~4Hz)
 *   R  -> JS  "rd:apply"       (restore a bookmarked state)
 */
(function () {
  "use strict";

  // --- constants ---------------------------------------------------------

  var GRID = 256;             // simulation is GRID x GRID on a torus
  var DU = 1.0;               // diffusion rate of chemical A
  var DV = 0.5;               // diffusion rate of chemical B
  var DT = 1.0;               // timestep (stable for the stencil below)
  var V_FULL = 0.4;           // concentration of B mapped to the top of the LUT

  // Parameter-map extent. Must match F_RANGE/K_RANGE in app.R.
  var F_MIN = 0.0, F_MAX = 0.1;
  var K_MIN = 0.03, K_MAX = 0.075;

  var STATE_DEBOUNCE_MS = 120;
  var ACTIVITY_INTERVAL_MS = 250;
  var SNAPSHOT_PX = 1024;

  // --- palettes ----------------------------------------------------------
  // Each palette is a list of [position, r, g, b] stops, expanded into a
  // 256-entry lookup table once at startup.

  var PALETTES = {
    "Inferno":  [[0.00,   4,   2,  18], [0.35,  70,  12,  80],
                 [0.62, 196,  60,  58], [0.84, 248, 165,  38], [1.00, 252, 253, 191]],
    "Abyss":    [[0.00,   3,   8,  26], [0.40,  10,  61, 110],
                 [0.72,  32, 170, 205], [1.00, 226, 251, 255]],
    "Jade":     [[0.00,   2,  14,  12], [0.42,  11,  84,  63],
                 [0.75,  52, 197, 130], [1.00, 227, 255, 235]],
    "Ember":    [[0.00,  10,   4,   4], [0.45, 110,  18,  20],
                 [0.75, 228,  96,  30], [1.00, 255, 234, 196]],
    "Orchid":   [[0.00,  14,   5,  26], [0.40, 104,  32, 140],
                 [0.74, 214,  95, 198], [1.00, 255, 226, 248]],
    "Bone":     [[0.00,   8,   8,  10], [0.55, 120, 122, 130], [1.00, 255, 255, 255]]
  };

  /** Expand colour stops into a flat RGBA lookup table of 256 entries. */
  function buildLut(stops) {
    var lut = new Uint8ClampedArray(256 * 4);
    var s = 0;
    for (var i = 0; i < 256; i++) {
      var t = i / 255;
      while (s < stops.length - 2 && t > stops[s + 1][0]) s++;
      var a = stops[s], b = stops[s + 1];
      var span = b[0] - a[0];
      var w = span > 0 ? (t - a[0]) / span : 0;
      if (w < 0) w = 0; else if (w > 1) w = 1;
      lut[i * 4 + 0] = a[1] + (b[1] - a[1]) * w;
      lut[i * 4 + 1] = a[2] + (b[2] - a[2]) * w;
      lut[i * 4 + 2] = a[3] + (b[3] - a[3]) * w;
      lut[i * 4 + 3] = 255;
    }
    return lut;
  }

  var LUTS = {};
  Object.keys(PALETTES).forEach(function (name) {
    LUTS[name] = buildLut(PALETTES[name]);
  });

  // --- model -------------------------------------------------------------

  /** Largest k that still admits a non-trivial steady state, for a given f. */
  function killLimit(f) {
    return Math.sqrt(f) / 2 - f;
  }

  function clamp(x, lo, hi) {
    return x < lo ? lo : (x > hi ? hi : x);
  }

  /**
   * Gray-Scott reactor on a periodic grid.
   *
   * Uses ping-pong Float32Array buffers and precomputed neighbour indices so
   * the inner loop is free of modulo arithmetic and allocation.
   */
  function Reactor(n) {
    this.n = n;
    var size = n * n;
    this.u = new Float32Array(size);
    this.v = new Float32Array(size);
    this.uNext = new Float32Array(size);
    this.vNext = new Float32Array(size);

    // Wrap-around neighbour tables: xm/xp are column offsets, ym/yp are
    // already multiplied by n so they can be added to a column index.
    this.xm = new Int32Array(n);
    this.xp = new Int32Array(n);
    this.ym = new Int32Array(n);
    this.yp = new Int32Array(n);
    for (var i = 0; i < n; i++) {
      this.xm[i] = (i - 1 + n) % n;
      this.xp[i] = (i + 1) % n;
      this.ym[i] = ((i - 1 + n) % n) * n;
      this.yp[i] = ((i + 1) % n) * n;
    }

    this.clear();
  }

  /** Uniform state: all A, no B. */
  Reactor.prototype.clear = function () {
    this.u.fill(1);
    this.v.fill(0);
  };

  /**
   * Inject chemical B in a disc.
   *
   * The cells are blended towards the half-and-half mixture (u = v = 0.5)
   * rather than having B piled on top of untouched A. That mixture is what
   * reliably ignites the reaction: seeds with a depleted-A core tend to burn
   * out before the autocatalytic front can organise.
   */
  Reactor.prototype.seed = function (cx, cy, radius, strength) {
    var n = this.n, r2 = radius * radius;
    var inner = radius * 0.75;
    var y0 = Math.floor(cy - radius), y1 = Math.ceil(cy + radius);
    var x0 = Math.floor(cx - radius), x1 = Math.ceil(cx + radius);
    for (var y = y0; y <= y1; y++) {
      var row = (((y % n) + n) % n) * n;
      for (var x = x0; x <= x1; x++) {
        var dx = x - cx, dy = y - cy;
        var d2 = dx * dx + dy * dy;
        if (d2 > r2) continue;
        var d = Math.sqrt(d2);
        // Flat core, short feathered rim, so the brush has a soft edge
        // without diluting the middle of the stroke.
        var edge = d <= inner ? 1 : (radius - d) / (radius - inner);
        var w = clamp(strength * edge, 0, 1);
        var c = row + (((x % n) + n) % n);
        this.u[c] += (0.5 - this.u[c]) * w;
        this.v[c] += (0.5 - this.v[c]) * w;
      }
    }
  };

  /** Scatter a handful of random blobs to kick the reaction off. */
  Reactor.prototype.randomize = function (blobs) {
    this.clear();
    var n = this.n;
    for (var i = 0; i < blobs; i++) {
      this.seed(Math.random() * n, Math.random() * n, 4 + Math.random() * 8, 1);
    }
  };

  /**
   * The canonical starting condition: one disc of B in the middle, plus a
   * whisper of noise.
   *
   * Without the noise a perfectly symmetric seed grows into a perfectly
   * symmetric mandala; the noise lets the patterns come out organic.
   */
  Reactor.prototype.seedCentre = function () {
    this.clear();
    this.seed(this.n / 2, this.n / 2, this.n / 16, 1);
    var v = this.v;
    for (var i = 0; i < v.length; i++) v[i] += Math.random() * 0.02;
  };

  /** Advance the system by one timestep. */
  Reactor.prototype.step = function (f, k) {
    var n = this.n;
    var u = this.u, v = this.v, uN = this.uNext, vN = this.vNext;
    var xm = this.xm, xp = this.xp, ym = this.ym, yp = this.yp;
    var fk = f + k;

    for (var y = 0; y < n; y++) {
      var row = y * n, up = ym[y], dn = yp[y];
      for (var x = 0; x < n; x++) {
        var c = row + x, xl = xm[x], xr = xp[x];
        var uc = u[c], vc = v[c];

        // 9-point Laplacian: 0.2 orthogonal, 0.05 diagonal, -1 centre.
        var lapU = 0.2 * (u[row + xl] + u[row + xr] + u[up + x] + u[dn + x]) +
                   0.05 * (u[up + xl] + u[up + xr] + u[dn + xl] + u[dn + xr]) - uc;
        var lapV = 0.2 * (v[row + xl] + v[row + xr] + v[up + x] + v[dn + x]) +
                   0.05 * (v[up + xl] + v[up + xr] + v[dn + xl] + v[dn + xr]) - vc;

        var reaction = uc * vc * vc;
        var un = uc + (DU * lapU - reaction + f * (1 - uc)) * DT;
        var vn = vc + (DV * lapV + reaction - fk * vc) * DT;

        uN[c] = un < 0 ? 0 : (un > 1 ? 1 : un);
        vN[c] = vn < 0 ? 0 : (vn > 1 ? 1 : vn);
      }
    }

    this.u = uN; this.uNext = u;
    this.v = vN; this.vNext = v;
  };

  /** Mean concentration of B, on a coarse sample of the grid. */
  Reactor.prototype.activity = function () {
    var v = this.v, total = 0, count = 0;
    for (var i = 0; i < v.length; i += 4) { total += v[i]; count++; }
    return total / count;
  };

  /** Paint the current state of B into an ImageData through a LUT. */
  Reactor.prototype.paint = function (image, lut) {
    var data = image.data, v = this.v, scale = 255 / V_FULL;
    for (var i = 0; i < v.length; i++) {
      var idx = (v[i] * scale) | 0;
      if (idx > 255) idx = 255;
      idx <<= 2;
      var o = i << 2;
      data[o] = lut[idx];
      data[o + 1] = lut[idx + 1];
      data[o + 2] = lut[idx + 2];
      data[o + 3] = 255;
    }
  };

  // --- canvas helpers ----------------------------------------------------

  /** Size a canvas to its CSS box at device pixel ratio. Returns {w, h}. */
  function fitCanvas(canvas, cssHeight) {
    var dpr = window.devicePixelRatio || 1;
    var rect = canvas.getBoundingClientRect();
    var w = Math.max(1, Math.round(rect.width));
    var h = Math.max(1, Math.round(cssHeight || rect.height));
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { w: w, h: h, ctx: ctx };
  }

  // --- parameter map -----------------------------------------------------

  /**
   * The f/k phase chart. The shaded region is not decorative: it is exactly
   * the set where k <= sqrt(f)/2 - f, i.e. where a non-trivial steady state
   * exists at all.
   */
  function ParameterMap(canvas, presets, onPick) {
    this.canvas = canvas;
    this.presets = presets;
    this.onPick = onPick;
    this.background = document.createElement("canvas");
    this.hover = null;
    this.f = 0; this.k = 0;
    this.dragging = false;

    var self = this;
    canvas.addEventListener("pointerdown", function (e) {
      if (e.button !== 0) return;
      self.dragging = true;
      canvas.setPointerCapture(e.pointerId);
      self.pick(e);
    });
    canvas.addEventListener("pointermove", function (e) {
      // e.buttons guards against a drag left stuck open by a pointerup we
      // never saw (window blur, capture stolen, devtools).
      if (self.dragging && e.buttons !== 0) self.pick(e);
      else { self.dragging = false; self.track(e); }
    });
    ["pointerup", "pointercancel", "lostpointercapture"].forEach(function (type) {
      canvas.addEventListener(type, function () { self.dragging = false; });
    });
    canvas.addEventListener("pointerleave", function () {
      if (self.hover !== null) { self.hover = null; self.draw(); }
    });
  }

  ParameterMap.prototype.toValue = function (e) {
    var rect = this.canvas.getBoundingClientRect();
    var nx = clamp((e.clientX - rect.left) / rect.width, 0, 1);
    var ny = clamp((e.clientY - rect.top) / rect.height, 0, 1);
    return {
      f: F_MIN + nx * (F_MAX - F_MIN),
      k: K_MAX - ny * (K_MAX - K_MIN)
    };
  };

  ParameterMap.prototype.pick = function (e) {
    var value = this.toValue(e);
    this.onPick(value.f, value.k);
  };

  /** Highlight the preset under the cursor, if any is close enough. */
  ParameterMap.prototype.track = function (e) {
    var rect = this.canvas.getBoundingClientRect();
    var px = e.clientX - rect.left, py = e.clientY - rect.top;
    var best = null, bestDist = 14 * 14;
    for (var i = 0; i < this.presets.length; i++) {
      var p = this.toPixels(this.presets[i].f, this.presets[i].k);
      var dx = p.x - px, dy = p.y - py;
      var d = dx * dx + dy * dy;
      if (d < bestDist) { bestDist = d; best = i; }
    }
    if (best !== this.hover) { this.hover = best; this.draw(); }
  };

  ParameterMap.prototype.toPixels = function (f, k) {
    return {
      x: (f - F_MIN) / (F_MAX - F_MIN) * this.size.w,
      y: (K_MAX - k) / (K_MAX - K_MIN) * this.size.h
    };
  };

  ParameterMap.prototype.resize = function () {
    this.size = fitCanvas(this.canvas);
    var dpr = window.devicePixelRatio || 1;
    this.background.width = this.canvas.width;
    this.background.height = this.canvas.height;
    var ctx = this.background.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.drawBackground(ctx);
    this.draw();
  };

  ParameterMap.prototype.drawBackground = function (ctx) {
    var w = this.size.w, h = this.size.h;
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = "#0b1119";
    ctx.fillRect(0, 0, w, h);

    // The self-sustaining region, filled column by column from the exact
    // boundary k = sqrt(f)/2 - f down to the bottom of the chart.
    var gradient = ctx.createLinearGradient(0, 0, w, h);
    gradient.addColorStop(0, "#15304a");
    gradient.addColorStop(0.5, "#1c4c63");
    gradient.addColorStop(1, "#123043");
    ctx.fillStyle = gradient;
    ctx.beginPath();
    ctx.moveTo(0, h);
    for (var px = 0; px <= w; px++) {
      var f = F_MIN + (px / w) * (F_MAX - F_MIN);
      var limit = clamp(killLimit(f), K_MIN, K_MAX);
      ctx.lineTo(px, (K_MAX - limit) / (K_MAX - K_MIN) * h);
    }
    ctx.lineTo(w, h);
    ctx.closePath();
    ctx.fill();

    ctx.strokeStyle = "rgba(125, 211, 252, 0.55)";
    ctx.lineWidth = 1.25;
    ctx.stroke();

    // Grid lines.
    ctx.strokeStyle = "rgba(148, 163, 184, 0.10)";
    ctx.lineWidth = 1;
    for (var i = 1; i < 5; i++) {
      ctx.beginPath();
      ctx.moveTo(w * i / 5, 0); ctx.lineTo(w * i / 5, h);
      ctx.moveTo(0, h * i / 5); ctx.lineTo(w, h * i / 5);
      ctx.stroke();
    }

    ctx.fillStyle = "rgba(148, 163, 184, 0.75)";
    ctx.font = "10px ui-monospace, SFMono-Regular, Menlo, monospace";
    ctx.textBaseline = "bottom";
    ctx.fillText("f →", w - 22, h - 4);
    ctx.textBaseline = "top";
    ctx.fillText("k ↑", 4, 4);
  };

  ParameterMap.prototype.setValue = function (f, k) {
    this.f = f; this.k = k;
    this.draw();
  };

  ParameterMap.prototype.draw = function () {
    if (!this.size) return;
    var ctx = this.size.ctx, w = this.size.w, h = this.size.h;
    var dpr = window.devicePixelRatio || 1;

    ctx.clearRect(0, 0, w, h);
    ctx.drawImage(this.background, 0, 0, this.background.width / dpr,
                  this.background.height / dpr);

    // Preset markers.
    for (var i = 0; i < this.presets.length; i++) {
      var p = this.toPixels(this.presets[i].f, this.presets[i].k);
      var active = i === this.hover;
      ctx.beginPath();
      ctx.arc(p.x, p.y, active ? 4.5 : 3, 0, Math.PI * 2);
      ctx.fillStyle = active ? "#fbbf24" : "rgba(226, 232, 240, 0.55)";
      ctx.fill();
    }

    // Crosshair for the current parameters.
    var c = this.toPixels(this.f, this.k);
    ctx.strokeStyle = "rgba(251, 191, 36, 0.35)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(c.x, 0); ctx.lineTo(c.x, h);
    ctx.moveTo(0, c.y); ctx.lineTo(w, c.y);
    ctx.stroke();

    ctx.beginPath();
    ctx.arc(c.x, c.y, 5.5, 0, Math.PI * 2);
    ctx.strokeStyle = "#fde68a";
    ctx.lineWidth = 2;
    ctx.stroke();

    // Hover label, nudged to stay inside the chart.
    if (this.hover !== null) {
      var preset = this.presets[this.hover];
      var pt = this.toPixels(preset.f, preset.k);
      ctx.font = "11px ui-monospace, SFMono-Regular, Menlo, monospace";
      var textWidth = ctx.measureText(preset.name).width;
      var bx = clamp(pt.x + 8, 2, w - textWidth - 12);
      var by = clamp(pt.y - 22, 2, h - 20);
      ctx.fillStyle = "rgba(2, 6, 12, 0.85)";
      ctx.fillRect(bx, by, textWidth + 10, 17);
      ctx.fillStyle = "#fde68a";
      ctx.textBaseline = "top";
      ctx.fillText(preset.name, bx + 5, by + 3);
    }
  };

  // --- application -------------------------------------------------------

  function init() {
    var canvas = document.getElementById("rd-canvas");
    var mapCanvas = document.getElementById("rd-map");
    if (!canvas || !mapCanvas) return;

    var presetsNode = document.getElementById("rd-presets");
    var presets = presetsNode ? JSON.parse(presetsNode.textContent) : [];

    var els = {
      feed: document.getElementById("rd-feed"),
      kill: document.getElementById("rd-kill"),
      speed: document.getElementById("rd-speed"),
      brush: document.getElementById("rd-brush"),
      palette: document.getElementById("rd-palette")
    };

    var reactor = new Reactor(GRID);
    reactor.seedCentre();

    // Off-screen buffer at grid resolution; the visible canvas just scales it.
    var buffer = document.createElement("canvas");
    buffer.width = GRID;
    buffer.height = GRID;
    var bufferCtx = buffer.getContext("2d");
    var image = bufferCtx.createImageData(GRID, GRID);

    var view = null;
    var state = {
      f: parseFloat(els.feed.value),
      k: parseFloat(els.kill.value),
      speed: parseInt(els.speed.value, 10),
      brush: parseInt(els.brush.value, 10),
      palette: Object.keys(PALETTES)[0],
      running: true
    };

    var map = new ParameterMap(mapCanvas, presets, function (f, k) {
      setParameters(f, k);
    });

    // --- Shiny bridge ----------------------------------------------------

    var connected = false;
    var pushTimer = null;

    function sendInput(name, value) {
      if (!connected || !window.Shiny) return;
      Shiny.setInputValue(name, value, { priority: "event" });
    }

    function pushState() {
      if (pushTimer) clearTimeout(pushTimer);
      pushTimer = setTimeout(function () {
        pushTimer = null;
        sendInput("sim_state", {
          f: state.f, k: state.k, speed: state.speed,
          palette: state.palette, running: state.running
        });
      }, STATE_DEBOUNCE_MS);
    }

    if (window.jQuery && window.Shiny) {
      jQuery(document).on("shiny:connected", function () {
        connected = true;
        pushState();
      });
      Shiny.addCustomMessageHandler("rd:apply", applyState);
      Shiny.addCustomMessageHandler("rd:url", receiveBookmarkUrl);
    }

    /** Show a short-lived message in the corner of the stage. */
    function toast(message) {
      var node = document.createElement("div");
      node.className = "rd-toast";
      node.textContent = message;
      document.body.appendChild(node);
      requestAnimationFrame(function () { node.classList.add("is-visible"); });
      setTimeout(function () {
        node.classList.remove("is-visible");
        setTimeout(function () { node.remove(); }, 400);
      }, 2400);
    }

    /** The server answered a bookmark request with a shareable URL. */
    function receiveBookmarkUrl(url) {
      try {
        window.history.replaceState(null, "", url);
      } catch (e) { /* cross-origin or file:// -- the clipboard still works */ }

      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(url).then(
          function () { toast("Link copied to clipboard"); },
          function () { toast("Link is in the address bar"); }
        );
      } else {
        toast("Link is in the address bar");
      }
    }

    /** Restore a state handed back by the server (bookmark restore). */
    function applyState(incoming) {
      if (!incoming) return;
      if (incoming.palette && LUTS[incoming.palette]) {
        state.palette = incoming.palette;
        els.palette.value = incoming.palette;
      }
      if (incoming.speed) {
        state.speed = clamp(parseInt(incoming.speed, 10), 1, 10);
        els.speed.value = state.speed;
      }
      if (typeof incoming.running === "boolean") setRunning(incoming.running);
      if (incoming.f != null && incoming.k != null) {
        setParameters(parseFloat(incoming.f), parseFloat(incoming.k));
      }
      syncLabels();
    }

    // --- parameter plumbing ----------------------------------------------

    function setParameters(f, k) {
      state.f = clamp(f, F_MIN, F_MAX);
      state.k = clamp(k, K_MIN, K_MAX);
      els.feed.value = state.f;
      els.kill.value = state.k;
      map.setValue(state.f, state.k);
      syncLabels();
      pushState();
    }

    function syncLabels() {
      setLabel("rd-feed", state.f.toFixed(4));
      setLabel("rd-kill", state.k.toFixed(4));
      setLabel("rd-speed", String(state.speed));
      setLabel("rd-brush", String(state.brush));
    }

    function setLabel(id, text) {
      var node = document.getElementById(id + "-value");
      if (node) node.textContent = text;
    }

    function setRunning(running) {
      state.running = running;
      var button = document.querySelector('[data-action="toggle"]');
      if (button) button.textContent = running ? "Pause" : "Play";
      pushState();
    }

    // --- painting on the canvas ------------------------------------------

    var painting = false;

    function paintAt(e) {
      var rect = canvas.getBoundingClientRect();
      var gx = (e.clientX - rect.left) / rect.width * GRID;
      var gy = (e.clientY - rect.top) / rect.height * GRID;
      reactor.seed(gx, gy, state.brush, 0.9);
    }

    canvas.addEventListener("pointerdown", function (e) {
      if (e.button !== 0) return;
      painting = true;
      canvas.setPointerCapture(e.pointerId);
      paintAt(e);
      e.preventDefault();
    });
    canvas.addEventListener("pointermove", function (e) {
      if (!painting) return;
      if (e.buttons === 0) { painting = false; return; }
      paintAt(e);
    });
    ["pointerup", "pointercancel", "lostpointercapture"].forEach(function (type) {
      canvas.addEventListener(type, function () { painting = false; });
    });

    // --- controls ---------------------------------------------------------

    els.feed.addEventListener("input", function () {
      setParameters(parseFloat(this.value), state.k);
    });
    els.kill.addEventListener("input", function () {
      setParameters(state.f, parseFloat(this.value));
    });
    els.speed.addEventListener("input", function () {
      state.speed = parseInt(this.value, 10);
      syncLabels();
      pushState();
    });
    els.brush.addEventListener("input", function () {
      state.brush = parseInt(this.value, 10);
      syncLabels();
    });

    Object.keys(PALETTES).forEach(function (name) {
      var option = document.createElement("option");
      option.value = name;
      option.textContent = name;
      els.palette.appendChild(option);
    });
    els.palette.value = state.palette;
    els.palette.addEventListener("change", function () {
      state.palette = this.value;
      render();
      pushState();
    });

    document.querySelectorAll(".rd-preset").forEach(function (button) {
      button.addEventListener("click", function () {
        setParameters(parseFloat(this.dataset.f), parseFloat(this.dataset.k));
        reactor.seedCentre();
      });
    });

    var actions = {
      toggle: function () { setRunning(!state.running); },
      reset: function () { reactor.seedCentre(); render(); },
      clear: function () { reactor.clear(); render(); },
      random: function () { randomParameters(); },
      snapshot: function () { saveSnapshot(); },
      share: function () {
        if (connected) sendInput("do_bookmark", Date.now());
        else toast("Not connected to the server");
      }
    };

    document.querySelectorAll("[data-action]").forEach(function (button) {
      var handler = actions[button.dataset.action];
      if (handler) button.addEventListener("click", handler);
    });

    /** Jump to a random point that is actually inside the pattern region. */
    function randomParameters() {
      var f, k, limit, guard = 0;
      do {
        f = 0.014 + Math.random() * 0.052;
        limit = killLimit(f);
        k = Math.max(K_MIN, limit - 0.012) + Math.random() * 0.012;
        guard++;
      } while (k > limit && guard < 50);
      setParameters(f, clamp(k, K_MIN, K_MAX));
      reactor.randomize(8);
    }

    function saveSnapshot() {
      var out = document.createElement("canvas");
      out.width = out.height = SNAPSHOT_PX;
      var ctx = out.getContext("2d");
      ctx.imageSmoothingEnabled = true;
      ctx.imageSmoothingQuality = "high";
      ctx.drawImage(buffer, 0, 0, SNAPSHOT_PX, SNAPSHOT_PX);

      var link = document.createElement("a");
      link.download = "reaction-diffusion-f" + state.f.toFixed(4) +
                      "-k" + state.k.toFixed(4) + ".png";
      link.href = out.toDataURL("image/png");
      link.click();
    }

    // --- render loop -------------------------------------------------------

    function render() {
      reactor.paint(image, LUTS[state.palette]);
      bufferCtx.putImageData(image, 0, 0);
      if (!view) return;
      view.ctx.imageSmoothingEnabled = true;
      view.ctx.imageSmoothingQuality = "high";
      view.ctx.clearRect(0, 0, view.w, view.h);
      view.ctx.drawImage(buffer, 0, 0, view.w, view.h);
    }

    var lastActivity = 0;

    function frame(now) {
      if (state.running) {
        for (var i = 0; i < state.speed; i++) reactor.step(state.f, state.k);
      }
      render();

      if (now - lastActivity > ACTIVITY_INTERVAL_MS) {
        lastActivity = now;
        sendInput("activity", reactor.activity());
      }
      requestAnimationFrame(frame);
    }

    // --- layout ------------------------------------------------------------

    function resize() {
      var stage = canvas.parentElement;
      var side = Math.max(120, Math.min(stage.clientWidth, stage.clientHeight));
      canvas.style.width = side + "px";
      canvas.style.height = side + "px";
      view = fitCanvas(canvas, side);
      map.resize();
      render();
    }

    window.addEventListener("resize", resize);
    document.addEventListener("visibilitychange", function () {
      // Nothing to do when hidden: rAF is already throttled by the browser,
      // but stop pretending we are measuring anything.
      if (!document.hidden) lastActivity = 0;
    });

    map.setValue(state.f, state.k);
    syncLabels();
    resize();
    requestAnimationFrame(frame);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
