# Reaction-Diffusion Lab -------------------------------------------------
#
# An interactive Gray-Scott reaction-diffusion sandbox.
#
#   du/dt = Du * lap(u) - u*v^2 + f*(1 - u)
#   dv/dt = Dv * lap(v) + u*v^2 - (f + k)*v
#
# The solver itself lives in www/rd.js and runs on the client, so dragging a
# slider never waits on the server. R owns the things R is good at: the preset
# table, the stability arithmetic, the live readout and bookmarking.
#
# Run with:  shiny::runApp()

library(shiny)

enableBookmarking("url")

# Parameter-map extent. Must match F_MIN/F_MAX/K_MIN/K_MAX in www/rd.js.
F_RANGE <- c(0.000, 0.100)
K_RANGE <- c(0.030, 0.075)

# Activity sparkline keeps this many samples (arriving at ~4 Hz).
HISTORY_LEN <- 180L

# Named regimes. This table is the single source of truth: it feeds both the
# preset buttons and the dots drawn on the parameter map.
#
# The values are not the ones usually copied around the web. Several of those
# (mitosis at k = 0.0649, solitons at k = 0.062) sit above the existence curve
# for this discretisation and simply fade out here, so every point below was
# chosen by sweeping the f/k plane, ranking by spatial structure rather than
# by brightness, and eyeballing the result. See in_pattern_zone() for why the
# interesting band hugs the curve from underneath.
PRESETS <- data.frame(
  name  = c("Amoeba", "Cells", "Mitosis", "Holes", "Fingerprint", "Lattice",
            "Maze", "Rings", "Worms", "Coral", "Loops"),
  f     = c(0.0140, 0.0180, 0.0180, 0.0220, 0.0220, 0.0300,
            0.0300, 0.0420, 0.0460, 0.0540, 0.0620),
  k     = c(0.0392, 0.0461, 0.0481, 0.0492, 0.0522, 0.0546,
            0.0566, 0.0595, 0.0612, 0.0622, 0.0615),
  blurb = c("Soft blobs that creep, merge and pull apart",
            "Bright clumps budding against a dark field",
            "Isolated spots that swell and pinch in two",
            "A solid field pitted with spreading pores",
            "Ridges that flow and fork like a thumbprint",
            "Dark dots settling into a regular grid",
            "Corridors that fill the plane and lock solid",
            "Concentric bands rippling out from the seed",
            "Wandering filaments that never quite touch",
            "Branching fronts crowding out into a reef",
            "Broad ribbons curling into closed loops"),
  stringsAsFactors = FALSE
)

# --- model arithmetic ----------------------------------------------------

#' Largest kill rate that still admits a non-trivial steady state.
#'
#' The uniform state (u = 1, v = 0) always exists. A second fixed point needs
#' f >= 4*(f + k)^2, i.e. k <= sqrt(f)/2 - f.
kill_limit <- function(f) {
  sqrt(f) / 2 - f
}

# How close to the curve counts as "on the edge" rather than clearly one side.
EDGE_MARGIN <- 0.002

#' Classify a point relative to the steady-state boundary.
#'
#' Deliberately three-valued. The boundary governs *steady states*, and the
#' most interesting patterns live in a narrow band hugging it -- some of them
#' just above it, persisting as travelling structures with no steady state to
#' settle into. Calling everything above the curve dead would be wrong.
#'
#' @return One of "stable", "edge" or "dying".
regime_class <- function(f, k) {
  slack <- kill_limit(f) - k
  if (slack > EDGE_MARGIN) "stable" else if (slack > -EDGE_MARGIN) "edge" else "dying"
}

VERDICTS <- list(
  stable = list(label = "Self-sustaining", class = "is-alive"),
  edge   = list(label = "On the knife edge", class = "is-edge"),
  dying  = list(label = "Fades to nothing", class = "is-dead")
)

#' Index of the preset closest to a point in (f, k) space.
nearest_preset <- function(f, k) {
  which.min((PRESETS$f - f)^2 + (PRESETS$k - k)^2)
}

# --- UI helpers ----------------------------------------------------------

#' A labelled range slider whose value is read directly by rd.js.
#'
#' Deliberately a plain <input type="range"> rather than sliderInput(): the
#' simulation reads it locally on every frame, with no server round-trip.
range_control <- function(id, label, min, max, step, value, hint = NULL) {
  div(
    class = "rd-control",
    tags$label(
      `for` = id,
      span(class = "rd-control-name", label),
      span(class = "rd-control-value", id = paste0(id, "-value"), value)
    ),
    tags$input(
      type = "range", class = "rd-range", id = id,
      # format() keeps small steps out of scientific notation, which not every
      # browser accepts in a step attribute.
      min = format(min, scientific = FALSE),
      max = format(max, scientific = FALSE),
      step = format(step, scientific = FALSE),
      value = format(value, scientific = FALSE)
    ),
    if (!is.null(hint)) div(class = "rd-hint", hint)
  )
}

#' One preset button, carrying its coordinates as data attributes.
preset_button <- function(row) {
  tags$button(
    type = "button", class = "rd-preset",
    `data-f` = row$f, `data-k` = row$k, title = row$blurb,
    row$name
  )
}

#' A plain toolbar button; rd.js dispatches on the data-action value.
tool_button <- function(action, label, primary = FALSE) {
  tags$button(
    type = "button",
    class = if (primary) "rd-button rd-button-primary" else "rd-button",
    `data-action` = action,
    label
  )
}

# --- UI ------------------------------------------------------------------

ui <- function(request) {
  tagList(
    tags$head(
      tags$title("Reaction-Diffusion Lab"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$link(rel = "stylesheet", href = "styles.css"),
      # Presets are embedded at render time so rd.js can draw the parameter
      # map on its first frame, with no message round-trip to wait for.
      tags$script(
        type = "application/json", id = "rd-presets",
        HTML(jsonlite::toJSON(PRESETS, dataframe = "rows", digits = 8))
      ),
      tags$script(src = "rd.js", defer = NA)
    ),

    div(
      class = "rd-app",

      tags$header(
        class = "rd-header",
        div(
          h1("Reaction-Diffusion Lab"),
          p(class = "rd-tagline",
            "Two chemicals, four numbers, and every pattern on an animal's back.")
        ),
        div(
          class = "rd-header-actions",
          tool_button("toggle", "Pause", primary = TRUE),
          tool_button("reset", "Reset"),
          tool_button("clear", "Clear"),
          tool_button("random", "Surprise me"),
          tool_button("snapshot", "Save PNG"),
          tool_button("share", "Copy link"),
          # Label and pressed state are owned by rd.js, which knows the theme
          # before the server does (it restores it from localStorage).
          tool_button("theme", "Day")
        )
      ),

      tags$main(
        class = "rd-main",

        tags$section(
          class = "rd-stage",
          tags$canvas(id = "rd-canvas", `aria-label` = "Reaction-diffusion simulation"),
          div(class = "rd-stage-hint", "Drag on the canvas to pour in chemical B")
        ),

        tags$aside(
          class = "rd-panel",

          div(
            class = "rd-card",
            h2("Parameter map"),
            p(class = "rd-hint",
              "Click anywhere. The lit region is where the reaction has a steady state to settle into. The best patterns hug its upper edge."),
            tags$canvas(id = "rd-map", `aria-label` = "Feed and kill rate map"),
            range_control("rd-feed", "Feed rate f", F_RANGE[1], F_RANGE[2], 0.0001, 0.0540),
            range_control("rd-kill", "Kill rate k", K_RANGE[1], K_RANGE[2], 0.0001, 0.0622)
          ),

          div(
            class = "rd-card",
            h2("Regime"),
            uiOutput("readout"),
            plotOutput("activity", height = "58px")
          ),

          div(
            class = "rd-card",
            h2("Presets"),
            div(
              class = "rd-presets",
              lapply(seq_len(nrow(PRESETS)), function(i) preset_button(PRESETS[i, ]))
            )
          ),

          div(
            class = "rd-card",
            h2("Rendering"),
            range_control("rd-speed", "Steps per frame", 1, 10, 1, 4,
                          "How much simulated time passes each frame."),
            range_control("rd-brush", "Brush size", 2, 30, 1, 8),
            div(
              class = "rd-control",
              tags$label(`for` = "rd-palette", span(class = "rd-control-name", "Palette")),
              tags$select(id = "rd-palette", class = "rd-select")
            )
          )
        )
      ),

      tags$footer(
        class = "rd-footer",
        "Gray-Scott model, solved on a 256x256 torus in your browser."
      )
    )
  )
}

# --- server --------------------------------------------------------------

server <- function(input, output, session) {

  # Latest simulation state, pushed from rd.js (debounced client-side).
  sim <- reactive({
    state <- input$sim_state
    if (is.null(state)) list(f = 0.0540, k = 0.0622) else state
  })

  # Rolling mean concentration of chemical B, sampled by the client at ~4 Hz.
  history <- reactiveVal(numeric(0))

  observeEvent(input$activity, {
    history(utils::tail(c(history(), as.numeric(input$activity)), HISTORY_LEN))
  })

  # The sparkline redraws far less often than samples arrive.
  history_throttled <- throttle(reactive(history()), 500)

  # The plot is drawn server-side, so it cannot pick up the CSS variables the
  # rest of the UI themes itself with; rd.js reports the active theme instead.
  SPARK_COLOURS <- list(
    dark  = list(fill = "#1d3a4d", line = "#38bdf8", tip = "#7dd3fc", label = "#5a6472"),
    light = list(fill = "#cbe2f4", line = "#0369a1", tip = "#0284c7", label = "#64748b")
  )

  spark_colours <- reactive({
    SPARK_COLOURS[[if (identical(input$theme, "light")) "light" else "dark"]]
  })

  output$readout <- renderUI({
    state <- sim()
    f <- state$f
    k <- state$k
    preset <- PRESETS[nearest_preset(f, k), ]
    verdict <- VERDICTS[[regime_class(f, k)]]
    limit <- kill_limit(f)

    tagList(
      div(
        class = "rd-readout-row",
        span(class = "rd-readout-key", "f"),
        span(class = "rd-readout-num", sprintf("%.4f", f)),
        span(class = "rd-readout-key", "k"),
        span(class = "rd-readout-num", sprintf("%.4f", k))
      ),
      div(
        class = paste("rd-verdict", verdict$class),
        verdict$label
      ),
      div(
        class = "rd-hint",
        sprintf("A steady state exists while k ≤ √f/2 − f = %.4f.", limit)
      ),
      div(
        class = "rd-nearest",
        strong(preset$name), span(class = "rd-hint", preset$blurb)
      )
    )
  })

  output$activity <- renderPlot(
    {
      h <- history_throttled()
      col <- spark_colours()
      op <- par(mar = c(0, 0, 0, 0), bg = NA)
      on.exit(par(op), add = TRUE)

      if (length(h) < 2) {
        plot.new()
        text(0.5, 0.5, "measuring…", col = col$label, cex = 1.1)
        return(invisible(NULL))
      }

      span <- range(h)
      pad <- max(diff(span) * 0.15, 1e-4)
      ylim <- c(span[1] - pad, span[2] + pad)
      x <- seq_along(h)

      plot(x, h, type = "n", axes = FALSE, xlab = "", ylab = "",
           xaxs = "i", yaxs = "i", ylim = ylim)
      polygon(c(x, rev(x)), c(h, rep(ylim[1], length(h))),
              col = col$fill, border = NA)
      lines(x, h, col = col$line, lwd = 2)
      points(length(h), h[length(h)], col = col$tip, pch = 19, cex = 0.9)
      invisible(NULL)
    },
    bg = "transparent",
    res = 96
  )

  # Bookmarking: the interesting state lives in the browser, so hand it over
  # on save and push it back on restore. The raw inputs are excluded -- only
  # the compact snapshot below belongs in a shareable URL.
  setBookmarkExclude(c("activity", "sim_state", "do_bookmark", "theme"))

  observeEvent(input$do_bookmark, session$doBookmark())

  onBookmark(function(state) {
    state$values$sim <- input$sim_state
  })

  # Defining this callback replaces Shiny's default modal, so the link can be
  # delivered straight to the clipboard instead.
  onBookmarked(function(url) {
    session$sendCustomMessage("rd:url", url)
  })

  onRestored(function(state) {
    if (!is.null(state$values$sim)) {
      session$sendCustomMessage("rd:apply", state$values$sim)
    }
  })
}

shinyApp(ui, server)
