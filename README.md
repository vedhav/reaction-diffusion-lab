# Reaction-Diffusion Lab

An interactive [Gray-Scott](https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system)
sandbox built with Shiny. Two chemicals diffuse and react on a 256x256 torus;
four numbers decide whether you get leopard spots, coral, mazes or nothing at
all.

```
du/dt = Du * lap(u) - u*v^2 + f*(1 - u)
dv/dt = Dv * lap(v) + u*v^2 - (f + k)*v
```

## Running

```r
shiny::runApp()
```

Requires only the `shiny` package.

## What to try

- **Drag on the canvas** to pour chemical B into the dish.
- **Click the parameter map.** The lit region is exactly where
  `k <= sqrt(f)/2 - f`, the condition for a non-trivial steady state to exist.
  The most interesting patterns hug its upper edge; well above the curve,
  whatever you paint fades away.
- **Preset chips** jump to named regimes. Their coordinates were picked by
  sweeping the plane and ranking by spatial structure, not copied from the
  values that circulate online -- several of those die in this discretisation.
- **Copy link** puts the current parameters in a shareable URL.

## Layout

| File             | Role                                                      |
| ---------------- | --------------------------------------------------------- |
| `app.R`          | UI, preset table, stability arithmetic, readout, bookmarks |
| `www/rd.js`      | The solver, canvas rendering, parameter map, Shiny bridge  |
| `www/styles.css` | Theme and layout                                           |

The simulation runs entirely in the browser so that interaction never waits on
the server; R is told about the state (parameters, and a mean-concentration
sample at 4 Hz) rather than driving it.
