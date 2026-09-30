library(shiny)
library(bslib)

# Board geometry --------------------------------------------------------------
# Everything is in mm, with the board centre at (0, 0).

R <- 170 # outer edge of the double ring; heatmap grid is -R:R in both directions
RINGS <- c(bull = 6.35, outer_bull = 15.9, triple_in = 99, triple_out = 107,
           double_in = 162, double_out = 170)

# segments counter-clockwise, starting at the right (6 is centred on 0 degrees)
SEGMENTS <- c(6, 13, 4, 18, 1, 20, 5, 12, 9, 14, 11, 8, 16, 7, 19, 3, 17, 2, 15, 10)
AIM_PRESETS <- list(T20 = c(0, 103), Bull = c(0, 0))

# Colours ---------------------------------------------------------------------

COL <- list(bg = "#1c1d22", sidebar = "#24252b", card = "#2a2b31", panel = "#33343b",
            border = "#43444d", fg = "#f1ede4", muted = "#a5a5ae", faint = "#7d7d87",
            cream = "#e6cda5", black = "#141414", red = "#df2623", green = "#11a551")
HEAT_PAL <- hcl.colors(64, "Inferno")

# Scoring ---------------------------------------------------------------------

# vectorised score of darts landing at (x, y)
scoreMM = function(x, y){
  r <- sqrt(x^2 + y^2)
  seg <- SEGMENTS[floor(((atan2(y, x) * 180 / pi + 9) %% 360) / 18) + 1]
  mult <- ifelse(r > RINGS[["double_out"]], 0,
          ifelse(r > RINGS[["double_in"]], 2,
          ifelse(r > RINGS[["triple_out"]], 1,
          ifelse(r > RINGS[["triple_in"]], 3, 1))))
  s <- seg * mult
  s[r <= RINGS[["outer_bull"]]] <- 25
  s[r <= RINGS[["bull"]]] <- 50
  s
}

# label of a single dart, e.g. "T20", "D16", "5", "Bull"
labelThrow = function(x, y){
  r <- sqrt(x^2 + y^2)
  s <- scoreMM(x, y)
  if(r <= RINGS[["bull"]]) return("Bull")
  if(r <= RINGS[["outer_bull"]] || s == 0) return(as.character(s))
  if(r > RINGS[["double_in"]]) return(paste0("D", s / 2))
  if(r > RINGS[["triple_in"]] && r <= RINGS[["triple_out"]]) return(paste0("T", s / 3))
  as.character(s)
}

# Expected scores -------------------------------------------------------------
# Darts aimed at a land at a + bias + N(0, sigma). The expected score of aiming
# at a is E(a) = sum_z score(z) phi(z - a - bias), i.e. a convolution of the
# score matrix with the (mirrored, shifted) throw density (Tibshirani et al.,
# 2010, plus a bias term). Done as a circular convolution via FFT on a
# 900 x 900 mm grid (900 = 2^2 3^2 5^2 keeps the FFT fast, and leaves enough
# room that the shifted density does not wrap around). The FFT of the score
# matrix is the same for everybody, so it is computed once at startup. The
# density is stored in wrapped order (offsets 0, 1, ..., -1), so the result
# needs no shifting.

N_GRID <- 900
MAX_BIAS <- 150                                     # mm per axis, keeps wrap-around negligible
GRID <- seq_len(N_GRID) - 1 - N_GRID / 2           # -450:449, index <-> coordinate
WRAP <- c(0:(N_GRID / 2 - 1), -(N_GRID / 2):-1)     # offsets in FFT order
GX <- matrix(WRAP, N_GRID, N_GRID)
GY <- t(GX)
SCORE_FFT <- fft(outer(GRID, GRID, scoreMM))
CROP <- match(-R:R, GRID)
INSIDE <- outer(-R:R, -R:R, function(x, y) x^2 + y^2 <= R^2)

# expected three-dart score for every aim point on the -R:R mm grid
expScores = function(sigma, bias = c(0, 0)){
  b <- pmin(pmax(bias, -MAX_BIAS), MAX_BIAS)
  d <- det(sigma)
  # phi is symmetric, so phi(z - a - b) = psi(a - z) with psi(u) = phi(u + b)
  X <- GX + b[1]
  Y <- GY + b[2]
  B <- exp(-(sigma[2,2] * X^2 - 2 * sigma[1,2] * X * Y + sigma[1,1] * Y^2) / (2 * d)) /
    (2 * pi * sqrt(d))
  E <- Re(fft(SCORE_FFT * fft(B), inverse = TRUE)) / N_GRID^2
  3 * E[CROP, CROP]
}

# value of an expected-score matrix at a point (mm)
atPoint = function(E, p) E[match(round(p[1]), -R:R), match(round(p[2]), -R:R)]

# covariance (mm^2) of the located throws, NULL if not (yet) estimable
estimateSigma = function(throws){
  xy <- as.matrix(throws[throws$located, c("x", "y")])
  if(nrow(xy) < 3) return(NULL)
  sigma <- cov(xy)
  if(det(sigma) < 1e-6) return(NULL)
  sigma
}

# Drawing ---------------------------------------------------------------------

wedge = function(a1, a2, r1, r2){
  t <- seq(a1, a2, length.out = 24) * pi / 180
  list(x = c(r2 * cos(t), r1 * cos(rev(t)), NA), y = c(r2 * sin(t), r1 * sin(rev(t)), NA))
}

circle = function(r, n = 200){
  t <- seq(0, 2 * pi, length.out = n)
  list(x = r * cos(t), y = r * sin(t))
}

# all board polygons, built once: 20 segments x 4 rings + bulls
BOARD <- local({
  rings <- list(c(RINGS[["outer_bull"]], RINGS[["triple_in"]]), c(RINGS[["triple_in"]], RINGS[["triple_out"]]),
                c(RINGS[["triple_out"]], RINGS[["double_in"]]), c(RINGS[["double_in"]], RINGS[["double_out"]]))
  x <- y <- col <- NULL
  for(i in seq_along(SEGMENTS)){
    a <- (i - 1) * 18
    single <- if(i %% 2 == 1) COL$cream else COL$black
    ring <- if(i %% 2 == 1) COL$green else COL$red
    for(k in 1:4){
      w <- wedge(a - 9, a + 9, rings[[k]][1], rings[[k]][2])
      x <- c(x, w$x); y <- c(y, w$y)
      col <- c(col, if(k %% 2 == 1) single else ring)
    }
  }
  list(x = x, y = y, col = col)
})

newCanvas = function(){
  par(mar = c(0, 0, 0, 0), bg = "transparent")
  plot.new()
  plot.window(xlim = c(-200, 200), ylim = c(-200, 200), asp = 1, xaxs = "i", yaxs = "i")
  polygon(circle(200), col = "#0b0b0c", border = "#55565f", lwd = 1.5)
}

drawNumbers = function(s){
  a <- (seq_along(SEGMENTS) - 1) * 18 * pi / 180
  text(185 * cos(a), 185 * sin(a), SEGMENTS, col = "#d8d3c8", font = 2, cex = 1.35 * s)
}

drawBoard = function(s){
  newCanvas()
  polygon(BOARD$x, BOARD$y, col = BOARD$col, border = "#9a9a9a", lwd = 0.6 * s)
  polygon(circle(RINGS[["outer_bull"]]), col = COL$green, border = "#9a9a9a", lwd = 0.6 * s)
  polygon(circle(RINGS[["bull"]]), col = COL$red, border = "#9a9a9a", lwd = 0.6 * s)
  drawNumbers(s)
}

# thin wire overlay for the heatmap
drawWires = function(s, alpha = 0.35){
  wire <- adjustcolor("white", alpha)
  for(r in RINGS) lines(circle(r), col = wire, lwd = s)
  a <- ((seq_along(SEGMENTS) - 1) * 18 + 9) * pi / 180
  segments(RINGS[["outer_bull"]] * cos(a), RINGS[["outer_bull"]] * sin(a),
           R * cos(a), R * sin(a), col = wire, lwd = s)
}

# white scope (ring + ticks, black halo) with a red centre dot
drawTarget = function(p, s){
  r <- 9    # ring radius (mm)
  tick <- 17 # outer end of the ticks (mm)
  ring <- circle(r, 60)
  for(k in list(list(col = "black", lwd = 6), list(col = "white", lwd = 2.6))){
    lines(p[1] + ring$x, p[2] + ring$y, col = k$col, lwd = k$lwd * s)
    segments(p[1] + c(-tick, r, 0, 0), p[2] + c(0, 0, -tick, r),
             p[1] + c(-r, tick, 0, 0), p[2] + c(0, 0, -r, tick), col = k$col, lwd = k$lwd * s, lend = 1)
  }
  points(p[1], p[2], pch = 21, cex = 1.1 * s, lwd = 1.5 * s, col = "white", bg = COL$red)
}

drawAim = function(p, s){
  points(p[1], p[2], pch = 1, cex = 2.8 * s, lwd = 7 * s, col = "black")
  points(p[1], p[2], pch = 1, cex = 2.8 * s, lwd = 3.5 * s, col = COL$cream)
  points(p[1], p[2], pch = 21, cex = 0.9 * s, lwd = 1.5 * s, col = "black", bg = COL$cream)
}

# arrow from the intended aim to where the darts land on average
drawBias = function(aim, mu, s){
  if(sqrt(sum((mu - aim)^2)) < 4) return()
  arrows(aim[1], aim[2], mu[1], mu[2], length = 0.12 * s, lwd = 5 * s, col = "black")
  arrows(aim[1], aim[2], mu[1], mu[2], length = 0.12 * s, lwd = 2.5 * s, col = COL$cream)
}

# e.g. "14 mm right and 9 mm below"
describeBias = function(b){
  parts <- c(if(abs(b[1]) >= 1) sprintf("%.0f mm %s", abs(b[1]), if(b[1] > 0) "right" else "left"),
             if(abs(b[2]) >= 1) sprintf("%.0f mm %s", abs(b[2]), if(b[2] > 0) "above" else "below"))
  if(length(parts) == 0) return("Your darts land centred on your aim.")
  paste0("Your darts land ", paste(parts, collapse = " and "), " your aim on average.")
}

drawThrows = function(thr, s){
  points(thr$x, thr$y, pch = 4, cex = 1.3 * s, lwd = 5 * s, col = "black")
  points(thr$x, thr$y, pch = 4, cex = 1.3 * s, lwd = 2.5 * s, col = "white")
}

drawSpread = function(sigma, mu, s){
  L <- t(chol(sigma))
  t <- seq(0, 2 * pi, length.out = 100)
  for(level in c(0.3, 0.6, 0.9)){
    e <- sqrt(qchisq(level, 2)) * L %*% rbind(cos(t), sin(t))
    lines(mu[1] + e[1, ], mu[2] + e[2, ], col = adjustcolor("white", 0.85), lwd = 1.8 * s)
  }
  points(mu[1], mu[2], pch = 16, cex = 1.1 * s, col = "white")
}

emptyThrows = function(){
  data.frame(round = integer(0), x = numeric(0), y = numeric(0),
             located = logical(0), score = numeric(0), label = character(0))
}

# UI --------------------------------------------------------------------------

theme <- bs_theme(
  version = 5, bg = COL$bg, fg = COL$fg, primary = COL$red, success = COL$green,
  secondary = COL$panel, border_color = COL$border,
  base_font = font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "sans-serif"),
  heading_font = font_collection("Barlow Condensed", "Arial Narrow", "sans-serif"),
  "card-bg" = COL$card, "card-border-color" = COL$border,
  "input-bg" = COL$panel, "input-border-color" = COL$border
)

css <- paste0(":root {", paste0("--", names(COL), ": ", COL, ";", collapse = " "), "}", "
  body { letter-spacing: 0.01em; }
  .bslib-sidebar-layout > .sidebar { background: var(--sidebar); border-right: 1px solid var(--border); }
  .navbar, .bslib-page-title { background: var(--sidebar) !important; border-bottom: 1px solid var(--border); }
  .app-title { font-family: 'Barlow Condensed', sans-serif; font-weight: 600; font-size: 1.6rem;
               text-transform: uppercase; letter-spacing: 0.08em; display: flex; align-items: center; gap: .6rem; }
  .app-title .dot { width: .8rem; height: .8rem; border-radius: 50%; background: var(--red);
                    box-shadow: 0 0 0 3px var(--green); }
  .app-title .light { color: var(--muted); font-weight: 400; }
  .num { font-family: 'Barlow Condensed', sans-serif; font-variant-numeric: tabular-nums; }
  .label { font-size: .72rem; text-transform: uppercase; letter-spacing: .12em; color: var(--muted); }

  .status { font-family: 'Barlow Condensed', sans-serif; font-size: 1.7rem; font-weight: 600; line-height: 1.1; }
  .scoreboard { display: flex; flex-direction: column; gap: .75rem; }
  #controls > * { margin-bottom: .75rem; }
  .btn-row { display: flex; gap: .5rem; }
  .btn-row .btn { flex: 1; }
  .btn-ghost { background: var(--panel); border: 1px solid var(--border); color: var(--fg); }
  .btn-ghost:hover { background: var(--border); color: var(--fg); border-color: var(--border); }
  .restart { color: var(--muted); text-decoration: none; }
  .restart:hover { color: var(--red); }
  .locked-n { display: flex; justify-content: space-between; align-items: baseline; }
  .locked-n .num { font-size: 1.4rem; }
  .locked-n + .locked-n { margin-top: .25rem; }

  .score-card { border: 1px solid var(--border); border-radius: .75rem; padding: .8rem 1rem; background: var(--panel);
                box-shadow: 0 1px 2px rgba(0,0,0,.25); }
  .score-head { display: flex; justify-content: space-between; align-items: center; }
  .score-main { font-size: 3.2rem; font-weight: 600; color: var(--cream); line-height: 1; margin-top: .3rem; }
  .score-main .unit { font-family: Inter, system-ui, sans-serif; font-size: .8rem; color: var(--muted);
                      font-weight: 400; margin-left: .4rem; letter-spacing: 0; }
  .score-main .empty { color: var(--faint); }
  .score-meta { font-size: .8rem; color: var(--muted); margin-top: .35rem; }
  .score-throws { font-size: .75rem; color: var(--faint); margin-top: .35rem; word-wrap: break-word; }
  .delta { font-family: 'Barlow Condensed', sans-serif; font-size: 1rem; font-weight: 600;
           padding: .05rem .55rem; border-radius: 999px; }
  .delta.pos { background: rgba(17,165,81,.22); color: #4fdb8c; }
  .delta.neg { background: rgba(223,38,35,.22); color: #ff7b78; }
  .exp-row { display: flex; justify-content: space-between; font-size: .85rem; margin-top: .3rem; }
  .exp-row .num { font-size: 1.15rem; }

  .stepper { display: flex; gap: .5rem; flex-wrap: wrap; margin-bottom: 1rem; }
  .step { border: 1px solid var(--border); background: var(--card); color: var(--muted); border-radius: 999px;
          padding: .3rem .9rem; font-size: .85rem; }
  .step b { font-family: 'Barlow Condensed', sans-serif; margin-right: .35rem; }
  .step.active { background: var(--red); border-color: var(--red); color: white; }
  .step.done { border-color: var(--green); color: #4fdb8c; }

  .card { box-shadow: 0 2px 8px rgba(0,0,0,.25); }
  .card-header { background: transparent; border-bottom: 1px solid var(--border); font-family: 'Barlow Condensed', sans-serif;
                 font-size: 1.2rem; font-weight: 600; letter-spacing: .03em; }
  .card-header .hint { font-family: Inter, system-ui, sans-serif; font-size: .8rem; font-weight: 400; color: var(--muted); }
  .legend { display: flex; align-items: center; gap: .6rem; font-size: .8rem; color: var(--muted); margin-top: .5rem; }
  .legend .bar { flex: 1; height: .5rem; border-radius: 999px; }
  #board { cursor: crosshair; }
")

ui <- page_sidebar(
  title = div(class = "app-title", span(class = "dot"), "Darts", span(class = "light", "Optimiser")),
  window_title = "DartsOptimizeR - Classroom",
  theme = theme,
  fillable = FALSE,
  tags$head(
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=Barlow+Condensed:wght@400;600&family=Inter:wght@400;500&display=swap"),
    tags$style(HTML(css))
  ),
  sidebar = sidebar(
    width = 320,
    uiOutput("settings"),
    div(class = "status", textOutput("status", inline = TRUE)),
    div(class = "actions",
        uiOutput("controls"),
        # always in the page; its label and enabled state follow the phase
        actionButton("next_step", "Start round 1", icon = icon("play"), class = "btn-primary w-100 mb-2"),
        div(class = "btn-row",
            actionButton("undo", "Undo", icon = icon("rotate-left"), class = "btn-ghost btn-sm"),
            actionButton("miss", "Miss (0)", icon = icon("ban"), class = "btn-ghost btn-sm"))),
    uiOutput("scoreboard"),
    actionButton("reset", "Restart", icon = icon("arrows-rotate"), class = "btn-link btn-sm restart")
  ),
  uiOutput("stepper"),
  layout_columns(
    col_widths = c(6, 6),
    card(fill = FALSE, card_header(uiOutput("board_title", inline = TRUE)),
         card_body(fill = FALSE, plotOutput("board", click = "plot_click", height = "auto"))),
    card(fill = FALSE, card_header("Expected points per 3 darts ", span(class = "hint", "when aiming here (brighter = better)")),
         card_body(fill = FALSE, plotOutput("expScore", height = "auto"), uiOutput("legend")))
  )
)

# Server ----------------------------------------------------------------------

server <- function(input, output, session) {
  # phase: "aim" (set intended aim), "r1" (throwing round 1), "r1_done",
  #        "r2" (throwing round 2), "done"
  val <- reactiveValues(phase = "aim", n = 9, use_bias = TRUE, aim = AIM_PRESETS$T20,
                        throws = emptyThrows(), frozen = NULL)
  firstPhase = function() if(val$use_bias) "aim" else "r1"

  currentRound = reactive(if(val$phase %in% c("aim", "r1", "r1_done")) 1 else 2)
  roundThrows = function(k) val$throws[val$throws$round == k, ]

  addThrow = function(x, y, located){
    if(!(val$phase %in% c("r1", "r2"))) return()
    k <- currentRound()
    score <- if(located) scoreMM(x, y) else 0
    label <- if(located) labelThrow(x, y) else "Miss"
    val$throws <- rbind(val$throws,
                        data.frame(round = k, x = x, y = y, located = located,
                                   score = score, label = label))
    if(nrow(roundThrows(k)) >= val$n) val$phase <- if(k == 1) "r1_done" else "done"
  }

  # settings: editable until the first throw, then locked
  output$settings <- renderUI({
    if(nrow(val$throws) == 0){
      tagList(
        numericInput("n", "Darts per round", value = isolate(val$n), min = 3, max = 30, step = 1),
        input_switch("use_bias", "Account for aiming bias", value = isolate(val$use_bias)))
    } else {
      tagList(
        div(class = "locked-n", span(class = "label", "Darts per round"), span(class = "num", val$n)),
        div(class = "locked-n", span(class = "label", "Aiming bias"),
            span(class = "num", if(val$use_bias) "on" else "off")))
    }
  })

  observeEvent(input$n, {
    if(nrow(val$throws) == 0 && !is.na(input$n)) val$n <- min(max(round(input$n), 3), 30)
  })

  # without bias correction there is no aim to set, so skip that step
  observeEvent(input$use_bias, {
    if(nrow(val$throws) > 0) return()
    val$use_bias <- input$use_bias
    if(val$phase %in% c("aim", "r1")) val$phase <- firstPhase()
  })

  observeEvent(input$plot_click, {
    x <- input$plot_click$x
    y <- input$plot_click$y
    if(val$phase == "aim"){
      if(x^2 + y^2 <= R^2) val$aim <- c(x, y)
    } else {
      addThrow(x, y, located = TRUE)
    }
  })

  observeEvent(input$aim_t20, val$aim <- AIM_PRESETS$T20)
  observeEvent(input$aim_bull, val$aim <- AIM_PRESETS$Bull)
  observeEvent(input$miss, addThrow(NA, NA, located = FALSE))

  observeEvent(input$undo, {
    idx <- which(val$throws$round == currentRound())
    if(length(idx) == 0){
      if(val$phase == "r1") val$phase <- firstPhase() # back to adjusting the aim
      return()
    }
    val$throws <- val$throws[-max(idx), ]
    if(val$phase == "r1_done") val$phase <- "r1"
    if(val$phase == "done") val$phase <- "r2"
  })

  restart = function(){
    val$phase <- firstPhase()
    val$aim <- AIM_PRESETS$T20
    val$throws <- emptyThrows()
    val$frozen <- NULL
  }

  observeEvent(input$reset, restart())

  observeEvent(input$next_step, {
    if(val$phase == "aim"){
      val$phase <- "r1"
    } else if(val$phase == "r1_done"){
      h <- heat()
      req(!is.null(h))
      val$frozen <- h
      val$phase <- "r2"
    } else if(val$phase == "done"){
      restart()
    }
  })

  observe({
    ready <- val$phase %in% c("aim", "done") || (val$phase == "r1_done" && !is.null(heat()))
    label <- switch(val$phase, aim = "Start round 1", done = "Restart", "Start round 2")
    updateActionButton(session, "next_step", label = label, disabled = !ready,
                       icon = icon(if(val$phase == "done") "arrows-rotate" else "play"))
  })

  # live heatmap from round 1, frozen once round 2 starts
  heat <- reactive({
    if(!is.null(val$frozen)) return(val$frozen)
    thr <- roundThrows(1)
    sigma <- estimateSigma(thr)
    if(is.null(sigma)) return(NULL)
    aim <- if(val$use_bias) val$aim
    bias <- if(val$use_bias) unname(colMeans(thr[thr$located, c("x", "y")]) - aim) else c(0, 0)
    E <- expScores(sigma, bias)
    best <- which(E == max(E), arr.ind = TRUE)[1, ]
    list(E = E, target = unname(best) - (R + 1), aim = aim, bias = bias)
  })

  # Sidebar -------------------------------------------------------------------

  output$status <- renderText({
    switch(val$phase,
           aim = "Set your aim",
           r1 = sprintf("Dart %d of %d", nrow(roundThrows(1)) + 1, val$n),
           r1_done = "Round 1 complete",
           r2 = sprintf("Dart %d of %d", nrow(roundThrows(2)) + 1, val$n),
           done = "Finished")
  })

  output$controls <- renderUI({
    if(val$phase == "aim"){
      return(tagList(
        p(class = "label", "Click the board where you will aim, or pick:"),
        div(class = "btn-row",
            actionButton("aim_t20", "T20", class = "btn-ghost btn-sm"),
            actionButton("aim_bull", "Bull", class = "btn-ghost btn-sm"))))
    }
    if(val$phase == "r1_done" && is.null(heat())){
      p(class = "label", "Not enough spread for a heatmap - undo and add throws.")
    }
  })

  roundBox = function(k, title, badge = NULL){
    thr <- roundThrows(k)
    n <- nrow(thr)
    div(class = "score-card",
        div(class = "score-head", span(class = "label", title), badge),
        div(class = "score-main num",
            if(n > 0) sprintf("%.1f", 3 * mean(thr$score)) else span(class = "empty", "\u2013"),
            span(class = "unit", "per 3 darts")),
        div(class = "score-meta",
            sprintf("Total %d  \u00b7  %s per dart  \u00b7  %d/%d darts", sum(thr$score),
                    if(n > 0) sprintf("%.1f", mean(thr$score)) else "\u2013", n, val$n)),
        if(n > 0) div(class = "score-throws", paste(thr$label, collapse = "  ")))
  }

  round2Box = function(){
    if(!(val$phase %in% c("r2", "done"))) return()
    badge <- NULL
    if(val$phase == "done"){
      diff <- 3 * (mean(roundThrows(2)$score) - mean(roundThrows(1)$score))
      badge <- span(class = paste("delta", if(diff >= 0) "pos" else "neg"),
                    sprintf("%+.1f", diff))
    }
    roundBox(2, "Round 2 \u00b7 optimal aim", badge)
  }

  expectedBox = function(){
    h <- heat()
    if(is.null(h)) return()
    if(is.null(h$aim)){
      return(div(class = "score-card",
                 span(class = "label", "Expected per 3 darts"),
                 div(class = "exp-row", "Aiming at T20",
                     span(class = "num", sprintf("%.1f", atPoint(h$E, AIM_PRESETS$T20)))),
                 div(class = "exp-row", sprintf("Optimal aim (%s)", labelThrow(h$target[1], h$target[2])),
                     span(class = "num", sprintf("%.1f", max(h$E))))))
    }
    div(class = "score-card",
        span(class = "label", "Expected per 3 darts"),
        div(class = "exp-row", sprintf("Your aim (%s)", labelThrow(h$aim[1], h$aim[2])),
            span(class = "num", sprintf("%.1f", atPoint(h$E, h$aim)))),
        div(class = "exp-row", sprintf("Optimal aim (%s)", labelThrow(h$target[1], h$target[2])),
            span(class = "num", sprintf("%.1f", max(h$E)))),
        div(class = "score-meta", describeBias(h$bias),
            if(sqrt(sum(h$bias^2)) >= 4){
              land <- h$target + pmin(pmax(h$bias, -MAX_BIAS), MAX_BIAS)
              sprintf(" So aim at %s to land around %s.", labelThrow(h$target[1], h$target[2]),
                      labelThrow(land[1], land[2]))
            }))
  }

  # one output for all cards, so hidden cards leave no gaps in the sidebar
  output$scoreboard <- renderUI({
    div(class = "scoreboard",
        if(val$phase != "aim") roundBox(1, "Round 1 \u00b7 usual aim"), round2Box(), expectedBox())
  })

  # Main panel ----------------------------------------------------------------

  output$stepper <- renderUI({
    step <- switch(val$phase, aim = 0, r1 = 1, r1_done = 1, r2 = 2, done = 3) + val$use_bias
    labels <- c(if(val$use_bias) "Set aim", "Round 1 \u00b7 usual aim", "Round 2 \u00b7 optimal aim", "Result")
    div(class = "stepper", lapply(seq_along(labels), function(i){
      cls <- if(i == step) "step active" else if(i < step || val$phase == "done") "step done" else "step"
      div(class = cls, tags$b(i), labels[i])
    }))
  })

  output$board_title <- renderUI({
    hint <- switch(val$phase,
                   aim = "click where you will aim",
                   r1 = , r1_done = if(val$use_bias) "aim at your marked aim" else "aim where you usually aim",
                   sprintf("aim at the marked target (%s)", labelThrow(val$frozen$target[1], val$frozen$target[2])))
    tagList("Board ", span(class = "hint", hint))
  })

  plotSize = function(id) function() max(session$clientData[[paste0("output_", id, "_width")]], 100)
  plotScale = function(id) plotSize(id)() / 520

  output$board <- renderPlot({
    s <- plotScale("board")
    drawBoard(s)
    k <- currentRound()
    thr <- roundThrows(k)
    thr <- thr[thr$located, ]
    if(k == 1){
      sigma <- estimateSigma(thr)
      if(!is.null(sigma)) drawSpread(sigma, colMeans(thr[, c("x", "y")]), s)
      if(val$use_bias){
        drawAim(val$aim, s)
        if(!is.null(sigma)) drawBias(val$aim, colMeans(thr[, c("x", "y")]), s)
      }
    } else {
      drawTarget(val$frozen$target, s)
    }
    drawThrows(thr, s)
  }, height = plotSize("board"), width = plotSize("board"), bg = "transparent")

  output$expScore <- renderPlot({
    s <- plotScale("expScore")
    newCanvas()
    h <- heat()
    if(is.null(h)){
      drawWires(s, alpha = 0.12)
      drawNumbers(s)
      text(0, 0, "Heatmap appears after 3 darts", col = COL$muted, cex = 1.3 * s)
      return()
    }
    E <- h$E
    E[!INSIDE] <- NA
    image(-R:R, -R:R, E, col = HEAT_PAL, add = TRUE, useRaster = TRUE)
    drawWires(s)
    drawNumbers(s)
    if(!is.null(h$aim)) drawAim(h$aim, s)
    drawTarget(h$target, s)
  }, height = plotSize("expScore"), width = plotSize("expScore"), bg = "transparent")

  output$legend <- renderUI({
    h <- heat()
    if(is.null(h)) return()
    rng <- range(h$E[INSIDE])
    stops <- substr(HEAT_PAL[round(seq(1, length(HEAT_PAL), length.out = 8))], 1, 7)
    div(class = "legend",
        span(class = "num", sprintf("%.0f", rng[1])),
        div(class = "bar", style = sprintf("background: linear-gradient(90deg, %s);", paste(stops, collapse = ", "))),
        span(class = "num", sprintf("%.0f", rng[2])))
  })
}

# Run the application
shinyApp(ui = ui, server = server)
