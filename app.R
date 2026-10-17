# STA2ASM Week 3: Discrete distributions
# Standalone Shiny application. Requires the shiny package only.
library(shiny)

COL_BIN <- "#E9B900"
COL_POIS <- "#C65D7B"
COL_NB <- "#8E70B5"
COL_OBS <- "#236A6E"

choices <- c("Binomial" = "bin", "Poisson" = "pois",
             "NB (traditional)" = "nbt", "NB (alternative parameterisation)" = "nba")

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)
math <- function(s) withMathJax(HTML(paste0("\\(", s, "\\)")))

# Each specification uses integer support, a probability vector, and an
# appropriate random generator. NB traditional counts TOTAL trials.
make_spec <- function(dist, n = 20, p = .30, lambda = 2.6 / 12,
                      exposure = 12, r = 3, mu = 2.6, theta = 2) {
  if (dist == "bin") {
    n <- as.integer(n); p <- as.numeric(p)
    mn <- n * p; vr <- n * p * (1 - p)
    maxx <- n; xx <- 0:maxx
    pr <- dbinom(xx, n, p)
    gen <- function(m) rbinom(m, n, p)
    return(list(x = xx, pr = pr, mean = mn, var = vr, gen = gen,
                color = COL_BIN, xlabel = "Number of patients achieving remission",
                label = "Binomial", variable = "X", mean_formula = "np", var_formula = "np(1-p)", note = "Number of successes among \\(n\\) independent trials",
                latex = sprintf("X \\sim \\operatorname{Binomial}(%d, %.2f)", n, p)))
  }
  if (dist == "pois") {
    mn <- lambda * exposure; vr <- mn
    maxx <- max(15L, as.integer(qpois(.9995, mn)))
    xx <- 0:maxx; pr <- dpois(xx, mn)
    gen <- function(m) rpois(m, mn)
    return(list(x = xx, pr = pr, mean = mn, var = vr, gen = gen,
                color = COL_POIS, xlabel = "Number of malaria cases in the selected period",
                label = "Poisson", variable = "X", mean_formula = "\\lambda t", var_formula = "\\lambda t", note = "Expected count is rate \\(\\lambda\\) multiplied by exposure \\(t\\)",
                latex = sprintf("X \\sim \\operatorname{Poisson}(\\mu),\\quad \\mu=\\lambda t=%.3f", mn)))
  }
  if (dist == "nbt") {
    r <- as.integer(r); p <- as.numeric(p)
    mn <- r / p; vr <- r * (1 - p) / p^2
    maxx <- max(r + 20L, r + as.integer(qnbinom(.9995, size = r, prob = p)))
    xx <- r:maxx; pr <- dnbinom(xx - r, size = r, prob = p)
    gen <- function(m) rnbinom(m, size = r, prob = p) + r
    return(list(x = xx, pr = pr, mean = mn, var = vr, gen = gen,
                color = COL_NB, xlabel = "Number of patients treated until the rth remission",
                label = "NB (traditional)", variable = "X", mean_formula = "\\frac{r}{p}", var_formula = "\\frac{r(1-p)}{p^2}", note = "\\(X\\) counts total trials; R's <code>dnbinom()</code> counts failures, so use <code>x - r</code>",
                latex = sprintf("X \\sim \\operatorname{NB}(r=%d,\\ p=%.2f)", r, p)))
  }
  mn <- mu; vr <- mu + mu^2 / theta
  maxx <- max(20L, as.integer(qnbinom(.9995, size = theta, mu = mu)))
  xx <- 0:maxx; pr <- dnbinom(xx, size = theta, mu = mu)
  gen <- function(m) rnbinom(m, size = theta, mu = mu)
  list(x = xx, pr = pr, mean = mn, var = vr, gen = gen,
       color = COL_NB, xlabel = "Number of malaria cases in one year",
       label = "NB (alternative parameterisation)", variable = "F", mean_formula = "\\mu", var_formula = "\\mu+\\frac{\\mu^2}{\\theta}", note = "Here, \\(F\\) represents the number of malaria cases in a year. The negative binomial distribution models these event counts; they are not failures before a specified number of successes.",
       latex = sprintf("F \\sim \\operatorname{NB}(\\mu=%.2f,\\ \\theta=%.2f)", mu, theta))
}

pmf_plot <- function(spec, selected = NULL, comparison = NULL, main = "Probability mass function") {
  ymax <- max(spec$pr, if (!is.null(comparison)) comparison$pr else 0) * 1.16
  plot(spec$x, spec$pr, pch = 21, bg = spec$color, col = "black", cex = 2.4,
       xlab = spec$xlabel, ylab = "Probability", main = main,
       ylim = c(0, ymax), xaxt = "n", las = 1)
  at <- pretty(range(spec$x), n = 12)
  at <- at[at >= min(spec$x) & at <= max(spec$x)]
  axis(1, at = at)
  if (!is.null(selected) && selected %in% spec$x)
    points(selected, spec$pr[match(selected, spec$x)], pch = 21,
           bg = spec$color, col = COL_OBS, lwd = 3, cex = 2.9)
  if (!is.null(comparison)) {
    points(comparison$x, comparison$pr, pch = 23, bg = comparison$color,
           col = "black", cex = 2.4)
    legend("topright", legend = c(spec$label, comparison$label),
           pch = c(21, 23), pt.bg = c(spec$color, comparison$color),
           pt.cex = c(2.4, 2.4), bty = "n", cex = .83)
  }
}

cdf_plot <- function(spec, selected = NULL) {
  cp <- cumsum(spec$pr)
  # Step function: F(x) is constant between successive integer values.
  plot(spec$x, cp, type = "s", col = spec$color, lwd = 2.3,
       xlab = spec$xlabel, ylab = "Cumulative probability",
       main = "Cumulative distribution function", ylim = c(0, 1.03),
       xaxt = "n", las = 1)
  at <- pretty(range(spec$x), n = 12)
  axis(1, at = at[at >= min(spec$x) & at <= max(spec$x)])
  points(spec$x, cp, pch = 21, bg = spec$color, col = "black", cex = 2.4)
  if (!is.null(selected) && selected %in% spec$x)
    points(selected, cp[match(selected, spec$x)], pch = 21,
           bg = spec$color, col = COL_OBS, lwd = 3, cex = 2.9)
}

# Limit the displayed x range without renormalising observed proportions.
observed_probs <- function(sample, support) {
  as.numeric(table(factor(sample, levels = support))) / length(sample)
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML('
      body {font-family:"Helvetica Neue",Arial,sans-serif; color:#333; background:#fff;}
      h1,h2,h3,h4 {color:#236A6E; font-weight:400;}
      .well {background:#eef5f5; border:1px solid #cfdfe0; border-left:4px solid #4d8386;}
      /* Compact, understated navigation matching the weekly readings */
      .navbar-default {background:#fff; border:0; border-bottom:1px solid #dce6e6; border-radius:0; box-shadow:none; min-height:43px; margin-bottom:20px;}
      .navbar-default .navbar-nav>li>a {color:#4a5b5c; background:transparent; padding:11px 16px 9px; border-bottom:3px solid transparent;}
      .navbar-default .navbar-nav>li>a:hover,
      .navbar-default .navbar-nav>li>a:focus {color:#236A6E; background:#f7fafa;}
      .navbar-default .navbar-nav>.active>a,
      .navbar-default .navbar-nav>.active>a:hover,
      .navbar-default .navbar-nav>.active>a:focus {color:#236A6E; font-weight:600; background:transparent; border-bottom:3px solid #236A6E;}
      .navbar-default .navbar-toggle {border-color:#cfdfe0;}
      .navbar-default .navbar-toggle .icon-bar {background:#236A6E;}
      .navbar-default .navbar-collapse {border-color:#dce6e6;}
      .statbox {background:#eef5f5; border-left:4px solid #4d8386; padding:10px 15px; margin:12px 0;}
      .note {color:#555; font-size:0.95em;}
      .plot-panel {border:1px solid #e0e8e8; border-radius:4px; padding:5px; margin-bottom:12px;}
      .table-striped>tbody>tr:hover {background:#e8eeee;}
      .form-group {margin-bottom:12px;}
      .math-note {padding:6px 0;}
    '))
  ),
  titlePanel("STA2ASM — Week 3: Discrete distributions"),
  p(class = "note", "Explore probability distributions, simulation, overdispersion and the Poisson approximation using the examples from the Week 3 readings."),
  navbarPage(title = NULL, id = "main_tab", collapsible = TRUE,
             tabPanel("Explore distributions", value = "explore",
                      fluidRow(
                        column(3, withMathJax(wellPanel(
                          radioButtons("dist", "Distribution", choices = choices, selected = "bin"),
                          uiOutput("params_explore"),
                          uiOutput("selected_x_ui"),
                          actionButton("reset_explore", "Restore reading defaults")
                        ))),
                        column(9,
                               uiOutput("explore_math"),
                               fluidRow(
                                 column(6, div(class="plot-panel", plotOutput("pmf", height="350px"))),
                                 column(6, div(class="plot-panel", plotOutput("cdf", height="350px")))
                               ),
                               uiOutput("explore_stats"),
                               p(class="note", withMathJax(HTML("The CDF at a selected integer \\(x\\) is the sum of the PMF probabilities up to and including \\(x\\). Both plots use the same horizontal scale.")))
                        )
                      )
             ),
             tabPanel("Simulation", value = "simulation",
                      fluidRow(
                        column(3, withMathJax(wellPanel(
                          radioButtons("sim_dist", "Distribution", choices = choices, selected = "bin"),
                          uiOutput("params_sim"),
                          radioButtons("m", "Number of repetitions", choices = c("10"=10,"100"=100,"1,000"=1000,"10,000"=10000), selected=1000, inline=TRUE),
                          actionButton("resim", "Resimulate", class="btn-primary"),
                          p(class="note", "The upper plot shows one sample of the selected size. The boxplots use 200 independent samples at each sample size. Changing parameters or clicking Resimulate generates new simulations.")
                        ))),
                        column(9,
                               uiOutput("sim_math"),
                               div(class="plot-panel", plotOutput("sim_plot", height="370px")),
                               h4("Theoretical and simulated summaries — one selected sample"),
                               tableOutput("sim_stats"),
                               h4("Sampling variability across 200 independent simulations"),
                               p(class="note", "Each boxplot summarises 200 estimates from independent samples of the indicated size. Dashed lines show the theoretical values."),
                               fluidRow(
                                 column(6, div(class="plot-panel", plotOutput("sim_mean_boxes", height="310px"))),
                                 column(6, div(class="plot-panel", plotOutput("sim_var_boxes", height="310px")))
                               )
                        )
                      )
             ),
             tabPanel("Overdispersion", value = "dispersion",
                      fluidRow(
                        column(3, withMathJax(wellPanel(
                          h4("Malaria surveillance"),
                          numericInput("disp_mu", HTML("Expected annual cases \\(\\mu\\)"), value=2.6, min=.1, max=30, step=.1),
                          sliderInput("disp_theta", HTML("Dispersion \\(\\theta\\)"), min=.5, max=50, value=2, step=.5),
                          numericInput("tail_cut", HTML("Threshold for \\(P(F>x)\\)"), value=8, min=0, max=100, step=1),
                          checkboxInput("show_conversion", HTML("Show equivalent \\(r\\) and \\(p\\) (when \\(\\theta\\) is an integer)"), FALSE)
                        ))),
                        column(9,
                               uiOutput("disp_math"),
                               div(class="plot-panel", plotOutput("disp_plot", height="410px")),
                               uiOutput("disp_stats"),
                               uiOutput("disp_conversion"),
                               p(class="note", HTML("Both distributions have the same expected count. Changing \\(\\theta\\) changes the negative binomial variance, not the mean."))
                        )
                      )
             ),
             tabPanel("Poisson approximation to binomial", value = "approx",
                      fluidRow(
                        column(3, withMathJax(wellPanel(
                          h4("Rare remission example"),
                          sliderInput("approx_n", HTML("Number of patients \\(n\\)"), min=20, max=500, value=100, step=10),
                          sliderInput("approx_p", HTML("Probability of remission \\(p\\)"), min=.001, max=.30, value=.02, step=.001),
                          checkboxInput("hold_np", HTML("Hold \\(np\\) fixed when changing \\(n\\)"), FALSE),
                          numericInput("approx_x", HTML("Compare \\(P(X=x)\\) at"), value=3, min=0, max=500, step=1),
                          actionButton("reset_approx", "Restore reading defaults")
                        ))),
                        column(9,
                               uiOutput("approx_math"),
                               div(class="plot-panel", plotOutput("approx_plot", height="410px")),
                               uiOutput("approx_stats"),
                               p(class="note", withMathJax(HTML("When \\(n\\) is large and \\(p\\) is small, the \\(\\operatorname{Poisson}(\\mu)\\) distribution with \\(\\mu=np\\) may be used to approximate the \\(\\operatorname{Binomial}(n,p)\\) distribution.")))
                        )
                      )
             )
  )
)

server <- function(input, output, session) {
  # Parameter controls are local to the first two tabs; each distribution
  # retains its own controls while the user changes the distribution.
  param_ui <- function(prefix, dist) {
    id <- function(s) paste0(prefix, "_", s)
    switch(dist,
           bin = tagList(numericInput(id("n"), HTML("Trials \\(n\\)"), 20, min=1, max=300, step=1),
                         sliderInput(id("p"), HTML("Success probability \\(p\\)"), 0, 1, .30, step=.01)),
           pois = tagList(numericInput(id("lambda"), HTML("Rate per month \\(\\lambda\\)"), 2.6/12, min=.001, step=.01),
                          numericInput(id("t"), HTML("Number of months \\(t\\)"), 12, min=.1, max=120, step=1)),
           nbt = tagList(numericInput(id("r"), HTML("Required remissions \\(r\\)"), 3, min=1, max=50, step=1),
                         sliderInput(id("p"), HTML("Remission probability \\(p\\)"), 0.01, 1, .30, step=.01)),
           nba = tagList(numericInput(id("mu"), HTML("Expected count \\(\\mu\\)"), 2.6, min=.01, max=100, step=.1),
                         numericInput(id("theta"), HTML("Dispersion \\(\\theta\\)"), 2, min=.1, max=100, step=.5))
    )
  }
  output$params_explore <- renderUI(withMathJax(param_ui("e", input$dist)))
  output$params_sim <- renderUI(withMathJax(param_ui("s", input$sim_dist)))
  
  spec_for <- function(prefix, dist) {
    val <- function(name, default) {
      v <- input[[paste0(prefix, "_", name)]]
      if (is.null(v) || !is.finite(v)) default else v
    }
    make_spec(dist, n=val("n",20), p=val("p",.30),
              lambda=val("lambda",2.6/12), exposure=val("t",12),
              r=val("r",3), mu=val("mu",2.6), theta=val("theta",2))
  }
  explore_spec <- reactive(spec_for("e", input$dist))
  sim_spec <- reactive(spec_for("s", input$sim_dist))
  
  observeEvent(input$reset_explore, {
    # Rebuild distribution-specific controls using their original defaults.
    updateRadioButtons(session, "dist", selected="bin")
    output$params_explore <- renderUI(withMathJax(param_ui("e", "bin")))
  })
  
  output$selected_x_ui <- renderUI({
    s <- explore_spec()
    withMathJax(numericInput("selected_x", HTML("Highlight integer \\(x\\)"), value=round(s$mean),
                             min=min(s$x), max=max(s$x), step=1))
  })
  output$explore_math <- renderUI({
    s <- explore_spec()
    div(class="statbox", math(s$latex), p(class="note", withMathJax(HTML(s$note))))
  })
  output$pmf <- renderPlot({
    s <- explore_spec()
    par(mar=c(5.3,4.5,3.2,1.0), cex=.82)
    pmf_plot(s, selected=input$selected_x)
  })
  output$cdf <- renderPlot({
    s <- explore_spec()
    par(mar=c(5.3,4.5,3.2,1.0), cex=.82)
    cdf_plot(s, selected=input$selected_x)
  })
  output$explore_stats <- renderUI({
    s <- explore_spec(); x <- input$selected_x
    if (is.null(x) || !is.finite(x)) x <- round(s$mean)
    x <- round(x)
    px <- if (x %in% s$x) s$pr[match(x,s$x)] else 0
    cpx <- sum(s$pr[s$x <= x])
    variable <- s$variable
    math_line <- function(lhs, formula, value) {
      div(math(sprintf("%s=%s=%s", lhs, formula, fmt(value, 4))))
    }
    div(class="statbox",
        fluidRow(
          column(6,
                 math_line(sprintf("\\operatorname{E}(%s)", variable), s$mean_formula, s$mean),
                 math_line(sprintf("\\operatorname{Var}(%s)", variable), s$var_formula, s$var)
          ),
          column(6,
                 div(math(sprintf("P(%s=%d)=%s", variable, x, fmt(px, 4)))),
                 div(math(sprintf("P(%s\\leq %d)=%s", variable, x, fmt(cpx, 4)))),
                 div(math(sprintf("P(%s>%d)=%s", variable, x, fmt(1-cpx, 4))))
          )
        )
    )
  })
  
  # eventReactive avoids re-running the sample on unrelated plot/UI updates.
  sample_data <- eventReactive(list(sim_spec(), input$m, input$resim), {
    s <- sim_spec()
    s$gen(as.integer(input$m))
  }, ignoreNULL=FALSE)
  output$sim_math <- renderUI({
    s <- sim_spec(); div(class="statbox", math(s$latex), p(class="note", withMathJax(HTML(s$note))))
  })
  output$sim_plot <- renderPlot({
    s <- sim_spec(); samp <- sample_data()
    observed <- observed_probs(samp, s$x)
    ymax <- max(s$pr, observed) * 1.18
    par(mar=c(5,4.5,3.3,1), cex=.86)
    plot(s$x, s$pr, pch=21, bg=s$color, col="black", cex=2.4,
         xlab=s$xlabel, ylab="Probability / observed proportion",
         main="Observed versus theoretical probabilities", las=1,
         xaxt="n", ylim=c(0,ymax))
    at <- pretty(range(s$x), n=15)
    axis(1, at=at[at>=min(s$x) & at<=max(s$x)])
    points(s$x, observed, pch=23, bg=COL_OBS, col="black", cex=1.35)
    legend("topright", c("Theoretical probability", "Observed proportion"),
           pch=c(21,23), pt.bg=c(s$color,COL_OBS),
           pt.cex=c(2.4,1.35), bty="n", cex=1.15)
  })
  output$sim_stats <- renderTable({
    s <- sim_spec(); samp <- sample_data()
    data.frame(Statistic=c("Mean", "Variance"),
               Theoretical=c(fmt(s$mean),fmt(s$var)),
               Simulated=c(fmt(mean(samp)),fmt(var(samp))),
               check.names=FALSE)
  }, striped=TRUE, hover=TRUE, bordered=FALSE, spacing="s")
  # Repeated independent simulations for the boxplots. Changing only the
  # selected sample size does not regenerate these 200-simulation results.
  repeated_sims <- eventReactive(list(sim_spec(), input$resim), {
    s <- sim_spec()
    sizes <- c(10L, 100L, 1000L, 10000L)
    means <- matrix(NA_real_, nrow=200L, ncol=length(sizes))
    variances <- matrix(NA_real_, nrow=200L, ncol=length(sizes))
    for (j in seq_along(sizes)) {
      for (i in seq_len(200L)) {
        x <- s$gen(sizes[j])
        means[i,j] <- mean(x)
        variances[i,j] <- var(x)
      }
    }
    list(sizes=sizes, means=means, variances=variances)
  }, ignoreNULL=FALSE)
  draw_sim_boxes <- function(estimates, sizes, target, title, ylabel, fill) {
    par(mar=c(4.7,4.6,3.2,1), cex=.88)
    boxplot(as.data.frame(estimates), names=format(sizes, big.mark=","),
            col=fill, border="black", outline=TRUE,
            main=title, xlab="Repetitions per sample (m)",
            ylab=ylabel, las=1, cex.axis=.9)
    abline(h=target, col="#444444", lwd=2, lty=2)
    legend("topright", legend=paste0("Theoretical = ",fmt(target,2)),
           col="#444444", lty=2, lwd=2, bty="n", cex=.9)
  }
  output$sim_mean_boxes <- renderPlot({
    s <- sim_spec(); d <- repeated_sims()
    draw_sim_boxes(d$means, d$sizes, s$mean,
                   "Sample means", "Sample mean", s$color)
  })
  output$sim_var_boxes <- renderPlot({
    s <- sim_spec(); d <- repeated_sims()
    draw_sim_boxes(d$variances, d$sizes, s$var,
                   "Sample variances", "Sample variance", s$color)
  })
  
  output$disp_math <- renderUI({
    mu <- input$disp_mu; theta <- input$disp_theta
    div(class="statbox",
        math(sprintf("\\mu=%.2f,\\quad \\theta=%.2f",mu,theta)),
        math(",\\quad \\operatorname{Var}(F)=\\mu+\\frac{\\mu^2}{\\theta}"))
  })
  output$disp_plot <- renderPlot({
    mu <- input$disp_mu; theta <- input$disp_theta
    xmax <- max(20L, as.integer(qnbinom(.999, size=theta, mu=mu)))
    xx <- 0:xmax; a <- dpois(xx,mu); b <- dnbinom(xx,size=theta,mu=mu)
    par(mar=c(5,4.5,3.3,1),cex=.85)
    plot(xx,a,pch=21,bg=COL_POIS,col="black",cex=2.4,las=1,
         xlab="Number of malaria cases in one year",ylab="Probability",
         main="Poisson and negative binomial distributions for the malaria example",xaxt="n",
         ylim=c(0,max(a,b)*1.17))
    at <- pretty(range(xx),n=15)
    axis(1,at=at[at>=0 & at<=xmax])
    points(xx,b,pch=23,bg=COL_NB,col="black",cex=2.4)
    legend("topright", legend=as.expression(list(
      bquote(paste("Poisson (", mu, " = ", .(round(mu,2)), ")")),
      bquote(paste("NB (", mu, " = ", .(round(mu,2)), ", ", theta, " = ", .(round(theta,2)), ")")))),
      pch=c(21,23),pt.bg=c(COL_POIS,COL_NB),pt.cex=1.6,bty="n",cex=1.15)
  })
  output$disp_stats <- renderUI({
    mu <- input$disp_mu; theta <- input$disp_theta; x <- round(input$tail_cut)
    rows <- list(
      c("Poisson", fmt(mu), fmt(mu), fmt(ppois(x, mu, lower.tail=FALSE), 5)),
      c("Negative binomial", fmt(mu), fmt(mu+mu^2/theta),
        fmt(pnbinom(x, size=theta, mu=mu, lower.tail=FALSE), 5))
    )
    withMathJax(tags$table(class="table table-striped table-hover",
                           tags$thead(tags$tr(
                             tags$th("Model"), tags$th("Mean"), tags$th("Variance"),
                             tags$th(HTML(sprintf("\\(P(F > %d)\\)", x)))
                           )),
                           tags$tbody(lapply(rows, function(row) tags$tr(lapply(row, tags$td))))
    ))
  })
  output$disp_conversion <- renderUI({
    if (!isTRUE(input$show_conversion)) return(NULL)
    theta <- input$disp_theta; mu <- input$disp_mu
    p <- theta/(theta+mu)
    div(class="statbox",
        math(sprintf("p=\\frac{\\theta}{\\theta+\\mu}=%.4f",p)),
        if (abs(theta-round(theta))<1e-8)
          withMathJax(tags$p(HTML(sprintf("Here \\(r=%d\\) and \\(p=%.4f\\): the same failures-count distribution can be constructed by waiting for \\(r\\) successes. Malaria cases are not literally failures.",as.integer(theta),p))))
        else withMathJax(tags$p(HTML("\\(\\theta\\) is not an integer, so there is no literal interpretation as waiting for a fractional number of successes."))))
  })
  
  observeEvent(input$reset_approx, {
    updateSliderInput(session,"approx_n",value=100)
    updateSliderInput(session,"approx_p",value=.02)
    updateCheckboxInput(session,"hold_np",value=FALSE)
    updateNumericInput(session,"approx_x",value=3)
  })
  # Capture np when the hold option is enabled; then change p when n changes.
  held_mean <- reactiveVal(NULL)
  observeEvent(input$hold_np, {
    if (isTRUE(input$hold_np)) held_mean(input$approx_n*input$approx_p)
    else held_mean(NULL)
  })
  observeEvent(input$approx_n, {
    if (isTRUE(input$hold_np) && !is.null(held_mean()))
      updateSliderInput(session,"approx_p",value=min(.30,max(.001,held_mean()/input$approx_n)))
  })
  output$approx_math <- renderUI({
    n <- input$approx_n; p <- input$approx_p
    div(class="statbox",math(sprintf("X\\sim\\operatorname{Binomial}(%d,%.3f),\\quad Y\\sim\\operatorname{Poisson}(\\mu=%.3f)",n,p,n*p)))
  })
  output$approx_plot <- renderPlot({
    n <- input$approx_n; p <- input$approx_p; mu <- n*p
    xmax <- max(12L,as.integer(qpois(.9995,mu)),as.integer(qbinom(.9995,n,p)))
    xx <- 0:xmax; a <- dbinom(xx,n,p); b <- dpois(xx,mu)
    par(mar=c(5,4.5,3.3,1),cex=.85)
    plot(xx,a,pch=21,bg=COL_BIN,col="black",cex=2.4,las=1,
         xlab="Number of patients achieving remission",ylab="Probability",
         main="Binomial and its Poisson approximation",xaxt="n",
         ylim=c(0,max(a,b)*1.18))
    at <- pretty(range(xx),n=15)
    axis(1,at=at[at>=0 & at<=xmax])
    points(xx,b,pch=23,bg=COL_POIS,col="black",cex=1.45)
    legend("topright",legend=c("Binomial (exact)","Poisson (approximation)"),
           pch=c(21,23),pt.bg=c(COL_BIN,COL_POIS),pt.cex=c(2.2,1.45),bty="n",cex=1.15)
  })
  output$approx_stats <- renderUI({
    n <- input$approx_n; p <- input$approx_p; mu <- n*p
    x <- max(0L,round(input$approx_x))
    rows <- list(
      c("Binomial (exact)", fmt(mu), fmt(mu*(1-p)), fmt(dbinom(x,n,p),4)),
      c("Poisson (approximation)", fmt(mu), fmt(mu), fmt(dpois(x,mu),4))
    )
    withMathJax(tags$table(class="table table-striped table-hover",
                           tags$thead(tags$tr(
                             tags$th("Model"), tags$th("Mean"), tags$th("Variance"),
                             tags$th(HTML("\\(P(X=x)\\)"))
                           )),
                           tags$tbody(lapply(rows, function(row) tags$tr(lapply(row, tags$td))))
    ))
  })
}

shinyApp(ui, server)
