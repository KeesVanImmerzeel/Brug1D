# Simulate time-dependent head and flow in a semi-infinite aquifer adjacent to 
# open water where the boundary condition changes from t = 0 according to a 
# specified course.
#
# The equation used describes one-dimensional flow for a situation with 
# spatially constant kD and S, without recharge (precipitation, evaporation, leakage).
# From: 'Van Edelman naar Bruggeman', Stromingen 12 (2006)

library(shiny)
library(shinyBS)
library(pracma)
library(ggplot2)
library(dplyr)
library(plotly)

#' Calculate the argument in the complementary error function in the Bruggeman formula
#'
#' @param x Distance from the river [L]
#' @param S Storage coefficient [-]
#' @param kD Hydraulic conductivity [L2/T]
#' @param t Time [T]
#' @example calc_u(x=100, S=0.15, kD=250, t=1)
#' @return The calculated value
calc_u <- function(x, S, kD, t) {
      sqrt(x^2 * S / (4 * kD * t))
}

#' Calculate the repeated integral of the complementary error function
#'
#' @param z Argument in the complementary error function [-]
#' @param n In the repeated integral [-]
#' @example ierfc(z=0.5, n=0)
#' @return The calculated value [-]
ierfc <- function(z, n) {
      if (n == -1) {
            return(2 / sqrt(pi) * exp(-z^2))
      } else if (n == 0) {
            return(erfc(z))
      } else {
            result <- -z / n * ierfc(z, n - 1) + 1 / (2 * n) * ierfc(z, n - 2)
            result <- max(result, 0)
            return(result)
      }
}

#' Calculate the change in head s [L] and the horizontal flux change q [L2/T] at distance x and time t due to a stress change a at t=0 at x=0
#'
#' @param x Distance from the river [L]
#' @param t Time [T]
#' @param S Storage coefficient [-]
#' @param kD Hydraulic conductivity [L2/T]
#' @param n Type of stress change (0, 1, 2 or 3): 
#' n=0: River stage changes suddenly with a fixed value a [L]; 
#' n=1: Constant infiltration a [L2/T] starts at x = 0; 
#' n=2: River stage rises at a constant rate a [L/T]; 
#' n=3: Infiltration at x = 0 increases at a constant rate a [L2/T]. 
#' @param a Constant in the definition of the stress [L] or [L2/T]
#' @example calc_stress_response(x=100, t=1, S=0.15, kD=250, n=0, a=5)
#' @return A dataframe with the columns x, t, s, q, n, and a
calc_stress_response <- function(x, t, S, kD, n, a) {
      expand.grid(x = x, t = t) |>
            dplyr::rowwise() |>
            dplyr::mutate(
                  u = calc_u(x, S, kD, t),
                  s = dplyr::case_when(
                        n == 0 ~ a * ierfc(u, 0) / ierfc(0, 0),
                        n == 1 ~ 2 * a * sqrt(t) * ierfc(u, 1) / (ierfc(0, 0) * sqrt(kD * S)),
                        n == 2 ~ a * t * ierfc(u, 2) / (ierfc(0, 2)),
                        n == 3 ~ 2 * a * t * sqrt(t) * ierfc(u, 3) / (ierfc(0, 2) * sqrt(kD * S))
                  ),
                  q = dplyr::case_when(
                        n == 0 ~ a * sqrt(kD * S) * ierfc(u, -1) / (ierfc(0, 0) * 2 * sqrt(t)),
                        n == 1 ~ a * ierfc(u, 0) / (ierfc(0, 0)),
                        n == 2 ~ a * sqrt(t) * sqrt(kD * S) * ierfc(u, 1) / (2 * ierfc(0, 2)),
                        n == 3 ~ a * t * ierfc(u, 2) / (ierfc(0, 2))
                  ),
                  n = n,
                  a = a
            ) |>
            dplyr::select(x, t, s, q, n, a)
}

#' Generate the title text for graph based on the type of stress change.
#' 
#' @param item: item description 
#' @param n: Type of stress change (0, 1, 2 or 3)
#' @param a: Constant in the definition of the stress [L] or [L2/T]
make_title_text <- function(item="Change in head s [L]:",n, a) {
      title_text <- switch(
            as.character(n),
            "0" = paste(item,
                  "River stage changes suddenly with a fixed value a=",
                  a,
                  "[L]"
            ),
            "1" = paste(item,
                  "Constant infiltration a=",
                  a,
                  "[L2/T] starts at x = 0"
            ),
            "2" = paste(item,
                  "River stage rises at a constant rate a=",
                  a,
                  "[L/T]"
            ),
            "3" = paste(item,
                  "Infiltration at x = 0 increases at a constant rate a=",
                  a,
                  "[L2/T]"
            ),
            paste("Unknown stress change")
      )
      return(title_text)
}

#' Plot the change in head s [L] as a function of distance x [L] for different time points t
#'
#' @param df Dataframe with the results of calc_stress_response
#' @return A ggplot object
plot_stress_response_x_s <- function(df) {
      # Generate the title text based on the type of stress change
      title_text <- make_title_text(n = df$n[1], a = df$a[1])
      p <- ggplot(df, aes(
            x = x,
            y = s,
            color = factor(t)
      )) +
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Distance from the river [L]",
                  y = "Change in head s [L]",
                  color = "Time [T]"
            ) +
            theme_minimal() +
            theme(
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, margin = margin(t = 30)),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 16, margin = margin(b = 20)),
                  panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        size = 0.5
                  ),
                  legend.box.background = element_rect(color = "black", size = 0.5)
            )
      # Convert ggplot to plotly
      plotly::ggplotly(p)
}

#' Plot the horizontal flux change q [L2/T] as a function of distance x [L] for different time points t
#'
#' @param df Dataframe with the results of calc_stress_response
#' @return A ggplot object
plot_stress_response_x_q <- function(df) {
      # Generate the title text based on the type of stress change
      title_text <- make_title_text(item = "Horizontal flux change q [L2/T]:",
                                    n = df$n[1],
                                    a = df$a[1])
      
      # Plot q versus x at different times t if there are multiple x values
      p <- ggplot(df, aes(
            x = x,
            y = q,
            color = factor(t)
      )) +
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Distance from the river [L]",
                  y = "Horizontal flux change q [L2/T]",
                  color = "Time [T]"
            ) +
            theme_minimal() +
            theme(
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, margin = margin(t = 30)),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 16, margin = margin(b = 20)),
                  panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        size = 0.5
                  ),
                  legend.box.background = element_rect(color = "black", size = 0.5)
            )
      # Convert ggplot to plotly
      plotly::ggplotly(p)
}

#' Plot the change in head s [L] as a function of time t [T] at different distances x [L]
#'
#' @param df Dataframe with the results of calc_stress_response
#' @return A ggplot object
plot_stress_response_t_s <- function(df) {
      # Generate the title text based on the type of stress change
      title_text <- make_title_text(n = df$n[1], a = df$a[1])
      p <- ggplot(df, aes(
            x = t,
            y = s,
            color = factor(x)
      )) +
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Time [T]",
                  y = "Change in head s [L]",
                  color = "Distance from the river [L]"
            ) +
            theme_minimal() +
            theme(
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, margin = margin(t = 30)),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 16, margin = margin(b = 20)),
                  panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        size = 0.5
                  ),
                  legend.box.background = element_rect(color = "black", size = 0.5)
            )
      # Convert ggplot to plotly
      plotly::ggplotly(p)
}

#' Plot the horizontal flux change q [L2/T] as a function of time t [T] at different distances x [L]
#'
#' @param df Dataframe with the results of calc_stress_response
#' @return A ggplot object
plot_stress_response_t_q <- function(df) {
      # Generate the title text based on the type of stress change
      title_text <- make_title_text(item = "Horizontal flux change q [L2/T]:",
                                    n = df$n[1],
                                    a = df$a[1])
      
      # Plot q versus x at different times t if there are multiple x values
      p <- ggplot(df, aes(
            x = t,
            y = q,
            color = factor(x)
      )) +
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Time [T]",
                  y = "Horizontal flux change q [L2/T]",
                  color = "Distance from the river [L]"
            ) +
            theme_minimal() +
            theme(
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, margin = margin(t = 30)),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 16, margin = margin(b = 20)),
                  panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        size = 0.5
                  ),
                  legend.box.background = element_rect(color = "black", size = 0.5)
            )
      # Convert ggplot to plotly
      plotly::ggplotly(p)
}

# Define UI for the Shiny app
ui <- fluidPage(
      titlePanel("Stress Response Plot"),
      tabsetPanel(
            tabPanel("System",
                     sidebarLayout(
                           sidebarPanel(
                                 radioButtons("n", "Type of stress change:",
                                              choices = list("0: River stage changes suddenly with a fixed value 'a' [L]" = 0,
                                                             "1: Constant infiltration 'a' starts at x = 0 [L2/T]" = 1,
                                                             "2: River stage rises at a constant rate 'a' [L/T]" = 2,
                                                             "3: Infiltration at x = 0 increases at a constant rate 'a' [L2/T]" = 3)),
                                 numericInput("a", "Value of a:", 5),
                                 numericInput("kD", "Hydraulic conductivity kD [L2/T]:", 250),
                                 numericInput("S", "Storage coefficient S [-]:", 0.15),
                                 bsTooltip("kD", "Must be greater than 1", "top", options = list(container = "body")),
                                 bsTooltip("S", "Must be greater than 0.0001 and less than 1", "top", options = list(container = "body")),
                                 downloadButton("downloadData", "Download input data"),
                                 fileInput("uploadData", "Upload input data", accept = c(".csv")),
                                 tags$a(href = "https://github.com/KeesVanImmerzeel/Brug1D/tree/master", "Documentation")
                           ),
                           mainPanel(
                                 
                           )
                     )
            ),
            tabPanel("System Plots", tabsetPanel(
                  tabPanel("x, result",
                           sidebarLayout(
                                 sidebarPanel(
                                       actionButton("plot_button", "Refresh", class = "btn-warning"),
                                       br(), br(),
                                       numericInput("num_points_x", "Number of points in x-array:", value = 100, min = 1),
                                       numericInput("min_x", "Minimum value of x:", value = 1, min = 0),
                                       numericInput("max_x", "Maximum value of x:", value = 1000, min = 0),
                                       numericInput("num_points_t", "Number of points in t-array:", value = 5, min = 1, max = 10),
                                       numericInput("min_t", "Minimum value of t:", value = 1, min = 0),
                                       numericInput("max_t", "Maximum value of t:", value = 1000, min = 0),
                                       bsTooltip("num_points_x", "value >= 1", "top", options = list(container = "body")),
                                       bsTooltip("num_points_t", "1 <= value <= 10", "top", options = list(container = "body")),
                                       bsTooltip("min_x", "value >= 0", "top", options = list(container = "body")),
                                       bsTooltip("max_x", "value > min_x", "top", options = list(container = "body")),
                                       bsTooltip("min_t", "value >= 0", "top", options = list(container = "body")),
                                       bsTooltip("max_t", "value >= min_t", "top", options = list(container = "body")),
                                       br(),
                                       downloadButton("downloadResults", "Download Results")
                                 ),
                                 mainPanel(
                                       plotly::plotlyOutput("stressPlotS"), br(), br(),
                                       plotly::plotlyOutput("stressPlotQ"), br(), br(),
                                       tableOutput("resultsTable")
                                 )
                           )
                  ),
                  tabPanel("t, result",
                           sidebarLayout(
                                 sidebarPanel(
                                       actionButton("plot_button_t", "Refresh", class = "btn-warning"),
                                       br(), br(),
                                       numericInput("num_points_x_t", "Number of points in x-array:", value = 6, min = 1, max=10),
                                       numericInput("min_x_t", "Minimum value of x:", value = 0, min = 0),
                                       numericInput("max_x_t", "Maximum value of x:", value = 1000, min = 0),
                                       numericInput("num_points_t_t", "Number of points in t-array:", value = 100, min = 1),
                                       numericInput("min_t_t", "Minimum value of t:", value = 1, min = 0),
                                       numericInput("max_t_t", "Maximum value of t:", value = 1000, min = 0),
                                       bsTooltip("num_points_x_t","1 <= value <= 10" , "top", options = list(container = "body")),
                                       bsTooltip("num_points_t_t", "value >= 1", "top", options = list(container = "body")),
                                       bsTooltip("min_x_t", "value >= 0", "top", options = list(container = "body")),
                                       bsTooltip("max_x_t", "value > min_x_t", "top", options = list(container = "body")),
                                       bsTooltip("min_t_t", "value >= 0", "top", options = list(container = "body")),
                                       bsTooltip("max_t_t", "value >= min_t_t", "top", options = list(container = "body")),
                                       br(),
                                       downloadButton("downloadResults_t", "Download Time Plot Results")
                                 ),
                                 mainPanel(
                                       plotly::plotlyOutput("stressPlotS_t"), br(), br(),
                                       plotly::plotlyOutput("stressPlotQ_t"), br(), br(),
                                       tableOutput("resultsTable_t")
                                 )
                           )
                  )
            ))
      )
)

# Define server logic for the Shiny app
server <- function(input, output, session) {
      result_x <- reactiveVal()
      result_t <- reactiveVal()
      
      observeEvent(input$plot_button, {
            x_values <- seq(input$min_x, input$max_x, length.out = input$num_points_x)
            t_values <- seq(input$min_t, input$max_t, length.out = input$num_points_t)
            result_x(calc_stress_response(x = x_values, t = t_values, S = input$S, kD = input$kD, n = input$n, a = input$a))
            output$stressPlotS <- renderPlotly({
                  plot_stress_response_x_s(result_x())
            })
            output$stressPlotQ <- renderPlotly({
                  plot_stress_response_x_q(result_x())
            })
            output$resultsTable <- renderTable({
                  result_x()
            })
      })
      
      observeEvent(input$plot_button_t, {
            x_values <- seq(input$min_x_t, input$max_x_t, length.out = input$num_points_x_t)
            t_values <- seq(input$min_t_t, input$max_t_t, length.out = input$num_points_t_t)
            result_t(calc_stress_response(x = x_values, t = t_values, S = input$S, kD = input$kD, n = input$n, a = input$a))
            output$stressPlotS_t <- renderPlotly({
                  plot_stress_response_t_s(result_t())
            })
            output$stressPlotQ_t <- renderPlotly({
                  plot_stress_response_t_q(result_t())
            })
            output$resultsTable_t <- renderTable({
                  result_t()
            })
      })
      
      output$downloadData <- downloadHandler(
            filename = function() {
                  paste("Brug1D_input_data_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  input_data <- data.frame(
                        n = input$n,
                        a = input$a,
                        kD = input$kD,
                        S = input$S,
                        num_points_x = input$num_points_x,
                        min_x = input$min_x,
                        max_x = input$max_x,
                        num_points_t = input$num_points_t,
                        min_t = input$min_t,
                        max_t = input$max_t,
                        num_points_x_t = input$num_points_x_t,
                        min_x_t = input$min_x_t,
                        max_x_t = input$max_x_t,
                        num_points_t_t = input$num_points_t_t,
                        min_t_t = input$min_t_t,
                        max_t_t = input$max_t_t
                  )
                  write.csv2(input_data, file, quote = FALSE, row.names = FALSE)
            }
      )
      
      output$downloadResults <- downloadHandler(
            filename = function() {
                  paste("Brug1D_x_vs_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(result_x(), file, quote = FALSE, row.names = FALSE)
            }
      )
      
      output$downloadResults_t <- downloadHandler(
            filename = function() {
                  paste("Brug1D_t_vs_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(result_t(), file, quote = FALSE, row.names = FALSE)
            }
      )
      
      observeEvent(input$uploadData, {
            req(input$uploadData)
            input_data <- read.csv(input$uploadData$datapath)
            updateRadioButtons(session, "n", selected = input_data$n)
            updateNumericInput(session, "a", value = input_data$a)
            updateNumericInput(session, "kD", value = input_data$kD)
            updateNumericInput(session, "S", value = input_data$S)
            updateNumericInput(session, "num_points_x", value = input_data$num_points_x)
            updateNumericInput(session, "min_x", value = input_data$min_x)
            updateNumericInput(session, "max_x", value = input_data$max_x)
            updateNumericInput(session, "num_points_t", value = input_data$num_points_t)
            updateNumericInput(session, "min_t", value = input_data$min_t)
            updateNumericInput(session, "max_t", value = input_data$max_t)
            updateNumericInput(session, "num_points_x_t", value = input_data$num_points_x_t)
            updateNumericInput(session, "min_x_t", value = input_data$min_x_t)
            updateNumericInput(session, "max_x_t", value = input_data$max_x_t)
            updateNumericInput(session, "num_points_t_t", value = input_data$num_points_t_t)
            updateNumericInput(session, "min_t_t", value = input_data$min_t_t)
            updateNumericInput(session, "max_t_t", value = input_data$max_t_t)
      })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)