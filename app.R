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

#' Plot the change in head s [L] for different time points t
#'
#' @param df Dataframe with the results of calc_stress_response
#' @example plot_stress_response_s(df = result_calc_stress_response_n0)
#' @return A ggplot object
plot_stress_response_s <- function(df) {
      title_text <- switch(as.character(df$n[1]),
                           "0" = paste("Change in head s: River stage changes suddenly with a fixed value a=", df$a[1], "[L]"),
                           "1" = paste("Change in head s: Constant infiltration a=", df$a[1], "[L2/T] starts at x = 0"),
                           "2" = paste("Change in head s: River stage rises at a constant rate a=", df$a[1], "[L/T]"),
                           "3" = paste("Change in head s: Infiltration at x = 0 increases at a constant rate a=", df$a[1], "[L2/T]"),
                           paste("Unknown stress change"))
      
      p <- ggplot(df, aes(x = x, y = s, color = factor(t))) +
            geom_line(linewidth = 1) +
            labs(title = title_text,
                 x = "Distance from the river [L]",
                 y = "Change in head s [L]",
                 color = "Time [T]") +
            theme_minimal() +
            theme(axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, margin = margin(t = 30)),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 16, margin = margin(b = 20)),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  legend.box.background = element_rect(color = "black", size = 0.5))
      plotly::ggplotly(p)
}

#' Plot the horizontal flux change q [L2/T] for different time points t
#'
#' @param df Dataframe with the results of calc_stress_response
#' @example plot_stress_response_q(df = result_calc_stress_response_n0)
#' @return A ggplot object
plot_stress_response_q <- function(df) {
      title_text <- switch(as.character(df$n[1]),
                           "0" = paste("Horizontal flux change q [L2/T]: River stage changes suddenly with a fixed value a=", df$a[1], "[L]"),
                           "1" = paste("Horizontal flux change q [L2/T]: Constant infiltration a=", df$a[1], "[L2/T] starts at x = 0"),
                           "2" = paste("Horizontal flux change q [L2/T]: River stage rises at a constant rate a=", df$a[1], "[L/T]"),
                           "3" = paste("Horizontal flux change q [L2/T]: Infiltration at x = 0 increases at a constant rate a=", df$a[1], "[L2/T]"),
                           paste("Unknown stress change"))
      
      p <- ggplot(df, aes(x = x, y = q, color = factor(t))) +
            geom_line(linewidth = 1) +
            labs(title = title_text,
                 x = "Distance from the river [L]",
                 y = "Horizontal flux change q [L2/T]",
                 color = "Time [T]") +
            theme_minimal() +
            theme(axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, margin = margin(t = 30)),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 16, margin = margin(b = 20)),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  legend.box.background = element_rect(color = "black", size = 0.5))
      plotly::ggplotly(p)
}

# Define UI for the Shiny app
ui <- fluidPage(
      titlePanel("Stress Response Plot"),
      tabsetPanel(
            tabPanel("Control",
                     sidebarLayout(
                           sidebarPanel(
                                 actionButton("plot_button", "Refresh results and plots", class = "btn-warning"),
                                 radioButtons("n", "Type of stress change:",
                                              choices = list("0: River stage changes suddenly with a fixed value 'a' [L]" = 0,
                                                             "1: Constant infiltration 'a' starts at x = 0 [L2/T]" = 1,
                                                             "2: River stage rises at a constant rate 'a' [L/T]" = 2,
                                                             "3: Infiltration at x = 0 increases at a constant rate 'a' [L2/T]" = 3)),
                                 numericInput("a", "Value of a:", 5),
                                 numericInput("kD", "Hydraulic conductivity kD [L2/T]:", 250),
                                 numericInput("S", "Storage coefficient S [-]:", 0.15),
                                 numericInput("num_points_x", "Number of points in x-array:", 100),
                                 numericInput("min_x", "Minimum value of x:", 0),
                                 numericInput("max_x", "Maximum value of x:", 1000),
                                 numericInput("num_points_t", "Number of points in t-array:", 10),
                                 numericInput("min_t", "Minimum value of t:", 1),
                                 numericInput("max_t", "Maximum value of t:", 100),
                                 bsTooltip("kD", "Must be greater than 1", "top", options = list(container = "body")),
                                 bsTooltip("S", "Must be greater than 0.0001 and less than 1", "top", options = list(container = "body")),
                                 bsTooltip("num_points_x", "Must be greater than 0", "top", options = list(container = "body")),
                                 bsTooltip("num_points_t", "Must be greater than 0", "top", options = list(container = "body")),
                                 bsTooltip("min_x", "Must be greater than or equal to 0", "top", options = list(container = "body")),
                                 bsTooltip("max_x", "Must be greater than min_x", "top", options = list(container = "body")),
                                 bsTooltip("min_t", "Must be greater than or equal to 1", "top", options = list(container = "body")),
                                 downloadButton("downloadData", "Download input data"),
                                 fileInput("uploadData", "Upload input data", accept = c(".csv")),
                                 tags$a(href = "https://github.com/KeesVanImmerzeel/Brug1D/tree/master", "Documentation")
                           ),
                           mainPanel(
                                 plotly::plotlyOutput("stressPlotS"),
                                 plotly::plotlyOutput("stressPlotQ")
                           )
                     )
            ),
            tabPanel("Results",
                     tableOutput("resultsTable"),
                     downloadButton("downloadResults", "Download Results")
            )
      )
)

# Define server logic for the Shiny app
server <- function(input, output, session) {
      result <- reactiveVal()
      
      observeEvent(input$plot_button, {
            x_values <- seq(input$min_x, input$max_x, length.out = input$num_points_x)
            t_values <- seq(input$min_t, input$max_t, length.out = input$num_points_t)
            
            result(calc_stress_response(x = x_values, t = t_values, S = input$S, kD = input$kD, n = input$n, a = input$a))
            
            output$stressPlotS <- renderPlotly({
                  plot_stress_response_s(result())
            })
            
            output$stressPlotQ <- renderPlotly({
                  plot_stress_response_q(result())
            })
            
            output$resultsTable <- renderTable({
                  result()
            })
      })
      
      output$downloadData <- downloadHandler(
            filename = function() {
                  paste("input_data-", Sys.Date(), ".csv", sep = "")
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
                        max_t = input$max_t
                  )
                  write.csv2(input_data, file, quote=FALSE, row.names = FALSE)
            }
      )
      
      output$downloadResults <- downloadHandler(
            filename = function() {
                  paste("results-", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(result(), file, quote=FALSE, row.names = FALSE)
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
      })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)