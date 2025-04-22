# Simulate time-dependent head and flow in a semi-infinite aquifer adjacent to 
# open water where the boundary condition changes from t = 0 according to a 
# specified course.
#
# The equation used describes one-dimensional flow for a situation with sim_min_t
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
#' @param x Distance [L]
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
#' @param x Distance [L]
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


############# Function related to simulation

#' Calculate the stress period index i, the stress a and the time t to the next output time for multiple values of t
#' 
#' @param t Vector of time values [T]
#' @param t_a Data frame with columns t [T] and a [L]
#' @example stress_period(t = c(1, 5, 10), t_a = data.frame(t = c(0, 2, 4), a = c(1, -1)))
#' @return A dataframe with columns t, sp, a, t0, and dt
stress_period <- function(t, t_a) {
      result_list <- lapply(t, function(single_t) {
            i <- which(t_a$t < single_t)
            data.frame(t = single_t, sp = i, a = t_a$a[i], t0 = t_a$t[i], dt = single_t - t_a$t[i])
      })
      result_df <- do.call(rbind, result_list)
      return(result_df)
}

#' Calculate the stress period data frame extended with a single value of x and multiple values of t
#' 
#' @param x Single value of distance [L]
#' @param t Vector of time values [T]
#' @param t_a Data frame with columns t [T] and a [L]
#' @example distance_stress_period_single(x = 5, t = c(1, 5, 10), t_a = data.frame(t = c(0, 2, 4), a = c(1, -1)))
#' @return A dataframe with columns x, t, sp, a, t0, and dt
distance_stress_period_single <- function(x, t, t_a) {
      data.frame(x = x, stress_period(t, t_a))
}

#' Calculate the stress period data frame extended with multiple values of x and multiple values of t
#' 
#' @param x Vector of distance values [L]
#' @param t Vector of time values [T]
#' @param t_a Data frame with columns t [T] and a [L]
#' @example distance_stress_period(x = c(5, 10, 15), t = c(1, 5, 10), t_a = data.frame(t = c(0, 2, 4), a = c(1, -1)))
#' @return A dataframe with columns x, t, sp, a, t0, and dt
distance_stress_period <- function(x, t, t_a) {
      result_list <- lapply(x, function(single_x) distance_stress_period_single(single_x, t, t_a))
      result_df <- do.call(rbind, result_list)
      return(result_df)
}

#' Simulate the stress response for given distances, output times, and stress-time series t_a
#' 
#' @param x Vector of distance values [L]
#' @param t Vector of output times [T]
#' @param S Storage coefficient [-]
#' @param kD Hydraulic conductivity [L2/T]
#' @param n Type of stress change (0, 1, 2 or 3): 
#' n=0: River stage changes suddenly with a fixed value a [L]; 
#' n=1: Constant infiltration a [L2/T] starts at x = 0; 
#' n=2: River stage rises at a constant rate a [L/T]; 
#' n=3: Infiltration at x = 0 increases at a constant rate a [L2/T]. 
#' @param t_a Data frame with columns t [T] and a [L]
#' @example simulate_stress_response(x = c(5, 10, 100), t = c(1, 10, 100), S = 0.15, kD = 250, n = 0, t_a = data.frame(t = c(0, 50), a = c(1, -1)))
#' @return A dataframe with summarized stress response
simulate_stress_response <- function(x, t, S, kD, n, t_a) {
      t_a <- t_a[order(t_a$t), ] # Order on time t
      t_a <- t_a[c(TRUE, diff(t_a$a) != 0), ] # Remove duplicates
      df <- distance_stress_period(x, t, t_a)
      df_stress_response <- data.frame(b = mapply(calc_stress_response, x=df$x, t=df$dt, a=df$a, MoreArgs = list(S = S, kD = kD, n = n))) |>
            t() |> as.data.frame() |> dplyr::rename(dt = t)
      row.names(df_stress_response) <- NULL
      df_stress_response$t <- df$t
      df_stress_response <- data.frame(x = unlist(df_stress_response$x), t = unlist(df_stress_response$t), s = unlist(df_stress_response$s), q = unlist(df_stress_response$q)) |>
            dplyr::group_by(x, t) |> dplyr::summarise(s = sum(s), q = sum(q))
      return(df_stress_response)
}

#############

#' Generate the text of the unit of a (constant in the definition of the stress [L] or [L2/T])
#' 
#' @param n: Type of stress change (0, 1, 2 or 3)
make_unit_text_of_a <- function(n) {
      unit_text <- switch(
            as.character(n),
            "0" = "[L]",
            "1" = "[L2/T]",
            "2" = "[L/T]",
            "3" = "[L2/T]",
            "[?]"
      )
}

#' Generate the title text for graph based on the type of stress change.
#' 
#' @param n: Type of stress change (0, 1, 2 or 3)
#' @param a: Constant in the definition of the stress [L] or [L2/T]
make_title_text <- function(n, a) {
      title_text <- switch(
            as.character(n),
            "0" = paste(
                  "River stage changes suddenly with a fixed value a=",
                  a, make_unit_text_of_a(n)
            ),
            "1" = paste(
                  "Constant infiltration a=",
                  a, make_unit_text_of_a(n),
                  "starts at x = 0"
            ),
            "2" = paste(
                  "River stage rises at a constant rate a=",
                  a, make_unit_text_of_a(n)
            ),
            "3" = paste(
                  "Infiltration at x = 0 increases at a constant rate a=",
                  a, make_unit_text_of_a(n)
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
      if (all(c("n", "a") %in% names(df))) {
            title_text <- make_title_text(n = df$n[1], a = df$a[1])
      }
      else {
            title_text <- ""
      }
      p <- ggplot(df, aes(
            x = x,
            y = s,
            color = factor(t)
      )) +
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Distance [L]",
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
      if (all(c("n", "a") %in% names(df))) {
            title_text <- make_title_text(n = df$n[1], a = df$a[1])
      }
      else {
            title_text <- ""
      }
      
      # Plot q versus x at different times t if there are multiple x values
      p <- ggplot(df, aes(
            x = x,
            y = q,
            color = factor(t)
      )) +
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Distance [L]",
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
plot_stress_response_t_s <- function(df, t_a=NULL) {
      # Generate the title text based on the type of stress change
      if (all(c("n", "a") %in% names(df))) {
            title_text <- make_title_text(n = df$n[1], a = df$a[1])
      }
      else {
            title_text <- ""
      }
      
      p <- ggplot(df, aes(
            x = t,
            y = s,
            color = factor(x)
      )) +
            
          #  {if (!is.null(t_a)) 
          #        ggplot2::geom_step(data = t_a, mapping = ggplot2::aes(x = t, y = a), color = "black", linetype = "dashed")} +

            
          #  {if (!is.null(t_a)) 
          #      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ ., name = "a"))} +
            
                        
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Time [T]",
                  y = "Change in head s [L]",
                  color = "Distance [L]"
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
      if (all(c("n", "a") %in% names(df))) {
            title_text <- make_title_text(n = df$n[1], a = df$a[1])
      }
      else {
            title_text <- ""
      }
      
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
                  color = "Distance [L]"
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
            tabPanel("System definition",
                     sidebarLayout(
                           sidebarPanel(
                                 radioButtons("n", "Type of stress:",
                                              choices = list("0: River stage changes suddenly with a fixed value 'a' [L]" = 0,
                                                             "1: Constant infiltration 'a' starts at x = 0 [L2/T]" = 1,
                                                             "2: River stage rises at a constant rate 'a' [L/T]" = 1,
                                                             "3: Infiltration at x = 0 increases at a constant rate 'a' [L2/T]" = 1)),
                                 numericInput("a", "Value of 'a' ([L] or [L2/T]):", 1),
                                 numericInput("kD", "Hydraulic conductivity kD [L2/T]:", value=250, min=0.001),
                                 numericInput("S", "Storage coefficient S [-]:", value=0.15, min=0.00001, max=1),
                                 bsTooltip("kD", "value > 0 [L2/T]", "top", options = list(container = "body")),
                                 bsTooltip("S", "0.00001 < value <= 1", "top", options = list(container = "body")),
                                 downloadButton("downloadData", "Download input data"), br(), br(),
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
                                       actionButton("plot_button_x", "Refresh", class = "btn-warning"),
                                       br(), br(),
                                       numericInput("num_points_x", "Number of points in x-array:", value = 100, min = 1),
                                       numericInput("min_x", "Minimum value of x:", value = 1, min = 0),
                                       numericInput("max_x", "Maximum value of x:", value = 1000, min = 0),
                                       numericInput("num_points_t", "Number of points in t-array:", value = 5, min = 1, max = 10),
                                       numericInput("min_t", "Minimum value of t:", value = 1, min = 0.001),
                                       numericInput("max_t", "Maximum value of t:", value = 1000, min = 0),
                                       bsTooltip("num_points_x", "value >= 1", "top", options = list(container = "body")),
                                       bsTooltip("num_points_t", "1 <= value <= 10", "top", options = list(container = "body")),
                                       bsTooltip("min_x", "value >= 0", "top", options = list(container = "body")),
                                       bsTooltip("max_x", "value > min_x", "top", options = list(container = "body")),
                                       bsTooltip("min_t", "value > 0", "top", options = list(container = "body")),
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
                                       numericInput("min_t_t", "Minimum value of t:", value = 1, min = 0.001),
                                       numericInput("max_t_t", "Maximum value of t:", value = 1000, min = 0),
                                       bsTooltip("num_points_x_t","1 <= value <= 10" , "top", options = list(container = "body")),
                                       bsTooltip("num_points_t_t", "value >= 1", "top", options = list(container = "body")),
                                       bsTooltip("min_x_t", "value > 0", "top", options = list(container = "body")),
                                       bsTooltip("max_x_t", "value > min_x_t", "top", options = list(container = "body")),
                                       bsTooltip("min_t_t", "value >= 0", "top", options = list(container = "body")),
                                       bsTooltip("max_t_t", "value >= min_t_t", "top", options = list(container = "body")),
                                       br(),
                                       downloadButton("downloadResults_t", "Download Results")
                                 ),
                                 mainPanel(
                                       plotly::plotlyOutput("stressPlotS_t"), br(), br(),
                                       plotly::plotlyOutput("stressPlotQ_t"), br(), br(),
                                       tableOutput("resultsTable_t")
                                 )
                           )
                  )
            )),
            tabPanel("Simulation definition",
                     sidebarLayout(
                           sidebarPanel(
                                 fileInput("file1", "Upload Spreadsheet",
                                           accept = c(".csv", ".xlsx"))
                           ),
                           mainPanel(
                                 plotly::plotlyOutput("stress_plot"),
                                 tableOutput("stress_sequence")
                           )
                     )
            ), 
            tabPanel("Simulation plots", tabsetPanel(
                  tabPanel("x, result",
                           sidebarLayout(
                                 sidebarPanel(
                                       actionButton("sim_plot_button_x", "Refresh", class = "btn-warning"),
                                       br(), br(),
                                       numericInput("sim_num_points_x", "Number of points in x-array:", value = 100, min = 1),
                                       numericInput("sim_min_x", "Minimum value of x:", value = 1, min = 0),
                                       numericInput("sim_max_x", "Maximum value of x:", value = 1000, min = 0),
                                       numericInput("sim_num_points_t", "Number of points in t-array:", value = 5, min = 1, max = 10),
                                       numericInput("sim_min_t", "Minimum value of t:", value = 1, min = 0.001),
                                       numericInput("sim_max_t", "Maximum value of t:", value = 1000, min = 0),
                                       bsTooltip("sim_num_points_x", "value >= 1", "top", options = list(container = "body")),
                                       bsTooltip("sim_num_points_t", "1 <= value <= 10", "top", options = list(container = "body")),
                                       bsTooltip("sim_min_x", "value >= 0", "top", options = list(container = "body")),
                                       bsTooltip("sim_max_x", "value > min_x", "top", options = list(container = "body")),
                                       bsTooltip("sim_min_t", "value > 0", "top", options = list(container = "body")),
                                       bsTooltip("sim_max_t", "value >= min_t", "top", options = list(container = "body")),
                                       br(),
                                       downloadButton("sim_downloadResults", "Download Results")
                                 ),
                                 mainPanel(
                                       plotly::plotlyOutput("sim_stressPlotS"), br(), br(),
                                       plotly::plotlyOutput("sim_stressPlotQ"), br(), br(),
                                       tableOutput("sim_resultsTable")
                                 )
                           )
                  ),
                  tabPanel("t, result",
                           sidebarLayout(
                                 sidebarPanel(
                                       actionButton("sim_plot_button_t", "Refresh", class = "btn-warning"),
                                       br(), br(),
                                       numericInput("sim_num_points_x_t", "Number of points in x-array:", value = 6, min = 1, max=10),
                                       numericInput("sim_min_x_t", "Minimum value of x:", value = 0, min = 0),
                                       numericInput("sim_max_x_t", "Maximum value of x:", value = 1000, min = 0),
                                       numericInput("sim_num_points_t_t", "Number of points in t-array:", value = 100, min = 1),
                                       numericInput("sim_min_t_t", "Minimum value of t:", value = 1, min = 0.001),
                                       numericInput("sim_max_t_t", "Maximum value of t:", value = 250, min = 0),
                                       bsTooltip("sim_num_points_x_t","1 <= value <= 10" , "top", options = list(container = "body")),
                                       bsTooltip("sim_num_points_t_t", "value >= 1", "top", options = list(container = "body")),
                                       bsTooltip("sim_min_x_t", "value > 0", "top", options = list(container = "body")),
                                       bsTooltip("sim_max_x_t", "value > min_x_t", "top", options = list(container = "body")),
                                       bsTooltip("sim_min_t_t", "value >= 0", "top", options = list(container = "body")),
                                       bsTooltip("sim_max_t_t", "value >= min_t_t", "top", options = list(container = "body")),
                                       br(),
                                       downloadButton("sim_downloadResults_t", "Download Results")
                                 ),
                                 mainPanel(
                                       plotly::plotlyOutput("sim_stressPlotS_t"), br(), br(),
                                       plotly::plotlyOutput("sim_stressPlotQ_t"), br(), br(),
                                       tableOutput("sim_resultsTable_t")
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
      sim_result_x <- reactiveVal()
      sim_result_t <- reactiveVal()
      
      observeEvent(input$plot_button_x, {
            x_values <- seq(input$min_x, input$max_x, length.out = input$num_points_x)
            #t_values <- seq(input$min_t, input$max_t, length.out = input$num_points_t)
            t_values <- round(pracma::logseq(input$min_t, input$max_t, input$num_points_t),1)
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
                        max_t_t = input$max_t_t,
                        
                        sim_num_points_x = input$sim_num_points_x,
                        sim_min_x = input$sim_min_x,
                        sim_max_x = input$sim_max_x,
                        sim_num_points_t = input$sim_num_points_t,
                        sim_min_t = input$sim_min_t,
                        sim_max_t = input$sim_max_t,
                        
                        sim_num_points_x_t = input$sim_num_points_x_t,
                        sim_min_x_t = input$sim_min_x_t,
                        sim_max_x_t = input$sim_max_x_t,
                        sim_num_points_t_t = input$sim_num_points_t_t,
                        sim_min_t_t = input$sim_min_t_t,
                        sim_max_t_t = input$sim_max_t_t
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
      
      output$sim_downloadResults <- downloadHandler(
            filename = function() {
                  paste("Brug1D_x_vs_sim_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(sim_result_x(), file, quote = FALSE, row.names = FALSE)
            }
      )
      
      output$sim_downloadResults_t <- downloadHandler(
            filename = function() {
                  paste("Brug1D_t_vs_sim_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(sim_result_t(), file, quote = FALSE, row.names = FALSE)
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
            
            updateNumericInput(session, "sim_num_points_x", value = input_data$sim_num_points_x)
            updateNumericInput(session, "sim_min_x", value = input_data$sim_min_x)
            updateNumericInput(session, "sim_max_x", value = input_data$sim_max_x)
            updateNumericInput(session, "sim_num_points_t", value = input_data$sim_num_points_t)
            updateNumericInput(session, "sim_min_t", value = input_data$sim_min_t)
            updateNumericInput(session, "sim_max_t", value = input_data$sim_max_t)
            
            updateNumericInput(session, "sim_num_points_x_t", value = input_data$sim_num_points_x_t)
            updateNumericInput(session, "sim_min_x_t", value = input_data$sim_min_x_t)
            updateNumericInput(session, "sim_max_x_t", value = input_data$sim_max_x_t)
            updateNumericInput(session, "sim_num_points_t_t", value = input_data$sim_num_points_t_t)
            updateNumericInput(session, "sim_min_t_t", value = input_data$sim_min_t_t)
            updateNumericInput(session, "sim_max_t_t", value = input_data$sim_max_t_t)
            
      })
      
      ######### Simulate in- and output
      
      stress_sequence <- reactive({
            req(input$file1)
            inFile <- input$file1
            
            if (grepl("\\.csv$", inFile$name)) {
                  df <- read.csv(inFile$datapath)
            } else if (grepl("\\.xlsx$", inFile$name)) {
                  df <- readxl::read_excel(inFile$datapath)
            }
            
            df <- df[, c("Time", "a")]
            df <- df[order(df$Time), ]
            #df <- df[c(TRUE, diff(df$a) != 0), ] # Remove duplicates
            return(df)
      })
      
      output$stress_plot <- renderPlotly({
            req(stress_sequence())
            df <- stress_sequence()
            
            # Ensure values remain the same until the next timestep
            df <- df %>%
                  mutate(Time_end = lead(Time, default = max(Time)),
                         a_end = a) %>%
                  tidyr::pivot_longer(cols = c("Time", "Time_end"), names_to = "Time_type", values_to = "Time_value") %>%
                  arrange(Time_value)
            
            ylab <- paste( "a", unit_str <- make_unit_text_of_a( input$n ))
            p <- ggplot(df, aes(x = Time_value, y = a)) +
                  geom_step(linewidth = 1) +
                  labs(title = "Stress Sequence Plot", x = "Time [T]", y = ylab) +
                  theme(
                        axis.text = element_text(size = 14),
                        axis.title = element_text(size = 16, margin = margin(t = 30)),
                        legend.text = element_text(size = 14),
                        legend.title = element_text(size = 16, margin = margin(b = 20)),
                        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                        legend.box.background = element_rect(color = "black", size = 0.5)
                  )
            # Convert ggplot to plotly
            plotly::ggplotly(p)
            
      })      
      
      output$stress_sequence <- renderTable({
            req(stress_sequence())
            stress_sequence()
      })
      
      observeEvent(input$sim_plot_button_x, {
            x_values <- seq(input$sim_min_x, input$sim_max_x, length.out = input$sim_num_points_x)
            t_values <- round(pracma::logseq(input$sim_min_t, input$sim_max_t, input$sim_num_points_t),1)
            
            t_a <- stress_sequence()
            names(t_a) <- c("t", "a")
            sim_result_x(simulate_stress_response(x = x_values, t = t_values, S = input$S, kD = input$kD, n = input$n, t_a = t_a ))
            output$sim_stressPlotS <- renderPlotly({
                  plot_stress_response_x_s(sim_result_x())
            })
            output$sim_stressPlotQ <- renderPlotly({
                  plot_stress_response_x_q(sim_result_x())
            })
            output$sim_resultsTable <- renderTable({
                  sim_result_x()
            })
      })
      
      observeEvent(input$sim_plot_button_t, {
            x_values <- seq(input$sim_min_x_t, input$sim_max_x_t, length.out = input$sim_num_points_x_t)
            t_values <- seq(input$sim_min_t_t, input$sim_max_t_t, length.out = input$sim_num_points_t_t)
            
            t_a <- stress_sequence()
            names(t_a) <- c("t", "a")
            sim_result_t(simulate_stress_response(x = x_values, t = t_values, S = input$S, kD = input$kD, n = input$n, t_a = t_a))
            output$sim_stressPlotS_t <- renderPlotly({
                  plot_stress_response_t_s(sim_result_t(), t_a)
            })
            output$sim_stressPlotQ_t <- renderPlotly({
                  plot_stress_response_t_q(sim_result_t())
            })
            output$sim_resultsTable_t <- renderTable({
                  sim_result_t()
            })
      })
      
}

# Run the Shiny app
shinyApp(ui = ui, server = server)