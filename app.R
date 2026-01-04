# =============================================================================
# Groningen Building Damage Explorer
#
# Interactive Shiny app for exploring spatial damage predictions
# =============================================================================

library(shiny)
library(tidyverse)
library(leaflet)
library(DT)
library(plotly)
library(sf)
library(here)

# =============================================================================
# CUSTOM UI HELPERS (must be defined before UI)
# =============================================================================

# Custom valueBox function (since shinydashboard not loaded)
valueBox <- function(value, subtitle, icon = NULL, color = "blue") {
  color_map <- c(
    blue = "#3498db", orange = "#f39c12",
    green = "#2ecc71", red = "#e74c3c"
  )
  bg_color <- color_map[color]

  div(
    class = "well",
    style = paste0(
      "background-color:", bg_color, "; color: white; ",
      "text-align: center; padding: 20px; border-radius: 5px;"
    ),
    h2(value, style = "margin: 0;"),
    p(subtitle, style = "margin: 5px 0 0 0;")
  )
}

renderValueBox <- function(expr, env = parent.frame(), quoted = FALSE) {
  func <- shiny::exprToFunction(expr, env, quoted)
  shiny::renderUI(func())
}

valueBoxOutput <- function(outputId, width = 12) {
  column(width, uiOutput(outputId))
}

# =============================================================================
# LOAD DATA
# =============================================================================

# Load model results
model_data <- list(
  additive = readRDS(here("outputs", "models", "spatial_damage_uncertainty_full.rds")),
  daf = readRDS(here("outputs", "models", "spatial_damage_daf.rds"))
)

# Model descriptions for UI
model_descriptions <- c(
  additive = "Additive (naive Markov)",
  daf = "Korswagen DAF (N-effect)"
)

# Default results (will be reactive)
results <- model_data$additive
earthquakes <- readRDS(here("outputs", "models", "groningen_earthquakes.rds"))

# Function to compute PC4 aggregates from building data
compute_pc4_damage <- function(data) {
  data |>
    group_by(pc4) |>
    summarise(
      n = n(),
      mean_delta_psi_virgin = mean(psi_virgin_mean, na.rm = TRUE),
      mean_delta_psi_predamage = mean(psi_predamage_mean, na.rm = TRUE),
      mean_psi_virgin = mean(psi_virgin_mean, na.rm = TRUE),
      mean_psi_predamage = mean(psi_predamage_mean, na.rm = TRUE),
      p10_delta_psi_virgin = mean(psi_virgin_p10, na.rm = TRUE),
      p50_delta_psi_virgin = mean(psi_virgin_p50, na.rm = TRUE),
      p90_delta_psi_virgin = mean(psi_virgin_p90, na.rm = TRUE),
      p10_delta_psi_predamage = mean(psi_predamage_p10, na.rm = TRUE),
      p50_delta_psi_predamage = mean(psi_predamage_p50, na.rm = TRUE),
      p90_delta_psi_predamage = mean(psi_predamage_p90, na.rm = TRUE),
      p10_psi_virgin = mean(psi_virgin_p10, na.rm = TRUE),
      p50_psi_virgin = mean(psi_virgin_p50, na.rm = TRUE),
      p90_psi_virgin = mean(psi_virgin_p90, na.rm = TRUE),
      p10_psi_predamage = mean(psi_predamage_p10, na.rm = TRUE),
      p50_psi_predamage = mean(psi_predamage_p50, na.rm = TRUE),
      p90_psi_predamage = mean(psi_predamage_p90, na.rm = TRUE),
      p_visible_virgin = mean(p_visible_virgin, na.rm = TRUE),
      p_visible_predamage = mean(p_visible_predamage, na.rm = TRUE),
      p_psi_1_virgin = mean(p_visible_virgin, na.rm = TRUE),
      p_psi_1_predamage = mean(p_visible_predamage, na.rm = TRUE),
      p_psi_2_5_virgin = mean(p_moderate_virgin, na.rm = TRUE),
      p_psi_2_5_predamage = mean(p_moderate_predamage, na.rm = TRUE),
      mean_max_pgv = mean(max_pgv_mean, na.rm = TRUE),
      .groups = "drop"
    )
}

# Function to compute zone/age summaries
compute_summaries <- function(data) {
  by_zone <- data |>
    group_by(zone) |>
    summarise(
      n = n(),
      psi_mean = mean(psi_virgin_mean, na.rm = TRUE),
      p_visible_virgin = mean(p_visible_virgin, na.rm = TRUE),
      p_visible_predamage = mean(p_visible_predamage, na.rm = TRUE),
      .groups = "drop"
    )

  by_age <- data |>
    group_by(age_category) |>
    summarise(
      n = n(),
      psi_mean = mean(psi_virgin_mean, na.rm = TRUE),
      p_visible_virgin = mean(p_visible_virgin, na.rm = TRUE),
      p_visible_predamage = mean(p_visible_predamage, na.rm = TRUE),
      .groups = "drop"
    )

  list(by_zone = by_zone, by_age = by_age)
}

# Load PC4 polygon boundaries
pc4_polygons <- st_read(here("datafiles", "pc4_groningen.gpkg"), quiet = TRUE)

# Pre-compute PC4 aggregates for all models (once at startup for performance)
pc4_damage_all <- lapply(model_data, compute_pc4_damage)
summaries_all <- lapply(model_data, compute_summaries)

# =============================================================================
# LOAD DAMAGE CLAIMS DATA
# =============================================================================

# Load weighted PC4 damage claims 2020-2024
claims_raw <- read_csv(here("datafiles", "img_physical_damage_pc4_2020_2024_weighted.csv"),
                       show_col_types = FALSE)

# Aggregate across years per PC4
claims_pc4 <- claims_raw |>
  mutate(pc4 = as.character(pc4)) |>
  group_by(pc4) |>
  summarise(
    total_claims = sum(claims_n, na.rm = TRUE),
    total_positive_decisions = sum(pos_besluiten_n, na.rm = TRUE),
    total_compensation_eur = sum(vergoeding_eur, na.rm = TRUE),
    n_addresses = max(adressen, na.rm = TRUE),
    .groups = "drop"
  )

# Function to join claims to PC4 predictions for a given model
join_claims_to_pc4 <- function(pc4_damage) {
  pc4_damage |>
    left_join(claims_pc4, by = "pc4") |>
    mutate(
      eur_per_building = ifelse(is.na(total_compensation_eur), 0, total_compensation_eur / n),
      claims_per_building = ifelse(is.na(total_claims), 0, total_claims / n),
      total_claims = replace_na(total_claims, 0),
      total_compensation_eur = replace_na(total_compensation_eur, 0)
    )
}

# Function to create map data for a given model
create_pc4_map <- function(pc4_damage) {
  pc4_polygons |>
    left_join(pc4_damage, by = "pc4") |>
    filter(!is.na(mean_delta_psi_virgin)) |>
    left_join(claims_pc4, by = "pc4") |>
    mutate(
      n_buildings = n,
      eur_per_building = ifelse(is.na(total_compensation_eur), 0, total_compensation_eur / n),
      total_claims = replace_na(total_claims, 0),
      total_compensation_eur = replace_na(total_compensation_eur, 0)
    )
}

# Pre-compute comparison data for all models
pc4_comparison_all <- lapply(pc4_damage_all, join_claims_to_pc4)
pc4_map_all <- lapply(pc4_damage_all, create_pc4_map)

# =============================================================================
# UI
# =============================================================================

ui <- navbarPage(
 "Groningen Damage Explorer",
  theme = bslib::bs_theme(bootswatch = "flatly"),

  # --- TAB 1: OVERVIEW ---
  tabPanel("Overview",
    fluidRow(
      column(8,
        h3("Groningen Building Damage Predictions"),
        p("Cumulative earthquake damage from 124 induced events (M ≥ 2.0, 1991-2024)")
      ),
      column(4,
        selectInput("model_select", "Accumulation Model",
          choices = c(
            "Additive (naive Markov)" = "additive",
            "Korswagen DAF (N-effect)" = "daf"
          ),
          selected = "daf",
          width = "100%"
        ),
        helpText(
          style = "font-size: 11px; color: #666;",
          HTML("<b>Additive:</b> State-dependent, each event adds independently<br>
                <b>DAF:</b> History-aware with N-effect for similar events")
        )
      )
    ),
    fluidRow(
      valueBoxOutput("n_buildings", width = 3),
      valueBoxOutput("n_earthquakes", width = 3),
      valueBoxOutput("pct_visible_virgin", width = 3),
      valueBoxOutput("pct_visible_predamage", width = 3)
    ),
    fluidRow(
      column(6,
        h4("Damage by Distance Zone"),
        plotlyOutput("zone_plot", height = "350px")
      ),
      column(6,
        h4("Damage by Building Age"),
        plotlyOutput("age_plot", height = "350px")
      )
    ),
    fluidRow(
      column(12,
        h4("Damage Distribution"),
        plotlyOutput("histogram", height = "300px")
      )
    )
  ),

  # --- TAB 2: MAP COMPARISON ---
  tabPanel("Map",
    fluidRow(
      column(12,
        h4("Predicted vs Observed Damage (Normalized 0-1 scale)"),
        fluidRow(
          column(4,
            selectInput("map_metric", "Predicted Metric",
              choices = c("Mean ΔΨ" = "mean_delta_psi",
                          "P(Ψ ≥ 1.0)" = "p_psi_1",
                          "P(Ψ ≥ 2.5)" = "p_psi_2_5"),
              selected = "mean_delta_psi")
          ),
          column(4,
            selectInput("map_scenario", "Scenario",
              choices = c("Virgin Walls" = "virgin", "Pre-Damage" = "predamage"))
          )
        )
      )
    ),
    fluidRow(
      column(6,
        h5("Predicted Damage", style = "text-align: center;"),
        leafletOutput("map_predicted", height = "600px")
      ),
      column(6,
        h5("Observed: € per Building", style = "text-align: center;"),
        leafletOutput("map_observed", height = "600px")
      )
    )
  ),

  # --- TAB 3: POSTCODE EXPLORER ---
  tabPanel("Postcode Explorer",
    sidebarLayout(
      sidebarPanel(width = 3,
        selectInput("pc4_select", "Select Postcode (PC4)",
          choices = sort(unique(pc4_damage_all$additive$pc4)),
          selected = "9651"),  # Loppersum
        hr(),
        h5("Postcode Summary"),
        verbatimTextOutput("pc4_summary")
      ),
      mainPanel(width = 9,
        fluidRow(
          column(6, plotlyOutput("pc4_histogram", height = "300px")),
          column(6, plotlyOutput("pc4_age_breakdown", height = "300px"))
        ),
        hr(),
        h4("Buildings in this Postcode"),
        DTOutput("pc4_table")
      )
    )
  ),

  # --- TAB 4: EARTHQUAKE CATALOG ---
  tabPanel("Earthquakes",
    fluidRow(
      column(8,
        h4("Earthquake Catalog (M ≥ 2.0)"),
        DTOutput("eq_table")
      ),
      column(4,
        h4("Magnitude Distribution"),
        plotlyOutput("eq_histogram", height = "250px"),
        h4("Timeline"),
        plotlyOutput("eq_timeline", height = "250px")
      )
    )
  ),

  # --- TAB 5: VALIDATION ---
  tabPanel("Validation",
    h3("Model Validation: Predicted Damage vs Observed Compensation"),
    p("Comparison of model predictions with IMG damage compensation data at PC4 level."),
    p("Uncertainty bands: 80% posterior credible intervals from 500 Monte Carlo samples (10 GMM × 50 fragility posterior draws)."),
    fluidRow(
      column(6,
        h4("Mean ΔΨ vs € per Building"),
        plotlyOutput("scatter_delta_eur", height = "400px")
      ),
      column(6,
        h4("P(Ψ ≥ 1.0) vs € per Building"),
        plotlyOutput("scatter_psi1_eur", height = "400px")
      )
    ),
    fluidRow(
      column(6,
        h4("P(Ψ ≥ 2.5) vs € per Building"),
        plotlyOutput("scatter_psi25_eur", height = "400px")
      ),
      column(6,
        h4("Mean Max PGV vs € per Building"),
        plotlyOutput("scatter_pgv_eur", height = "400px")
      )
    ),
    hr(),
    fluidRow(
      column(12,
        h4("Correlation Summary"),
        verbatimTextOutput("correlation_summary")
      )
    )
  ),

  # --- TAB 6: MODEL COMPARISON ---
  tabPanel("Model Comparison",
    h3("Comparison: Additive vs DAF"),
    p("Side-by-side comparison of damage accumulation models across all 542,957 buildings."),
    fluidRow(
      column(12,
        h4("Model Summary Statistics"),
        tableOutput("model_comparison_table")
      )
    ),
    fluidRow(
      column(6,
        h4("Mean Ψ by Zone"),
        plotlyOutput("compare_zone_plot", height = "350px")
      ),
      column(6,
        h4("P(Visible) by Zone"),
        plotlyOutput("compare_visible_plot", height = "350px")
      )
    ),
    fluidRow(
      column(12,
        h4("Damage Distribution Comparison"),
        plotlyOutput("compare_histogram", height = "350px")
      )
    )
  ),

  # --- TAB 7: DATA DOWNLOAD ---
  tabPanel("Download",
    h3("Download Results"),
    p("Download the full building-level results for further analysis."),
    fluidRow(
      column(4,
        h4("Building Results"),
        downloadButton("download_buildings", "Download CSV (542k rows)")
      ),
      column(4,
        h4("Postcode Summary"),
        downloadButton("download_postcodes", "Download CSV")
      ),
      column(4,
        h4("Earthquake Catalog"),
        downloadButton("download_earthquakes", "Download CSV")
      )
    ),
    hr(),
    h4("Data Dictionary"),
    tags$ul(
      tags$li(tags$b("psi_virgin"), " - Cumulative damage (Ψ) assuming virgin walls (Ψ₀=0)"),
      tags$li(tags$b("psi_predamage"), " - Cumulative damage with age-based initial damage"),
      tags$li(tags$b("initial_psi"), " - Assumed initial damage based on construction year"),
      tags$li(tags$b("max_pgv"), " - Maximum PGV experienced (mm/s)"),
      tags$li(tags$b("visible_*"), " - Binary: Ψ ≥ 1 (visible damage threshold)"),
      tags$li(tags$b("moderate_*"), " - Binary: Ψ ≥ 2 (moderate damage threshold)"),
      tags$li(tags$b("severe_*"), " - Binary: Ψ ≥ 3 (severe damage threshold)")
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # --- REACTIVE DATA BASED ON MODEL SELECTION ---
  selected_results <- reactive({
    model_data[[input$model_select]]
  })

  selected_summaries <- reactive({
    summaries_all[[input$model_select]]
  })

  selected_pc4_damage <- reactive({
    pc4_damage_all[[input$model_select]]
  })

  selected_pc4_map <- reactive({
    pc4_map_all[[input$model_select]]
  })

  selected_pc4_comparison <- reactive({
    pc4_comparison_all[[input$model_select]]
  })

  # --- VALUE BOXES ---
  output$n_buildings <- renderValueBox({
    valueBox(
      format(nrow(selected_results()), big.mark = ","),
      "Buildings",
      icon = icon("building"),
      color = "blue"
    )
  })

  output$n_earthquakes <- renderValueBox({
    valueBox(
      nrow(earthquakes),
      "Earthquakes (M≥2)",
      icon = icon("house-crack"),
      color = "orange"
    )
  })

  output$pct_visible_virgin <- renderValueBox({
    data <- selected_results()
    valueBox(
      paste0(round(100 * mean(data$p_visible_virgin), 1), "%"),
      "Visible Damage (Virgin)",
      icon = icon("chart-line"),
      color = "green"
    )
  })

  output$pct_visible_predamage <- renderValueBox({
    data <- selected_results()
    valueBox(
      paste0(round(100 * mean(data$p_visible_predamage), 1), "%"),
      "Visible Damage (Pre-Damage)",
      icon = icon("chart-line"),
      color = "red"
    )
  })

  # --- OVERVIEW PLOTS ---
  output$zone_plot <- renderPlotly({
    summaries <- selected_summaries()
    plot_data <- summaries$by_zone |>
      select(zone, Virgin = p_visible_virgin, `Pre-Damage` = p_visible_predamage) |>
      pivot_longer(-zone, names_to = "Scenario", values_to = "p_visible") |>
      mutate(zone = factor(zone, levels = c("Near (<10km)", "Mid (10-25km)", "Far (>25km)")))

    p <- ggplot(plot_data, aes(x = zone, y = p_visible, fill = Scenario)) +
      geom_col(position = "dodge") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = c("Virgin" = "#3498db", "Pre-Damage" = "#e74c3c")) +
      labs(x = NULL, y = "P(Visible Damage)") +
      theme_minimal()

    ggplotly(p, tooltip = c("x", "y", "fill"))
  })

  output$age_plot <- renderPlotly({
    summaries <- selected_summaries()
    plot_data <- summaries$by_age |>
      select(age_category, Virgin = p_visible_virgin, `Pre-Damage` = p_visible_predamage) |>
      pivot_longer(-age_category, names_to = "Scenario", values_to = "p_visible") |>
      mutate(age_category = factor(age_category,
        levels = c("Pre-1920", "1920-1949", "1950-1969", "1970-1990", "Post-1991")))

    p <- ggplot(plot_data, aes(x = age_category, y = p_visible, fill = Scenario)) +
      geom_col(position = "dodge") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = c("Virgin" = "#3498db", "Pre-Damage" = "#e74c3c")) +
      labs(x = NULL, y = "P(Visible Damage)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggplotly(p, tooltip = c("x", "y", "fill"))
  })

  output$histogram <- renderPlotly({
    data <- selected_results()
    plot_data <- data |>
      select(psi_virgin_mean, psi_predamage_mean) |>
      sample_n(min(50000, n())) |>
      pivot_longer(everything(), names_to = "Scenario", values_to = "psi") |>
      mutate(Scenario = ifelse(Scenario == "psi_virgin_mean", "Virgin", "Pre-Damage"))

    p <- ggplot(plot_data, aes(x = psi, fill = Scenario)) +
      geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("Virgin" = "#3498db", "Pre-Damage" = "#e74c3c")) +
      labs(x = "Cumulative Ψ (posterior mean)", y = "Count") +
      xlim(0, 4) +
      theme_minimal()

    ggplotly(p)
  })

  # --- SIDE-BY-SIDE MAPS ---

  # Initialize predicted map
  output$map_predicted <- renderLeaflet({
    leaflet() |>
      addProviderTiles(providers$CartoDB.Positron) |>
      setView(lng = 6.72, lat = 53.35, zoom = 10)
  })

  # Initialize observed map
  output$map_observed <- renderLeaflet({
    leaflet() |>
      addProviderTiles(providers$CartoDB.Positron) |>
      setView(lng = 6.72, lat = 53.35, zoom = 10)
  })

  # Update predicted map
  observe({
    metric_col <- paste0(input$map_metric, "_", input$map_scenario)
    map_data <- selected_pc4_map()
    map_data$val <- map_data[[metric_col]]

    # Use colorBin with quantile-based breaks (handles zeros/duplicates gracefully)
    vals <- map_data$val[!is.na(map_data$val) & map_data$val > 0]
    if (length(vals) > 10) {
      breaks <- unique(c(0, quantile(vals, probs = seq(0.2, 1, by = 0.2), na.rm = TRUE)))
    } else {
      breaks <- c(0, max(map_data$val, na.rm = TRUE))
    }
    pal <- colorBin("YlOrRd", domain = map_data$val, bins = breaks, na.color = "#cccccc")

    # Create popup with model info
    model_name <- names(model_descriptions)[match(input$model_select, c("additive", "daf"))]
    map_data$popup <- paste0(
      "<b>PC4: ", map_data$pc4, "</b><br>",
      "Model: ", input$model_select, "<br>",
      "Buildings: ", map_data$n_buildings, "<br>",
      input$map_metric, ": ", round(map_data$val, 3)
    )

    leafletProxy("map_predicted") |>
      clearShapes() |>
      clearControls() |>
      addPolygons(
        data = map_data,
        fillColor = ~pal(val),
        fillOpacity = 0.7,
        color = "#333333",
        weight = 1,
        popup = ~popup,
        highlightOptions = highlightOptions(weight = 3, color = "#000", fillOpacity = 0.9, bringToFront = TRUE)
      ) |>
      addLegend(position = "bottomright", pal = pal, values = map_data$val,
                title = input$map_metric, opacity = 0.8)
  })

  # Update observed map (always € per building)
  observe({
    # Trigger when model or map inputs change
    input$model_select
    input$map_metric
    input$map_scenario

    map_data <- selected_pc4_map()

    # Use colorBin with quantile-based breaks (handles zeros/duplicates gracefully)
    vals <- map_data$eur_per_building[!is.na(map_data$eur_per_building) & map_data$eur_per_building > 0]
    if (length(vals) > 10) {
      breaks <- unique(c(0, quantile(vals, probs = seq(0.2, 1, by = 0.2), na.rm = TRUE)))
    } else {
      breaks <- c(0, max(map_data$eur_per_building, na.rm = TRUE))
    }
    pal <- colorBin("YlOrRd", domain = map_data$eur_per_building, bins = breaks, na.color = "#cccccc")

    # Create popup
    map_data$popup <- paste0(
      "<b>PC4: ", map_data$pc4, "</b><br>",
      "Buildings: ", map_data$n_buildings, "<br>",
      "€/Building: €", format(round(map_data$eur_per_building), big.mark = ",")
    )

    leafletProxy("map_observed") |>
      clearShapes() |>
      clearControls() |>
      addPolygons(
        data = map_data,
        fillColor = ~pal(eur_per_building),
        fillOpacity = 0.7,
        color = "#333333",
        weight = 1,
        popup = ~popup,
        highlightOptions = highlightOptions(weight = 3, color = "#000", fillOpacity = 0.9, bringToFront = TRUE)
      ) |>
      addLegend(position = "bottomright", pal = pal, values = map_data$eur_per_building,
                title = "€/Building", opacity = 0.8)
  })

  # --- POSTCODE EXPLORER ---
  pc4_buildings <- reactive({
    selected_results() |> filter(pc4 == input$pc4_select)
  })

  output$pc4_summary <- renderPrint({
    data <- pc4_buildings()
    cat("Buildings:", nrow(data), "\n")
    cat("Mean Ψ (virgin):", round(mean(data$psi_virgin_mean), 2), "\n")
    cat("  [p10-p90]:", round(mean(data$psi_virgin_p10), 2), "-", round(mean(data$psi_virgin_p90), 2), "\n")
    cat("Mean Ψ (pre-damage):", round(mean(data$psi_predamage_mean), 2), "\n")
    cat("P(Visible) virgin:", round(100 * mean(data$p_visible_virgin), 1), "%\n")
    cat("P(Visible) pre-damage:", round(100 * mean(data$p_visible_predamage), 1), "%\n")
    cat("Mean Max PGV:", round(mean(data$max_pgv_mean), 1), "mm/s\n")
  })

  output$pc4_histogram <- renderPlotly({
    data <- pc4_buildings()

    p <- ggplot(data, aes(x = psi_virgin_mean)) +
      geom_histogram(binwidth = 0.2, fill = "#3498db", alpha = 0.7) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      labs(title = "Damage Distribution (Virgin)", x = "Ψ (posterior mean)", y = "Count") +
      theme_minimal()

    ggplotly(p)
  })

  output$pc4_age_breakdown <- renderPlotly({
    data <- pc4_buildings() |>
      group_by(age_category) |>
      summarise(
        n = n(),
        p_visible = mean(p_visible_virgin),
        .groups = "drop"
      ) |>
      mutate(age_category = factor(age_category,
        levels = c("Pre-1920", "1920-1949", "1950-1969", "1970-1990", "Post-1991")))

    p <- ggplot(data, aes(x = age_category, y = n, fill = p_visible)) +
      geom_col() +
      scale_fill_gradient(low = "#3498db", high = "#e74c3c", limits = c(0, 1)) +
      labs(title = "Buildings by Age", x = NULL, y = "Count", fill = "P(Visible)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggplotly(p)
  })

  output$pc4_table <- renderDT({
    pc4_buildings() |>
      select(identificatie, bouwjaar, psi_virgin_mean, psi_virgin_p10, psi_virgin_p90,
             psi_predamage_mean, max_pgv_mean, p_visible_virgin, age_category) |>
      mutate(across(c(psi_virgin_mean, psi_virgin_p10, psi_virgin_p90,
                      psi_predamage_mean, max_pgv_mean, p_visible_virgin), ~round(., 2))) |>
      datatable(
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE,
        colnames = c("ID", "Year", "Ψ Mean", "Ψ p10", "Ψ p90", "Ψ Pre-Damage", "Max PGV", "P(Visible)", "Age")
      )
  })

  # --- EARTHQUAKE CATALOG ---
  output$eq_table <- renderDT({
    earthquakes |>
      arrange(desc(mag)) |>
      select(time_utc, mag, lat, lon, depth_km) |>
      mutate(
        time_utc = format(time_utc, "%Y-%m-%d %H:%M"),
        mag = round(mag, 2),
        lat = round(lat, 3),
        lon = round(lon, 3),
        depth_km = round(depth_km, 1)
      ) |>
      datatable(
        options = list(pageLength = 20),
        rownames = FALSE,
        colnames = c("Date/Time", "Magnitude", "Latitude", "Longitude", "Depth (km)")
      )
  })

  output$eq_histogram <- renderPlotly({
    p <- ggplot(earthquakes, aes(x = mag)) +
      geom_histogram(binwidth = 0.1, fill = "#f39c12", color = "white") +
      labs(x = "Magnitude", y = "Count") +
      theme_minimal()

    ggplotly(p)
  })

  output$eq_timeline <- renderPlotly({
    eq_year <- earthquakes |>
      mutate(year = year(time_utc)) |>
      count(year)

    p <- ggplot(eq_year, aes(x = year, y = n)) +
      geom_col(fill = "#f39c12") +
      labs(x = "Year", y = "Events") +
      theme_minimal()

    ggplotly(p)
  })

  # --- VALIDATION PLOTS ---
  # Credible intervals are TRUE posterior CIs from 500 MC samples (10 GMM × 50 posterior)
  output$scatter_delta_eur <- renderPlotly({
    pc4_comp <- selected_pc4_comparison()
    p <- ggplot(pc4_comp, aes(x = mean_delta_psi_virgin, y = eur_per_building)) +
      geom_segment(aes(x = p10_delta_psi_virgin, xend = p90_delta_psi_virgin,
                       y = eur_per_building, yend = eur_per_building),
                   alpha = 0.3, color = "#2ecc71", linewidth = 0.5) +
      geom_point(aes(size = n, text = paste0("PC4: ", pc4, "\nBuildings: ", n,
                                              "\nMean Ψ: ", round(mean_delta_psi_virgin, 3),
                                              "\n80% CI: [", round(p10_delta_psi_virgin, 3), "-", round(p90_delta_psi_virgin, 3), "]",
                                              "\n€/bldg: €", format(round(eur_per_building), big.mark = ","))),
                 alpha = 0.6, color = "#2ecc71") +
      geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", linetype = "dashed") +
      scale_y_continuous(labels = scales::dollar_format(prefix = "€", big.mark = ",")) +
      labs(x = "Predicted Mean Ψ (bars: 80% posterior CI)", y = "Observed € per Building") +
      theme_minimal() +
      theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })

  output$scatter_psi1_eur <- renderPlotly({
    pc4_comp <- selected_pc4_comparison()
    p <- ggplot(pc4_comp, aes(x = p_psi_1_virgin, y = eur_per_building,
                                     text = paste0("PC4: ", pc4, "\nBuildings: ", n,
                                                   "\nP(Ψ≥1): ", round(100*p_psi_1_virgin, 1), "%",
                                                   "\nΨ [p10-p90]: ", round(p10_psi_virgin, 2), "-", round(p90_psi_virgin, 2),
                                                   "\n€/bldg: €", format(round(eur_per_building), big.mark = ",")))) +
      geom_point(aes(size = n), alpha = 0.6, color = "#2ecc71") +
      geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", linetype = "dashed") +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::dollar_format(prefix = "€", big.mark = ",")) +
      labs(x = "Predicted P(Ψ ≥ 1.0)", y = "Observed € per Building") +
      theme_minimal() +
      theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })

  output$scatter_psi25_eur <- renderPlotly({
    pc4_comp <- selected_pc4_comparison()
    p <- ggplot(pc4_comp, aes(x = p_psi_2_5_virgin, y = eur_per_building,
                                     text = paste0("PC4: ", pc4, "\nBuildings: ", n,
                                                   "\nP(Ψ≥2.5): ", round(100*p_psi_2_5_virgin, 1), "%",
                                                   "\nΨ [p10-p90]: ", round(p10_psi_virgin, 2), "-", round(p90_psi_virgin, 2),
                                                   "\n€/bldg: €", format(round(eur_per_building), big.mark = ",")))) +
      geom_point(aes(size = n), alpha = 0.6, color = "#2ecc71") +
      geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", linetype = "dashed") +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::dollar_format(prefix = "€", big.mark = ",")) +
      labs(x = "Predicted P(Ψ ≥ 2.5)", y = "Observed € per Building") +
      theme_minimal() +
      theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })

  output$scatter_pgv_eur <- renderPlotly({
    pc4_comp <- selected_pc4_comparison()
    p <- ggplot(pc4_comp, aes(x = mean_max_pgv, y = eur_per_building)) +
      geom_point(aes(size = n, text = paste0("PC4: ", pc4, "\nBuildings: ", n,
                                              "\nMean PGV: ", round(mean_max_pgv, 1), " mm/s",
                                              "\n€/bldg: €", format(round(eur_per_building), big.mark = ","))),
                 alpha = 0.6, color = "#2ecc71") +
      geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", linetype = "dashed") +
      scale_y_continuous(labels = scales::dollar_format(prefix = "€", big.mark = ",")) +
      labs(x = "Mean Max PGV (mm/s)", y = "Observed € per Building") +
      theme_minimal() +
      theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })

  output$correlation_summary <- renderPrint({
    pc4_comp <- selected_pc4_comparison()
    valid_data <- pc4_comp |> filter(eur_per_building > 0)

    cat("=== Correlation Analysis (", input$model_select, " model) ===\n")
    cat("PC4s with compensation > 0: n =", nrow(valid_data), "\n\n")

    cor_delta <- cor(valid_data$mean_delta_psi_virgin, valid_data$eur_per_building, use = "complete.obs")
    cor_psi1 <- cor(valid_data$p_psi_1_virgin, valid_data$eur_per_building, use = "complete.obs")
    cor_psi25 <- cor(valid_data$p_psi_2_5_virgin, valid_data$eur_per_building, use = "complete.obs")
    cor_pgv <- cor(valid_data$mean_max_pgv, valid_data$eur_per_building, use = "complete.obs")

    cat("Pearson Correlations with € per Building:\n")
    cat("  Mean ΔΨ:       r =", round(cor_delta, 3), "\n")
    cat("  P(Ψ ≥ 1.0):    r =", round(cor_psi1, 3), "\n")
    cat("  P(Ψ ≥ 2.5):    r =", round(cor_psi25, 3), "\n")
    cat("  Mean Max PGV:  r =", round(cor_pgv, 3), "\n\n")

    spear_delta <- cor(valid_data$mean_delta_psi_virgin, valid_data$eur_per_building, method = "spearman", use = "complete.obs")
    spear_psi1 <- cor(valid_data$p_psi_1_virgin, valid_data$eur_per_building, method = "spearman", use = "complete.obs")

    cat("Spearman (Rank) Correlations:\n")
    cat("  Mean ΔΨ:       ρ =", round(spear_delta, 3), "\n")
    cat("  P(Ψ ≥ 1.0):    ρ =", round(spear_psi1, 3), "\n")
  })

  # --- MODEL COMPARISON ---
  output$model_comparison_table <- renderTable({
    tibble(
      Model = c("Additive (naive Markov)", "DAF (Korswagen)"),
      `Mean Ψ` = c(
        mean(model_data$additive$psi_virgin_mean),
        mean(model_data$daf$psi_virgin_mean)
      ),
      `P(Visible)` = c(
        mean(model_data$additive$p_visible_virgin),
        mean(model_data$daf$p_visible_virgin)
      ),
      `P(Moderate)` = c(
        mean(model_data$additive$p_moderate_virgin),
        mean(model_data$daf$p_moderate_virgin)
      ),
      `DAF/Additive` = c(
        1.0,
        mean(model_data$daf$psi_virgin_mean) / mean(model_data$additive$psi_virgin_mean)
      )
    ) |>
      mutate(
        `Mean Ψ` = round(`Mean Ψ`, 4),
        `P(Visible)` = paste0(round(100 * `P(Visible)`, 2), "%"),
        `P(Moderate)` = paste0(round(100 * `P(Moderate)`, 2), "%"),
        `DAF/Additive` = paste0(round(`DAF/Additive`, 2), "×")
      )
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$compare_zone_plot <- renderPlotly({
    plot_data <- bind_rows(
      summaries_all$additive$by_zone |> mutate(Model = "Additive"),
      summaries_all$daf$by_zone |> mutate(Model = "DAF")
    ) |>
      mutate(
        zone = factor(zone, levels = c("Near (<10km)", "Mid (10-25km)", "Far (>25km)")),
        Model = factor(Model, levels = c("Additive", "DAF"))
      )

    p <- ggplot(plot_data, aes(x = zone, y = psi_mean, fill = Model)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("Additive" = "#e74c3c", "DAF" = "#3498db")) +
      labs(x = NULL, y = "Mean Ψ (virgin)") +
      theme_minimal()

    ggplotly(p)
  })

  output$compare_visible_plot <- renderPlotly({
    plot_data <- bind_rows(
      summaries_all$additive$by_zone |> mutate(Model = "Additive"),
      summaries_all$daf$by_zone |> mutate(Model = "DAF")
    ) |>
      mutate(
        zone = factor(zone, levels = c("Near (<10km)", "Mid (10-25km)", "Far (>25km)")),
        Model = factor(Model, levels = c("Additive", "DAF"))
      )

    p <- ggplot(plot_data, aes(x = zone, y = p_visible_virgin, fill = Model)) +
      geom_col(position = "dodge") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = c("Additive" = "#e74c3c", "DAF" = "#3498db")) +
      labs(x = NULL, y = "P(Visible Damage)") +
      theme_minimal()

    ggplotly(p)
  })

  output$compare_histogram <- renderPlotly({
    set.seed(42)
    n_sample <- 30000
    plot_data <- bind_rows(
      model_data$additive |> sample_n(n_sample) |> select(psi_virgin_mean) |> mutate(Model = "Additive"),
      model_data$daf |> sample_n(n_sample) |> select(psi_virgin_mean) |> mutate(Model = "DAF")
    ) |>
      mutate(Model = factor(Model, levels = c("Additive", "DAF")))

    p <- ggplot(plot_data, aes(x = psi_virgin_mean, fill = Model)) +
      geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
      scale_fill_manual(values = c("Additive" = "#e74c3c", "DAF" = "#3498db")) +
      labs(x = "Cumulative Ψ (virgin)", y = "Count") +
      xlim(0, 4) +
      theme_minimal()

    ggplotly(p)
  })

  # --- DOWNLOADS ---
  output$download_buildings <- downloadHandler(
    filename = function() {
      paste0("groningen_building_damage_", input$model_select, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(selected_results(), file)
    }
  )

  output$download_postcodes <- downloadHandler(
    filename = function() {
      paste0("groningen_postcode_summary_", input$model_select, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(selected_pc4_damage(), file)
    }
  )

  output$download_earthquakes <- downloadHandler(
    filename = function() {
      paste0("groningen_earthquakes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(earthquakes, file)
    }
  )
}

# =============================================================================
# RUN APP
# =============================================================================

shinyApp(ui = ui, server = server)
