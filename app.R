# =============================================================================
# Groningen Building Damage Explorer
#
# Interactive Shiny app matching Groningen_analyse style
# Shows claims vs predicted damage with side-by-side maps and scatter plots
# =============================================================================

library(shiny)
library(tidyverse)
library(leaflet)
library(DT)
library(plotly)
library(sf)
library(here)
library(scales)

# =============================================================================
# GRONINGEN COLOR PALETTE & THEME
# =============================================================================

groningen_colors <- list(
  # Primary data colors (colorblind-safe)
  claims      = "#7A4070",   # Warm purple/magenta - observed/claims
  predicted   = "#207068",   # Dark teal/cyan - predicted/model
  accent      = "#C96B4C",   # Terracotta - emphasis
  highlight   = "#D4A84B",   # Warm gold - highlights

  # UI colors
  background  = "#F7F5F2",   # Warm off-white
  panel       = "#FFFFFF",   # Pure white
  panel_alt   = "#F0EDE8",   # Alternating panel background
  border      = "#E5E2DE",   # Warm gray border
  grid        = "#E8E8E8",   # Grid lines (subtle)

  # Text colors
  text_dark   = "#2D2D2D",   # Near-black for headers
  text_body   = "#525252",   # Dark gray for body
  text_muted  = "#6B6B6B",   # Medium gray for captions

  # Extended palette
  claims_light    = "#A87BA0",   # Lighter purple for ribbons
  predicted_light = "#5A9A92"    # Lighter teal for ribbons
)

# ggplot2 theme matching Groningen_analyse style
corriethema <- theme(
  legend.title = element_blank(),
  legend.position = "top",
  legend.text = element_text(family = "sans", size = 9, color = groningen_colors$text_muted),
  legend.justification = 'left',
  legend.direction = 'horizontal',
  legend.key.size = unit(0.8, "lines"),
  strip.text = element_text(size = 10, family = "serif", color = groningen_colors$text_body),
  plot.background = element_rect(fill = groningen_colors$background, linetype = "blank"),
  panel.background = element_rect(fill = groningen_colors$panel),
  panel.grid.major = element_line(colour = groningen_colors$grid, linetype = "dashed", linewidth = 0.3),
  panel.grid.minor = element_blank(),
  plot.title = element_text(
    margin = margin(t = 5, r = 10, b = 5, l = 0),
    size = 14,
    family = "sans",
    face = "bold",
    color = groningen_colors$text_dark,
    lineheight = 1.1
  ),
  strip.background = element_blank(),
  plot.title.position = "plot",
  plot.subtitle = element_text(
    margin = margin(t = 2, r = 10, b = 15, l = 0),
    size = 11,
    family = "sans",
    colour = groningen_colors$text_body,
    lineheight = 1.1
  ),
  plot.caption = element_text(
    size = 8,
    family = "serif",
    face = "italic",
    colour = groningen_colors$text_muted
  ),
  axis.title.x = element_text(
    margin = margin(t = 10, r = 0, b = 5, l = 0),
    size = 10,
    family = "serif",
    colour = groningen_colors$text_body
  ),
  axis.title.y = element_text(
    margin = margin(t = 0, r = 10, b = 0, l = 5),
    size = 10,
    family = "serif",
    colour = groningen_colors$text_body
  ),
  axis.text = element_text(size = 9, colour = groningen_colors$text_muted),
  legend.background = element_rect(fill = groningen_colors$background),
  legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
  plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
)

# CSS for Shiny styling
shiny_css <- paste0('
body { background-color: ', groningen_colors$background, '; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }
h4 { margin-top: 20px; margin-bottom: 10px; color: ', groningen_colors$text_dark, '; }
h2 { color: ', groningen_colors$text_dark, '; font-family: sans-serif; }
.well { background-color: #FFFFFF; border: 1px solid #E5E2DE; border-radius: 8px; }
.info-box { background: #F7F5F2; padding: 12px; border-radius: 6px; margin-top: 12px; font-size: 12px; color: #525252; border-left: 3px solid #207068; }
')

# =============================================================================
# LOAD DATA
# =============================================================================

# Load model results
model_data <- list(
  additive = readRDS(here("outputs", "models", "spatial_damage_uncertainty_full.rds")),
  daf = readRDS(here("outputs", "models", "spatial_damage_daf.rds"))
)

earthquakes <- readRDS(here("outputs", "models", "groningen_earthquakes.rds"))

# Function to compute PC4 aggregates from building data
compute_pc4_damage <- function(data) {
  data |>
    group_by(pc4) |>
    summarise(
      n_buildings = n(),
      mean_psi = mean(psi_virgin_mean, na.rm = TRUE),
      p10_psi = mean(psi_virgin_p10, na.rm = TRUE),
      p90_psi = mean(psi_virgin_p90, na.rm = TRUE),
      p_visible = mean(p_visible_virgin, na.rm = TRUE),
      p_moderate = mean(p_moderate_virgin, na.rm = TRUE),
      mean_max_pgv = mean(max_pgv_mean, na.rm = TRUE),
      .groups = "drop"
    )
}

# Load PC4 polygon boundaries
pc4_polygons <- st_read(here("datafiles", "pc4_groningen.gpkg"), quiet = TRUE)

# Pre-compute PC4 aggregates for all models
pc4_damage_all <- lapply(model_data, compute_pc4_damage)

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

# Function to join claims to PC4 predictions and create map data
create_map_data <- function(pc4_damage) {
  pc4_polygons |>
    left_join(pc4_damage, by = "pc4") |>
    filter(!is.na(n_buildings), n_buildings > 0) |>
    left_join(claims_pc4, by = "pc4") |>
    mutate(
      total_claims = replace_na(total_claims, 0),
      total_compensation_eur = replace_na(total_compensation_eur, 0),
      # Per 1000 buildings
      claims_per_1000 = (total_claims / n_buildings) * 1000,
      predicted_per_1000 = p_visible * 1000,  # P(visible) * 1000
      eur_per_building = total_compensation_eur / n_buildings
    )
}

# Pre-compute map data for all models
pc4_map_all <- lapply(pc4_damage_all, create_map_data)

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML(shiny_css))),

  titlePanel("Groningen Schadeclaims vs Voorspelde Schade"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Model Parameters"),

      selectInput("model_select", "Accumulatiemodel:",
        choices = c(
          "Additive (naive Markov)" = "additive",
          "Korswagen DAF (N-effect)" = "daf"
        ),
        selected = "daf"
      ),
      helpText(
        style = "font-size: 11px; color: #666;",
        HTML("<b>Additive:</b> State-dependent, elk event telt onafhankelijk<br>
              <b>DAF:</b> History-aware met N-effect voor vergelijkbare events")
      ),

      hr(),

      h4("Kaart Opties"),

      selectInput("damage_metric", "Schademaat:",
        choices = c(
          "P(Zichtbaar) Ψ ≥ 1" = "p_visible",
          "Gemiddelde Ψ" = "mean_psi"
        ),
        selected = "p_visible"
      ),

      hr(),

      h4("Samenvatting"),
      verbatimTextOutput("stats_summary"),

      hr(),

      div(
        class = "info-box",
        p(strong("Data bronnen:"), br(),
          em("Claims:"), " IMG fysieke schade 2020-2024", br(),
          em("Gebouwen:"), " BAG 2024 (542.957 panden)", br(),
          em("Aardbevingen:"), " KNMI M≥2.0 (1991-2024)")
      )
    ),

    mainPanel(
      width = 9,

      tabsetPanel(
        id = "main_tabs",

        # TAB 1: Maps
        tabPanel(
          "Kaarten",
          value = "maps",

          fluidRow(
            column(6,
              h5("Schadeclaims per 1.000 gebouwen", style = "text-align: center; color: #7A4070;"),
              leafletOutput("map_claims", height = "450px")
            ),
            column(6,
              h5("Voorspeld per 1.000 gebouwen", style = "text-align: center; color: #207068;"),
              leafletOutput("map_predicted", height = "450px")
            )
          ),

          fluidRow(
            column(12,
              plotOutput("scatter_plot", height = "450px", width = "100%")
            )
          )
        ),

        # TAB 2: Euro comparison
        tabPanel(
          "Euro's per PC4",
          value = "euros",

          fluidRow(
            column(6,
              h5("Uitgekeerd per gebouw (€)", style = "text-align: center; color: #7A4070;"),
              leafletOutput("euro_map_claims", height = "450px")
            ),
            column(6,
              h5("Voorspelde schade (Ψ)", style = "text-align: center; color: #207068;"),
              leafletOutput("euro_map_predicted", height = "450px")
            )
          ),

          fluidRow(
            column(12,
              plotOutput("euro_scatter", height = "450px", width = "100%")
            )
          )
        ),

        # TAB 3: Model Comparison
        tabPanel(
          "Model Vergelijking",
          value = "comparison",

          h4("Vergelijking: Additive vs DAF"),
          p("Alle 542.957 gebouwen, virgin walls scenario."),

          fluidRow(
            column(12,
              tableOutput("model_comparison_table")
            )
          ),

          fluidRow(
            column(12,
              plotOutput("compare_histogram", height = "350px")
            )
          ),

          fluidRow(
            column(6,
              plotOutput("compare_scatter", height = "400px")
            ),
            column(6,
              verbatimTextOutput("correlation_summary")
            )
          )
        ),

        # TAB 4: Earthquakes
        tabPanel(
          "Aardbevingen",
          value = "earthquakes",

          fluidRow(
            column(8,
              h4("Catalogus (M ≥ 2.0)"),
              DTOutput("eq_table")
            ),
            column(4,
              h4("Magnitude Verdeling"),
              plotlyOutput("eq_histogram", height = "250px"),
              h4("Tijdlijn"),
              plotlyOutput("eq_timeline", height = "250px")
            )
          )
        ),

        # TAB 5: Download
        tabPanel(
          "Download",
          value = "download",

          h4("Download Resultaten"),
          p("Download de volledige resultaten voor verdere analyse."),

          fluidRow(
            column(4,
              downloadButton("download_buildings", "Gebouwen CSV (542k rijen)")
            ),
            column(4,
              downloadButton("download_postcodes", "Postcode Samenvatting CSV")
            ),
            column(4,
              downloadButton("download_earthquakes", "Aardbevingen CSV")
            )
          )
        )
      ),

      fluidRow(
        column(12,
          div(
            style = "margin-top: 20px; font-size: 12px; color: #666;",
            p(em("brms hurdle-gamma model | Korswagen fragility curves | Bommer GMM"))
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # --- REACTIVE DATA ---
  selected_results <- reactive({
    model_data[[input$model_select]]
  })

  selected_map_data <- reactive({
    pc4_map_all[[input$model_select]]
  })

  # Common max value for both maps (so they use same scale)
  map_max_val <- reactive({
    data <- selected_map_data()
    max(
      quantile(data$claims_per_1000, 0.95, na.rm = TRUE),
      quantile(data$predicted_per_1000, 0.95, na.rm = TRUE),
      na.rm = TRUE
    )
  })

  # Color palette function
  get_color_pal <- function(max_val) {
    colorNumeric(
      palette = "magma",
      domain = c(0, max_val),
      reverse = TRUE,
      na.color = "#f0f0f0"
    )
  }

  # --- STATISTICS SUMMARY ---
  output$stats_summary <- renderText({
    data <- selected_map_data()
    results <- selected_results()

    total_claims <- sum(data$total_claims, na.rm = TRUE)
    total_predicted <- sum(data$n_buildings * data$p_visible, na.rm = TRUE)
    total_buildings <- sum(data$n_buildings, na.rm = TRUE)

    ratio <- if (total_predicted > 0) total_claims / total_predicted else NA

    # Correlation
    cor_val <- cor(data$predicted_per_1000, data$claims_per_1000, use = "complete.obs")

    paste0(
      "=== ", toupper(input$model_select), " MODEL ===\n",
      "\nClaims (2020-2024): ", format(total_claims, big.mark = "."), "\n",
      "Voorspeld (P visible): ", format(round(total_predicted), big.mark = "."), "\n",
      "Ratio claims/voorspeld: ", if (!is.na(ratio)) sprintf("%.2fx", ratio) else "N/A", "\n",
      "\nCorrelatie (r): ", sprintf("%.3f", cor_val), "\n",
      "\nGebouwen: ", format(total_buildings, big.mark = "."), "\n",
      "PC4 gebieden: ", nrow(data), "\n",
      "Aardbevingen: ", nrow(earthquakes)
    )
  })

  # --- CLAIMS MAP ---
  output$map_claims <- renderLeaflet({
    data <- selected_map_data()
    data_sf <- st_transform(data, 4326)

    max_val <- map_max_val()
    pal <- get_color_pal(max_val)

    # Popup
    data_sf$popup <- sprintf(
      "<strong>PC4: %s</strong><br/>
       Claims: %.1f per 1.000 geb.<br/>
       Totaal claims: %s<br/>
       Gebouwen: %s",
      data_sf$pc4,
      round(data_sf$claims_per_1000, 1),
      format(data_sf$total_claims, big.mark = "."),
      format(data_sf$n_buildings, big.mark = ".")
    )

    leaflet(data_sf) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addPolygons(
        fillColor = ~pal(pmin(claims_per_1000, max_val)),
        fillOpacity = 0.7,
        weight = 0.5,
        color = "grey",
        popup = ~popup,
        highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.9, bringToFront = TRUE)
      ) |>
      addLegend(
        position = "bottomleft",
        pal = pal,
        values = c(0, max_val),
        title = "Per 1.000 geb.",
        opacity = 0.7
      ) |>
      addControl(
        html = "<strong>Schadeclaims (2020-2024)</strong>",
        position = "bottomright"
      )
  })

  # --- PREDICTED MAP ---
  output$map_predicted <- renderLeaflet({
    data <- selected_map_data()
    data_sf <- st_transform(data, 4326)

    max_val <- map_max_val()
    pal <- get_color_pal(max_val)

    model_label <- if (input$model_select == "daf") "DAF" else "Additive"

    # Popup
    data_sf$popup <- sprintf(
      "<strong>PC4: %s</strong><br/>
       Voorspeld: %.1f per 1.000 geb.<br/>
       P(visible): %.1f%%<br/>
       Gemiddelde Psi: %.3f<br/>
       Gebouwen: %s",
      data_sf$pc4,
      round(data_sf$predicted_per_1000, 1),
      round(data_sf$p_visible * 100, 1),
      round(data_sf$mean_psi, 3),
      format(data_sf$n_buildings, big.mark = ".")
    )

    leaflet(data_sf) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addPolygons(
        fillColor = ~pal(pmin(predicted_per_1000, max_val)),
        fillOpacity = 0.7,
        weight = 0.5,
        color = "grey",
        popup = ~popup,
        highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.9, bringToFront = TRUE)
      ) |>
      addLegend(
        position = "bottomleft",
        pal = pal,
        values = c(0, max_val),
        title = "Per 1.000 geb.",
        opacity = 0.7
      ) |>
      addControl(
        html = sprintf("<strong>%s Model (1991-2024)</strong>", model_label),
        position = "bottomright"
      )
  })

  # --- SCATTER PLOT ---
  output$scatter_plot <- renderPlot({
    data <- as.data.frame(selected_map_data())
    data <- data[!is.na(data$claims_per_1000) & !is.na(data$predicted_per_1000), ]

    if (nrow(data) == 0) return(NULL)

    cor_val <- cor(data$predicted_per_1000, data$claims_per_1000, use = "complete.obs")
    max_val <- max(c(data$claims_per_1000, data$predicted_per_1000), na.rm = TRUE)

    model_label <- if (input$model_select == "daf") "DAF" else "Additive"

    ggplot(data, aes(x = predicted_per_1000, y = claims_per_1000)) +
      geom_point(aes(size = n_buildings), alpha = 0.6, color = groningen_colors$claims) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = groningen_colors$text_dark, linewidth = 0.8) +
      geom_smooth(method = "lm", se = FALSE, color = groningen_colors$text_dark, linewidth = 0.8) +
      scale_size_continuous(range = c(3, 14), guide = "none") +
      scale_x_continuous(limits = c(0, max_val * 1.05), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, max_val * 1.05), expand = c(0, 0)) +
      labs(
        title = "Claims vs Voorspelling",
        subtitle = sprintf("%s model | r = %.3f", model_label, cor_val),
        x = "Voorspeld per 1.000 gebouwen",
        y = "Claims per 1.000 gebouwen",
        caption = "Elke punt = PC4 | Grootte ~ gebouwen | Streepjeslijn = perfecte match"
      ) +
      corriethema
  })

  # --- EURO MAPS ---
  euro_max_val <- reactive({
    data <- selected_map_data()
    quantile(data$eur_per_building[data$eur_per_building > 0], 0.95, na.rm = TRUE)
  })

  output$euro_map_claims <- renderLeaflet({
    data <- selected_map_data()
    data_sf <- st_transform(data, 4326)

    max_val <- euro_max_val()
    pal <- colorNumeric(palette = "plasma", domain = c(0, max_val), reverse = TRUE, na.color = "#f0f0f0")

    data_sf$popup <- sprintf(
      "<strong>PC4: %s</strong><br/>
       Uitgekeerd: EUR %s per gebouw<br/>
       Totaal: EUR %s<br/>
       Gebouwen: %s",
      data_sf$pc4,
      format(round(data_sf$eur_per_building), big.mark = "."),
      format(round(data_sf$total_compensation_eur), big.mark = "."),
      format(data_sf$n_buildings, big.mark = ".")
    )

    leaflet(data_sf) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addPolygons(
        fillColor = ~pal(pmin(eur_per_building, max_val)),
        fillOpacity = 0.7,
        weight = 0.5,
        color = "grey",
        popup = ~popup,
        highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.9, bringToFront = TRUE)
      ) |>
      addLegend(
        position = "bottomleft",
        pal = pal,
        values = c(0, max_val),
        title = "EUR/gebouw",
        opacity = 0.7
      )
  })

  output$euro_map_predicted <- renderLeaflet({
    data <- selected_map_data()
    data_sf <- st_transform(data, 4326)

    # Use mean_psi for predicted damage
    max_val <- quantile(data_sf$mean_psi, 0.95, na.rm = TRUE)
    pal <- colorNumeric(palette = "YlOrRd", domain = c(0, max_val), na.color = "#f0f0f0")

    model_label <- if (input$model_select == "daf") "DAF" else "Additive"

    data_sf$popup <- sprintf(
      "<strong>PC4: %s</strong><br/>
       Gemiddelde Psi: %.3f<br/>
       P(visible): %.1f%%<br/>
       Gebouwen: %s",
      data_sf$pc4,
      round(data_sf$mean_psi, 3),
      round(data_sf$p_visible * 100, 1),
      format(data_sf$n_buildings, big.mark = ".")
    )

    leaflet(data_sf) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addPolygons(
        fillColor = ~pal(pmin(mean_psi, max_val)),
        fillOpacity = 0.7,
        weight = 0.5,
        color = "grey",
        popup = ~popup,
        highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.9, bringToFront = TRUE)
      ) |>
      addLegend(
        position = "bottomleft",
        pal = pal,
        values = c(0, max_val),
        title = "Mean Psi",
        opacity = 0.7
      ) |>
      addControl(
        html = sprintf("<strong>%s Model</strong>", model_label),
        position = "bottomright"
      )
  })

  output$euro_scatter <- renderPlot({
    data <- as.data.frame(selected_map_data())
    data <- data[data$eur_per_building > 0, ]

    if (nrow(data) == 0) return(NULL)

    cor_val <- cor(data$mean_psi, data$eur_per_building, use = "complete.obs")
    model_label <- if (input$model_select == "daf") "DAF" else "Additive"

    ggplot(data, aes(x = mean_psi, y = eur_per_building)) +
      geom_point(aes(size = n_buildings), alpha = 0.6, color = groningen_colors$predicted) +
      geom_smooth(method = "lm", se = TRUE, color = groningen_colors$text_dark, linetype = "dashed") +
      scale_size_continuous(range = c(3, 14), guide = "none") +
      scale_y_continuous(labels = dollar_format(prefix = "EUR ", big.mark = ".")) +
      labs(
        title = "Uitgekeerde Schade vs Voorspelde Schade",
        subtitle = sprintf("%s model | r = %.3f", model_label, cor_val),
        x = "Voorspelde gemiddelde Psi",
        y = "Uitgekeerd per gebouw (EUR)",
        caption = "Elke punt = PC4 | Grootte ~ gebouwen"
      ) +
      corriethema
  })

  # --- MODEL COMPARISON ---
  output$model_comparison_table <- renderTable({
    tibble(
      Model = c("Additive (naive Markov)", "DAF (Korswagen)"),
      `Mean Psi` = c(
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
        `Mean Psi` = round(`Mean Psi`, 4),
        `P(Visible)` = paste0(round(100 * `P(Visible)`, 2), "%"),
        `P(Moderate)` = paste0(round(100 * `P(Moderate)`, 2), "%"),
        `DAF/Additive` = paste0(round(`DAF/Additive`, 2), "x")
      )
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$compare_histogram <- renderPlot({
    set.seed(42)
    n_sample <- 30000
    plot_data <- bind_rows(
      model_data$additive |> sample_n(n_sample) |> select(psi_virgin_mean) |> mutate(Model = "Additive"),
      model_data$daf |> sample_n(n_sample) |> select(psi_virgin_mean) |> mutate(Model = "DAF")
    ) |>
      mutate(Model = factor(Model, levels = c("Additive", "DAF")))

    ggplot(plot_data, aes(x = psi_virgin_mean, fill = Model)) +
      geom_histogram(binwidth = 0.05, alpha = 0.6, position = "identity") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
      scale_fill_manual(values = c("Additive" = groningen_colors$claims, "DAF" = groningen_colors$predicted)) +
      labs(
        title = "Verdeling Cumulatieve Schade",
        x = "Cumulatieve Psi (virgin walls)",
        y = "Aantal gebouwen",
        caption = "Verticale lijn = zichtbare schade drempel (Psi = 1)"
      ) +
      xlim(0, 2) +
      corriethema
  })

  output$compare_scatter <- renderPlot({
    # Merge model predictions at PC4 level
    additive_pc4 <- pc4_damage_all$additive |> select(pc4, psi_additive = mean_psi, p_visible_additive = p_visible)
    daf_pc4 <- pc4_damage_all$daf |> select(pc4, psi_daf = mean_psi, p_visible_daf = p_visible)

    compare_data <- inner_join(additive_pc4, daf_pc4, by = "pc4")

    cor_val <- cor(compare_data$psi_additive, compare_data$psi_daf, use = "complete.obs")

    ggplot(compare_data, aes(x = psi_additive, y = psi_daf)) +
      geom_point(alpha = 0.5, color = groningen_colors$accent) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = groningen_colors$text_dark) +
      labs(
        title = "DAF vs Additive per PC4",
        subtitle = sprintf("r = %.3f", cor_val),
        x = "Additive Mean Psi",
        y = "DAF Mean Psi"
      ) +
      corriethema
  })

  output$correlation_summary <- renderPrint({
    additive_map <- pc4_map_all$additive
    daf_map <- pc4_map_all$daf

    cat("=== CORRELATIE MET CLAIMS ===\n\n")

    # Additive
    cor_add_claims <- cor(additive_map$predicted_per_1000, additive_map$claims_per_1000, use = "complete.obs")
    cor_add_eur <- cor(additive_map$mean_psi, additive_map$eur_per_building, use = "complete.obs")

    cat("ADDITIVE MODEL:\n")
    cat("  P(visible) vs claims/1000: r =", round(cor_add_claims, 3), "\n")
    cat("  Mean Psi vs EUR/gebouw:    r =", round(cor_add_eur, 3), "\n\n")

    # DAF
    cor_daf_claims <- cor(daf_map$predicted_per_1000, daf_map$claims_per_1000, use = "complete.obs")
    cor_daf_eur <- cor(daf_map$mean_psi, daf_map$eur_per_building, use = "complete.obs")

    cat("DAF MODEL:\n")
    cat("  P(visible) vs claims/1000: r =", round(cor_daf_claims, 3), "\n")
    cat("  Mean Psi vs EUR/gebouw:    r =", round(cor_daf_eur, 3), "\n\n")

    # Model comparison
    cat("=== MODEL VERGELIJKING ===\n\n")
    cat("Mean Psi (Additive):  ", round(mean(model_data$additive$psi_virgin_mean), 4), "\n")
    cat("Mean Psi (DAF):       ", round(mean(model_data$daf$psi_virgin_mean), 4), "\n")
    cat("DAF/Additive ratio:   ", round(mean(model_data$daf$psi_virgin_mean) / mean(model_data$additive$psi_virgin_mean), 2), "x\n")
  })

  # --- EARTHQUAKES ---
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
        colnames = c("Datum/Tijd", "Magnitude", "Lat", "Lon", "Diepte (km)")
      )
  })

  output$eq_histogram <- renderPlotly({
    p <- ggplot(earthquakes, aes(x = mag)) +
      geom_histogram(binwidth = 0.1, fill = groningen_colors$accent, color = "white") +
      labs(x = "Magnitude", y = "Aantal") +
      theme_minimal()
    ggplotly(p)
  })

  output$eq_timeline <- renderPlotly({
    eq_year <- earthquakes |>
      mutate(year = year(time_utc)) |>
      count(year)

    p <- ggplot(eq_year, aes(x = year, y = n)) +
      geom_col(fill = groningen_colors$accent) +
      labs(x = "Jaar", y = "Events") +
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
      write_csv(pc4_damage_all[[input$model_select]], file)
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
