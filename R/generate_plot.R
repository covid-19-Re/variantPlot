library(plotly)
library(tidyverse)
library(lubridate)


generate_plot <- function(output_dir) {
  system(paste0("mkdir -p ", output_dir))

  data_raw <- read_csv(
    "https://raw.githubusercontent.com/covid-19-Re/variantPlot/master/data/data.csv",
    col_types = cols(
      year = col_integer(),
      week = col_integer(),
      lab = col_factor(
        levels = c("Viollier", "Risch", "Total")
      ),
      n = col_integer(),
      b117 = col_integer(),
      s501yv2 = col_integer()
    )
  )

  data <- data_raw %>%
    pivot_longer(c(b117, s501yv2), names_to = "variant", values_to = "k") %>%
    drop_na(k)

  # Prepare data
  total <- data %>%
    group_by(year, week, variant) %>%
    summarize(
      lab = "Total",
      n = sum(n),
      k = sum(k),
      .groups = "drop"
    ) %>%
    # TODO Currently, only Viollier's sequences will be used for SA variant detection
    # so that it does not make sense to plot the total.
    filter(variant != "s501yv2")
  data <- rbind(data, total)

  ci <- binom::binom.confint(data$k, data$n, methods = "wilson")
  data <- data %>%
    mutate(
      p = round(ci$mean, digits = 4) * 100,
      p_lower = round(ci$lower, digits = 4) * 100,
      p_upper = round(ci$upper, digits = 4) * 100,
      year_week = paste0(year, "-", week)
    )

  gen <- function(variant_name, plot_title) {
    fig <- plot_ly(
      data %>% filter(variant == variant_name), type = "box", x = ~year_week, color = ~lab,
      median = ~p, q1 = ~p, q3 = ~p,
      lowerfence = ~p_lower, upperfence = ~p_upper,
      hovertemplate = "x",
      text = ~p
    ) %>% layout(
      title = plot_title,
      boxmode = "group",
      xaxis = list(
        title = "Calendar Week"
      ),
      yaxis = list(
        title = "Estimated Percentage of New Variant",
        ticksuffix = "%"
      ),
      showlegend = TRUE,
      legend = list(
        title = list(
          text = "<b>Diagnostic Lab</b>"
        )
      )
    ) %>% config(
      displaylogo = FALSE,
      modeBarButtons = list(list("zoom2d", "toImage", "resetScale2d", "pan2d")),
      toImageButtonOptions = list(format = "png", width = 1200, height = 800, scale = 1)
    )
    return(fig)
  }

  fig_b117 <- gen("b117", "B.1.1.7 Variant in Switzerland")
  fig_s501yv2 <- gen("s501yv2", "S.501Y.V2 Variant in Switzerland")

  # Save to file
  htmlwidgets::saveWidget(fig_b117, paste0(output_dir, "/variantPlot_b117.html"), selfcontained = FALSE,
                          title = "B.1.1.7 Variant in Switzerland", libdir = paste0(output_dir, "/lib_b117"))
  htmlwidgets::saveWidget(fig_s501yv2, paste0(output_dir, "/variantPlot_s501yv2.html"), selfcontained = FALSE,
                          title = "S.501Y.V2 Variant in Switzerland", libdir = paste0(output_dir, "/lib_s501yv2"))
  # save data for Lagebericht
  write_csv(data, paste0(output_dir, "/variantPlot_data.csv"))
}


args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("Please provide an output path.")
}
generate_plot(args[1])
