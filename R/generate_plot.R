library(plotly)
library(tidyverse)
library(lubridate)


db_connection <- NULL
get_database_connection <- function() {
  if (is.null(db_connection)) {
    db_config <- config::get(value = "database", file = here::here("config.yml"), use_parent = FALSE)
    db_connection <<- DBI::dbConnect(
      RPostgres::Postgres(),
      host = db_config$host,
      port = db_config$port,
      user = db_config$username,
      password = db_config$password,
      dbname = db_config$dbname
    )
  }
  return(db_connection)
}


generate_switzerland_plot <- function(output_dir, internal = FALSE) {
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
    )
  data <- rbind(data, total)

  ci <- binom::binom.confint(data$k, data$n, methods = "wilson")
  data <- data %>%
    mutate(
      p = round(ci$mean, digits = 4) * 100,
      p_lower = round(ci$lower, digits = 4) * 100,
      p_upper = round(ci$upper, digits = 4) * 100,
      year_week = paste0(year, "-", week)
    )

  if (internal) {
    db_connection <- get_database_connection()
    data_case_numbers <- DBI::dbGetQuery(db_connection, "
    select
      extract(isoyear from fall_dt)::int as year,
      extract(week from fall_dt)::int as week,
      count(*)::int as total_cases
    from bag_dashboard_meldeformular
    group by year, week;
  ") %>%
      mutate(year_week = paste0(year, "-", week)) %>%
      inner_join(data %>% filter(lab == "Total"), by = "year_week") %>%
      mutate(
        cases = as.integer(round(p * total_cases / 100)),
        cases_lower = as.integer(floor(p_lower * total_cases / 100)),
        cases_upper = as.integer(ceiling(p_upper * total_cases / 100))
      )
  }


  gen <- function(variant_name, plot_title) {
    fig <- plot_ly() %>%
      add_trace(
        data = data %>% filter(variant == variant_name), type = "box", x = ~year_week, color = ~lab,
        median = ~p, q1 = ~p, q3 = ~p,
        lowerfence = ~p_lower, upperfence = ~p_upper,
        hovertemplate = "x",
        text = ~p,
        alpha = 1
      )

    if (internal) {
      data_case_numbers <- data_case_numbers %>% filter(variant == variant_name)
      fig <- fig %>%
        add_trace(
          data = data_case_numbers, type = "scatter", mode = "lines",
          x = ~year_week, y = ~cases, yaxis = "y2",
          fill = 'tozeroy',
          fillcolor = "rgba(230,215,176,0.15)",
          line = list(
            color = 'rgba(230,215,176,0.3)'
          ),
          name = 'Estimated Case Number'
        ) %>%
        add_trace(
          data = data_case_numbers, type = "scatter", mode = "lines",
          x = ~year_week, y = ~cases_upper, yaxis = "y2",
          line = list(
            color = 'rgba(230,215,176,1)',
            dash = "dash"
          ),
          name = ''
          # , showlegend = FALSE
        ) %>%
        add_trace(
          data = data_case_numbers, type = "scatter", mode = "lines",
          x = ~year_week, y = ~cases_lower, yaxis = "y2",
          line = list(
            color = 'rgba(230,215,176,1)',
            dash = "dash"
          ),
          name = ''
          # , showlegend = FALSE
        )
    }

    yaxis2 <- list()
    if (internal) {
      yaxis2 <- list(
        overlaying = "y",
        side = "right",
        title = "Estimated Case Number of New Variant",
        range = c(0, max(data_case_numbers$cases_upper))
      )
    }

    fig <- fig %>%
      layout(
        title = plot_title,
        boxmode = "group",
        xaxis = list(
          title = "Calendar Week"
        ),
        yaxis = list(
          title = "Estimated Percentage of New Variant",
          ticksuffix = "%",
          range = c(0, max(data$p_upper))
        ),
        yaxis2 = yaxis2,
        showlegend = TRUE,
        legend = list(
          title = list(
            text = "<b>Diagnostic Lab</b>"
          )
        )
      ) %>%
      config(
        displaylogo = FALSE,
        modeBarButtons = list(list("zoom2d", "toImage", "resetScale2d", "pan2d")),
        toImageButtonOptions = list(format = "svg", width = 1600, height = 800, filename = plot_title)
      )
    return(fig)
  }

  fig_b117 <- gen("b117", "B.1.1.7 Variant in Switzerland")
  fig_s501yv2 <- gen("s501yv2", "S.501Y.V2 Variant in Switzerland")


  # Save to file
  suffix <- ""
  if (internal) {
    suffix <- "_internal"
  }

  htmlwidgets::saveWidget(fig_b117, paste0(output_dir, "/variantPlot_b117", suffix, ".html"), selfcontained = FALSE,
                          title = "B.1.1.7 Variant in Switzerland", libdir = paste0(output_dir, "/lib_b117", suffix))
  htmlwidgets::saveWidget(fig_s501yv2, paste0(output_dir, "/variantPlot_s501yv2", suffix, ".html"), selfcontained = FALSE,
                          title = "S.501Y.V2 Variant in Switzerland", libdir = paste0(output_dir, "/lib_s501yv2", suffix))

  # save data for Lagebericht
  write_csv(data, paste0(output_dir, "/variantPlot_data.csv"))
}


generate_comparison_plot <- function(output_dir, logscale = FALSE) {
  # Denmark https://www.ssi.dk/-/media/cdn/files/notat_engelsk_virusvariant090121.pdf?la=da
  # or https://www.covid19genomics.dk/statistics

  # UK S-dropout
  # https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
  # https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/950841/Variant_of_Concern_VOC_202012_01_Technical_Briefing_3_England_Data.ods
  # Fig 4 table

  data_raw <- read_csv(
    "https://raw.githubusercontent.com/covid-19-Re/variantPlot/master/data/data_comparison.csv",
    col_types = cols(
      country = col_character(),
      year = col_integer(),
      week = col_integer(),
      n = col_integer(),
      b117 = col_integer()
    )
  )

  data <- data_raw %>%
    mutate(
      date = ISOweek::ISOweek2date(paste0(year, "-W", str_pad(week, 2, pad = "0"), "-", 1))
    )

  ci <- binom::binom.confint(data$b117, data$n, methods = "wilson")
  data <- data %>%
    mutate(
      p = round(ci$mean, digits = 4) * 100,
      p_lower = round(ci$lower, digits = 4) * 100,
      p_upper = round(ci$upper, digits = 4) * 100
    ) %>%
    mutate(
      text = paste0("Week ", week, "\n", round(p, digits = 2), "% (", round(p_lower, digits = 2),
                    "%-", round(p_upper, digits = 2), "%)")
    )

  min_date <- min(data$date)
  max_date <- max(data$date)

  if (logscale) {
    data <- data %>%
      mutate(b117 = na_if(b117, 0)) %>%
      drop_na(b117)
  }
  fig <- plot_ly(data, x = ~date, y = ~p, color = ~country, colors = c("blue", "red", "black", "gray"),
                 type = "scatter", mode = "lines+markers", line = list(width = 2),
                 text = ~text,
                 hoverinfo = "text"
  ) %>%
    add_trace(x = ~date, y = ~p_lower, color = ~country,
              type = "scatter", mode = "lines", line = list(width = 1, dash = "dash"),
              showlegend = FALSE
    ) %>%
    add_trace(x = ~date, y = ~p_upper, color = ~country,
              type = "scatter", mode = "lines", line = list(width = 1, dash = "dash"),
              showlegend = FALSE
    ) %>%
    layout(
      title = "B.1.1.7 Variant - International Comparison",
      xaxis = list(
        title = "",
        type = "date",
        # Nice! Just in time. %V became possible only a few months ago!
        # https://github.com/plotly/plotly.js/issues/3052#issuecomment-688834265
        tickformat = "Week %V<br>%Y",
        range = c(min_date - 7, max_date + 7)
      ),
      yaxis = list(
        title = "Estimated Percentage of B.1.1.7",
        ticksuffix = "%"
      )
    )

  if (logscale) {
    fig <- fig %>% layout(yaxis = list(type = "log"))
  }
  fig <- fig %>%
    config(
      displaylogo = FALSE,
      modeBarButtons = list(list("zoom2d", "toImage", "resetScale2d", "pan2d")),
      toImageButtonOptions = list(format = "svg", width = 1600, height = 800, filename = "variantPlot_b117_international")
    )

  # Save to file
  suffix <- ""
  if (logscale) {
    suffix <- "_log"
  }
  htmlwidgets::saveWidget(fig, paste0(output_dir, "/variantPlot_b117_international", suffix, ".html"), selfcontained = FALSE,
                          title = "B.1.1.7 Variant - International Comparison",
                          libdir = paste0(output_dir, "/lib_b117_international", suffix))
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please provide an output path.")
}
print(paste0("Here is \"here\": ", here::here()))
generate_switzerland_plot(args[1])
generate_switzerland_plot(args[1], TRUE)
generate_comparison_plot(args[1])
generate_comparison_plot(args[1], TRUE)
