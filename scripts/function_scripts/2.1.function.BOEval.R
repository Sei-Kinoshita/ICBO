###### Function to plot metrics from BO result ######
SaveBoPlotPdf <- function(df,
                          pdf_file,
                          GenCor_value = NULL,
                          Gamma_value = NULL,
                          h2_value = NULL,
                          width = 8,
                          height = 11) {
  
  d <- df
  
  if (!is.null(GenCor_value)) {
    d <- d %>% filter(GenCor == GenCor_value)
  }
  if (!is.null(Gamma_value)) {
    d <- d %>% filter(Gamma == Gamma_value)
  }
  if (!is.null(h2_value)) {
    d <- d %>% filter(h2 == h2_value)
  }
  
  ### change to long format
  d_long <- d %>%
    pivot_longer(
      cols = all_of(metricOrder),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = factor(Metric, levels = metricOrder)
    ) %>%
    mutate(
      Value = as.numeric(Value)
    )
  
  ## page keys for boNo
  page_keys <- d_long %>%
    distinct(Scenario, BoNo) %>%
    arrange(Scenario, BoNo)
  
  pdf(pdf_file, width = width, height = height, onefile = TRUE)
  
  for (i in seq_len(nrow(page_keys))) {
    sc_i <- page_keys$Scenario[i]
    bo_i <- page_keys$BoNo[i]
    
    dd <- d_long %>%
      filter(Scenario == sc_i, BoNo == bo_i)
    
    p <- ggplot(dd, aes(x = Cycle, y = Value, color = Index, group = Index)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.1) +
      facet_grid(Metric ~ ., scales = "free_y", switch = "y") +
      theme_bw(base_size = 12) +
      theme(
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey95", color = "grey70"),
        panel.spacing.y = unit(0.5, "lines"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
      ) +
      labs(
        title = paste0(sc_i, " | BoNo=", bo_i),
        x = "Cycle",
        y = NULL,
        color = "Method"
      )
    
    print(p)
  }
  
  dev.off()
}


###### Function to plot metrics from BO result (summary all rep) ######
SaveBoMeanPlotPdf <- function(df,
                              pdf_file,
                              GenCor_value = NULL,
                              Gamma_value = NULL,
                              h2_value = NULL,
                              width = 8,
                              height = 11,
                              minBoForPlot = 20,
                              ciLevel = 0.95,
                              metricYRangeDf = NULL) {
  
  d <- df
  
  if (!is.null(GenCor_value)) {
    d <- d %>% filter(GenCor == GenCor_value)
  }
  if (!is.null(Gamma_value)) {
    d <- d %>% filter(Gamma == Gamma_value)
  }
  if (!is.null(h2_value)) {
    d <- d %>% filter(h2 == h2_value)
  }
  
  ### change to long format
  d_long <- d %>%
    pivot_longer(
      cols = all_of(metricOrder),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = factor(Metric, levels = metricOrder),
      Value = as.numeric(Value)
    )
  
  ### mena & 95%CI for all BoNo in each cycle
  alphaNow <- 1 - ciLevel
  
  d_sum <- d_long %>%
    group_by(Scenario, Index, Metric, Cycle) %>%
    summarise(
      nBo = n_distinct(BoNo),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(nBo),
      tcrit = qt(1 - alphaNow / 2, df = pmax(nBo - 1, 1)),
      Lower = Mean - tcrit * SE,
      Upper = Mean + tcrit * SE,
      .groups = "drop"
    ) %>%
    mutate(
      Lower = ifelse(is.na(Lower), Mean, Lower),
      Upper = ifelse(is.na(Upper), Mean, Upper)
    )
  
  ### only plot cycle with nBo >= minBoForPlot
  d_sum_plot <- d_sum %>%
    filter(nBo >= minBoForPlot)
  
  ### no plot for Mu with metric using postSd
  d_sum_plot <- d_sum_plot %>%
    filter(!(Metric %in% c("postSd_gTrue", "postSd_error") & Index == "Mu"))
  
  ### 1 scenario in 1 pdf page
  page_keys <- d_sum_plot %>%
    distinct(Scenario) %>%
    arrange(Scenario)
  
  pdf(pdf_file, width = width, height = height, onefile = TRUE)
  
  for (i in seq_len(nrow(page_keys))) {
    sc_i <- page_keys$Scenario[i]
    
    dd <- d_sum_plot %>%
      filter(Scenario == sc_i)
    yRangeNow <- metricYRangeDf %>%
      mutate(Cycle = min(dd$Cycle))
    
    p <- ggplot(dd, aes(x = Cycle, y = Mean, color = Index, fill = Index, group = Index)) +
      geom_blank(
        data = yRangeNow,
        aes(x = Cycle, y = ymin),
        inherit.aes = FALSE
      ) +
      geom_blank(
        data = yRangeNow,
        aes(x = Cycle, y = ymax),
        inherit.aes = FALSE
      ) +
      geom_ribbon(
        data = dd %>% filter(Metric != "nTop_Match"),
        aes(ymin = Lower, ymax = Upper),
        alpha = 0.18, linewidth = 0
      ) +
      geom_line(linewidth = 0.9) + 
      geom_point(size = 1.2) + 
      facet_grid(Metric ~ ., scales = "free_y", switch = "y") +
      theme_bw(base_size = 12) +
      theme(
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey95", color = "grey70"),
        panel.spacing.y = unit(0.5, "lines"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)
      ) +
      labs(
        title = paste0(sc_i, "  (only cycles with nBo ≥ ", minBoForPlot, ")"), 
        x = "Cycle",
        y = NULL,
        color = "Method",
        fill = "Method"
      )
    
    print(p)
  }
  
  dev.off()
}

###### Function to plot metrics from BO result (summary all rep, 9 cell) ######
SaveBoMeanPlotPdf_9Panel <- function(df,
                                     pdf_file,
                                     h2_value = NULL,
                                     width = 14,
                                     height = 12,
                                     minBoForPlot = 20,
                                     ciLevel = 0.95,
                                     metricYRangeDf = NULL) {
  
  d <- df
  
  if (!is.null(h2_value)) {
    d <- d %>% filter(h2 == h2_value)
  }
  
  ### long format
  d_long <- d %>%
    pivot_longer(
      cols = all_of(metricOrder),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = factor(Metric, levels = metricOrder),
      Value = as.numeric(Value)
    )
  
  ### summary
  alphaNow <- 1 - ciLevel
  
  d_sum <- d_long %>%
    group_by(GenCor, Gamma, h2, Scenario, Index, Metric, Cycle) %>%
    summarise(
      nBo = n_distinct(BoNo),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(nBo),
      tcrit = qt(1 - alphaNow / 2, df = pmax(nBo - 1, 1)),
      Lower = Mean - tcrit * SE,
      Upper = Mean + tcrit * SE,
      .groups = "drop"
    ) %>%
    mutate(
      Lower = ifelse(is.na(Lower), Mean, Lower),
      Upper = ifelse(is.na(Upper), Mean, Upper)
    )
  
  ### same filtering as before
  d_sum_plot <- d_sum %>%
    filter(nBo >= minBoForPlot) %>%
    filter(!(Metric %in% c("postSd_gTrue", "postSd_error") & Index == "Mu")) %>%
    mutate(
      GammaF = factor(Gamma, levels = sort(unique(Gamma))),
      GenCorF = factor(GenCor, levels = sort(unique(GenCor)))
    )
  
  ### x-axis max actually used in this plot
  maxCycleForPlot <- max(d_sum_plot$Cycle, na.rm = TRUE)
  
  ### y range for each metric
  if (is.null(metricYRangeDf)) {
    metricYRangeDf <- d_sum_plot %>%
      pivot_longer(
        cols = c(Mean, Lower, Upper),
        names_to = "Stat",
        values_to = "AxisValue"
      ) %>%
      group_by(Metric) %>%
      summarise(
        ymin = min(AxisValue, na.rm = TRUE),
        ymax = max(AxisValue, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Metric = factor(Metric, levels = metricOrder)
      )
  }
  
  yRangeNow <- expand.grid(
    Metric = levels(factor(d_sum_plot$Metric, levels = metricOrder)),
    GammaF = levels(d_sum_plot$GammaF),
    GenCorF = levels(d_sum_plot$GenCorF)
  ) %>%
    left_join(metricYRangeDf, by = "Metric") %>%
    mutate(
      Metric = factor(Metric, levels = metricOrder),
      GammaF = factor(GammaF, levels = levels(d_sum_plot$GammaF)),
      GenCorF = factor(GenCorF, levels = levels(d_sum_plot$GenCorF)),
      Cycle = 1
    )
  
  pdf(pdf_file, width = width, height = height, onefile = TRUE)
  
  for (metricNow in metricOrder) {
    
    dd <- d_sum_plot %>%
      filter(Metric == metricNow)
    
    yRangeMetricNow <- yRangeNow %>%
      filter(Metric == metricNow)
    
    p <- ggplot(dd, aes(x = Cycle, y = Mean, color = Index, fill = Index, group = Index)) +
      geom_blank(
        data = yRangeMetricNow,
        aes(x = Cycle, y = ymin),
        inherit.aes = FALSE
      ) +
      geom_blank(
        data = yRangeMetricNow,
        aes(x = Cycle, y = ymax),
        inherit.aes = FALSE
      ) +
      geom_ribbon(
        data = dd %>% filter(Metric != "nTop_Match"),
        aes(ymin = Lower, ymax = Upper),
        alpha = 0.18, linewidth = 0
      ) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.0) +
      facet_grid(GenCorF ~ GammaF) +
      scale_x_continuous(
        limits = c(1, maxCycleForPlot)
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey70"),
        panel.spacing = unit(0.6, "lines"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
      ) +
      labs(
        title = paste0("h2 = ", h2_value, " | Metric = ", metricNow,
                       "  (only cycles with nBo ≥ ", minBoForPlot, ")"),
        x = "Cycle",
        y = NULL,
        color = "Method",
        fill = "Method"
      )
    
    print(p)
  }
  
  dev.off()
}
