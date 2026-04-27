#' Create significance letters for model comparisons
#'
#' Perform ANOVA and Tukey's HSD test within each combination of
#' `genCor`, `gamma`, and `CVMethod`, and create compact letter
#' displays for annotation.
#'
#' @param dat A data frame containing `genCor`, `gamma`, `CVMethod`,
#'   `Model`, and `Accuracy`.
#'
#' @return A data frame containing significance letters and y-positions
#'   for plotting.
MakeLettersDf <- function(dat) {
  dat %>%
    group_by(genCor, gamma, CVMethod) %>%
    group_modify(~{
      fit <- aov(Accuracy ~ Model, data = .x)
      tk <- TukeyHSD(fit, "Model")
      
      let <- multcompView::multcompLetters(tk$Model[, "p adj"])
      
      cldTbl <- data.frame(
        Model = names(let$Letters),
        label = unname(let$Letters),
        stringsAsFactors = FALSE
      )
      
      ymaxTbl <- .x %>%
        group_by(Model) %>%
        summarise(
          y = max(Accuracy, na.rm = TRUE) + 0.02,
          .groups = "drop"
        )
      
      left_join(cldTbl, ymaxTbl, by = "Model")
    }) %>%
    ungroup()
}

#' Create a faceted boxplot of accuracy
#'
#' Generate a boxplot of `Accuracy` by `CVMethod`, grouped by `Model`,
#' faceted by `genCor` and `gamma`, and annotated with significance letters.
#'
#' @param dat A data frame containing `genCor`, `gamma`, `CVMethod`,
#'   `Model`, and `Accuracy`.
#' @param plotTitle A character string for the plot title. Default is `NULL`.
#'
#' @return A ggplot object.
MakePlot <- function(dat, plotTitle = NULL) {
  
  dat$CVMethod <- factor(dat$CVMethod, 
                         levels = c("All", "Once"), labels = c("CV1", "CV2"))
  dat$Model <- factor(dat$Model)
  dat$genCor <- factor(dat$genCor)
  dat$gamma <- factor(dat$gamma)
  
  lettersDf <- MakeLettersDf(dat)
  
  ggplot(dat, aes(x = CVMethod, y = Accuracy, fill = Model)) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.7,
      outlier.shape = 16,
      outlier.size = 1.5
    ) +
    geom_text(
      data = lettersDf,
      aes(x = CVMethod, y = y, label = label, group = Model),
      position = position_dodge(width = 0.8),
      vjust = 0,
      size = 3.5,
      show.legend = FALSE
    ) +
    facet_grid(genCor ~ gamma) +
    scale_y_continuous(
      limits = c(0, NA)
    ) +
    labs(
      title = plotTitle,
      x = "CVMethod",
      y = "Accuracy",
      fill = "Model"
    ) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 15, face = "bold"),
          legend.position = "none")
}
