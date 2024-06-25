

#' plot motif count or fraction as bar plot
#' @import  ggplot2 tidyr dplyr
#' @param x
#' @param ylim
#' @param x_title
#' @param plot_type
#' @param integrate_mut
#' @param bar_color
#' @param motif_levels
#' @param ...
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
plotMotif <- function(x,
                      ylim,
                      x_title = "Motif",
                      plot_type = c("Fraction"),
                      integrate_mut = FALSE,
                      bar_color = c("A" = "cornflowerblue",
                                    "C" = "darksalmon",
                                    "G" = "palevioletred2",
                                    "T" = "forestgreen"),
                      motif_levels = c("C", "G", "A", "T"),
                      ...) {

  if (integrate_mut == FALSE) {
    y_use <- NULL
    x <- tibble::as_tibble(x)
    plot_type <-  stringr::str_to_lower(plot_type)
    max_count <- max(dplyr::pull(x, .data$n))
    max_fraction <- max(dplyr::pull(x, .data$fraction))

    if (plot_type == c("fraction")) {

      y_use <- "fraction"

      if (missing(ylim)) {
        ylim <- c(0, max_fraction * 1.1)
        message("ylim", " was set as: ", ylim[1], " - ", ylim[2])
      }

    } else if (plot_type == c("n")) {
      y_use <- "n"

      if(missing(ylim)) {
        ylim <- c(0, max_count * 1.2)

        message("ylim", " was set as: ", ylim[1], " - ", ylim[2])

      }

    }

    x <- x %>%
      dplyr::mutate(group = stringr::str_extract(.data$motif,
                                                 pattern = "^[ATCG]")) %>%
      dplyr::filter(.data$group %in% motif_levels)

    x$group <- factor(x$group, levels = motif_levels)

    x <- x %>%
      dplyr::arrange(.data$group)

    x$motif <-  factor(x$motif, levels = x$motif)

    p <- ggplot(data = x,
                mapping = aes(x = .data$motif,
                              y = .data[[y_use]],
                              fill = .data$group)) +
      geom_col() +
      labs(x = x_title, y = stringr::str_to_title(plot_type)) +
      scale_fill_manual(values = bar_color) +
      scale_y_continuous(limits = ylim) +
      theme_motif_plot()

  } else if (integrate_mut == TRUE) {

    x <- as_tibble(x)
    plot_type <- tolower(plot_type)

    # Create color mapping for REF and MUT using adjustcolor
    ref_colors <- setNames(adjustcolor(bar_color, alpha.f = 0.5),
                           paste0(names(bar_color), "_REF"))
    mut_colors <- setNames(adjustcolor(bar_color, alpha.f = 1.5),
                           paste0(names(bar_color), "_MUT"))
    all_colors <- c(ref_colors, mut_colors)

    # Reshape the data for stacking
    if (plot_type == "fraction") {
      x <- x %>%
        pivot_longer(cols = starts_with("fraction"),
                     names_to = "type", values_to = "value")
    } else {
      x <- x %>%
        pivot_longer(cols = starts_with("n_"),
                     names_to = "type", values_to = "value")
    }

    # Set ylim if not provided
    if (is.null(ylim)) {
      max_val <- max(x$value, na.rm = TRUE)
      ylim <- c(0, max_val * 1.2)
    }

    # Extract the first letter of motif for color grouping and adjust type
    x <- x %>%
      mutate(group = stringr::str_extract(motif, "^[ATCG]"),
             type = ifelse(type == "fraction_ref" | type == "n_ref",
                           paste0(group, "_REF"), paste0(group, "_MUT"))) %>%
      filter(group %in% motif_levels) %>%
      arrange(match(group, motif_levels))  # Ensure correct ordering

    x$group <- factor(x$group, levels = motif_levels)
    x$motif <- factor(x$motif, levels = unique(x$motif))

    # Plotting with adjusted theme
    p <- ggplot(x, aes(x = motif, y = value, fill = type)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(x = x_title, y = stringr::str_to_title(plot_type), fill = "Type") +
      scale_fill_manual(values = all_colors) +
      scale_y_continuous(limits = ylim) +
      theme_classic() %+replace%
        theme(
          axis.text.x = element_text(angle = 45,
                                     size = 5),
          axis.text.y = element_text(size = 5),
          axis.title.x = element_text(size = 5, face = "bold",
                                      margin = margin(t = 0,
                                                      r = 0,
                                                      b = 0,
                                                      l = 0)),
          axis.title.y = element_text(size = 5,
                                      face = "bold"),
          legend.title = element_text(size = 4),  # Adjust legend title size
          legend.text = element_text(size = 4),  # Adjust legend text size
          panel.background = element_rect(fill = "white",
                                          colour = NA),
          legend.key.size = unit(0.15, "cm"),
          plot.background = element_rect(fill = "white",
                                         colour = NA),
          legend.margin = margin(t = 0, r = 0.1, b = 0, l = -0.3, unit = "cm")
          # Reduce the size of the legend keys
        )
  }
  return(p)
}
