#' Plot Trinucleotide Data
#'
#' This function processes and plots trinucleotide data. It first applies specified
#' filters and transformations to the data and then generates a visual representation
#' of the results. The function handles data normalization, exclusion, and retention
#' based on provided column names, and it creates detailed plots with options for
#' customization of plot aesthetics.
#'
#' @importFrom ggplot2 ggplot geom_bar aes_string scale_fill_manual theme
#' @import ggpattern
#' @import patchwork
#'
#' @param trinuc_df DataFrame containing trinucleotide data.
#' @param exclude_if_type_present Vector of strings specifying column names whose non-zero presence should trigger exclusion of rows.
#' @param retain_if_type_present Vector of strings specifying column names whose non-zero presence is necessary to retain rows.
#' @param remove_type Vector of strings specifying column names where all values should be set to zero.
#' @param overlap_type Logical; whether to stratify SBS trinucleotides by read pair overlap type.
#' @param normalize_counts Logical; if TRUE, normalizes the 'value' column to sum to 1.
#' @param show_overlap_type Logical; if TRUE, show read-pair overlap types.
#' @param min_y_val Numeric; minimum y-axis value for the plots.
#' @param max_y_val Numeric; maximum y-axis value for the plots.
#' @param plot_title String; title for the plots.
#' @param y_axis_title String; title for the y-axis.
#' @param draw_x_axis_labels Logical; whether to draw x-axis labels.
#' @param draw_y_axis_labels Logical; whether to draw y-axis labels.
#' @param draw_y_axis_title Logical; whether to display a title for the y-axis.
#' @param save_plot Logical; whether to save the object as a pdf file.
#' @param output_file String; name and path of output pdf file.
#'
#' @return A trinucleotide SBS plot object and an optional pdf file.
#' @examples
#' \dontrun{
#'  plotTrinucleotide(trinuc_df)
#' }
#'
plotTrinucleotide <- function(trinuc_df,
                              exclude_if_type_present = NULL,
                              retain_if_type_present = NULL,
                              remove_type = NULL,
                              normalize_counts = TRUE,
                              show_overlap_type = TRUE,
                              min_y_val = 0.0,
                              max_y_val = 0.5,
                              plot_title = "Trinucleotide Profile",
                              y_axis_title = "Percentage of Single Base Substitutions",
                              draw_x_axis_labels = TRUE,
                              draw_y_axis_labels = TRUE,
                              draw_y_axis_title = TRUE,
                              save_plot = TRUE,
                              output_file = "./trinucleotide_profile.pdf") {

  # Step 1: Process trinucleotide data
  count_df_processed <- processTrinucleotideData(
    trinuc_df = trinuc_df,
    exclude_if_type_present = exclude_if_type_present,
    retain_if_type_present = retain_if_type_present,
    remove_type = remove_type,
    normalize_counts = normalize_counts
  )

  # Step 2: Plot the processed data
  plotTrinucData(
    count_df = count_df_processed,
    min_y_val = min_y_val,
    max_y_val = max_y_val,
    show_overlap_type = show_overlap_type,
    plot_title = plot_title,
    y_axis_title = y_axis_title,
    draw_x_axis_labels = draw_x_axis_labels,
    draw_y_axis_labels = draw_y_axis_labels,
    draw_y_axis_title = draw_y_axis_title,
    save_plot = TRUE,
    output_file = output_file
  )
}

# Internal functions
processTrinucleotideData <- function(trinuc_df,
                                     exclude_if_type_present = NULL,
                                     retain_if_type_present = NULL,
                                     remove_type = NULL,
                                     normalize_counts = TRUE) {

  # Exclude rows based on exclude_if_type_present
  if (!is.null(exclude_if_type_present) &&
        length(exclude_if_type_present) > 0) {
          trinuc_df <- trinuc_df[!apply(trinuc_df[exclude_if_type_present] > 0,
                                  1, any), ]
  }

  # Retain rows based on retain_if_type_present
  if (!is.null(retain_if_type_present) && length(retain_if_type_present) > 0) {
    trinuc_df <- trinuc_df[apply(trinuc_df[retain_if_type_present] > 0,
                                 1,
                                 any), ]
  }

  # Convert all counts of specified types to 0 in remove_type columns
  if (!is.null(remove_type) && length(remove_type) > 0) {
    trinuc_df[remove_type] <- lapply(trinuc_df[remove_type],
                                     function(x) ifelse(x > 0, 0, x))
  }

  # Derive consensus_mismatch_type based on specified logic
  trinuc_df$consensus_mismatch_type <- ifelse(
    grepl("MUT", trinuc_df$consensus_mismatch),
    ifelse(trinuc_df$CO_MUT >= trinuc_df$SO_MUT, "CO_MUT", "SO_MUT"),
    ifelse(
      grepl("discordant",
            trinuc_df$consensus_mismatch), "DO",
      ifelse(
        grepl("other_base_concordant",
              trinuc_df$consensus_mismatch), "CO_OTHER",
        ifelse(
          grepl("other_base_single_read",
                trinuc_df$consensus_mismatch), "SO_OTHER",
          NA
        )
      )
    )
  )

  # Create a 96 channel trinucleotide matrix from the 6 SBS types
  mutation.types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  bases <- c("A", "T", "C", "G")
  combinations <- expand.grid(start = bases, mut = mutation.types, end = bases)
  sub.types.96 <- with(combinations, paste0(start, "[", mut, "]", end))

  # Count occurrences of each subtype in trinuc_df$SBS96
  subtype_counts <- table(factor(trinuc_df$SBS96, levels = sub.types.96))
  count_df <- data.frame(sample = as.integer(subtype_counts))
  rownames(count_df) <- sub.types.96

  # Count the consensus_mismatch_type occurrences for each SBS96
  for (subtype in sub.types.96) {
    subtype_rows <- trinuc_df[trinuc_df$SBS96 == subtype, ]
    count_df[subtype, "CO_MUT"] <- sum(
      subtype_rows$consensus_mismatch_type == "CO_MUT")
    count_df[subtype, "SO_MUT"] <- sum(
      subtype_rows$consensus_mismatch_type == "SO_MUT")
    count_df[subtype, "DO"] <- sum(
      subtype_rows$consensus_mismatch_type == "DO")
    count_df[subtype, "SO_OTHER"] <- sum(
      subtype_rows$consensus_mismatch_type == "SO_OTHER")
    count_df[subtype, "CO_OTHER"] <- sum(
      subtype_rows$consensus_mismatch_type == "CO_OTHER")
  }

  #Add the counts of SO_OTHER_sample to SO_MUT_sample
  count_df$SO_MUT <- count_df$SO_MUT + count_df$SO_OTHER

  # Add the counts of CO_OTHER_sample to CO_MUT_sample
  count_df$CO_MUT <- count_df$CO_MUT + count_df$CO_OTHER

  # Remove the SO_OTHER_sample and CO_OTHER_sample columns
  count_df <- count_df[, !colnames(count_df) %in% c("SO_OTHER", "CO_OTHER")]

  # Normalize counts if normalize_counts is TRUE
  if (normalize_counts) {
    total_sample_sum <- sum(count_df$sample)
    count_df <- count_df %>%
      mutate(across(c(CO_MUT, SO_MUT, DO), ~ .x / total_sample_sum))
  }

  # Remove the sample column
  count_df <- count_df[, !colnames(count_df) %in% c("sample")]

  # Add the SBS names as a column and convert to long format
  count_df$SBS <- rownames(count_df)
  count_df <- pivot_longer(count_df, cols = c(CO_MUT, SO_MUT, DO),
                           names_to = "overlap_type",
                           values_to = "value")

  # Mutation types in the desired order
  mutation.types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

  # Function to extract the mutation type from the SBS column
  extract_mutation_type <- function(s) {
    sub(".*\\[([A-Z]>[A-Z])\\].*", "\\1", s)
  }

  # Add a column for the extracted mutation type
  count_df <- count_df %>%
    mutate(mutation_type = extract_mutation_type(SBS))

  # Convert the mutation_type column to a factor with the specified levels
  count_df <- count_df %>%
    mutate(mutation_type = factor(mutation_type, levels = mutation.types))

  # Arrange the dataframe based on the mutation_type factor
  count_df <- count_df %>%
    arrange(mutation_type) %>%
    select(-mutation_type)  # Remove the helper column

  count_df <- as.data.frame(count_df)

  return(count_df)
}



plotTrinucData <- function(
    count_df,
    min_y_val = 0.0,
    max_y_val = 0.5,
    show_overlap_type = TRUE,
    plot_title = "Trinucleotide Profile",
    y_axis_title = "Proportion of Single Base Substitutions",
    draw_x_axis_labels = TRUE,
    draw_y_axis_labels = TRUE,
    draw_y_axis_title = TRUE,
    save_plot = TRUE,
    output_file = output_file) {

  # Define the mutation types
  mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

  # Define Plot colours
  plot.colors <- c("#232f7c", "#000000",
                   "#612370", "#474242",
                   "#436e2d", "#ad5c74")

  strip.text.x.colors <- c("white", "white", "white", "white", "white", "white")

  # Initialize a list to hold data frames for each mutation type
  df.plot.data <- list()

  # Populate the list with segmented data frames for each mutation type
  for (i in 1:6) {
    start.idx <- ((i - 1) * 48) + 1
    end.idx <- i * 48

    # Create a temporary data frame with sliced data
    df.temp <- data.frame(
      mutation_subtype = count_df$SBS[start.idx:end.idx],
      value = count_df$value[start.idx:end.idx],
      overlap_type = count_df$overlap_type[start.idx:end.idx],
      title = rep(mutation_types[i], 48),
      stringsAsFactors = FALSE
    )

    # Add the temporary data frame to the list
    df.plot.data[[i]] <- df.temp
  }

  plots <- list()

  # Define axis text settings based on input parameters
  axis.text.x <- if (draw_x_axis_labels) {
    element_text(size = 4,
                 family = "Helvetica",
                 colour = "#888888",
                 angle = 90,
                 vjust = 0.5,
                 hjust = 0.5)
  } else {
    element_blank()
  }

  axis.text.y <- if (draw_y_axis_labels) {
    element_text(size = 5,
                 family = "Helvetica",
                 colour = "#888888")
  } else {
    element_blank()
  }

  # Apply factor levels for 'overlap_type'
  overlap_levels <- c("DO", "SO_MUT", "CO_MUT")

  # Loop through each data frame in the list and create plots
  for (i in seq_along(df.plot.data)) {
    df.plot.data[[i]]$overlap_type <- factor(df.plot.data[[i]]$overlap_type,
                                             levels = overlap_levels)
  }

  # Initialise the plot
  plots <- list()

  # Compile the plot
  for (i in 1:6) {
    if (show_overlap_type == TRUE) {
      p <- ggplot(df.plot.data[[i]], aes_string(x = 'mutation_subtype',
                                                y = 'value',
                                                pattern = 'overlap_type')) +

        geom_bar_pattern(stat = "identity", width = 0.9,
                        fill = plot.colors[i], pattern_density = 0.005,
                        pattern_spacing = 0.05, pattern_colour = "#87cfe3",
                        pattern_fill = "#87cfe3") +

        scale_pattern_manual(
          values = c("CO_MUT" = "none",
                    "SO_MUT" = "stripe",
                    "DO" = "crosshatch"),
          labels = c("CO_MUT" = "Concordant",
                    "SO_MUT" = "Single-Read",
                    "DO" = "Discordant")) +

        ylab(NULL) + xlab(NULL) +
        scale_y_continuous(limits = c(min_y_val,
                                      max_y_val),
                          expand = c(0, 0)) +
        theme(axis.text.x = axis.text.x,
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major.y = element_line(colour = "#DCDCDC",
                                                linetype = "solid"),
              panel.grid.major.x = element_blank(),
              panel.background = element_blank())

        # Conditional logic to hide legend except for the 4th plot
        if (i != 4) {
          p <- p + theme(legend.position = "none")
        } else {
          # Customize the legend appearance for the 4th plot
          p <- p +
            guides(pattern = guide_legend(
              title = "Overlap Type",
              title.position = "top", label.position = "right",
              override.aes = list(fill = "#f3f6f4"))) +
            theme(
              legend.title = element_text(size = 6, face = "bold"),
              legend.text = element_text(size = 5),
              legend.key.size = unit(0.63, "cm"),
              legend.spacing = unit(0.3, "cm"),
              legend.margin = margin(0.1, 0.1, 0.1, 0.1),
              legend.box.margin = margin(0.1, 0.1, 0.1, 0.1)
            )
        }
    } else if (show_overlap_type == FALSE) {
      # Compile the plot
      p <- ggplot(df.plot.data[[i]],
                    aes_string(x = 'mutation_subtype', y = 'value')) +
        geom_col(fill = plot.colors[i], width = 0.9)  # Simple colored bars

      # Conditional logic to hide legend for all plots
      p <- p + theme(legend.position = "none")

      p <- p + ylab(NULL) + xlab(NULL) +
      scale_y_continuous(limits = c(min_y_val, max_y_val), expand = c(0, 0)) +
      theme(axis.text.x = axis.text.x,
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_line(colour = "#DCDCDC",
                                              linetype = "solid"),
            panel.grid.major.x = element_blank(),
            panel.background = element_blank())
    }

    # Draw top strip
    p <- p +
      theme(strip.background = element_rect(fill = plot.colors[i],
                                            colour = plot.colors[i]),
            strip.text.x = element_text(size = 7,
                                        family = "Helvetica",
                                        face = "bold",
                                        colour = strip.text.x.colors[i])) +
      facet_grid(~title)

    # Margins for the left-most mutation type sub-plot
    if (i == 1) {
      p <- p + theme(axis.text.y = axis.text.y,
                     plot.margin = unit(c(0.3,
                                          0,
                                          0.3,
                                          0.3), "cm"))

      # Margins for the right-most mutation type sub-plot
    } else if(i == 6) {
      p <- p + theme(axis.text.y = element_text(
                                                size = 5,
                                                family = "Helvetica",
                                                colour = "#FFFFFF00"),
      plot.margin = unit(c(0.3,
                           0.3,
                           0.3,
                           0), "cm"))

      # Margins for the middle mutation types
    } else {
      p <- p + theme(axis.text.y = element_text(size = 5,
                                                family = "Helvetica",
                                                colour = "#FFFFFF00"),
                     plot.margin = unit(c(0.3,
                                          0,
                                          0.3,
                                          0), "cm"))
    }
    plots[[i]] <- p
  }

  # Patchwork combine the plots
  plot_combined <- plots[[1]] + plots[[2]] +
    plots[[3]] + plots[[4]] +
    plots[[5]] + plots[[6]] +
    plot_layout(guides = 'collect', ncol = 6) &
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))

  # Use the tag label as a y-axis label
  p_output <- wrap_elements(plot_combined) +
    labs(tag = y_axis_title) +
    theme(
      plot.tag = element_text(size = rel(0.5), angle = 90),
      plot.tag.position = "left"
    )

  # Save the plot if required
  if (save_plot) {
    ggsave(output_file,
           plot = p_output, width = 17,
           height = 6, units = "cm",
           device = "pdf")
  }

  #return(p_output)

}
