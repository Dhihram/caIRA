#' Extract Clusters from a BEAST Phylogenetic Tree
#'
#' @param beast_tree A BEAST tree object of class `phylo4d` or similar.
#' @param metadata A metadata dataframe containing columns `label`, `location`, and `date`.
#' @param post_threshold Numeric value specifying the posterior probability threshold for including a cluster (default = 0.70).
#' @param date_range Numeric value specifying the maximum allowable date range (in days) within clusters (default = 90).
#' @param samearea Logical value indicating whether clusters must originate from the same geographic area (default = FALSE).
#'
#' @return A dataframe summarizing the clusters, including labels, posterior probabilities, areas, date ranges, and the number of tips.
#' @importFrom dplyr group_by summarize ungroup filter select arrange rename mutate inner_join
#' @importFrom tidyr separate_rows
#' @importFrom ggtree ggtree theme_tree2 geom_range
#' @importFrom progress progress_bar
#' @export
#' @examples
#' # Example usage:
#' # Assuming `beast_tree` is a BEAST tree object and `metadata` is a dataframe
#' data_csv <- system.file("extdata", "metadata_samp.csv", package = "caIRA")
#' metadata<-read.csv(data_csv)
#' data_beast <- system.file("extdata", "pox_strict_comb.tree", package = "caIRA")
#' beast_tree <- treeio::read.beast(data_beast)
#' # with the required columns:
#' beastclus(beast_tree, metadata, post_threshold = 0.50, date_range = 90, samearea = TRUE)
beastclus <- function(beast_tree, metadata, post_threshold = 0.70, date_range = 90, samearea = FALSE) {
  library(treeio)
  library(ape)
  library(ggtree)
  library(dplyr)
  library(tidyr)

  # Calculate the total number of iterations for progress tracking
  tip_count <- length(beast_tree@phylo$tip.label)
  total_iterations <- (tip_count * (tip_count - 1)) / 2

  # Set up the progress bar
  pb <- progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = total_iterations,
    complete = "=",   # Completion bar character
    incomplete = "-", # Incomplete bar character
    current = ">",    # Current bar character
    clear = FALSE,    # If TRUE, clears the bar when finish
    width = 100       # Width of the progress bar
  )
  # Number of tips in the tree (leaf nodes)
  tree <- beast_tree@phylo
  num_tips <- length(tree$tip.label)

  # Extract internal node numbers using the edge matrix (parent-child relationships)
  internal_nodes <- unique(tree$edge[, 1][tree$edge[, 1] > num_tips])

  # Create an empty list to store clades and their leaf nodes
  clade_leaf_nodes <- list()

  # Loop through each internal node and extract the clade
  for (i in internal_nodes) {
    # Extract the clade rooted at the current internal node
    clade <- extract.clade(tree, node = i)

    # Get the leaf (tip) names for the clade
    leaf_nodes <- clade$tip.label

    # Store the clade and its leaf nodes in the list
    clade_leaf_nodes[[as.character(i)]] <- leaf_nodes  # Assign node number as the name
  }

  # Convert the list to a data frame for easier visualization
  clade_df <- data.frame(
    nodes = rep(names(clade_leaf_nodes), sapply(clade_leaf_nodes, length)),
    label = unlist(clade_leaf_nodes)
  )

  # Merging posterior probability from tree data
  tre <- ggtree(beast_tree, mrsd = min(metadata$date)) +
    theme_tree2() + geom_range(range = 'length_0.95_HPD', color = 'red', alpha = .6, size = 2)
  post_dat <- tre$data
  post_dat <- select(post_dat, node, posterior, x)

  clade_df <- merge(clade_df, post_dat, by.x = "nodes", by.y = "node", all.x = TRUE)

  # Merging metadata
  metadata_2 <- select(metadata, label, location, date)
  clade_df <- merge(clade_df, metadata_2, by.x = "label", by.y = "label", all.x = TRUE)

  # Summarize data
  summary_df <- clade_df %>%
    group_by(nodes) %>%
    summarize(
      label = paste(label, collapse = ", "),                     # Merge all labels for the group
      Posterior = round(mean(posterior, na.rm = TRUE), 2),       # Use mean posterior (or adjust as needed)
      AreaName = ifelse(length(unique(location)) > 1,            # Check the number of unique locations
                        paste(unique(location), collapse = ", "),
                        unique(location)),                       # Use single location if only one
      min_date = min(as.Date(date)),                             # Minimum date in the group
      max_date = max(as.Date(date)),                             # Maximum date in the group
      NumTips = n(),                                             # Count of label entries
      .groups = "drop"                                           # Avoid grouping in the output
    )

  summary_df$dif <- as.numeric(summary_df$max_date - summary_df$min_date)
  summary_df = rename(summary_df, c('ParentNode' = 'nodes'))

  # Filtering the data based on the posterior and date range criteria
  filtered_df <- summary_df %>%
    filter(Posterior > post_threshold,
           dif <= date_range,
           (samearea == TRUE & !grepl(",", AreaName)) | samearea == FALSE)

  # Separate the label into different rows based on ',' separator
  filtered_df <- filtered_df %>%
    separate_rows(label, sep = ", ")  # Use ", " as separator to split the label

  # Remove duplicate labels, keeping the one with the maximum NumTips
  filtered_df <- filtered_df %>%
    group_by(label) %>%
    filter(NumTips == max(NumTips)) %>%
    ungroup()

  # Merge the `label` values with a comma separator
  final_df <- filtered_df %>%
    group_by(ParentNode) %>%
    summarize(
      label = paste(label, collapse = ", "),   # Merge labels with a comma separator
      Posterior = first(Posterior),            # Keep the first value of `Posterior`
      AreaName = first(AreaName),              # Keep the first value of `AreaName`
      min_date = min(min_date),                # Get the minimum date
      max_date = max(max_date),                # Get the maximum date
      NumTips = first(NumTips),                  # Sum NumTips (if needed)
      dif = max(dif),                          # Get the maximum difference
      .groups = "drop"                         # Remove the grouping after summarizing
    )

  # Convert ParentNode to numeric
  final_df$ParentNode <- as.numeric(final_df$ParentNode)

  #order date
  final_df <- final_df[order(final_df$min_date),]

  return(final_df)
}

utils::globalVariables(c("node", "posterior", "label", "location", "n"))
