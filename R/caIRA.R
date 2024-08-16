#' Collecting the Cluster Based Tree and Metadata
#'
#' This function identifies monophyletic groups within a phylogenetic tree based on provided criteria, such as bootstrap threshold, date range, and geographic location.
#'
#' @param tree A phylogenetic tree object (likely of class `phylo` from the `ape` package).
#' @param metat A metadata dataframe containing at least the columns `label`, `location`, and `date`.
#' @param bootstrap_threshold Numeric value specifying the minimum bootstrap support required for a clade to be considered.
#' @param date_range Numeric value specifying the range of dates (in days) within which tips must fall to be considered as part of the same group.
#' @param samearea Logical flag indicating whether only groups from the same geographic area should be considered.
#'
#' @return A dataframe containing the monophyletic groups that meet the specified criteria.
#' @export
#' @importFrom ape getMRCA extract.clade
#' @importFrom dplyr %>% filter select mutate rowwise ungroup inner_join arrange desc group_by slice bind_rows
#' @importFrom stringr str_split
#' @importFrom progress progress_bar
#'
#' @examples
#' # Load necessary packages
#' library(ape)
#' library(dplyr)
#' library(stringr)
#' library(progress)
#'
#' # Generate a random tree with 20 tips
#' tree <- rtree(n = 20)
#' # Generate random bootstrap values for each node
#' set.seed(666)
#' bootstrap_values <- sample(50:100, size = tree$Nnode, replace = TRUE)
#' tree$node.label <- bootstrap_values
#' # Generate metadata
#' set.seed(666)
#' areas <- sample(c("Area1", "Area2", "Area3"), size = length(tree$tip.label), replace = TRUE)
#' start_date <- as.Date("2020-01-01")
#' end_date <- as.Date("2020-01-20")
#' random_dates <- sample(seq(start_date, end_date, by = "day"), size = length(tree$tip.label), replace = TRUE)
#' # Create a data frame with the metadata
#' metadata <- data.frame(label = tree$tip.label, location = areas, date = random_dates)
#' # Run the genclus function
#' genclus(tree, metadata, bootstrap_threshold = 90, date_range = 30, samearea = TRUE)
genclus <- function(tree, metat, bootstrap_threshold, date_range, samearea) {

  bootstrap_threshold <- as.numeric(bootstrap_threshold)
  date_range <- as.numeric(date_range)
  samearea <- as.logical(samearea)

  # Check if the tree has bootstrap values
  if (is.null(tree$node.label) || all(is.na(tree$node.label))) {
    stop("The tree does not contain bootstrap values.")
  }

  # Identify the root node (usually it's the first node after the tips)
  root_node <- length(tree$tip.label) + 1

  # Calculate the total number of iterations for progress tracking
  tip_count <- length(tree$tip.label)
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

  # Extract monophyletic groups with bootstrap value filtering first
  groups <- list()

  for (i in 1:(tip_count - 1)) {
    for (j in (i + 1):tip_count) {
      # Update progress
      pb$tick()

      common_ancestor <- getMRCA(tree, c(tree$tip.label[i], tree$tip.label[j]))

      # Determine bootstrap value, set to NA if the node is the root
      if (common_ancestor == root_node) {
        bootstrap_value <- NA
      } else {
        bootstrap_value <- tree$node.label[common_ancestor - length(tree$tip.label)]
      }

      # Proceed only if the bootstrap value is NA (root) or meets the threshold
      if (is.na(bootstrap_value) || bootstrap_value >= bootstrap_threshold) {
        clade <- extract.clade(tree, common_ancestor)

        if (!is.null(clade) && all(clade$tip.label %in% tree$tip.label)) {
          group <- sort(clade$tip.label)
          group_name <- paste(group, collapse = ", ")

          if (!group_name %in% names(groups)) {
            groups[[group_name]] <- list(
              tips = group,
              bootstrap_value = bootstrap_value,
              parent_node = common_ancestor
            )
          }
        }
      }
    }
  }

  # Convert groups to data frame
  monophyletic_df <- do.call(rbind, lapply(names(groups), function(name) {
    data.frame(
      Group = name,
      Tips = paste(groups[[name]]$tips, collapse = ", "),
      Bootstrap = groups[[name]]$bootstrap_value,
      ParentNode = groups[[name]]$parent_node,
      stringsAsFactors = FALSE
    )
  }))

  # Split the tips and check their areas
  monophyletic_df <- monophyletic_df %>%
    rowwise() %>%
    mutate(
      TipList = strsplit(Tips, ", "),
      Areas = list(metat$location[match(unlist(TipList), metat$label)]),
      Dates = list(metat$date[match(unlist(TipList), metat$label)])
    ) %>%
    ungroup()

  # Add the AllSameArea column and extract the area name
  monophyletic_df <- monophyletic_df %>%
    rowwise() %>%
    mutate(
      AllSameArea = length(unique(unlist(Areas))) == 1,
      AreaName = ifelse(AllSameArea, unique(unlist(Areas))[1], paste(unique(unlist(Areas)), collapse = ", ")),
      MinDate = min(unlist(Dates)),
      MaxDate = max(unlist(Dates)),
      diff = as.numeric(difftime(MaxDate, MinDate, units = "days"))
    ) %>%
    ungroup()

  # Convert Bootstrap to numeric (though already filtered above, it's good to ensure this)
  monophyletic_df <- monophyletic_df %>%
    mutate(Bootstrap = as.numeric(Bootstrap))

  # Filter the DataFrame based on the same area if required and date range within given days
  if (samearea) {
    filtered_df <- monophyletic_df %>%
      filter(AllSameArea == TRUE, diff <= date_range) %>%
      select(Group, Tips, Bootstrap, ParentNode, AreaName, MinDate, MaxDate, diff)
  } else {
    filtered_df <- monophyletic_df %>%
      filter(diff <= date_range) %>%
      select(Group, Tips, Bootstrap, ParentNode, AreaName, MinDate, MaxDate, diff)
  }

  # Add the column to count the number of tips
  filtered_df <- filtered_df %>%
    mutate(NumTips = sapply(strsplit(Tips, ", "), length))

  # Handle duplicates within the same function
  # Split the Group column into individual IDs
  split_ids <- str_split(filtered_df$Group, ",\\s*")

  # Flatten the list of IDs into a single vector and check for duplicates
  all_ids <- unlist(split_ids)
  duplicated_ids <- all_ids[duplicated(all_ids)]

  if (length(duplicated_ids) > 0) {
    # Create a data frame where each row is an ID from the groups
    id_df <- data.frame(
      Group = rep(filtered_df$Group, sapply(split_ids, length)),
      ID = all_ids,
      stringsAsFactors = FALSE
    )

    # Filter out rows where ID is duplicated, keeping the Group with the highest NumTips for each ID
    id_df <- id_df %>%
      inner_join(filtered_df, by = "Group") %>%
      filter(ID %in% duplicated_ids) %>%
      arrange(ID, desc(NumTips)) %>%
      group_by(ID) %>%
      slice(1) %>%
      ungroup() %>%
      select(Group)

    # Create a set of groups to keep based on the ID filtering
    groups_to_keep <- unique(id_df$Group)

    # Filter the original data frame to keep only the selected groups
    final_df <- filtered_df %>%
      filter(!sapply(Group, function(group) any(str_split(group, ",\\s*")[[1]] %in% duplicated_ids))) %>%
      bind_rows(filtered_df %>% filter(Group %in% groups_to_keep))
  } else {
    final_df <- filtered_df
  }

  return(final_df)
}

# Declare global variables to avoid R CMD check warnings
utils::globalVariables(c("Tips", "TipList", "Areas", "AllSameArea", "Dates", "MaxDate", "MinDate", "Bootstrap", "Group", "ParentNode", "AreaName", "ID", "NumTips"))

