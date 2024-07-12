# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

huebreaker <- function(tree, metat, bootstrap_threshold, date_range, samearea) {
  bootstrap_threshold <- as.numeric(bootstrap_threshold)
  date_range <- as.numeric(date_range)
  samearea <- as.logical(samearea)

  # Extract monophyletic groups
  groups <- list()
  tip_count <- length(tree$tip.label)
  for (i in 1:(tip_count - 1)) {
    for (j in (i + 1):tip_count) {
      common_ancestor <- getMRCA(tree, c(tree$tip.label[i], tree$tip.label[j]))
      clade <- extract.clade(tree, common_ancestor)
      if (!is.null(clade) && all(clade$tip.label %in% tree$tip.label)) {
        group <- sort(clade$tip.label)
        group_name <- paste(group, collapse = ", ")
        if (!group_name %in% names(groups)) {
          groups[[group_name]] <- list(
            tips = group,
            bootstrap_value = tree$node.label[common_ancestor - length(tree$tip.label)],
            parent_node = common_ancestor
          )
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

  # Convert Bootstrap to numeric
  monophyletic_df <- monophyletic_df %>%
    mutate(Bootstrap = as.numeric(Bootstrap))

  # Filter the DataFrame based on bootstrap values, same area if required, and date range within given days
  if (samearea) {
    filtered_df <- monophyletic_df %>%
      filter(Bootstrap >= bootstrap_threshold, AllSameArea == TRUE, diff <= date_range) %>%
      select(Group, Tips, Bootstrap, ParentNode, AreaName, MinDate, MaxDate, diff)
  } else {
    filtered_df <- monophyletic_df %>%
      filter(Bootstrap >= bootstrap_threshold, diff <= date_range) %>%
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
