test_that("beastclus returns expected output", {
  library(treeio)
  library(ggtree)
  library(dplyr)
  library(tidyverse)
  # Example inputs
  metadata<-read.csv('https://raw.githubusercontent.com/Dhihram/caIRA/refs/heads/master/data/metadata_samp.csv')
  beast_tree <- read.beast('https://raw.githubusercontent.com/Dhihram/caIRA/refs/heads/master/data/pox_strict_comb.tree')

  # Run the function
  result <- beastclus(beast_tree, metadata, post_threshold = 0.50, date_range = 90, samearea = TRUE)

  # Expectations
  expect_s3_class(result, "data.frame")
  expect_named(result, c("ParentNode", "label", "Posterior", "AreaName", "min_date", "max_date", "NumTips", "dif"))
})
