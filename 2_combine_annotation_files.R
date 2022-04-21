library(tidyverse)

# Globals
parsed_csv_dir <- "./output/parsed_annotations"
combined_fn <- "./output/combined_annotations.csv"
test_data_fn <- "./annotations_tests/Annotation_tests_CRLM_cohort.csv"
is_test = FALSE # TRUE #  Test consistency of parsed test annotations with test dataset (csv > ndpa > parse)

# Remove previous combined file
if (file.exists(combined_fn)) {
  file.remove(combined_fn)
}

# Combined dataframe
df <- data.frame()
# List of previously parsed cvs files with annotation data, one for each slide
file_list <- list.files(path = parsed_csv_dir, pattern = "*.csv")
# Iterate csv files and combine in one large dataframe
for (fn in file_list){
    temp_df <-read.csv(paste(parsed_csv_dir, "/",fn, sep=""), row.names=NULL, colClasses = "character") # colClasses needed to avoid conversion of blocks "F" to boolean ("FALSE")
    
    df <- rbind(df, temp_df)
}
# Write combined dataframe as csv
write.csv(df, combined_fn, row.names=FALSE)

# Test
if(is_test) {
  
  annot_levels <- c("D", "R", "R2", "P", "%")
  
  print("Testing parsed GP and regression annotation data")
  combined_data <- read.csv(combined_fn)
  combined_data <- combined_data %>% unite(slide, ids, tumors, blocks, sep = "-") %>% rename(label = annotation_types, value = lengths_um) %>% mutate(label = recode(label, 'Tumor' = '%'))
  combined_data <- combined_data %>% mutate(value = ifelse( label == '%', percents, value)) %>% select(-percents)
  combined_data$slide <- as.factor(combined_data$slide)
  combined_data$label <- fct_relevel(combined_data$label, annot_levels)

  test_data <- read.csv(test_data_fn) 
  test_data <- test_data %>% select(slide:value)
  test_data$label <- fct_relevel(test_data$label, annot_levels)
  # Change a value in a dataset to provoke for test failure
  # test_data[1, 3] <- test_data[1, 3] + as.integer(1000)
  print(all_equal(combined_data, test_data))
  stopifnot(all_equal(combined_data, test_data))
  print("Test PASSED, parsed dataset is identical till test dataset")
}
