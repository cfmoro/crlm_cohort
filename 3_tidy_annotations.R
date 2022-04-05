library(tidyverse)

# Globals
combined_fn <- "./output/combined_annotations.csv"

combined_slide_fn <- "./output/gp_annotations_by_slide.csv"
combined_tumor_fn <- "./output/gp_annotations_by_tumor.csv"
combined_probe_fn <- "./output/gp_annotations_by_probe.csv"

regression_slide_fn <- "./output/regression_by_slide.csv"
regression_tumor_fn <- "./output/regression_by_tumor.csv"
regression_probe_fn <- "./output/regression_by_probe.csv"

test_data_fn <- "./annotations/Annotation_tests_CRLM_cohort.csv"
is_test = FALSE #TRUE   # Test consistency of parsed annotations with test dataset (csv > ndpa > parse)

# Read all annotations
df <- read.csv(combined_fn, row.names=NULL)


# Tumor regression by slide
regression_slide <- df %>% select(-lengths_um) %>% filter(annotation_types == "Tumor")

write.csv(regression_slide, regression_slide_fn, row.names=FALSE)

# Tumor regression by tumor
regression_tumor <- regression_slide %>% group_by(ids, tumors) %>% summarise(avg_percent = round(mean(percents), 2))

write.csv(regression_tumor, regression_tumor_fn, row.names=FALSE)

# Tumor regression by probe - summarizing the regression in every slide (not by tumors, to diminish bias by tumor size)
regression_probe <- regression_slide %>% group_by(ids) %>% summarise(avg_percent = round(mean(percents), 2))
# Tumor regression by probe - summarizing the regression by tumors, as with regression by report
#regression_probe <- regression_tumor %>% group_by(ids) %>% summarise(avg_percent = round(mean(avg_percent), 2))

write.csv(regression_probe, regression_probe_fn, row.names=FALSE)


# Inv front annotations
inv_front <- df %>% select(-percents) %>% filter(annotation_types != "Tumor")

# Remove GP annotations from tumors with complete regression
tumors_complete_regression <- filter(regression_tumor, avg_percent == 0) %>% mutate(id_tumor = paste(ids, tumors, sep = "-"))

inv_front <- inv_front %>% mutate(id_tumor = paste(ids, tumors, sep = "-"))
inv_front <- inv_front %>% filter(!(id_tumor %in% tumors_complete_regression$id_tumor))

# Inv front annotations, sum and % by GP and slide
sum_inv_front_slide <-  inv_front %>% group_by(ids, tumors, blocks, annotation_types) %>% summarise(length_um = sum(lengths_um))

percent_inv_front_slide <- sum_inv_front_slide %>% group_by(ids, tumors, blocks) %>% mutate(percent_gp = round(100 / sum(length_um) * length_um, 2))

write.csv(percent_inv_front_slide, combined_slide_fn, row.names=FALSE)

# Inv front annotations, sum and % by GP and tumor
sum_inv_front_tumor <- inv_front %>% group_by(ids, tumors, annotation_types) %>% summarise(length_um = sum(lengths_um))

percent_inv_front_tumor <- sum_inv_front_tumor %>% group_by(ids, tumors) %>% mutate(percent_gp = round(100 / sum(length_um) * length_um, 2))

write.csv(percent_inv_front_tumor, combined_tumor_fn, row.names=FALSE)

# Inv front annotations, sum and % by GP and probe - summarizing the GPs in every slide (not by tumors, to diminish bias by tumor size)
sum_inv_front_probe <- inv_front %>% group_by(ids, annotation_types) %>% summarise(length_um = sum(lengths_um))

percent_inv_front_probe <- sum_inv_front_probe %>% group_by(ids) %>% mutate(percent_gp = round(100 / sum(length_um) * length_um, 2))

write.csv(percent_inv_front_probe, combined_probe_fn, row.names=FALSE)


# Tests
# TODO Add tests for tumor regression
if(is_test) {
  
  gp_levels <- c("D", "R", "R2", "P")
  # Test GP by slide
  parsed_by_slide <- read.csv(combined_slide_fn)
  parsed_by_slide <- parsed_by_slide %>% unite(slide, ids, tumors, blocks, sep = "-") %>% rename(label = annotation_types, by.slide = length_um) %>% select(-percent_gp)
  parsed_by_slide$slide <- as.factor(parsed_by_slide$slide)
  parsed_by_slide <- parsed_by_slide %>% mutate(label = fct_relevel(label, gp_levels))

  test_data <- read.csv(test_data_fn) 
  test_data <- test_data %>% select(slide, label, by.slide) %>% filter(!is.na(by.slide), label != '%')
  test_data <- test_data %>% mutate(label = fct_relevel(label, gp_levels)) %>% mutate(label = fct_drop(label))
  print("Testing GPs by slide")
  print(all_equal(parsed_by_slide, test_data))
  stopifnot(all_equal(parsed_by_slide, test_data))
  print("Test PASSED, parsed dataset is identical till test datase")  
  
  # Test GP by tumor
  parsed_by_tumor <- read.csv(combined_tumor_fn)
  parsed_by_tumor <- parsed_by_tumor %>% unite(tumor, ids, tumors, sep = "-") %>% rename(label = annotation_types, by.tumor = length_um) %>% select(-percent_gp)
  parsed_by_tumor$tumor <- as.factor(parsed_by_tumor$tumor)
  parsed_by_tumor <- parsed_by_tumor %>% mutate(label = fct_relevel(label, gp_levels))
  
  print("")
  print("Testing GPs by tumor")
  test_data <- read.csv(test_data_fn) 
  test_data <- test_data %>% select(slide, label, by.tumor) %>% filter(!is.na(by.tumor), label != '%')
  test_data <- test_data %>% mutate(label = fct_relevel(label, gp_levels)) %>% mutate(label = fct_drop(label))
  test_data <- test_data %>% mutate(tumor = str_extract(slide, "\\d+-[a-p]")) %>% select(-slide)
  test_data <- test_data %>% relocate(tumor)
  test_data$tumor <- as.factor(test_data$tumor)
  print(all_equal(parsed_by_tumor, test_data))
  stopifnot(all_equal(parsed_by_tumor, test_data))
  print("Test PASSED, parsed dataset is identical till test datase")
  
  # Test GP by probe
  parsed_by_probe <- read.csv(combined_probe_fn)
  parsed_by_probe <- parsed_by_probe %>% rename(label = annotation_types, by.probe = length_um) %>% rename(probe = ids) %>% select(-percent_gp)
  parsed_by_probe$probe <- as.factor(parsed_by_probe$probe)
  parsed_by_probe <- parsed_by_probe %>% mutate(label = fct_relevel(label, gp_levels))
  
  print("")
  print("Testing GPs by probe")
  test_data <- read.csv(test_data_fn) 
  test_data <- test_data %>% select(slide, label, by.probe) %>% filter(!is.na(by.probe), label != '%')
  test_data <- test_data %>% mutate(label = fct_relevel(label, gp_levels)) %>% mutate(label = fct_drop(label))
  test_data <- test_data %>% mutate(probe = str_extract(slide, "\\d+"))  %>% select(-slide)
  test_data <- test_data %>% relocate(probe)
  test_data$probe <- as.factor(test_data$probe)
  print(all_equal(parsed_by_probe, test_data))
  stopifnot(all_equal(parsed_by_probe, test_data))
  print("Test PASSED, parsed dataset is identical till test datase")      
  
}
