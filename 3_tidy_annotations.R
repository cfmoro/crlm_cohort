library(tidyverse)

# Globals
base_dir <-"/home/bibu/Workspace/crlm_cohort"
combined_fn <- paste(base_dir, "/output/combined_annotations.csv", sep="")

combined_slide_fn <- paste(base_dir, "/output/gp_annotations_by_slide.csv", sep="")
combined_tumor_fn <- paste(base_dir, "/output/gp_annotations_by_tumor.csv", sep="")
combined_probe_fn <- paste(base_dir, "/output/gp_annotations_by_probe.csv", sep="")

regression_slide_fn <- paste(base_dir, "/output/regression_by_slide.csv", sep="")
regression_tumor_fn <- paste(base_dir, "/output/regression_by_tumor.csv", sep="")
regression_probe_fn <- paste(base_dir, "/output/regression_by_probe.csv", sep="")

test_data_fn <- paste(base_dir, "/annotations/Annotation_tests_CRLM_cohort.csv", sep="")
is_test = TRUE #FALSE   # Test consistency of parsed annotations with test dataset (csv > ndpa > parse)

# Read all annotations
df <- read.csv(combined_fn, row.names=NULL)

# Inv front annotations
inv_front <- df %>% select(-percents) %>% filter(annotation_types != "Tumor")

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


# Tumor regression by slide
regression_slide <- df %>% select(-lengths_um) %>% filter(annotation_types == "Tumor")

write.csv(regression_slide, regression_slide_fn, row.names=FALSE)

# Tumor regression by tumor
regression_tumor <- regression_slide %>% group_by(ids, tumors) %>% summarise(avg_percent = round(mean(percents), 2))

write.csv(regression_tumor, regression_tumor_fn, row.names=FALSE)

# Tumor regression by probe - summarizing the regression in every slide (not by tumors, to diminish bias by tumor size)
regression_probe <- regression_slide %>% group_by(ids) %>% summarise(avg_percent = round(mean(percents), 2))

write.csv(regression_probe, regression_probe_fn, row.names=FALSE)

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
