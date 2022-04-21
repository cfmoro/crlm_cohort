library(xml2)
library(rgeos)
library(stringr)

# Set globals
is_test = FALSE # TRUE #

if(is_test) {
  ndpa_dir <- "./annotations_tests/annotations/"  
} else {
  ndpa_dir <- "./input/annotations_220421/"  
}

out_dir <- "./output"

# Parse annotation XML files and build an initial df of form: ids-tumors-blocks-annotation_type-annotation_value
# where annotation_type = R(eplacement):length | D(esmoplastic):length | T(umor residual):%
parse_annotation_file <- function (ndpa_file) {
  
  # Split filename to obtain the different fields: id-tumor-block
  fn_tokens <- unlist(str_split(ndpa_file, "-"))
  id <- fn_tokens[1]
  tumor <- fn_tokens[2]
  block <- str_replace(fn_tokens[3], ".ndpi.ndpa", "")

  print(paste("parse_annotation_file: Parsing annotation file:", ndpa_file))
  print(paste("id:", id, "; tumor:", tumor, "block:", block))
  
  
  # Parse annotations in xml format
  annotations_df <- parse_geometry_data(paste(ndpa_dir, "/", ndpa_file, sep = ""))
  
  # Build df parsed annotations
  ids <- c(rep(id, nrow(annotations_df)))
  tumors <- c(rep(tumor, nrow(annotations_df)))
  blocks <- c(rep(block, nrow(annotations_df)))
  
  annotation_slide_data <- cbind(ids, tumors, blocks, annotations_df)
  
  print(paste("End of parsing annotation file: ", ndpa_file))
  print("")
  
  return(annotation_slide_data)
}

# Parse invasion front annotation and return its type and length
get_annotation_length <- function(point_list) {
    point_list_val <- as_list(point_list)
    point_list_val <- unlist(point_list_val)
    x_points <- point_list_val[seq(1, length(point_list_val), by=2)]
    y_points <- point_list_val[seq(2, length(point_list_val), by=2)]
    
    # Create a well-known text representation of geometry (https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) 
    # from point coordinates of annotation
    my_wkt <- "LINESTRING("
    for(cont in 1:length(x_points)) {   
      if(cont == 1) {
        my_wkt <- paste(my_wkt, x_points[[cont]], " ", y_points[[cont]], sep="")
      } else {
        my_wkt <- paste(my_wkt, ",", x_points[[cont]], " ", y_points[[cont]], sep="")
      }
    }
    my_wkt <- paste(my_wkt,")", sep="")
    print("Polygon points read")
    
    my_polygon <- readWKT(my_wkt)
    print("Polygon created")
    
    # Obtain and return length in um from polygon
    return(round(gLength(my_polygon) / 10e+2)) # um
}

# Parse annotation XML files and build an initial df of form: ids-tumors-blocks-annotation_type-annotation_value
# where annotation_type = R(eplacement):length | D(esmoplastic):length | T(umor residual):%
# Invasion front in xml as annotation type="freehand" and title for the type: R(eplacemnt type 1), R2 (replacement type 2), D(esmoplastic)
# T(umor residual) in xml as annotation type="pointer" and % value in title 
parse_geometry_data <- function(ndpa_file) {
  
  print(paste("parse_geometry_data: Parsing annotation file:", ndpa_file))

  # Values to obtain for each annotation the xml
  annotation_types <- c() 
  lengths_um <- c()
  percents <- c()
  
  # Read xml file and annotations
  ndpa_xml <- read_xml(ndpa_file)
  annotations <- xml_children(ndpa_xml)
  
  # iterate annotations and obtain values
  has_tumor_percent_annot <- FALSE
  for(annotation in annotations) {
    
    # Initialize as empty annotation values to be obtained
    id_val <- NA
    title_val <- ""
    annotation_type_val <- ""

    # Parse atomic annotation values
    id_val <- as.integer(xml_attr(annotation, "id"))
    title_val <- xml_text(xml_find_first(annotation, ".//title"))
    annotation_type_val <- xml_attr(xml_find_first(annotation, ".//annotation"), "displayname")
    
    # Sanity check
    if (!(annotation_type_val %in% c("AnnotateFreehandLine", "AnnotatePointer"))) {
      warning(paste("Wrong annotation type: In file: ", ndpa_file, "; id: ", id_val, "; annotation_type_val:", annotation_type_val))
      next
    }
    
    if(annotation_type_val == "AnnotateFreehandLine") {
      # Sanity check
      if (!(title_val %in% c("R", "r", "R1", "r1", "R2", "r2", "D", "d", "P", "p"))) {
        warning(paste("Wrong inv front annotation label: In file: ", ndpa_file, "; id: ", id_val, "; title_val", title_val))
        next
      }
      # Obtain point coordinates for inv front annotation
      ann_length <- get_annotation_length(xml_find_first(annotation, ".//pointlist")) 
      print(paste("Inv front annotation length obtained:", title_val, ":", ann_length, " um"))
      
      # Complete setting of final values
      ann_type <- str_to_upper(title_val)
      if(ann_type == 'R1') {
        ann_type <- 'R'
      }
      ann_percent <- NA
    } else if(annotation_type_val == "AnnotatePointer") {
      # Sanity check
      if(!str_detect(title_val, "%")) {
        warning(paste("Wrong tumour % label, missing '%' sign: In file: ", ndpa_file, "; id: ", id_val, "; title_val", title_val))
        next
      }
      # Remove % sign from label
      ann_percent <- str_replace(title_val, "%", "")
      # Sanity check
      if(!grepl("^[0-9]+$", ann_percent)) {
        warning(paste("Wrong tumour % label, wrong number: In file: ", ndpa_file, "; id: ", id_val, "; title_val", title_val))
        next
      }
      has_tumor_percent_annot <- TRUE
      print(paste("**** Tumour content (%) value obtained:", title_val))
      
      # Complete setting of final values
      ann_type <- "Tumor"
      ann_length <- NA 
    }
    
    # Add values to df  
    annotation_types <- c(annotation_types, ann_type)
    lengths_um <- c(lengths_um, ann_length)
    percents <- c(percents, ann_percent)
  }
  
  if(!has_tumor_percent_annot) {
    warning(paste("Missing tumour % annotation: In file: ", ndpa_file))
  }
  
  annotation_df <- cbind(annotation_types, lengths_um, percents)
  
  return(annotation_df)
}

# Main execution block

# Clean output directory of previously generated files
out_files <- list.files(path = paste(out_dir, "/parsed_annotations", sep=""), pattern = "*.csv")
for(out_file in out_files) {
  unlink(paste(out_dir, "/parsed_annotations/", out_file, sep=""))
}

# Set up logging - Obs for this to work execute the script by selecting all + CrlEnter, and NOT by CtrAlt-R
log_file_name <- paste(out_dir, "/parsed_annotations/errors.txt", sep="")
log_file <- file(log_file_name, open="wt")
sink(log_file, type="message")
  
# Iterate and parse ndpa files
ndpa_files <- list.files(path = ndpa_dir, pattern = "*.ndpa")
for(ndpa_file in ndpa_files) {
  
  print(paste("File: ", ndpa_file))
  annotation_slide_data <- parse_annotation_file(ndpa_file)
  
  output_fn <- paste(out_dir, "/parsed_annotations/", str_replace(ndpa_file, ".ndpi.ndpa", ""), ".csv", sep="")
  write.csv(annotation_slide_data, output_fn, row.names = FALSE)
  
} 

## reset message sink and close the file connection
sink(type="message")
close(log_file)

## Display the log file
readLines(log_file_name)