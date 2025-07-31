###############

## Extract ptm info from CPTAC pathways

###############
## Setup
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(xml2)

## Change working dir to local folder containing gpmls.
setwd("~/Downloads/wikipathways-20250510-gpml-Homo_sapiens/")

## Define output dir.
dir.path <- "~/github/network-ptm-integration/pathways/ptm_info/"

#files <- c("Hs_Melanoma_WP4685_20250303.gpml", "Hs_MAPK_signaling_WP382_20250303.gpml") ##manually defined list
files <- list.files()

## Define combined data frame
wp.cptac.ptm.all <- data.frame() ##for combined table

for (f in files) {
print(f)
wpid <- str_extract(f, "WP\\d+")
gpml <- read_xml(f)
ns <- xml_ns(gpml)
states <- xml_find_all(gpml, ".//d1:State", ns)

if (length(states) > 0){
  state_info <- map_dfr(states, function(node) {
    comment_node <- xml_find_first(node, ".//d1:Comment", ns)
    comment_text <- xml_text(comment_node)

    # Skip this <State> if the comment doesn't contain 'parentid='
    if (is.na(comment_text) || !str_detect(comment_text, "parentid=")) {
      return(NULL)
    }
    
    # Parse key=value pairs from the comment
    kv <- str_split(comment_text, ";\\s*")[[1]]
    kv_split <- str_split_fixed(kv, "=", 2)
    comment_df <- as.data.frame(kv_split, stringsAsFactors = FALSE)
    
    colnames(comment_df) <- c("key", "value")
    parsed <- pivot_wider(comment_df, names_from = key, values_from = value)
    # Add GraphId and TextLabel from the <State> attributes
    parsed$GraphId <- xml_attr(node, "GraphId")
    parsed$TextLabel <- xml_attr(node, "TextLabel")
   
    parsed
  })

if ("parentid" %in% colnames(state_info)){
state_info <- state_info %>%
  filter(parentid != "") %>%
  select(parentid, parentsymbol, position) %>%
  distinct() %>%
  mutate(WPID=wpid)

## Save to combined data frame
wp.cptac.ptm.all <- rbind(wp.cptac.ptm.all, state_info)

## Export individual files
file.name <- paste0(wpid, '-ptm.txt')
write.table(state_info,
            paste0(dir.path, file.name),
            sep = "\t", row.names = FALSE, quote = FALSE)
}
}
}

## Export combined file
write.table(wp.cptac.ptm.all, 
            paste0(dir.path, "wp-cptac-all-ptm.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
