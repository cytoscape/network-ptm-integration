###############

## Extract ptm info from CPTAC pathways

###############

## Setup
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

BiocManager::install("biomaRt")
BiocManager::install("rWikiPathways")
BiocManager::install("RCy3")
install.packages("dplyr")
BiocManager::install("ndexr")

library(biomaRt)
library(rWikiPathways)
library(RCy3)
library(dplyr)
library(stringr)

## Ensure Cytoscape is running and connected
cytoscapePing()

## Change working dir. Adjust path based on setup.
setwd("~/github/network-ptm-integration/scripts")
dir.path <- "../pathways/ptm_info/"

###############
## Read in list of pathways with manually curated ptms
## This list was manually constructed
wp.cptac <- read.csv("../pathways/wp-cptac-list.txt", stringsAsFactors = FALSE, sep = "\t")

###############

## Define combined data frame
wp.cptac.ptm.all <- data.frame() ##for combined table

## Loop through pathways with manually curated ptms
for (wp in wp.cptac$wpid) {

## Open the relevant WP in Cytoscape using rWikiPathways
openwp.cmd <- paste0('wikipathways import-as-pathway id="', wp, '"')
RCy3::commandsRun(openwp.cmd)

## Get the full node table for the pathway
node.table <- RCy3::getTableColumns(table = "node")

## Select the ptm node info
node.table.ptm <- node.table %>% 
  dplyr::select(parentid, parentsymbol, position) %>%
  filter(parentid != "") %>%
  distinct()  %>%
  mutate(WPID=wp)
  
## Save to combined data frame
cptac.wp.ptm.all <- rbind(cptac.wp.ptm.all, node.table.ptm)
  
## Export individual files
file.name <- paste0(wp, '-ptm.txt')
write.table(node.table.ptm,
            paste0(dir.path, file.name),
            sep = "\t", row.names = FALSE, quote = FALSE)
}

## Export combined file
write.table(cptac.wp.ptm.all, 
            paste0(dir.path, "wp-cptac-all-ptm.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
