###############

## PSP site annotation helper

###############

## Setup
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

BiocManager::install("biomaRt")
BiocManager::install("rWikiPathways")
BiocManager::install("RCy3")
install.packages("dplyr")

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

## Read in PSP data
psp.data <- read.csv("../annotations/kinase/Kinase_Substrate_Dataset.txt", stringsAsFactors = FALSE, sep = "\t")

## Read in local BioMART file
biomart <- read.csv("../annotations/mapping/ensembl_mappings.txt", stringsAsFactors = FALSE, sep = "\t") %>%
  dplyr::select(ensembl_gene_id,uniprotswissprot) %>% 
  filter(uniprotswissprot != "") %>%
  distinct()

psp.data.human <- psp.data %>%
  filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human" & IN_VIVO_RXN == "X") %>%
  dplyr::select(GENE, KINASE, KIN_ACC_ID, SUBSTRATE, SUB_ACC_ID, SUB_GENE_ID, SUB_MOD_RSD)

psp.data.human.full <- psp.data %>%
 filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human")

## Open the relevant WP in Cytoscape using rWikiPathways
wp = "WP4223"
openwp.cmd <- paste0('wikipathways import-as-pathway id="', wp, '"')
RCy3::commandsRun(openwp.cmd)

## Get the full node table for the pathway
node.table <- RCy3::getTableColumns(table = "node")

## Select the ptm node info
node.table.gp <- node.table %>% 
  filter(Type %in% c("Protein", "GeneProduct")) %>%
  dplyr::select(SUID, name, XrefId, Ensembl)

## Map the node table to uniprot
node.table.gp.mapped <- merge(node.table.gp, biomart, by.x ="Ensembl", by.y = "ensembl_gene_id")

curation.nodes <- psp.data.human %>% 
  filter (SUB_ACC_ID %in% node.table.gp.mapped$uniprotswissprot) %>%
  filter (KIN_ACC_ID %in% node.table.gp.mapped$uniprotswissprot)

curation.nodes.suid <- merge(curation.nodes, node.table.gp.mapped, by.x ="SUB_ACC_ID", by.y ="uniprotswissprot") %>%
  dplyr::select(SUID) %>% 
  pull()

setNodeColorBypass(curation.nodes.suid, '#FF0088')

print(curation.nodes)

############
