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
  dplyr::select(GENE, KINASE, KIN_ACC_ID, SUBSTRATE, SUB_ACC_ID, SUB_GENE_ID, SUB_GENE, SUB_MOD_RSD, SITE_...7_AA) %>%
  mutate(SUB_MOD_RSD_2 = case_when(
    str_starts(SUB_MOD_RSD, "S") ~ str_replace(SUB_MOD_RSD, "^S", "ser"),
    str_starts(SUB_MOD_RSD, "T") ~ str_replace(SUB_MOD_RSD, "^T", "thr"),
    str_starts(SUB_MOD_RSD, "Y") ~ str_replace(SUB_MOD_RSD, "^Y", "tyr"),
    TRUE ~ SUB_MOD_RSD  # leave unchanged if it doesn't match S/T/Y
  ))

psp.data.human.full <- psp.data %>%
 filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human")

## Open the relevant WP in Cytoscape using rWikiPathways
wp = "WP179"
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
  filter (KIN_ACC_ID %in% node.table.gp.mapped$uniprotswissprot) %>%
  mutate (comment = paste0("parentid=",SUB_ACC_ID,"; parentsymbol=",SUBSTRATE,"; site=",SITE_...7_AA,"; position=",
                           SUB_MOD_RSD_2,"; ptm=p; direction="))

curation.nodes.suid <- merge(curation.nodes, node.table.gp.mapped, by.x ="SUB_ACC_ID", by.y ="uniprotswissprot") %>%
  dplyr::select(SUID) %>% 
  pull()

setNodeColorBypass(curation.nodes.suid, '#FF7799')

## Manually enter a substrate to print only those rows
substrate = "RBL2"
substrate_subset <- curation.nodes %>%
  filter(SUB_GENE == substrate) %>%
  select(GENE, SUB_GENE, SUB_MOD_RSD_2, comment)
print(substrate_subset)

############
