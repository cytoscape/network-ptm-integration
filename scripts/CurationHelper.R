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

## Read in PSP data
psp.data <- read.csv("../annotations/kinase/Kinase_Substrate_Dataset.txt", stringsAsFactors = FALSE, sep = "\t")

## Read in local BioMART file
biomart <- read.csv("../annotations/mapping/ensembl_mappings.txt", stringsAsFactors = FALSE, sep = "\t") %>%
  dplyr::select(ensembl_gene_id,uniprotswissprot) %>% 
  filter(uniprotswissprot != "") %>%
  distinct()

psp.data.human <- psp.data %>%
  filter(KIN_ORGANISM == "human" & IN_VIVO_RXN == "X") %>% 
  dplyr::select(GENE, KINASE, KIN_ACC_ID, SUB_ACC_ID, SUB_MOD_RSD)

## Data mapping PSP data
psp.data.human.mapped <- merge(psp.data.human, biomart, by.x = "KIN_ACC_ID", by.y = "uniprotswissprot") %>%
  mutate(kinase_ensembl_id = ensembl_gene_id) %>%
  mutate(prot_site = paste0(ensembl_peptide_id, "_", SUB_MOD_RSD)) %>%
  dplyr::select(GENE, KINASE, KIN_ACC_ID, kinase_ensembl_id, SUB_ACC_ID, ensembl_peptide_id, prot_site, SUB_MOD_RSD)

## Open the relevant WP in Cytoscape using rWikiPathways
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

## Get list of proteins on the pathway with ptms according to PSP data, and filter for those with kinases in the pathway
curation.nodes <- psp.data.human %>% 
  filter (SUB_ACC_ID %in% node.table.gp.mapped$uniprotswissprot) %>%
  filter (KIN_ACC_ID %in% node.table.gp.mapped$uniprotswissprot) %>%
  dplyr::select(SUB_ACC_ID) %>%
  distinct()

curation.nodes.suid <- merge(curation.nodes, node.table.gp.mapped, by.x ="SUB_ACC_ID", by.y ="uniprotswissprot") %>%
  dplyr::select(SUID) %>% 
  pull()

setNodeColorBypass(curation.nodes.suid, '#FF0088')

############
