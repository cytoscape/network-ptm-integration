###############

## This script demonstrates how to perform basic data visualization on manually cureated ptms on WikiPathways models.
## The example data files are from CPTAC (https://proteomics.cancer.gov/programs/cptac) and use Ensembl Protein (peptide) identifiers.
## Adjustments to the script will need to be made for data with other identifier types, like Uniprot.
## The strategy is to visualize the same data type, "val" for both proteomics data and phosphoproteomics data on the main data node 
## and the ptm "state" node respectively.

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
cytoscapePing()

## Change working dir for local dev. Adjust path based on setup.
setwd("~/github/network-ptm-integration/scripts")

###############
## Read in data files. The example data files are from CPTAC (https://proteomics.cancer.gov/programs/cptac). 

## Phosphoprotein example from CPTAC, for Clear cell renal cell carcinoma (CCRCC) vs normal
## The "pval" and "val" statistics are from Wilcoxon rank sum test. Positive values means abundance is higher in tumor.
cptac.phospho <- read.csv("../datasets/CPTAC_phospho_tn.txt", stringsAsFactors = F, sep = "\t")

## Protein data from CPTAC, for Clear cell renal cell carcinoma (CCRCC) vs normal. 
## The "pval" and "val" statistics are from Wilcoxon rank sum test. Positive values means abundance is higher in tumor.
cptac.protein <- read.csv("../datasets/CPTAC_protein_tn.txt", stringsAsFactors = F, sep = "\t")

## Get BioMART mapping info for Ensembl protein id to Uniprot-Swissprot, to enable data mapping between node table and CPTAC data using EnsemblProt id.
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensemblprot_swissprot <- getBM(attributes=c('ensembl_peptide_id','hgnc_symbol','uniprotswissprot'), mart = ensembl)
ensemblprot_swissprot_sub <- ensemblprot_swissprot %>%
  filter(!is.na(uniprotswissprot), uniprotswissprot != "")

ensemblprot_swissprot_sub <- ensemblprot_swissprot_sub[, c("ensembl_peptide_id", "uniprotswissprot")]

###############
## Open and process WP, get relevant phospho data nodes
## Open the relevant WP. The below example uses the EGFR pathway.
RCy3::commandsRun('wikipathways import-as-pathway id=WP4806') 

## Get all relevant protein/gene product data nodes in pathway
node.table <- RCy3::getTableColumns(table = "node")
## Get the ptm nodes. This will be used later for manipulating the visual style of ptm nodes
ptm.nodes <- node.table %>%
  filter(ptm == 'p')

## Add new column, data_mapping_id, for data mapping.
node.table.ptm <- node.table %>% 
  filter(ptm == 'p') %>%
  mutate(data_mapping_id=paste0(parentid, "_", position)) %>% 
  mutate(name=position) %>% 
  select(SUID, data_mapping_id, name)

loadTableData(node.table.ptm, data.key.column="SUID", "node", table.key.column = 'SUID') ## load the new column into Cytoscape

## Add Ensembl ids to data_mapping_id column
node.table.ensembl <- node.table %>% 
  filter(Type == 'GeneProduct') %>%
  mutate(data_mapping_id=Ensembl) %>% 
  select(SUID, data_mapping_id)

loadTableData(node.table.ensembl, data.key.column="SUID", "node", table.key.column = 'SUID') ## load the new column into Cytoscape

## Map phosphoproteomics data from Ensembl Prot ids to Swissprot using Biomart data mapping
cptac.phospho.mapped <- merge(cptac.phospho, ensemblprot_swissprot_sub, by.x = "protein", by.y = "ensembl_peptide_id") %>%
  mutate(site = str_replace(site, "^S", "ser"),
         site = str_replace(site, "^T", "thr"),
         site = str_replace(site, "^Y", "tyr")) %>%
  mutate(prot_site=paste0(uniprotswissprot, "_", site))

## Load data
loadTableData(cptac.phospho.mapped, data.key.column="prot_site", "node", table.key.column = 'data_mapping_id') ##load phospho data
loadTableData(cptac.protein, data.key.column="ensembl", "node", table.key.column = 'data_mapping_id') ##load phospho data

## Set visual style
style.name = "WikiPathways"
setNodeColorDefault('#FFFFFF', style.name = style.name)
setNodeBorderColorDefault("#737373", style.name = style.name)

## Remove bypasses for ptm nodes that will interfere with data visualization
setNodeFontSizeBypass(ptm.nodes$SUID, 9)
##clearNodePropertyBypass(node.names = ptm.nodes$SUID, visual.property = "NODE_FILL_COLOR") ## this doesnt work
##clearNodePropertyBypass(node.names = ptm.nodes$SUID, visual.property = "NODE_LABEL") ## this doesnt work

RCy3::setNodeColorMapping('CCRCC.val', colors=paletteColorBrewerRdBu, style.name = style.name) 

###############