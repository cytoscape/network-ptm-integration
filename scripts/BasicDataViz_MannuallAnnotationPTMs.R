###############

## This script demonstrates how to perform basic data visualization in Cytoscape on manually curated ptms in WikiPathways models.
## The example data files are from CPTAC (https://proteomics.cancer.gov/programs/cptac) and use Ensembl Protein (peptide) identifiers.
## Adjustments to the script will need to be made for data with other identifier types, like Uniprot.
## The strategy is to visualize the same data type, "val" in the example data, for both proteomics data and phosphoproteomics data on the main data node 
## and the ptm "state" node respectively. "val" and "pval" are from Wilcoxon rank sum test. Positive values means abundance is higher in tumor.

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

# Ensure Cytoscape is running and connected
cytoscapePing()

## Change working dir. Adjust path based on setup.
setwd("~/github/network-ptm-integration/scripts")

###############
## Read in data files. The example data files are from CPTAC (https://proteomics.cancer.gov/programs/cptac). 

## Phosphoproteomics data from CPTAC, for multiple cancer types.
## The "pval" and "val" statistics are from Wilcoxon rank sum test. Positive values means abundance is higher in tumor.
cptac.phospho <- read.csv("../datasets/CPTAC_phospho_tn.txt", stringsAsFactors = F, sep = "\t")

## Proteomics data from CPTAC, for multiple cancer types.
## The "pval" and "val" statistics are from Wilcoxon rank sum test. Positive values means abundance is higher in tumor.
cptac.protein <- read.csv("../datasets/CPTAC_protein_tn.txt", stringsAsFactors = F, sep = "\t")

## Get BioMART mapping info for Ensembl protein id to Uniprot-Swissprot, to enable data mapping between node table and CPTAC data using EnsemblProt id.
## Manually annotated ptms on pathways are annotated with Uniprot-Swissprot identifiers.
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensemblprot_swissprot <- getBM(attributes=c('ensembl_peptide_id','hgnc_symbol','uniprotswissprot', 'ensembl_gene_id'), mart = ensembl)
ensemblprot_swissprot_sub <- ensemblprot_swissprot %>%
  filter(!is.na(uniprotswissprot), uniprotswissprot != "")

ensemblprot_swissprot_sub <- ensemblprot_swissprot_sub[, c("ensembl_peptide_id", "uniprotswissprot")]
ensembl_swissprot_sub <- ensemblprot_swissprot[, c("ensembl_gene_id", "uniprotswissprot")] %>%
  filter(!is.na(uniprotswissprot), uniprotswissprot != "")
ensembl_swissprot_sub <- unique(ensembl_swissprot_sub)

###############
## Open the relevant WP in Cytoscape. The below example uses the Head and neck squamous cell carcinoma pathway.
RCy3::commandsRun('wikipathways import-as-pathway id=WP4674') 

## Get the full node table for the pathway.
node.table <- RCy3::getTableColumns(table = "node")
## Get the ptm nodes. This will be used later for manipulating the visual style of ptm nodes.
ptm.nodes <- node.table %>%
  filter(ptm == 'p')

## Add a new column, data_mapping_id, for data mapping. A new data frame is created with the new column and then read into the Cytoscape node table.
## The new data mapping column will contain a new id which is a composite of the parent id (Uniprot) and ptm site information.
## The site information is also switched from the format from a three-letter amino acid code to one-letter code, for example ser233 to S233, since the
## proteomics data has the ptm site information in this format. The resulting composite id is in the form P27361_T202.

node.table.ptm <- ptm.nodes %>% 
  mutate(data_mapping_id=paste0(parentid, "_", position)) %>% 
  mutate(data_mapping_id = str_replace(data_mapping_id, "ser", "S"),
         data_mapping_id = str_replace(data_mapping_id, "thr", "T"),
         data_mapping_id = str_replace(data_mapping_id, "tyr", "Y")) %>%
  mutate(name=position) %>%
  mutate(name = str_replace(name, "^ser", "S"),
         name = str_replace(name, "^thr", "T"),
         name = str_replace(name, "^tyr", "Y")) %>%
  dplyr::select(SUID, data_mapping_id, name)

loadTableData(node.table.ptm, data.key.column="SUID", "node", table.key.column = 'SUID') ## load the new column into Cytoscape

## Copy Ensembl ids from "Ensembl" column to data_mapping_id column. For the proteomics data, this will be the mapping id.
node.table.ensembl <- node.table %>% 
  filter(Type == 'GeneProduct') %>%
  mutate(data_mapping_id=Ensembl) %>% 
  dplyr::select(SUID, data_mapping_id)

loadTableData(node.table.ensembl, data.key.column="SUID", "node", table.key.column = 'SUID') ## load the new column into Cytoscape

## The ptm node label "P" is controlled by a node label bypass. We will remove this bypass to allow for the ptm nodes to be labeled by the amino acid site.
for (s in ptm.nodes$SUID) {
  clearNodePropertyBypass(s, visual.property = "NODE_FILL_COLOR")
  clearNodePropertyBypass(s, visual.property = "NODE_LABEL")
}

## The above for loop is a workaround due to a bug in RCy3. Once the bug is fixed, it should be updated to:
##clearNodePropertyBypass(node.names = ptm.nodes$SUID, visual.property = "NODE_FILL_COLOR") ## this doesnt work
##clearNodePropertyBypass(node.names = ptm.nodes$SUID, visual.property = "NODE_LABEL") ## this doesnt work

setNodeFontSizeBypass(ptm.nodes$SUID, 9)

## Map phosphoproteomics data from Ensembl Prot ids to Swissprot using Biomart data mapping.
cptac.phospho.mapped <- merge(cptac.phospho, ensemblprot_swissprot_sub, by.x = "protein", by.y = "ensembl_peptide_id") %>%
  mutate(prot_site=paste0(uniprotswissprot, "_", site))

## Load phosphoproteomics and proteomics data into Cytoscape, using the newly added data_mapping_id column.
loadTableData(cptac.phospho.mapped, data.key.column="prot_site", "node", table.key.column = 'data_mapping_id') 
loadTableData(cptac.protein, data.key.column="ensembl", "node", table.key.column = 'data_mapping_id') 

## Set visual style
style.name = "WikiPathways"
setNodeColorDefault('#FFFFFF', style.name = style.name)
setNodeBorderColorDefault("#737373", style.name = style.name)

## Map node color to HNSCC.val (Head and neck squamous cell carcinoma)
RCy3::setNodeColorMapping('HNSCC.val', colors=paletteColorBrewerRdBu, style.name = style.name) 

###############