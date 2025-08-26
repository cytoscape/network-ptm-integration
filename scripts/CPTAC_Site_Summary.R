###############

## Data viz of CPTAC site summary information

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
library(ndexr)

# Ensure Cytoscape is running and connected
cytoscapePing()

## Change working dir. Adjust path based on setup.
setwd("~/github/network-ptm-integration/scripts")

## Define Cytoscape version variable dynamically
cytoscapeVersion <- cytoscapeVersionInfo()[2]
cytoscapeSampleDataPath <- paste0("/Applications/Cytoscape_v", cytoscapeVersion, "/sampleData/")

###############
## Read in data file

## Read in local BioMART file
biomart <- read.csv("../annotations/id-mapping/ensembl_mappings.txt", stringsAsFactors = FALSE, sep = "\t") %>%
  dplyr::select(ensembl_gene_id,ensembl_peptide_id)

## CPTAC site summary data
## Note: The file is not in the repo. Adjust path for locally stored file.
cptac.site.summary <- read.csv("site_summary_for_pathway.csv", stringsAsFactors = F, sep = ",")
cptac.site.summary <- merge(cptac.site.summary, biomart, by.x = "harm_protein", by.y = "ensembl_peptide_id") %>%
  mutate(newnodename = paste0(ensembl_gene_id, "_ptm")) ## add Ensembl gene id and newnodename id for data mapping.

################
##Setup individual data frames for each count type
## Prep data frames to keep top 8 for each count

cptac.site.summary.text <- cptac.site.summary %>% 
  dplyr::select(substrate, harm_protein, harm_site, text_evidence_count, ensembl_gene_id, newnodename) %>%
  filter(text_evidence_count != "NA")
cptac.site.summary.text <- cptac.site.summary.text %>%
  group_by(harm_protein) %>%
  mutate(site_count = n()) %>%                 # count how many times harm_protein appears
  arrange(desc(text_evidence_count)) %>%       # sort by text_evidence_count descending
  slice_head(n = 8) %>%                        # keep top 8 per harm_protein
  ungroup()

cptac.site.summary.figure <- cptac.site.summary %>% 
  dplyr::select(substrate, harm_protein, harm_site, figure_evidence_count, ensembl_gene_id, newnodename) %>%
  filter(figure_evidence_count != "NA") 
cptac.site.summary.figure <- cptac.site.summary.figure %>%
  group_by(harm_protein) %>%
  mutate(site_count = n()) %>%                 # count how many times harm_protein appears
  arrange(desc(figure_evidence_count)) %>%     # sort by figure_evidence_count descending
  slice_head(n = 8) %>%                        # keep top 8 per harm_protein
  ungroup()

cptac.site.summary.exp <- cptac.site.summary %>% 
  dplyr::select(substrate, harm_protein, harm_site, experiment_evidence_count, ensembl_gene_id, newnodename) %>%
  filter(experiment_evidence_count != "NA")
cptac.site.summary.exp <- cptac.site.summary.exp %>%
  group_by(harm_protein) %>%
  mutate(site_count = n()) %>%                 # count how many times harm_protein appears
  arrange(desc(experiment_evidence_count)) %>%     # sort by figure_evidence_count descending
  slice_head(n = 8) %>%                        # keep top 8 per harm_protein
  ungroup()

cptac.site.summary.cohort <- cptac.site.summary %>% 
  dplyr::select(substrate, harm_protein, harm_site, cohort_evidence_count, ensembl_gene_id, newnodename) %>%
  filter(cohort_evidence_count != "NA") 
cptac.site.summary.cohort <- cptac.site.summary.cohort %>%
  group_by(harm_protein) %>%
  mutate(site_count = n()) %>%                 # count how many times harm_protein appears
  arrange(desc(cohort_evidence_count)) %>%     # sort by cohort_evidence_count descending
  slice_head(n = 8) %>%                        # keep top 8 per harm_protein
  ungroup()

cptac.site.summary.total <- cptac.site.summary %>% 
  dplyr::select(substrate, harm_protein, harm_site, total_count, ensembl_gene_id, newnodename) %>%
  filter(total_count != "NA") 
cptac.site.summary.total <- cptac.site.summary.total %>%
  group_by(harm_protein) %>%
  mutate(site_count = n()) %>%                 # count how many times harm_protein appears
  arrange(desc(total_count)) %>%               # sort by total_count descending
  slice_head(n = 8) %>%                        # keep top 8 per harm_protein
  ungroup()

################

##Export files. This is necessary for omicsvisualizer
write.table(cptac.site.summary.text, 
            paste0(cytoscapeSampleDataPath, "cptac.site.summary.text.txt"),
            sep = "\t", row.names = FALSE)
write.table(cptac.site.summary.figure, 
            paste0(cytoscapeSampleDataPath, "cptac.site.summary.figure.txt"),
            sep = "\t", row.names = FALSE)
write.table(cptac.site.summary.exp, 
            paste0(cytoscapeSampleDataPath, "cptac.site.summary.exp.txt"),
            sep = "\t", row.names = FALSE)
write.table(cptac.site.summary.cohort, 
            paste0(cytoscapeSampleDataPath, "cptac.site.summary.cohort.txt"),
            sep = "\t", row.names = FALSE)
write.table(cptac.site.summary.total, 
            paste0(cytoscapeSampleDataPath, "cptac.site.summary.total.txt"),
            sep = "\t", row.names = FALSE)

###############
# ## Open the relevant WP in Cytoscape using rWikiPathways. Below are some exmaples.
#RCy3::commandsRun('wikipathways import-as-pathway id=WP4806') ## EGFR pathway
#RCy3::commandsRun('wikipathways import-as-pathway id=WP2828') #bladder cancer
#RCy3::commandsRun('wikipathways import-as-pathway id=WP4155') #endometrial cancer. issue with adding ptm nodes

RCy3::commandsRun('wikipathways import-as-pathway id=WP4018') # Clear cell renal cell carcinoma

## Get the full node table for the pathway.
node.table <- RCy3::getTableColumns(table = "node")
## Find and delete the ptm nodes. 
selectNodes(nodes = "p", by.col = "name", preserve.current.selection = FALSE) 
deleteSelectedNodes()

## Add EnsemblProt ids for data mapping
node.table.prot <- node.table %>%
  filter(Type %in% c("Protein", "GeneProduct"))

## Get node positions for all protein nodes
node.positions <- getNodePosition(node.names = node.table.prot$SUID)
node.positions <- cbind(suid = rownames(node.positions), node.positions)
rownames(node.positions) <- 1:nrow(node.positions)
node.width <- getNodeWidth(node.names = node.table.prot$SUID)
node.width <- as.data.frame(node.width)
node.width <- cbind(suid = rownames(node.width), node.width)
rownames(node.width) <- 1:nrow(node.width)
node.height <- getNodeHeight(node.names = node.table.prot$SUID)
node.height <- as.data.frame(node.height)
node.height <- cbind(suid = rownames(node.height), node.height)
rownames(node.height) <- 1:nrow(node.height)
node.layout <- merge(node.width, node.height, by = "suid")
node.layout <- merge(node.layout, node.positions, by = "suid")
node.layout.pie <- node.layout  # For alternative (pie chart) viz

##Calculate layout for new nodes
node.layout.pie <- node.layout.pie %>%
  mutate(x.ptm = (x_location + node.width / 2 + 20)) %>%
  mutate(y.ptm = y_location)

matching.nodes.prot.pie <- node.table.prot %>% 
  filter(Ensembl %in% cptac.site.summary.cohort$ensembl_gene_id) %>% 
  mutate(name = paste0(Ensembl, "_ptm")) %>% 
  dplyr::select(SUID, name)

ptmnodes <- data.frame()

for (p in matching.nodes.prot.pie$SUID) {
  ptm.name <- matching.nodes.prot.pie$name[matching.nodes.prot.pie$SUID == p]
  print(ptm.name)
  suid.list <- addCyNodes(node.names = ptm.name, skip.duplicate.names = FALSE)
  suid.ptm <- suid.list[[1]]$SUID
  x <- node.layout.pie$x.ptm[node.layout.pie$suid == p]
  y <- node.layout.pie$y.ptm[node.layout.pie$suid == p]
  setNodePositionBypass(node.names = suid.ptm, x, y)
  # setNodeWidthBypass(node.names = suid.ptm, 40) ##doesnt work
  # setNodeHeightBypass(node.names = suid.ptm, 40) ##doesnt work
  ptmnodes <- rbind(ptmnodes, suid.ptm)
}

ptmnodeslist <- ptmnodes %>% pull() ## dataframe to list

ovload.cmd <- paste('ov load',
                    'dataTypeList="string,string,string,int,string,int,string"', 
                    paste0('file="', cytoscapeSampleDataPath, '"cptac.site.summary.total.txt"'),
                    'newTableName="cptac.site.summary.total"', 
                    'startLoadRow="2"')
commandsRun(ovload.cmd)
ovconnect.cmd <- paste('ov connect',
                       'mappingColNet="Ensembl"', 
                       'mappingColTable="ensembl_gene_id"')
commandsRun(ovconnect.cmd)
ovviz.cmd <- paste('ov viz apply inner continuous',
                   paste0('attributes="total_count"'), 
                   'paletteName="Orange shades"', 
                   'labels="harm_site"')
commandsRun(ovviz.cmd)

#Select all ptm nodes, update the node label to space (empty) by setting bypass
setNodeLabelBypass(ptmnodeslist, ' ')

## workaround: select ptm nodes to manually update node size in UI
selectNodes(nodes = ptmnodeslist, preserve.current.selection = FALSE) 

