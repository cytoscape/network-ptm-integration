## Setup
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

BiocManager::install("rWikiPathways")
BiocManager::install("RCy3")
#BiocManager::install("BridgeDbR")
#install.packages("biomaRt")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")

library(rWikiPathways)
library(RCy3)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(purrr)
library(ggplot2)
#library('biomaRt')
cytoscapePing()

##Change working dir for local dev. Adjust path based on setup.
setwd("~/Dropbox (Gladstone)/Work/github/network-ptm-integration/scripts")

###############
## Read in data files

## Phosphoprotein data from CPTAC. Note that this includes NA values. 
## The "pval" and "val" statistics are from Wilcoxon rank sum test. Positive values means abundance is higher in tumor.
## Add unique id symbol_site for each protein/site combination
cptac.phospho <- read.csv("../datasets/CPTAC_phospho_tn.txt", stringsAsFactors = F, sep = "\t")
cptac.phospho <- cptac.phospho %>% 
  mutate(symbol_site=paste0(symbol, "_", site))

## Extracting phospho data for one cancer type: CCRCC
## Select significant rows based on significant p-value in CCRCC vs control
cptac.phospho.ccrcc.sig <- cptac.phospho %>% 
  filter(CCRCC.pval > 0 & CCRCC.pval < 0.05) %>%
  select(symbol, site, protein, symbol_site, CCRCC.val, CCRCC.pval)

##Testing using omicsvisualizer
cptac.phospho.ccrcc.ov <- cptac.phospho.ccrcc.sig %>%
  filter(symbol_site %in% cptac.progeny.egfr.ccrcc.pos$symbol_site) %>%
  mutate(symbol_ptm=paste0(symbol, "_ptm"))

# for testing purposes
#write.table(cptac.phospho.ccrcc.ov,"cptac.phospho.ccrcc.sig.txt",sep="\t",row.names=FALSE)

## PROGENy EGFR pathway correlated phospho data. 
## Sites are either positively correlated with the EGFR pathway (contributing to pathway activity), or negatively correlated (inhibitory) 
## cohort.val: Spearman correlation coefficient;  cohort.pval: p-value
cptac.progeny.egfr <- read.csv("../datasets/CPTAC_PROGENy__EGFR.txt", stringsAsFactors = F, sep = "\t")

## Positively regulated PROGENy EGFR pathway and CCRCC
## To-Do: Not sure how to filter this data using val and pval
cptac.progeny.egfr.ccrcc.pos <- cptac.progeny.egfr %>%
  filter(CCRCC.pval > 0 & CCRCC.pval < 0.05) %>%
  mutate(symbol_site=paste0(symbol, "_", site)) %>% 
  select(symbol, protein, site, symbol_site, CCRCC.val, CCRCC.pval) 

## Negatively regulated PROGENy EGFR pathway and CCRCC
cptac.progeny.egfr.ccrcc.neg <- cptac.progeny.egfr %>%
  filter(CCRCC.pval > -0.05 & CCRCC.pval < 0) %>%
  mutate(symbol_site=paste0(symbol, "_", site)) %>%
  select(symbol, protein, site, symbol_site, CCRCC.val, CCRCC.pval)

## Protein data from CPTAC. Note that this includes NA values.
## The "pval" and "val" statistics are from Wilcoxon rank sum test. Positive values means abundance is higher in tumor.
cptac.protein <- read.csv("../datasets/CPTAC_protein_tn.txt", stringsAsFactors = F, sep = "\t")

## Extracting protein data for one cancer type: CCRCC
cptac.protein.ccrcc <- cptac.protein %>% 
  select(ensembl, symbol, CCRCC.val, CCRCC.pval) %>%
  filter(!is.na(CCRCC.pval) & !is.na(CCRCC.pval))

## Phosphosite kinase-substrate data
psp.data <- read.csv("../annotations/Kinase_Substrate_Dataset.txt", stringsAsFactors = F, sep = "\t")
psp.data.human<- psp.data %>%
  filter(KIN_ORGANISM == "human" & IN_VIVO_RXN == "X")  %>% 
  mutate(symbol_site=paste0(SUBSTRATE, "_", SUB_MOD_RSD)) %>%
  select(GENE, symbol_site, SUBSTRATE, SUB_MOD_RSD) #human only, in vivo evidence only

## Import BioMART mapping data and rename columns. To enable data mapping between node table and CPTAC data using EnsemblProt id.
biomart <- read.csv("../annotations/mart_export.txt", stringsAsFactors = F, sep = ",")

biomart <- biomart %>%
  rename(Ensembl = Gene.stable.ID, EnsemblProt = Protein.stable.ID) %>%
  select(Ensembl, EnsemblProt)

###############
## Open and process WP, get relevant phospho data nodes and cross-check against kinase-substrate data
RCy3::commandsRun('wikipathways import-as-pathway id=WP4806') 

## Remove existing ptm states from nodes. These are encoded as nodes with the label P.
selectNodes(nodes = "p", by.col = "name", preserve.current.selection = FALSE) 
deleteSelectedNodes()

## Default mode is "all", meaning add all ptms from the data regardless if the kinase is on the pathway or not.
mode <- "all"

## Get all relevant protein/gene product data nodes in pathway
node.table <- RCy3::getTableColumns(table = "node")
node.table.prot <- subset(node.table, Type == 'Protein' | Type == 'GeneProduct') %>%
  select(SUID, name, XrefId, Ensembl)

## Add Ensembl prot column to use for all data mapping
node.table.prot.mapped <- merge(node.table.prot, biomart, by="Ensembl") %>%
  filter(!is.na(EnsemblProt)) ##this contains rows where EnsemblProt is empty

## Add all matching phospho nodes: Intersection between pathway protein/gene nodes and significant phospho data
## Using CPTAC data
# matching.nodes.prot <- node.table.prot %>% 
#   filter(HGNC %in% cptac.phospho.ccrcc.sig$symbol) %>% 
#   select(SUID, name) 

## Find pathway protein nodes to add phospho nodes to. Intersection between pathway protein/gene nodes and positively regulated PROGENy sites.
## Using PROGENy data, positive regulation
matching.nodes.prot <- node.table.prot.mapped %>% 
  filter(EnsemblProt %in% cptac.progeny.egfr.ccrcc.pos$protein) %>% 
  select(SUID, name) 

##Get position for all relevant pathway protein nodes, store as data frame
node.positions <- getNodePosition(node.names=matching.nodes.prot$SUID)
node.positions <- cbind(suid = rownames(node.positions), node.positions)
rownames(node.positions) <- 1:nrow(node.positions)

node.width <- getNodeWidth(node.names=matching.nodes.prot$SUID)
node.width <- as.data.frame(node.width)
node.width  <- cbind(suid = rownames(node.width), node.width)
rownames(node.width) <- 1:nrow(node.width)

node.height <- getNodeHeight(node.names=matching.nodes.prot$SUID)
node.height <- as.data.frame(node.height)
node.height <- cbind(suid = rownames(node.height), node.height)
rownames(node.height) <- 1:nrow(node.height)

node.layout <- merge(node.width, node.height, by="suid")
node.layout <- merge(node.layout, node.positions, by="suid")
node.layout.pie <- node.layout #copy

## Traditional ptm viz using 4 sites per node.
## Calculate and add x and y pos.
node.layout <- node.layout %>%
  mutate(x.ptm = pmap(list(x_location, node.width), 
                            function(x_val, width_val) {
                              list(x_val + (width_val / 2) + 16,x_val + (width_val / 2) + 16,x_val - (width_val / 2) - 16,x_val - (width_val / 2) - 16)}))
node.layout <- node.layout %>%
  mutate(y.ptm = pmap(list(y_location, node.height), 
                      function(y_val, height_val) {
                        list(y_val + (height_val / 2),y_val - (height_val / 2),y_val - (height_val / 2),y_val + (height_val / 2))}))

## Alternative viz with one ptm node with pie chart
## Calculate and add x and y pos.
node.layout.pie <- node.layout.pie %>%
  mutate(x.ptm = (x_location + node.width / 2 + 20)) %>%
  mutate(y.ptm = y_location )

## Need this for matching with kinase substrate data
# matching.nodes.phospho <- cptac.phospho.ccrcc.sig %>%
#   filter(symbol %in% node.table.prot$name) %>%
#   mutate(symbol_site=paste0(symbol, "_", site)) %>%
#   select(symbol, site, symbol_site)

## Use PROGENy data to select ptms to add
matching.nodes.phospho <- cptac.progeny.egfr.ccrcc.pos %>%
  filter(symbol %in% node.table.prot$name) %>%
  mutate(symbol_site=paste0(symbol, "_", site)) %>%
  select(symbol, site, symbol_site)

###############
## Kinase-restricted mode: Only those sites with kinases on pathway: Subset of matching nodes with relevant kinases in the pathway
## To-Do: redo for SUID
kinase.pw <- intersect(psp.data.human$GENE, node.table.prot$name) #Kinases on pathway. GENE column contains kinase.

## Overlap between phospho nodes and kinase-substrate data, where the kinase (GENE) is on the pathway
matching.nodes.phospho.kin <- inner_join(matching.nodes.phospho, psp.data.human) %>%
  filter(GENE %in% kinase.pw) %>%
  select(symbol, symbol_site, site)

matching.nodes.phospho.kin.sites <- matching.nodes.phospho.kin %>%
  select(symbol, symbol_site) 

#mode <- "kinase"

###############

## Traditional viz with max 4 ptms per protein node
## Adding phospho nodes to pathway
## Add phospho nodes by looping through relevant protein nodes

site.table <- matching.nodes.phospho 

## If mode is "kinase only", update inputs
if (mode == "kinase"){
  site.table <- matching.nodes.phospho.kin.sites
}

ptms.all <- data.frame()

for (p in matching.nodes.prot$SUID){
  prot <- matching.nodes.prot$name[matching.nodes.prot$SUID == p]
  #find relevant subset of cptac phospho data, these are the ptm nodes to add
  phospho.nodes <- site.table %>% 
    filter(symbol == prot) %>% 
    select(symbol_site)
  phospho.nodes <- head(phospho.nodes, 4) ##take first 4
  
  ##Add a node for each ptm, for the particular protein. Collect the SUID for the added ptm, and use it for updating position below.
  suids <- addCyNodes(node.names = phospho.nodes$symbol_site, skip.duplicate.names = FALSE)
  ptms <- data.frame()

  for (l in suids){
    ptms <- rbind(ptms, l)
  }
  
  ##Set node position bypass for ptm nodes.
  pos <- 1
  for (ptm in ptms$SUID){
    #addCyEdges(c(p, ptm)) #need the SUID for this
    x <- node.layout$x.ptm[node.layout$suid == p][[1]][[pos]]
    y <- node.layout$y.ptm[node.layout$suid == p][[1]][[pos]]
      setNodePositionBypass(ptm, x, y)
    pos <- pos+1
  }
  ptms.all <- rbind(ptms.all, ptms)
}

###############

## Alt viz mockup: pie charts using OmicsVisualizer

#make new data frame with correct name
matching.nodes.prot.pie <- node.table.prot.mapped %>% 
  filter(EnsemblProt %in% cptac.progeny.egfr.ccrcc.pos$protein) %>% 
  mutate(name=paste0(name,"_ptm")) %>% 
  select(SUID, name) 

#Add new ptm node for all relevant proteins
for (p in matching.nodes.prot.pie$SUID){
  #add node
  ptm.name <- matching.nodes.prot.pie$name[matching.nodes.prot.pie$SUID == p]
  print(ptm.name)
  suid.list <- addCyNodes(node.names = ptm.name, skip.duplicate.names = FALSE)
  suid.ptm <- suid.list[[1]]$SUID
  print(suid.ptm)
  
  #move node
  x <- node.layout.pie$x.ptm[node.layout.pie$suid == p]
  y <- node.layout.pie$y.ptm[node.layout.pie$suid == p]
  
  print(x)
  print(y)
  setNodePositionBypass(node.names = suid.ptm, x, y)
}

#OmicsVisualizer
#Import table from data frame (cptac.phospho.ccrcc.ov). Not sure this is possible via ov commands.....
#Current solution for tesing: export table to csv (line 48) and import using UI

#Link table
#commandsRun('ov connect symbol_ptm symbol_site')

#Add chart visualization
#commandsRun(paste0('ov ov viz apply inner continuous idAttribute="Ensembl" linkSetFiles="', drugbank, '"') )

##############
## Data viz
style.name = "WikiPathways"

## Update defaults
setNodeColorDefault('#FFFFFF', style.name = style.name)
setNodeBorderColorDefault("#737373", style.name = style.name)

## Update style for ptm nodes
setNodeFontSizeBypass(ptms.all$SUID, 9)
setNodeWidthBypass(ptms.all$SUID, 35) 
setNodeHeightBypass(ptms.all$SUID, 20)
#to-do: bring phospho nodes to the front. not sure if this is possible with RCy3

## Load protein data and visualize as node color.
## Alt 1: Visualize CPTAC protein and CPTAC phospho data on node fill for protein and ptm nodes.
## For this we can use the same visualization since the value is the same, Wilcoxon rank sum test

loadTableData(cptac.protein.ccrcc, data.key.column="ensembl", "node", table.key.column = 'Ensembl') ##load protein data
loadTableData(cptac.phospho.ccrcc.sig, data.key.column="symbol_site", "node", table.key.column = 'shared name') ##load phospho data
RCy3::setNodeColorMapping('CCRCC.val', colors=paletteColorBrewerRdBu, style.name = style.name) 
ptms.all.data <- inner_join(ptms.all, cptac.phospho, by = join_by(name == symbol_site)) ##add back site info etc

# ## Alt 2: Visualize CPTAC protein on protein node, and PROGENy data on ptm nodes.
# ## Since the data is different (Wilcoxon vs Spearman correlation), the ptm node color will be done via bypass.
# loadTableData(cptac.egfr.ccrcc.pos, data.key.column="symbol_site", "node", table.key.column = 'shared name') ##load phospho data
# 
# ptms.all.data <- inner_join(ptms.all, cptac.egfr.ccrcc.pos, by = join_by(name == symbol_site)) %>%
#   mutate(color = case_when(CCRCC.val > 0.5 ~ "#f06262", CCRCC.val < 0.5 ~ "#f0a3a3", CCRCC.val == 0.5 ~ "#f0a3a3"))
# setNodeColorBypass(node.names = ptms.all.data$SUID, new.colors = ptms.all.data$color)


##Create new df to map to network to update ptm node label
ptms.all.data.site <- ptms.all.data %>%
  select(SUID, site) %>%
  rename(name = site)
loadTableData(ptms.all.data.site,data.key.column="SUID", "node", table.key.column = "SUID")
