## Setup
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
BiocManager::install("rWikiPathways")
}

if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(rWikiPathways)
library(RCy3)
install.packages("dplyr")
library(dplyr)
library(stringr)
cytoscapePing()

setwd("~/Dropbox (Gladstone)/Work/github/network-ptm-integration/scripts")

###############
## Read in data

##Phosphoprotein data
cptac.data <- read.csv("~/Dropbox (Gladstone)/Work/CPTAC/PROGENy__EGFR.txt", stringsAsFactors = F, sep = "\t")
cptac.data.brca.pval <- cptac.data %>% filter(BRCA.pval > 0 & BRCA.pval < 0.05) ##dataset specific: keep only rows with data. Change cancer type if necessary.
# phospho.data <- read.csv("../datasets/Francavilla_PMID28355574_mmc4_phospho.csv", stringsAsFactors = F)

##Protein data. Rename columns.
# protein.data <- read.csv("../datasets/Francavilla_PMID28355574_mmc3_protein.csv", stringsAsFactors = F)
# protein.data.avgexp <- protein.data %>%
#   dplyr::mutate(AvgExpProtein = AvgExp) %>%
#   select(Gene.names, AvgExpProtein)

##Phosphosite kinase-substrate data
psp.data <- read.csv("~/Dropbox (Gladstone)/Work/CPTAC/network-ptm-integration/CPTAC data downloads/Kinase_Substrate_Dataset.txt", stringsAsFactors = F, sep = "\t")
psp.data.human<- psp.data %>% filter(KIN_ORGANISM == "human" & IN_VIVO_RXN == "X") #human only, in vivo evidence only

###############
##Open and process WP, get relevant phospho data nodes and cross-check against kinase-substrate data
RCy3::commandsRun('wikipathways import-as-pathway id=WP4806') 
mode <- "all"

##Map all nodes to Uniprot. Skip this if mapping via symbol.
#mapped.cols <- mapTableColumn('XrefId','Human','Ensembl','Uniprot', force.single=FALSE) ##doesn't work like expected

##Get all relevant data nodes in pathway
node.table <- getTableColumns()
node.table.prot <- subset(node.table, Type == 'Protein' | Type == 'GeneProduct') %>% select(SUID, name, XrefId, Ensembl) ##get node table entries for relevant nodes

##All matching nodes: Intersection between pathway protein/gene nodes and phospho data
matching.nodes.prot <- node.table.prot %>% filter(name %in% cptac.data.brca.pval$symbol)
matching.nodes.phospho <- cptac.data.brca.pval %>% filter(symbol %in% node.table.prot$name)

#cptac.data.ptm <- cptac.data.pval[cptac.data.pval$symbol %in% ptm.nodes, ] #Relevant subset of phospho data
#phospho.nodes <- cptac.data.ptm$site 

source.target <- matching.nodes.phospho %>% select(symbol, site)
##Convert to list format that RCy3 expects
source.target.list <- list() 
for (i in 1:nrow(source.target)) {
  row <- source.target[i, ]  
  row_list <- as.list(row)  
  source.target.list[[i]] <- row_list  
}

##Only those with kinases on pathway: Subset of matching nodes with relevant kinases in the pathway
kinase.pw <- intersect(psp.data$GENE, node.table.prot$name) #Read in kinase-substrate data

#create temp data frames for comparisons
psp.data.sub <- psp.data %>% mutate(comp_id=paste0(SUBSTRATE, "_", SUB_MOD_RSD)) %>% select(GENE, comp_id)
cptac.data.sub <- matching.nodes.phospho %>% mutate(comp_id=paste0(symbol, "_", site)) %>% select(symbol, comp_id, site)
psp.cptac.sub <- inner_join(psp.data.sub, cptac.data.sub) #combined
psp.cptac.kin <- psp.cptac.sub[psp.cptac.sub$GENE %in% kinase.pw, ] #subset for only sites with kinases on the pathway

#define phospho site nodes and edges
phospho.nodes.kin <- psp.cptac.kin$site 
source.target.kin <- psp.cptac.kin %>% select(symbol, site)
ptm.nodes.kin <- psp.cptac.kin$symbol

###############
##Adding phospho nodes to pathway

##To-Do: Select mode
if (mode == "kinase"){
  phospho.nodes <- phospho.nodes.kin
  source.target <- source.target.kin
  ptm.nodes <- ptm.nodes.kin
}

##To-Do: Change this to use composite id to make it be always unique
phospho.suid <- addCyNodes(node.names=matching.nodes.phospho$site, skip.duplicate.names = FALSE) ##phospho sites are not unique
#addCyNodes(phospho.nodes, skip.duplicate.names = FALSE) ##phospho sites are not unique
addCyEdges(source.target.list)
node.table.new <- getTableColumns() #get updated node table
phospho.suid <- as.data.frame(do.call(rbind, phospho.suid))


## Update ptm node position to be slightly offset from the parent node
## To-Do: better logic for new node position, multiple nodes etc

#Plan: Loop through all regular nodes. For each, lookup which phospho nodes belong, then adjust node position accordingly for each
node.positions <- getNodePosition(node.names = matching.nodes.prot$SUID)
node.positions <- cbind(SUID=rownames(node.positions), node.positions)

#for (n in cptac.data.sub$symbol){
for (n in matching.nodes.prot$SUID){
  #define positions for phospho nodes based on x,y
  pos.x <- node.positions$x_location[node.positions$symbol == n]
  pos.y <- node.positions$y_location[node.positions$symbol == n]
  positions.phospho <- data.frame(position=1:4, x=c(pos.x+30, pos.x+30, pos.x-30, pos.x-30), y=c(pos.y+10, pos.y-10, pos.y-10, pos.y+10))
  # 
  ##To-Do: complete this section using SUIDs
  # #get phospho sites for this node and check how many. If more than 4, delete extras
  # p.nodes <- matching.data.phospho %>% subset(symbol == n)
  # if (length(p.nodes) > 4){
  #   p.nodes <- p.nodes[1:4, ]
  # }
  # 
  # pos <- 1
  # #test <- subset(positions.phospho, position == "1")$x
  # for (p in p.nodes){
  #   x <- subset(positions.phospho, position == pos )$x
  #   y <- subset(positions.phospho, position == pos)$y
  #   setNodePositionBypass(p, x, y)
  #   pos <- pos+1
  # }
}

##############
## Data viz
style.name = "WikiPathways"

#To-Do: use cancer-type specific data (pval) to create mapping
#setNodeColorBypass

#setNodeShapeBypass(phospho.node.names, 'ellipse')
#setNodeColorMapping('AvgExp', colors=paletteColorBrewerRdBu, style.name=style.name)

## Update label for phospho nodes
# phospho.nodes.position <- phospho.data.sub %>% pull(Aa.Position) %>% as.list()
# selectNodes(phospho.node.names, by='name')
# setNodeLabelBypass(phospho.node.names, phospho.data.sub[,c("Aa.Position")])

## Load protein data
#loadTableData(protein.data.avgexp,data.key.column="Gene.names", "node", table.key.column = 'name')
#setNodeCustomHeatMapChart(c("AvgExpProtein"), style.name=style.name)

## Load data as table
#loadTableData(combined.data.column,data.key.column="node.name", "node", table.key.column = 'name')

## Apply viz
#setNodeColorMapping('new.avg.exp', colors=paletteColorBrewerRdBu, style.name=style.name)
