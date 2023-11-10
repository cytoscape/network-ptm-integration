## Setup

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rWikiPathways")

if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(RCy3)

install.packages("dplyr")
library(dplyr)
library(stringr)

cytoscapePing()
setwd("~/Dropbox (Gladstone)/Work/github/network-ptm-integration/scripts")

###############
## Read in data

##Phosphoprotein data
phospho.data <- read.csv("../datasets/Francavilla_PMID28355574_mmc4_phospho.csv", stringsAsFactors = F)
phospho.data.pval <- phospho.data[phospho.data$EOC.vs.FTE.p.value > 0, ] ##dataset specific: keep only rows with data

##Protein data. Rename columns.
protein.data <- read.csv("../datasets/Francavilla_PMID28355574_mmc3_protein.csv", stringsAsFactors = F)
protein.data.avgexp <- protein.data %>%
  dplyr::mutate(AvgExpProtein = AvgExp) %>%
  select(Gene.names, AvgExpProtein)

###############
## Open and process WP.
RCy3::commandsRun('wikipathways import-as-pathway id=WP4205') 

## Map all nodes to Uniprot
mapped.cols <- mapTableColumn('XrefId','Human','Ensembl','Uniprot', force.single=FALSE) ##doesn't work like expected

## Get all relevant data nodes
node.table <- getTableColumns()
node.table <- subset(node.table, Type == 'Protein' | Type == 'GeneProduct') ##get node table entries for relevant nodes

## Intersection between pathway and phosphoprotein data (ALL mode)
intersection <- intersect(phospho.data.pval$Gene.names, node.table$name)

## Subset of phosphoprotein data relevant to intersection
phospho.data.sub <- phospho.data.pval[phospho.data.pval$Gene.names %in% intersection, ]
phospho.nodes <- phospho.data.sub$UniProt_AAPosition 
source.target <- phospho.data.sub %>% select(Gene.names, UniProt_AAPosition)

## Convert to list format that RCy3 expects
source.target.list <- list() 
for (i in 1:nrow(source.target)) {
  row <- source.target[i, ]  
  row_list <- as.list(row)  
  source.target.list[[i]] <- row_list  
}

addCyNodes(phospho.nodes)
addCyEdges(source.target.list)

## Update ptm node position to be slightly offset from the parent node
## To-Do: better logic for new node position, multiple nodes etc
## This part doesn't work
node.positions <- getNodePosition(node.names = intersection)
node.positions <- cbind(Gene.names=rownames(node.positions), node.positions)
#new.node.positions <- node.positions %>% mutate(~ . + 30)
source.target.positions <- left_join(source.target, node.positions, by = "Gene.names")

new.x.positions <- as.list(source.target.positions$x_location)
new.y.positions <- as.list(source.target.positions$y_location)

phospho.node.names <- as.list(phospho.nodes)

setNodePositionBypass(phospho.node.names, new.x.positions, new.y.positions) ##this doesn't work

##############
## Data viz
loadTableData(phospho.data.pval,data.key.column="UniProt_AAPosition", "node", table.key.column = 'name')
style.name = "WikiPathways"
setNodeShapeBypass(phospho.node.names, 'ellipse')
setNodeColorMapping('AvgExp', colors=paletteColorBrewerRdBu, style.name=style.name)

## Update label for phospho nodes
phospho.nodes.position <- phospho.data.sub %>% pull(Aa.Position) %>% as.list()
selectNodes(phospho.node.names, by='name')
setNodeLabelBypass(phospho.node.names, phospho.data.sub[,c("Aa.Position")])

## Load protein data
loadTableData(protein.data.avgexp,data.key.column="Gene.names", "node", table.key.column = 'name')
#setNodeCustomHeatMapChart(c("AvgExpProtein"), style.name=style.name)

## To-Do: consolidate data in one new column to facilitate using the same data viz for both types of data

## Get two columns from phospho data and rename columns
phospho.data.sub2 <- phospho.data %>% 
  mutate(node.name = UniProt_AAPosition) %>%
  mutate(new.avg.exp = AvgExp) %>%
  select(node.name, new.avg.exp)

## Get two columns from protein data and rename columns
protein.data.sub2 <-  protein.data %>%
  mutate(node.name = Gene.names) %>%
  mutate(new.avg.exp = AvgExp) %>%
  select(node.name,new.avg.exp)

## Combine into one
combined.data.column <- rbind(phospho.data.sub2, protein.data.sub2)

## Load into network
loadTableData(combined.data.column,data.key.column="node.name", "node", table.key.column = 'name')

## Apply viz
setNodeColorMapping('new.avg.exp', colors=paletteColorBrewerRdBu, style.name=style.name)
