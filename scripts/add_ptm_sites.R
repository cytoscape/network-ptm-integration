## Setup
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
BiocManager::install("rWikiPathways")
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

##Phosphoprotein data, add unique id comp_id for each protein/site combination
cptac.data <- read.csv("../datasets/PROGENy__EGFR.txt", stringsAsFactors = F, sep = "\t")
cptac.data.brca.pval <- cptac.data %>% 
  mutate(comp_id=paste0(symbol, "_", site)) %>%
  filter(BRCA.pval > 0 & BRCA.pval < 0.05)

##BRCA dataset without NAs
cptac.data.brca <- cptac.data %>% 
  mutate(comp_id=paste0(symbol, "_", site)) %>%
  select(BRCA.val, BRCA.pval, symbol, site) %>%
  filter(!is.na(BRCA.pval) & !is.na(BRCA.pval))

##Positively regulated with PROGENy
cptac.data.brca.pos <- cptac.data.brca %>%
  filter(BRCA.pval > 0 & BRCA.pval < 0.05)

##Negatively regulated with PROGENy
cptac.data.brca.neg <- cptac.data.brca %>%
  filter(BRCA.pval > -0.05 & BRCA.pval < 0)

##Phosphosite kinase-substrate data
psp.data <- read.csv("../annotations/Kinase_Substrate_Dataset.txt", stringsAsFactors = F, sep = "\t")
psp.data.human<- psp.data %>% filter(KIN_ORGANISM == "human" & IN_VIVO_RXN == "X") #human only, in vivo evidence only

###############
##Open and process WP, get relevant phospho data nodes and cross-check against kinase-substrate data
RCy3::commandsRun('wikipathways import-as-pathway id=WP4806') 
mode <- "all"

##Get all relevant data nodes in pathway
node.table <- getTableColumns()
node.table.prot <- subset(node.table, Type == 'Protein' | Type == 'GeneProduct') %>%
  select(SUID, name, XrefId, Ensembl) ##get node table entries for relevant nodes

##All matching nodes: Intersection between pathway protein/gene nodes and phospho data
matching.nodes.prot <- node.table.prot %>% 
  filter(name %in% cptac.data.brca.pval$symbol) %>% 
  rename(SUID_prot = SUID) %>% 
  rename(symbol = name) %>%
  select(SUID_prot, symbol) #these are the nodes we are going to add sites to

#need this for matching with kinase substrate data
matching.nodes.phospho <- cptac.data.brca.pval %>%
  filter(symbol %in% node.table.prot$name) %>%
  mutate(comp_id=paste0(symbol, "_", site)) %>%
  select(symbol, site, comp_id, BRCA.val, BRCA.pval)
# 
# matching.nodes.final <- merge(matching.nodes.prot, matching.nodes.phospho)

# summarized <- matching.nodes.final %>% group_by(SUID_prot) %>% 
#   summarize(site = paste(sort(unique(site)),collapse=","))

matching.nodes.prot.final <- merge(matching.nodes.prot, summarized)

###############
##Only those with kinases on pathway: Subset of matching nodes with relevant kinases in the pathway
##To-Do: need to refactor for alt strategy
kinase.pw <- intersect(psp.data$GENE, node.table.prot$name) #Read in kinase-substrate data

#create temp data frames for comparisons
psp.data.sub <- psp.data %>% mutate(comp_id=paste0(SUBSTRATE, "_", SUB_MOD_RSD)) %>% select(GENE, SUB_MOD_RSD, comp_id)
cptac.data.sub <- matching.nodes.phospho %>% mutate(comp_id=paste0(symbol, "_", site)) %>% select(symbol, site, comp_id)
psp.cptac.sub <- inner_join(psp.data.sub, cptac.data.sub) #combined
psp.cptac.kin <- psp.cptac.sub[psp.cptac.sub$GENE %in% kinase.pw, ] #subset for only sites with kinases on the pathway

#define phospho site nodes and edges
source.target.kin <- psp.cptac.kin %>% select(symbol, comp_id)
ptm.nodes.kin <- psp.cptac.kin$symbol

#If mode is "kinase only", update inputs
if (mode == "kinase"){
  phospho.nodes <- phospho.nodes.kin
  source.target <- source.target.kin
  ptm.nodes <- ptm.nodes.kin
}
###############
##Adding phospho nodes to pathway

##Alt strategy: Add phospho nodes by looping through protein nodes
#p <- "38082" #testing
for (p in matching.nodes.prot$SUID_prot){
  prot <- matching.nodes.prot$symbol[matching.nodes.prot$SUID_prot == p]
  #get node position for protein node
  p.position <- getNodePosition(node.names = p) #requires node name
  #calculate possible phospho site positions, limit to 4 sites
  phospho.positions <- data.frame()
  pos.x <- p.position$x_location
  pos.y <- p.position$y_location
  phospho.positions <- data.frame(position=1:4, x=c(pos.x+30, pos.x+30, pos.x-30, pos.x-30), y=c(pos.y+10, pos.y-10, pos.y-10, pos.y+10))
  
  #find relevant subset of cptac data
  phospho.nodes <- cptac.data.brca.pval %>% 
    filter(symbol == prot)
  
  #add a node for each
  suids <- addCyNodes(node.names = phospho.nodes$comp_id)
  #suids_2 <- as.data.frame(do.call(rbind, suids)) ##doesn't work
  #addCyEdges() #need the SUID for this
  
  #set node position bypass for phospho nodes
  pos <- 1
  for (ptm in phospho.nodes$comp_id){
    print(pos)
    x <- subset(phospho.positions, position == pos)$x
    y <- subset(phospho.positions, position == pos)$y
    setNodePositionBypass(ptm, x, y)
    pos <- pos+1
  }
}

##############
## Data viz
#Copy current visual style to make new style 
copyVisualStyle("WikiPathways", "WikiPathways PTM")
style.name = "WikiPathways PTM"
setVisualStyle(style.name)

setNodeColorDefault('#FFFFFF', style.name = style.name)
setNodeBorderColorDefault("#737373", style.name = style.name)

## Update label for phospho nodes
matching.nodes.phospho.names <- matching.nodes.phospho %>% 
  select(comp_id, site) %>% 
  rename(name = site) 
loadTableData(matching.nodes.phospho.names, data.key.column = "comp_id", table = "node", table.key.column = "shared name") #read in simple site as name for node label

setNodeFontSizeBypass(matching.nodes.phospho$comp_id, 9)
setNodeWidthBypass(matching.nodes.phospho$site, 40) ##this should use SUID
setNodeHeightBypass(matching.nodes.phospho$site, 20) ##this should use SUID

#To-Do: use cancer-type specific data (pval) to create mapping via comp_id
#Define colors based on cutoffs.
matching.nodes.phospho <- matching.nodes.phospho %>%
  mutate(color = case_when(BRCA.pval < 0.01 ~ "#FF0000", BRCA.pval > 0.01 ~ "#F496AF"))
setNodeColorBypass(node.names = matching.nodes.phospho$site, new.colors = matching.nodes.phospho$color)

## Load protein data
#loadTableData(protein.data.avgexp,data.key.column="Gene.names", "node", table.key.column = 'name')
#setNodeCustomHeatMapChart(c("AvgExpProtein"), style.name=style.name)
