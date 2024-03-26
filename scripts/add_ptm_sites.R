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
library(RColorBrewer)
cytoscapePing()

###############
## Read in data

## Phosphoprotein data. Note that this includes NA values. 
## Add unique id comp_id for each protein/site combination
cptac.phospho <- read.csv("../datasets/CPTAC_phospho_tn.txt", stringsAsFactors = F, sep = "\t")
cptac.phospho <- cptac.phospho %>% 
  mutate(comp_id=paste0(symbol, "_", site))

## Extracting phospho data for one cancer type - CCRC example
## CCRC phosphoprotein abundance dataset without NAs
cptac.phospho.ccrcc <- cptac.phospho %>% 
  select(symbol, site, protein, comp_id, CCRCC.val, CCRCC.pval) %>%
  filter(!is.na(CCRCC.pval) & !is.na(CCRCC.pval))

cptac.phospho.ccrcc.sig <- cptac.phospho.ccrcc %>% 
  filter(CCRCC.pval > 0 & CCRCC.pval < 0.05)

## PROGENy EGFR pathway correlated data. Not sure we need this.
cptac.progeny.egfr <- read.csv("../datasets/CPTAC_PROGENy__EGFR.txt", stringsAsFactors = F, sep = "\t")
cptac.progeny.egfr <- cptac.progeny.egfr %>% 
  mutate(comp_id=paste0(symbol, "_", site))

## Positively regulated PROGEny EGFR pathway and CCRC
cptac.egfr.ccrcc.pos <- cptac.progeny.egfr %>%
  filter(CCRCC.pval > 0 & CCRCC.pval < 0.05) ##Automatically gets rid of NAs

## Negatively regulated PROGEny EGFR pathway and CCRC
cptac.egfr.ccrcc.neg <- cptac.progeny.egfr %>%
  filter(CCRCC.pval > -0.05 & CCRCC.pval < 0)

## Protein data. Note that this includes NA values.
cptac.protein <- read.csv("../datasets/CPTAC_protein_tn.txt", stringsAsFactors = F, sep = "\t")

## Extracting protein data for one cancer type - CCRC example
## CCRC protein dataset without NAs
cptac.protein.ccrcc <- cptac.protein %>% 
  select(ensembl, symbol, CCRCC.val, CCRCC.pval) %>%
  filter(!is.na(CCRCC.pval) & !is.na(CCRCC.pval))

## Phosphosite kinase-substrate data
psp.data <- read.csv("../annotations/Kinase_Substrate_Dataset.txt", stringsAsFactors = F, sep = "\t")
psp.data.human<- psp.data %>%
  filter(KIN_ORGANISM == "human" & IN_VIVO_RXN == "X")  %>% 
  mutate(comp_id=paste0(SUBSTRATE, "_", SUB_MOD_RSD)) %>%
  select(GENE, comp_id, SUBSTRATE, SUB_MOD_RSD) #human only, in vivo evidence only

###############
## Open and process WP, get relevant phospho data nodes and cross-check against kinase-substrate data
RCy3::commandsRun('wikipathways import-as-pathway id=WP4806') 
## Remove existing ptm states from nodes. These are encoded as nodes with the label P.
selectNodes(nodes = "p", by.col = "name", preserve.current.selection = FALSE) 
deleteSelectedNodes()

mode <- "all"

## Get all relevant protein/gene product data nodes in pathway
node.table <- RCy3::getTableColumns(table = "node")
node.table.prot <- subset(node.table, Type == 'Protein' | Type == 'GeneProduct') %>%
  select(SUID, name, XrefId, Ensembl) ##get node table entries for relevant nodes. not sure we need this

## All matching nodes: Intersection between pathway protein/gene nodes and significant phospho data
## These are the nodes we are going to add phospho nodes to.
matching.nodes.prot <- node.table.prot %>% 
  filter(name %in% cptac.phospho.ccrcc.sig$symbol) %>% 
  select(SUID, name) 

## Need this for matching with kinase substrate data
matching.nodes.phospho <- cptac.phospho.ccrcc.sig %>%
  filter(symbol %in% node.table.prot$name) %>%
  mutate(comp_id=paste0(symbol, "_", site)) %>%
  select(symbol, site, comp_id)

###############
## Only those sites with kinases on pathway: Subset of matching nodes with relevant kinases in the pathway
## To-Do: redo for SUID
kinase.pw <- intersect(psp.data.human$GENE, node.table.prot$name) #Kinases on pathway. GENE column contains kinase.

## Overlap between phospho nodes and kinase-substrate data, where the kinase (GENE) is on the pathway
matching.nodes.phospho.kin <- inner_join(matching.nodes.phospho, psp.data.human) %>%
  filter(GENE %in% kinase.pw) %>%
  select(symbol, comp_id, site)

matching.nodes.phospho.kin.sites <- matching.nodes.phospho.kin %>%
  select(symbol, comp_id) 

mode <- "kinase"
###############
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
  #Get node position for protein node
  p.position <- getNodePosition(node.names = p)
  p.width <- getNodeWidth(node.names = p)
  p.height <- getNodeHeight(node.names = p)
  #Calculate possible phospho site positions, limit to 4 sites for now.
  phospho.positions <- data.frame()
  pos.x <- p.position$x_location 
  pos.y <- p.position$y_location
  #phospho.positions <- data.frame(position=1:4, x=c(pos.x+40, pos.x+40, pos.x-40, pos.x-40), y=c(pos.y+10, pos.y-10, pos.y-10, pos.y+10))
  phospho.positions <- data.frame(position=1:4, x=c(pos.x+(p.width/2)+18, pos.x+(p.width/2)+18, pos.x-+(p.width/2)-18, pos.x-(p.width/2)-18), y=c(pos.y+(p.height/2), pos.y-(p.height/2), pos.y-(p.height/2), pos.y+(p.height/2)))
  
  #find relevant subset of cptac phospho data, these are the ptm nodes to add
  phospho.nodes <- site.table %>% 
    filter(symbol == prot) %>% 
    select(comp_id)
  phospho.nodes <- head(phospho.nodes, 4) ##take first 4
  
  ##Add a node for each ptm, for the particular protein. Collect the SUID for the added ptm, and use it for updating position below.
  suids <- addCyNodes(node.names = phospho.nodes$comp_id, skip.duplicate.names = FALSE)
  ptms <- data.frame()
  # ptms <- do.call(rbind, suids) #doesnt work
  
  for (l in suids){
    ptms <- rbind(ptms, l)
  }
  
  ##Set node position bypass for ptm nodes. Use 4 sites on each node.
  pos <- 1
  #for (ptm in phospho.nodes$comp_id){
  for (ptm in ptms$SUID){
    addCyEdges(c(p, ptm)) #need the SUID for this
    x <- subset(phospho.positions, position == pos)$x
    y <- subset(phospho.positions, position == pos)$y
    print(x)
    print(y)
    setNodePositionBypass(ptm, x, y)
    pos <- pos+1
  }
  ptms.all <- rbind(ptms.all, ptms)
}

##############
## Data viz
style.name = "WikiPathways"

## Load protein data and visualize as node color. Details of this will depend on what the data is, current test data is pval.

loadTableData(cptac.protein, data.key.column="symbol", "node", table.key.column = 'shared name')
RCy3::setNodeColorMapping('CCRCC.pval', colors=paletteColorBrewerReds, style.name = style.name) 

setNodeColorDefault('#FFFFFF', style.name = style.name)
setNodeBorderColorDefault("#737373", style.name = style.name)
 
## Update style for ptm nodes.
setNodeFontSizeBypass(ptms.all$SUID, 9)
setNodeWidthBypass(ptms.all$SUID, 40) 
setNodeHeightBypass(ptms.all$SUID, 20)
#to-do: bring phospho nodes to the front. not sure if this is possible with RCy3

## Define colors based on cutoffs, set as defalut node fill, update phospho node label
ptms.all.data <- inner_join(ptms.all, matching.nodes.phospho, by = join_by(name == comp_id)) %>% 
  mutate(color = case_when(CCRCC.pval < 0.01 ~ "#FF0000", CCRCC.pval > 0.01 ~ "#F496AF"))

ptms.all.data.site <- ptms.all.data %>%
  select(SUID, site) %>%
  rename(name = site)

loadTableData(ptms.all.data.site,data.key.column="SUID", "node", table.key.column = "SUID")
setNodeColorBypass(node.names = ptms.all.data$SUID, new.colors = ptms.all.data$color)



