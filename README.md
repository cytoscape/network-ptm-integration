# network-ptm-integration
This repo tracks the work of integrating networks (pathway models) with ptm data. This includes data-driven visualization of phosphorylation sites on proteins in networks, visualization of phosphoproteomics data on manually curated phosphorylation sites, and utility scripts for related tasks. 

The RShiny tool is designed as an interactive tool where users can visualize data on relevant phospho sites on pathways. Phospho sites can be visualized in multiple ways: an oval PTM node offset from the parent node for each phospho site, or as a pie chart offset from the node visualizing all phospho sites. The user can choose to visualize all sites in the data, or restrict the sites added based on various methods/data (relevant kinases/phosphatase on the pathway, pan-cancer PROGENy data, inhibitor data). 

# How to use add_ptm_sites.R
This app is designed to facilitate exploration of proteomics and phosphoproteomics on phospho sites in pathway models. Phospho sites can be added to the pathways in multiple data-driven modes, or manually annotated phospho sites already in the pathway models can be used.

- Import and preprocess data
- Connect to Cytoscape
- Visualize PTM data, two visualization modes:
   - Traditional PTM: Places up to four PTM nodes around each parent protein node.
   - Pie Chart: Creates one pie chart visualization of multiple phospho sites per parent protein node.

## Prerequisites

- Install Cytoscape. (See https://cytoscape.org/download.html)
- Install R Studio. (See https://posit.co/download/rstudio-desktop/)

## Run App
1. Launch Cytoscape on your local machine.
2. Then launch the Shiny App by running the R script containing the app code in R or RStudio.

## Using custom data
Custom data files must be located in the dirs corresponding to the type of data:
- Phosphoproteomics data: in <code>datasets/phospho</code>
- Proteomics data: in <code>datasets/protein</code>
- Phosphosite data: in <code>datasets/phosphosite</code>

### Formatting requirements for custom data - UNDER CONSTRUCTION
You can use custom data to define the phosphosites to add to the pathway, as well as custom proteomics and phosphoproteomics data. 

Phosphosite and phosphoproteomics data: 
- The data *must* contain a column containing protein identifiers for substrate nodes, meaning the nodes to add phospho sites to
- The data *must* contain a column containing the site specific information in the following format: <code>S16</code> for serine at residue 16.

Proteomics data
- The data *must* contain a column containing protein identifiers
---

# Other scripts

## BasicDataViz_ManualAnnotationPTM.R

Example script showing how to do data visualization of proteomics and phosphoproteomics on WikiPathways with manually annotated PTM sites.

## Extract_WP-PTM_info.R

Extract all manually annotated PTMs in WikiPathways.
