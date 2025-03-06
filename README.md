# network-ptm-integration
This repo tracks the work of phosphoproteome data-driven addition of phosphorylation sites to proteins in networks (and pathways). The idea is to create an interactive tool where users can visualize data on relevant phospho sites on pathways. Phospho sites can be visualized in multiple ways: an oval offset from the parent node for each phospho site, or as a pie chart offset from the node visualizing all phospho sites (idea from CPTAC group). The user could choose to visualize all sites in the data, or restrict the sites added based on various methods/data (relevant kinases/phosphatase on the pathway, pan-cancer PROGENy data, inhibitor data). 
# How to Use add_ptm_sites.R
This app is designed to:

Import and preprocess data: Load phosphoproteomics, proteomics, PROGENy, kinase–substrate, and BioMart mapping files.
 
Connect to Cytoscape: Ensure Cytoscape is running, and use RCy3 to import pathways and update network layouts.
 
Visualize PTM data: Offer two visualization modes:
 
Traditional PTM: Places up to four PTM nodes around each protein node.
 
Pie Chart: Creates a pie chart visualization by linking phospho information with network nodes.
 
## Prerequisites

Install Cytoscape on your PC. (See https://cytoscape.org/download.html)

Install R Studio on your PC. (See https://posit.co/download/rstudio-desktop/)

## Run App

Launch Cytoscape on your local machine.

Then launch the Shiny App:

Run the R script containing the app code in R or RStudio.
 
The app will open in your default web browser.
