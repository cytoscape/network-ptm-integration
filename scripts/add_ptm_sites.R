# Load required libraries
library(shiny)
library(shinythemes)
library(RCy3)
library(rWikiPathways)
library(dplyr)
library(purrr)
# Uncomment additional libraries if needed:
# library(stringr)
# library(RColorBrewer)

# Ensure Cytoscape is running and connected
cytoscapePing()
installApp("Omics Visualizer")

## Define Cytoscape version variable dynamically
cytoscapeVersion <- cytoscapeVersionInfo()[2]
cytoscapeSampleDataPath <- paste0("/Applications/Cytoscape_v", cytoscapeVersion, "/sampleData/")

## Helper function: List all .txt files in a given folder.
list_txt_files <- function(path) {
  files <- list.files(path = path, pattern = ".*\\.txt$", full.names = TRUE)
  # Use the base file names as the display names:
  setNames(files, basename(files))
}

## Dynamically generate file choices:
# For the datasets folder (all .txt files)
datasetsChoices <- list_txt_files("../datasets")
# For the PROGENy subfolder (all .txt files)
progenyChoices  <- list_txt_files("../datasets/PROGENy")
# For the annotations folder (all .txt files)
annotationsChoices <- list_txt_files("../annotations")

## -----------------------
## Shiny UI
## -----------------------
ui <- navbarPage(
  "Network PTM Integration",
  id = "mainTabs",                      # ID to allow tab switching
  theme = shinytheme("united"),         # Orange-inspired theme
  
  ## Tab 1: Data File Selection
  tabPanel("Data Files",
           sidebarLayout(
             sidebarPanel(
               h4("Select Data Files"),
               selectInput("phosphoFile", "Phospho File (datasets):",
                           choices = datasetsChoices,
                           selected = names(datasetsChoices)[grep("phospho", names(datasetsChoices), ignore.case = TRUE)[1]]),
               selectInput("proteinFile", "Protein File (datasets):",
                           choices = datasetsChoices,
                           selected = names(datasetsChoices)[grep("protein", names(datasetsChoices), ignore.case = TRUE)[1]]),
               selectInput("progenyFile", "PROGENy File (datasets/PROGENy):",
                           choices = progenyChoices,
                           selected = names(progenyChoices)[1]),
               selectInput("kinaseFile", "Kinase-Substrate File (annotations):",
                           choices = annotationsChoices,
                           selected = names(annotationsChoices)[grep("Kinase", names(annotationsChoices), ignore.case = TRUE)[1]]),
               selectInput("biomartFile", "BioMart Mapping File (annotations):",
                           choices = annotationsChoices,
                           selected = names(annotationsChoices)[grep("mart", names(annotationsChoices), ignore.case = TRUE)[1]]),
               br(),
               actionButton("goToAnalysis", "Proceed to Analysis")
             ),
             mainPanel(
               h4("Current File Selections"),
               verbatimTextOutput("selectedFiles")
             )
           )
  ),
  
  ## Tab 2: Analysis Controls & Output
  tabPanel("Analysis",
           sidebarLayout(
             sidebarPanel(
               h4("Analysis Parameters"),
               textInput("wpid", "WikiPathway ID:", value = "WP4806"),
               selectInput("vizMode", "Visualization Mode:",
                           choices = c("Traditional PTM", "Pie Chart")),
               selectInput("analysisMode", "Analysis Mode:",
                           choices = c("all", "kinase")),
               br(),
               actionButton("run", "Run Analysis")
             ),
             mainPanel(
               h4("Status"),
               verbatimTextOutput("status")
             )
           )
  )
)

## -----------------------
## Shiny Server
## -----------------------
server <- function(input, output, session) {
  
  ## Display selected file names in the "Data Files" tab.
  output$selectedFiles <- renderPrint({
    list(
      Phospho = input$phosphoFile,
      Protein = input$proteinFile,
      PROGENy = input$progenyFile,
      Kinase_Substrate = input$kinaseFile,
      BioMart = input$biomartFile
    )
  })
  
  ## Automatically switch to Analysis tab when "Proceed to Analysis" is clicked.
  observeEvent(input$goToAnalysis, {
    updateNavbarPage(session, "mainTabs", selected = "Analysis")
  })
  
  observeEvent(input$run, {
    req(input$wpid, input$vizMode, input$analysisMode,
        input$phosphoFile, input$proteinFile, input$progenyFile,
        input$kinaseFile, input$biomartFile)
    
    # Retrieve user selections
    wpid         <- input$wpid
    vizMode      <- input$vizMode        # "Traditional PTM" or "Pie Chart"
    analysisMode <- input$analysisMode     # "all" or "kinase"
    
    output$status <- renderText({
      paste("Starting analysis with WikiPathway ID:", wpid, 
            "\nVisualization Mode:", vizMode, 
            "\nAnalysis Mode:", analysisMode)
    })
    
    ## -----------------------
    ## Import Files Based on User Selection
    ## -----------------------
    cptac.phospho <- read.csv(input$phosphoFile, stringsAsFactors = FALSE, sep = "\t") %>% 
      mutate(prot_site = paste0(protein, "_", site))
    cptac.protein <- read.csv(input$proteinFile, stringsAsFactors = FALSE, sep = "\t")
    cptac.progeny.egfr <- read.csv(input$progenyFile, stringsAsFactors = FALSE, sep = "\t")
    psp.data <- read.csv(input$kinaseFile, stringsAsFactors = FALSE, sep = "\t")
    biomart <- read.csv(input$biomartFile, stringsAsFactors = FALSE, sep = ",") %>%
      rename(Ensembl = Gene.stable.ID, EnsemblProt = Protein.stable.ID) %>%
      select(Ensembl, EnsemblProt)
    
    ## -----------------------
    ## Preprocess Data from Datasets
    ## -----------------------
    cptac.progeny.egfr.ccrcc.pos <- cptac.progeny.egfr %>%
      filter(CCRCC.pval > 0 & CCRCC.pval < 0.05) %>%
      mutate(prot_site = paste0(protein, "_", site)) %>% 
      select(symbol, protein, site, prot_site, CCRCC.val, CCRCC.pval)
    cptac.phospho.ccrcc <- cptac.phospho %>% 
      select(symbol, site, protein, prot_site, CCRCC.val, CCRCC.pval)
    cptac.protein.ccrcc <- cptac.protein %>% 
      select(ensembl, symbol, CCRCC.val, CCRCC.pval) %>%
      filter(!is.na(CCRCC.pval) & !is.na(CCRCC.pval))
    
    ## -----------------------
    ## Import Pathway into Cytoscape
    ## -----------------------
    import_cmd <- paste0('wikipathways import-as-pathway id=', wpid)
    RCy3::commandsRun(import_cmd)
    
    # Remove existing PTM nodes (nodes labeled "p")
    selectNodes(nodes = "p", by.col = "name", preserve.current.selection = FALSE) 
    deleteSelectedNodes()
    
    ## -----------------------------
    ## Common Steps: Node & Data Setup
    ## -----------------------------
    node.table <- RCy3::getTableColumns(table = "node")
    node.table.prot <- node.table %>% 
      filter(Type %in% c("Protein", "GeneProduct")) %>%
      select(SUID, name, XrefId, Ensembl)
    node.table.prot.mapped <- merge(node.table.prot, biomart, by = "Ensembl") %>%
      filter(EnsemblProt != "")
    matching.nodes.prot <- node.table.prot.mapped %>% 
      filter(EnsemblProt %in% cptac.progeny.egfr.ccrcc.pos$protein) %>% 
      select(SUID, name, EnsemblProt)
    node.positions <- getNodePosition(node.names = matching.nodes.prot$SUID)
    node.positions <- cbind(suid = rownames(node.positions), node.positions)
    rownames(node.positions) <- 1:nrow(node.positions)
    node.width <- getNodeWidth(node.names = matching.nodes.prot$SUID)
    node.width <- as.data.frame(node.width)
    node.width <- cbind(suid = rownames(node.width), node.width)
    rownames(node.width) <- 1:nrow(node.width)
    node.height <- getNodeHeight(node.names = matching.nodes.prot$SUID)
    node.height <- as.data.frame(node.height)
    node.height <- cbind(suid = rownames(node.height), node.height)
    rownames(node.height) <- 1:nrow(node.height)
    node.layout <- merge(node.width, node.height, by = "suid")
    node.layout <- merge(node.layout, node.positions, by = "suid")
    node.layout.pie <- node.layout  # For alternative (pie chart) viz
    matching.nodes.phospho <- cptac.progeny.egfr.ccrcc.pos %>%
      filter(protein %in% node.table.prot.mapped$EnsemblProt) %>%
      mutate(prot_site = paste0(protein, "_", site)) %>%
      select(symbol, protein, site, prot_site)
    
    if (analysisMode == "kinase") {
      psp.data.human <- psp.data %>%
        filter(KIN_ORGANISM == "human" & IN_VIVO_RXN == "X") %>% 
        mutate(prot_site = paste0(SUBSTRATE, "_", SUB_MOD_RSD)) %>% 
        select(GENE, prot_site, SUBSTRATE, SUB_MOD_RSD)
      kinase.pw <- intersect(psp.data.human$GENE, node.table.prot$name)
      matching.nodes.phospho <- inner_join(matching.nodes.phospho, psp.data.human, 
                                           by = c("symbol" = "GENE")) %>%
        filter(symbol %in% kinase.pw) %>%
        select(symbol, prot_site, site)
    }
    
    ## -----------------------------
    ## Visualization Branching
    ## -----------------------------
    if (vizMode == "Traditional PTM") {
      node.layout <- node.layout %>%
        mutate(x.ptm = pmap(list(x_location, node.width), 
                            function(x_val, width_val) {
                              list(x_val + (width_val / 2) + 16,
                                   x_val + (width_val / 2) + 16,
                                   x_val - (width_val / 2) - 16,
                                   x_val - (width_val / 2) - 16)
                            }))
      node.layout <- node.layout %>%
        mutate(y.ptm = pmap(list(y_location, node.height), 
                            function(y_val, height_val) {
                              list(y_val + (height_val / 2),
                                   y_val - (height_val / 2),
                                   y_val - (height_val / 2),
                                   y_val + (height_val / 2))
                            }))
      site.table <- matching.nodes.phospho
      ptms.all <- data.frame()
      for (p in matching.nodes.prot$SUID) {
        prot <- matching.nodes.prot$EnsemblProt[matching.nodes.prot$SUID == p]
        message("Processing protein (Traditional): ", prot)
        phospho.nodes <- site.table %>% 
          filter(protein == prot) %>% 
          select(prot_site)
        phospho.nodes <- head(phospho.nodes, 4)
        suids <- addCyNodes(node.names = phospho.nodes$prot_site, skip.duplicate.names = FALSE)
        ptms <- data.frame()
        for (l in suids) {
          ptms <- rbind(ptms, l)
        }
        pos <- 1
        for (ptm in ptms$SUID) {
          x <- node.layout$x.ptm[node.layout$suid == p][[1]][[pos]]
          y <- node.layout$y.ptm[node.layout$suid == p][[1]][[pos]]
          setNodePositionBypass(ptm, x, y)
          pos <- pos + 1
        }
        ptms.all <- rbind(ptms.all, ptms)
      }
      style.name <- "WikiPathways"
      setNodeColorDefault('#FFFFFF', style.name = style.name)
      setNodeBorderColorDefault("#737373", style.name = style.name)
      setNodeFontSizeBypass(ptms.all$SUID, 9)
      setNodeWidthBypass(ptms.all$SUID, 35)
      setNodeHeightBypass(ptms.all$SUID, 20)
      loadTableData(cptac.protein.ccrcc, data.key.column = "ensembl", 
                    table = "node", table.key.column = 'Ensembl')
      loadTableData(cptac.phospho.ccrcc, data.key.column = "prot_site", 
                    table = "node", table.key.column = 'shared name')
      RCy3::setNodeColorMapping('CCRCC.val', colors = paletteColorBrewerRdBu, 
                                style.name = style.name)
      ptms.all.data <- inner_join(ptms.all, cptac.phospho, by = join_by(name == prot_site))
      ptms.all.data.site <- ptms.all.data %>%
        select(SUID, site) %>%
        rename(name = site)
      loadTableData(ptms.all.data.site, data.key.column = "SUID", 
                    table = "node", table.key.column = "SUID")
      setNodeLabelMapping('name', style.name = style.name)
      output$status <- renderText({
        paste("Traditional PTM visualization complete for WikiPathway", wpid)
      })
    } else if (vizMode == "Pie Chart") {
      node.layout.pie <- node.layout.pie %>%
        mutate(x.ptm = (x_location + node.width / 2 + 20)) %>%
        mutate(y.ptm = y_location)
      cptac.phospho.ccrcc.full.ov <- cptac.phospho.ccrcc %>%
        filter(prot_site %in% cptac.progeny.egfr.ccrcc.pos$prot_site) %>%
        mutate(symbol_ptm = paste0(symbol, "_ptm"))
      write.table(cptac.phospho.ccrcc.full.ov, 
                  paste0(cytoscapeSampleDataPath, "cptac.phospho.ccrcc.full.txt"),
                  sep = "\t", row.names = FALSE)
      matching.nodes.prot.pie <- node.table.prot.mapped %>% 
        filter(EnsemblProt %in% cptac.progeny.egfr.ccrcc.pos$protein) %>% 
        mutate(name = paste0(name, "_ptm")) %>% 
        select(SUID, name)
      for (p in matching.nodes.prot.pie$SUID) {
        ptm.name <- matching.nodes.prot.pie$name[matching.nodes.prot.pie$SUID == p]
        suid.list <- addCyNodes(node.names = ptm.name, skip.duplicate.names = FALSE)
        suid.ptm <- suid.list[[1]]$SUID
        x <- node.layout.pie$x.ptm[node.layout.pie$suid == p]
        y <- node.layout.pie$y.ptm[node.layout.pie$suid == p]
        setNodePositionBypass(node.names = suid.ptm, x, y)
        setNodeWidthBypass(node.names = suid.ptm, 40)
        setNodeHeightBypass(node.names = suid.ptm, 40)
      }
      ovload.cmd <- paste('ov load',
                          'dataTypeList="string,string,string,string,double,double,string"', 
                          paste0('file="', cytoscapeSampleDataPath, "cptac.phospho.ccrcc.full.txt\""),
                          'newTableName="cptac.phospho.ccrcc.full"', 
                          'startLoadRow="2"')
      commandsRun(ovload.cmd)
      ovconnect.cmd <- paste('ov connect',
                             'mappingColNet="shared name"', 
                             'mappingColTable="symbol_ptm"')
      commandsRun(ovconnect.cmd)
      ovviz.cmd <- paste('ov viz apply inner continuous',
                         'attributes="CCRCC.val"', 
                         'paletteName="Red-Blue"', 
                         'labels="site"')
      commandsRun(ovviz.cmd)
      style.name <- "WikiPathways"
      setNodeColorDefault('#FFFFFF', style.name = style.name)
      setNodeBorderColorDefault("#737373", style.name = style.name)
      setNodeFontSizeBypass(matching.nodes.prot.pie$SUID, 9)
      loadTableData(cptac.protein.ccrcc, data.key.column = "ensembl", 
                    table = "node", table.key.column = 'Ensembl')
      RCy3::setNodeColorMapping('CCRCC.val', colors = paletteColorBrewerRdBu, 
                                style.name = style.name)
      output$status <- renderText({
        paste("Pie Chart visualization complete for WikiPathway", wpid)
      })
    }
    
  })
}

## Run the Shiny App
shinyApp(ui, server)