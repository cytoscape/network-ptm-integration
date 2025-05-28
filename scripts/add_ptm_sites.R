if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RCy3")
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
# For the kinase-substrate mapping folder (all .txt files)
kinasemappingChoices <- list_txt_files("../annotations/kinase")
# For the identifier mapping folder (all .txt files)
identifiermappingChoices <- list_txt_files("../annotations/mapping")

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
               selectInput("phosphoFile", "Phosphoproteomics Data File:",
                           choices = datasetsChoices,
                           selected = names(datasetsChoices)[grep("phospho", names(datasetsChoices), ignore.case = TRUE)[1]]),
               selectInput("proteinFile", "Proteomics Data File (datasets):",
                           choices = datasetsChoices,
                           selected = names(datasetsChoices)[grep("protein", names(datasetsChoices), ignore.case = TRUE)[1]]),
               selectInput("progenyFile", "PROGENy Data File:",
                           choices = progenyChoices,
                           selected = names(progenyChoices)[1]),
               selectInput("kinaseFile", "Kinase-Substrate Mapping File:",
                           choices = kinasemappingChoices,
                           selected = names(kinasemappingChoices)[grep("Kinase", names(kinasemappingChoices), ignore.case = TRUE)[1]]),
               selectInput("biomartFile", "Identifier Mapping File:",
                           choices = identifiermappingChoices,
                           selected = names(identifiermappingChoices)[grep("mart", names(identifiermappingChoices), ignore.case = TRUE)[1]]),
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
               selectInput("wpid", "WikiPathway ID:",
                           choices = c(
                             "4-hydroxytamoxifen, dexamethasone, and retinoic acids regulation of p27 expression - WP3879 (Homo sapiens)" = "WP3879",
                             "ATM signaling - WP2516 (Homo sapiens)" = "WP2516",
                             "ATM signaling in development and disease - WP3878 (Homo sapiens)" = "WP3878",
                             "ATR signaling - WP3875 (Homo sapiens)" = "WP3875",
                             "Altered glycosylation of MUC1 in tumor microenvironment - WP4480 (Homo sapiens)" = "WP4480",
                             "Amino acid metabolism in triple-negative breast cancer cells - WP5213 (Homo sapiens)" = "WP5213",
                             "Androgen receptor network in prostate cancer - WP2263 (Homo sapiens)" = "WP2263",
                             "Angiogenesis - WP1539 (Homo sapiens)" = "WP1539",
                             "Apoptosis - WP254 (Homo sapiens)" = "WP254",
                             "Autophagy in pancreatic ductal adenocarcinoma - WP5331 (Homo sapiens)" = "WP5331",
                             "Base excision repair - WP4752 (Homo sapiens)" = "WP4752",
                             "Bladder cancer - WP2828 (Homo sapiens)" = "WP2828",
                             "Breast cancer pathway - WP4262 (Homo sapiens)" = "WP4262",
                             "CCL18 signaling - WP5097 (Homo sapiens)" = "WP5097",
                             "Cancer immunotherapy by CTLA4 blockade - WP4582 (Homo sapiens)" = "WP4582",
                             "Cancer immunotherapy by PD-1 blockade - WP4585 (Homo sapiens)" = "WP4585",
                             "Cancer pathways - WP5434 (Homo sapiens)" = "WP5434",
                             "Cell cycle - WP179 (Homo sapiens)" = "WP179",
                             "Chemokine signaling - WP3929 (Homo sapiens)" = "WP3929",
                             "Chromosomal and microsatellite instability in colorectal cancer - WP4216 (Homo sapiens)" = "WP4216",
                             "Clear cell renal cell carcinoma pathways - WP4018 (Homo sapiens)" = "WP4018",
                             "Cytokines and inflammatory response - WP530 (Homo sapiens)" = "WP530",
                             "Cytosolic DNA-sensing pathway - WP4655 (Homo sapiens)" = "WP4655",
                             "DNA IR-damage and cellular response via ATR - WP4016 (Homo sapiens)" = "WP4016",
                             "DNA IR-double strand breaks and cellular response via ATM - WP3959 (Homo sapiens)" = "WP3959",
                             "DNA damage response - WP707 (Homo sapiens)" = "WP707",
                             "DNA damage response (only ATM dependent) - WP710 (Homo sapiens)" = "WP710",
                             "DNA mismatch repair - WP531 (Homo sapiens)" = "WP531",
                             "Direct reversal repair - WP4931 (Homo sapiens)" = "WP4931",
                             "EGFR tyrosine kinase inhibitor resistance - WP4806 (Homo sapiens)" = "WP4806",
                             "EPO receptor signaling - WP581 (Homo sapiens)" = "WP581",
                             "Endometrial cancer - WP4155 (Homo sapiens)" = "WP4155",
                             "Epithelial to mesenchymal transition in colorectal cancer - WP4239 (Homo sapiens)" = "WP4239",
                             "ErbB signaling - WP673 (Homo sapiens)" = "WP673",
                             "FBXL10 enhancement of MAP/ERK signaling in diffuse large B-cell lymphoma - WP4553 (Homo sapiens)" = "WP4553",
                             "Fatty acid beta-oxidation - WP143 (Homo sapiens)" = "WP143",
                             "Fluoropyrimidine activity - WP1601 (Homo sapiens)" = "WP1601",
                             "Focal adhesion: PI3K-Akt-mTOR-signaling - WP3932 (Homo sapiens)" = "WP3932",
                             "G1 to S cell cycle control - WP45 (Homo sapiens)" = "WP45",
                             "Glycolysis and gluconeogenesis - WP534 (Homo sapiens)" = "WP534",
                             "H19, Rb-E2F1, and CDK-beta-catenin in colorectal cancer - WP3969 (Homo sapiens)" = "WP3969",
                             "Head and neck squamous cell carcinoma - WP4674 (Homo sapiens)" = "WP4674",
                             "Hedgehog signaling - WP4249 (Homo sapiens)" = "WP4249",
                             "Hepatitis B infection - WP4666 (Homo sapiens)" = "WP4666",
                             "Hepatitis C and hepatocellular carcinoma - WP3646 (Homo sapiens)" = "WP3646",
                             "Hepatocyte growth factor receptor signaling - WP313 (Homo sapiens)" = "WP313",
                             "Hereditary leiomyomatosis and renal cell carcinoma pathway - WP4206 (Homo sapiens)" = "WP4206",
                             "Hippo-Merlin signaling dysregulation - WP4541 (Homo sapiens)" = "WP4541",
                             "IL1 signaling - WP195 (Homo sapiens)" = "WP195",
                             "IL6 signaling - WP364 (Homo sapiens)" = "WP364",
                             "Imatinib and chronic myeloid leukemia - WP3640 (Homo sapiens)" = "WP3640",
                             "Immune infiltration in pancreatic cancer - WP5285 (Homo sapiens)" = "WP5285",
                             "Inhibition of exosome biogenesis and secretion by Manumycin A in CRPC cells - WP4301 (Homo sapiens)" = "WP4301",
                             "Integrin-mediated cell adhesion - WP185 (Homo sapiens)" = "WP185",
                             "Interactions between immune cells and microRNAs in tumor microenvironment - WP4559 (Homo sapiens)" = "WP4559",
                             "Interferon type I signaling - WP585 (Homo sapiens)" = "WP585",
                             "Irinotecan pathway - WP229 (Homo sapiens)" = "WP229",
                             "MAPK signaling - WP382 (Homo sapiens)" = "WP382",
                             "MET in type 1 papillary renal cell carcinoma - WP4205 (Homo sapiens)" = "WP4205",
                             "MFAP5 effect on permeability and motility of endothelial cells via cytoskeleton rearrangement - WP4560 (Homo sapiens)" = "WP4560",
                             "MFAP5-mediated ovarian cancer cell motility and invasiveness - WP3301 (Homo sapiens)" = "WP3301",
                             "MSMP expression regulation in cancer cells and its pro-angiogenic role in ovarian tumors - WP4397 (Homo sapiens)" = "WP4397",
                             "Measles virus infection - WP4630 (Homo sapiens)" = "WP4630",
                             "Melanoma - WP4685 (Homo sapiens)" = "WP4685",
                             "Metabolic reprogramming in colon cancer - WP4290 (Homo sapiens)" = "WP4290",
                             "Metabolic reprogramming in pancreatic cancer - WP5220 (Homo sapiens)" = "WP5220",
                             "Methylation pathways - WP704 (Homo sapiens)" = "WP704",
                             "MicroRNA for targeting cancer growth and vascularization in glioblastoma - WP3593 (Homo sapiens)" = "WP3593",
                             "NRF2-ARE regulation - WP4357 (Homo sapiens)" = "WP4357",
                             "NRP1-triggered signaling in pancreatic cancer - WP5144 (Homo sapiens)" = "WP5144",
                             "Neural crest cell migration in cancer - WP4565 (Homo sapiens)" = "WP4565",
                             "Non-homologous end joining - WP438 (Homo sapiens)" = "WP438",
                             "Non-small cell lung cancer - WP4255 (Homo sapiens)" = "WP4255",
                             "Notch signaling - WP268 (Homo sapiens)" = "WP268",
                             "Notch signaling - WP61 (Homo sapiens)" = "WP61",
                             "Nucleotide excision repair - WP4753 (Homo sapiens)" = "WP4753",
                             "One-carbon metabolism - WP241 (Homo sapiens)" = "WP241",
                             "PDGF pathway - WP2526 (Homo sapiens)" = "WP2526",
                             "PDGFR-beta pathway - WP3972 (Homo sapiens)" = "WP3972",
                             "PI3K-AKT-mTOR signaling and therapeutic opportunities in prostate cancer - WP3844 (Homo sapiens)" = "WP3844",
                             "Pancreatic adenocarcinoma pathway - WP4263 (Homo sapiens)" = "WP4263",
                             "Photodynamic therapy-induced HIF-1 survival signaling - WP3614 (Homo sapiens)" = "WP3614",
                             "Photodynamic therapy-induced NF-kB survival signaling - WP3617 (Homo sapiens)" = "WP3617",
                             "Pleural mesothelioma - WP5087 (Homo sapiens)" = "WP5087",
                             "Ras signaling - WP4223 (Homo sapiens)" = "WP4223",
                             "Regulation of Wnt / B-catenin signaling by small molecule compounds - WP3664 (Homo sapiens)" = "WP3664",
                             "Regulation of sister chromatid separation at the metaphase-anaphase transition - WP4240 (Homo sapiens)" = "WP4240",
                             "Retinoblastoma gene in cancer - WP2446 (Homo sapiens)" = "WP2446",
                             "Small cell lung cancer - WP4658 (Homo sapiens)" = "WP4658",
                             "TCA cycle nutrient use and invasiveness of ovarian cancer - WP2868 (Homo sapiens)" = "WP2868",
                             "TGF-beta receptor signaling - WP560 (Homo sapiens)" = "WP560",
                             "TGF-beta signaling in thyroid cells for epithelial-mesenchymal transition - WP3859 (Homo sapiens)" = "WP3859",
                             "TGF-beta signaling pathway - WP366 (Homo sapiens)" = "WP366",
                             "TP53 network - WP1742 (Homo sapiens)" = "WP1742",
                             "Target of rapamycin signaling - WP1471 (Homo sapiens)" = "WP1471",
                             "Transcriptional activation by NRF2 in response to phytochemicals - WP3 (Homo sapiens)" = "WP3",
                             "Translation inhibitors in chronically activated PDGFRA cells - WP4566 (Homo sapiens)" = "WP4566",
                             "Tumor suppressor activity of SMARCB1 - WP4204 (Homo sapiens)" = "WP4204",
                             "Type 2 papillary renal cell carcinoma - WP4241 (Homo sapiens)" = "WP4241",
                             "Type II interferon signaling - WP619 (Homo sapiens)" = "WP619",
                             "Ultraconserved region 339 modulation of tumor suppressor microRNAs in cancer - WP4284 (Homo sapiens)" = "WP4284",
                             "VEGFA-VEGFR2 signaling - WP3888 (Homo sapiens)" = "WP3888",
                             "Warburg effect modulated by deubiquitinating enzymes and their substrates - WP5216 (Homo sapiens)" = "WP5216",
                             "Wnt signaling - WP363 (Homo sapiens)" = "WP363",
                             "Wnt signaling - WP428 (Homo sapiens)" = "WP428",
                             "Wnt/beta-catenin signaling in leukemia - WP3658 (Homo sapiens)" = "WP3658",
                             "lncRNA in canonical Wnt signaling and colorectal cancer - WP4258 (Homo sapiens)" = "WP4258",
                             "lncRNA-mediated mechanisms of therapeutic resistance - WP3672 (Homo sapiens)" = "WP3672",
                             "miRNA regulation of p53 pathway in prostate cancer - WP3982 (Homo sapiens)" = "WP3982",
                             "miRNA regulation of prostate cancer signaling - WP3981 (Homo sapiens)" = "WP3981",
                             "ncRNAs in Wnt signaling in hepatocellular carcinoma - WP4336 (Homo sapiens)" = "WP4336",
                             "ncRNAs involved in STAT3 signaling in hepatocellular carcinoma - WP4337 (Homo sapiens)" = "WP4337"
                           ),
                           selected = "WP4806"),
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
      
      ##Get the new node table to get the pie chart nodes
      final.node.table <- RCy3::getTableColumns(table = "node")
      ptm.nodes <- final.node.table[grepl("_", final.node.table$name), ]$SUID
      
      style.name <- "WikiPathways"
      setNodeColorDefault('#FFFFFF', style.name = style.name)
      setNodeBorderColorDefault("#737373", style.name = style.name)
      ##The next line will remove the main node label on pie chart nodes, leaving the omicsvisualizer label for each site
      ##setNodeLabelBypass(ptm.nodes, '')
      setNodeFontSizeBypass(ptm.nodes, 9)
      
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