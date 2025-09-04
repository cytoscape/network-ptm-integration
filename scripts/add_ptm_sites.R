if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RCy3")
# Load required libraries
library(shiny)
library(shinyjs)
library(shinythemes)
library(RCy3)
library(dplyr)
library(purrr)
# Uncomment additional libraries if needed:
library(stringr)
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
phosphodatasetsChoices <- list_txt_files("../datasets/phospho")
proteindatasetsChoices <- list_txt_files("../datasets/protein")
# For the PROGENy subfolder (all .txt files)
progenyChoices  <- list_txt_files("../datasets/PROGENy")
# For the CPTAC cancer type
cptacChoices <- c("CCRCC", "COAD", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")
# For the kinase-substrate mapping folder (all .txt files)
#kinasemappingChoices <- list_txt_files("../annotations/kinase-substrate")
# For the identifier mapping folder (all .txt files)
identifiermappingChoices <- list_txt_files("../annotations/id-mapping")

## -----------------------
## Shiny UI
## -----------------------
ui <- navbarPage(
  "Phosphosite - Pathway Integration",
  id = "mainTabs",                      # ID to allow tab switching
  theme = shinytheme("cerulean"),         # Orange-inspired theme
  
  ## Tab 1: Data File Selection
  tabPanel("Data Files",
           useShinyjs(),
           wellPanel(
             h3("Phosphosite - Pathway Integration"),
             p("This tool allows you to add phosphorylation sites to proteins in WikiPathways models in Cytoscape, and to visualize data on ptms and parent nodes."),
             p("Phosphosites can be added from data (phosphoproteomics data, PROGENy data), or manually curated phosphosites already present in WikiPathways models can be used directly for data visualization. 
             The tool includes pre-loaded example data (CPTAC pan-cancer phosphoproteomics and proteomics data, PROGENy pathway activity data).")
           ),
           sidebarLayout(
             sidebarPanel(
               h4("Select Method for Adding Phospho Sites"),
               selectInput("phosphoMode", "Phospho Site Mode:",
                           choices = c("Data-driven PROGENy", "Data-driven phosphoproteomics data", "Manually curated")),
               h4("Select Data Files"),
               selectInput("phosphoFile", "Phosphoproteomics Data File:",
                           choices = phosphodatasetsChoices,
                           selected = names(phosphodatasetsChoices)[grep("phospho", names(phosphodatasetsChoices), ignore.case = TRUE)[1]]),
               selectInput("proteinFile", "Proteomics Data File:",
                           choices = proteindatasetsChoices,
                           selected = names(proteindatasetsChoices)[grep("protein", names(proteindatasetsChoices), ignore.case = TRUE)[1]]),
               selectInput("progenyFile", "PROGENy Data File:",
                           choices = progenyChoices,
                           selected = names(progenyChoices)[1]),
               selectInput("cptacType", "CPTAC Cancer Type:",
                           choices = cptacChoices,
                           selected = names(cptacChoices)[1]),
               # selectInput("kinaseFile", "Kinase-Substrate Mapping File:",
               #             choices = kinasemappingChoices,
               #             selected = names(kinasemappingChoices)[grep("Kinase", names(kinasemappingChoices), ignore.case = TRUE)[1]]),
               # selectInput("biomartFile", "Identifier Mapping File:",
               #             choices = identifiermappingChoices,
               #             selected = names(identifiermappingChoices)[grep("mart", names(identifiermappingChoices), ignore.case = TRUE)[1]]),
               br(),
               actionButton("goToAnalysis", "Proceed to Analysis")
             ),
             mainPanel(
               h4("Current Selections"),
               uiOutput("selectedFiles")
             )
           )
  ),
  
  ## Tab 2: Analysis Controls & Output
  tabPanel("Analysis",
           useShinyjs(),
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
               # selectInput("analysisMode", "Analysis Mode:",
               #             choices = c("all", "kinase")),
               br(),
               actionButton("run", "Run Analysis"),
               textOutput("status")
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
  
  output$selectedFiles <- renderUI({
    tags$ul(style = "list-style-type: none; padding-left: 0;",
            tags$li(style="margin-bottom: 8px;", tags$b("Mode:"), input$phosphoMode),
            tags$li(style="margin-bottom: 8px;", tags$b("Phospho:"), input$phosphoFile),
            tags$li(style="margin-bottom: 8px;", tags$b("Protein:"), input$proteinFile),
            tags$li(style="margin-bottom: 8px;", tags$b("PROGENy:"), input$progenyFile),
            tags$li(style="margin-bottom: 8px;", tags$b("CPTAC Cancer Type:"), input$cptacType)
            # tags$li(style="margin-bottom: 8px;", tags$b("Kinase Substrate:"), input$kinaseFile),
            # tags$li(style="margin-bottom: 8px;", tags$b("BioMart:"), input$biomartFile)
    )
  })
  
  ## Automatically switch to Analysis tab when "Proceed to Analysis" is clicked.
  observeEvent(input$goToAnalysis, {
    updateNavbarPage(session, "mainTabs", selected = "Analysis")
  })
  
  observeEvent(input$run, {
    req(input$wpid, input$vizMode, input$phosphoMode,
        input$phosphoFile, input$proteinFile, input$progenyFile,
        input$cptacType)
    
    # Retrieve user selections
    wpid         <- input$wpid
    vizMode      <- input$vizMode        # "Traditional PTM" or "Pie Chart"
    #analysisMode <- input$analysisMode     # "all" or "kinase"
    mode <- input$phosphoMode  #which phospho mode
    
    #hard-coded mode
    analysisMode <- "all"
    
    # #This doesn't work, it is displayed AFTER completed analysis
    # output$status <- renderText({
    #   paste("Starting analysis with WikiPathway ID:", wpid,
    #         "\nPhospho Mode:", mode,
    #         "\nVisualization Mode:", vizMode,
    #         "\nAnalysis Mode:", analysisMode)
    # })

    ## -----------------------
    ## Import Files Based on User Selection
    ## -----------------------
    
    ##Import id mapping file from Biomart and kinase-substrate mapping from PSP
    #psp.data <- read.csv("../annotations/kinase-substrate/PSP_Kinase_Substrate_Dataset.txt", stringsAsFactors = FALSE, sep = "\t")
    biomart <- read.csv("../annotations/id-mapping/ensembl_mappings.txt", stringsAsFactors = FALSE, sep = "\t")

    cptac.phospho <- read.csv(input$phosphoFile, stringsAsFactors = FALSE, sep = "\t") %>% 
      mutate(prot_site = paste0(protein, "_", site))
    
    cptac.protein <- read.csv(input$proteinFile, stringsAsFactors = FALSE, sep = "\t")
    cptac.progeny <- read.csv(input$progenyFile, stringsAsFactors = FALSE, sep = "\t")
    # psp.data <- read.csv(input$kinaseFile, stringsAsFactors = FALSE, sep = "\t")
    # biomart <- read.csv(input$biomartFile, stringsAsFactors = FALSE, sep = "\t")
    
    ## -----------------------
    ## Preprocess Data from Datasets
    ## -----------------------
    
    type.pval <- paste0(input$cptacType, ".pval") ##get relevant column header based on cancer type selected
    type.val <- paste0(input$cptacType, ".val")
  
    if (mode == "Data-driven PROGENy"){
      ##Get relevant sites based on positive PROGENy scores
      cptac.progeny.pos <- cptac.progeny %>%
        filter(.data[[type.pval]] > 0 & .data[[type.pval]] < 0.05) %>%
        mutate(prot_site = paste0(protein, "_", site)) %>%
        dplyr::select(symbol, protein, site, prot_site, all_of(type.val), all_of(type.pval))
    }
    
    cptac.phospho <- cptac.phospho %>% 
      dplyr::select(symbol, site, protein, prot_site, all_of(type.val), all_of(type.pval))
    cptac.protein <- cptac.protein %>% 
      dplyr::select(ensembl, symbol, all_of(type.val), all_of(type.pval)) %>%
      filter(!is.na(.data[[type.pval]]) & !is.na(.data[[type.pval]]))
    
    if (mode == "Data-driven phosphoproteomics data"){
      cptac.phospho.threshold <- cptac.phospho %>%
      filter(.data[[type.pval]] > 0 & .data[[type.pval]] < 0.05)
    }
    
    ## -----------------------
    ## Import Pathway into Cytoscape
    ## -----------------------
    import_cmd <- paste0('wikipathways import-as-pathway id=', wpid)
    RCy3::commandsRun(import_cmd)
    
    ## Get nodes in the pathway, specifically protein nodes, and map them to Ensembl protein ids
    node.table <- RCy3::getTableColumns(table = "node")
    
    if (mode %in% c("Data-driven PROGENy", "Data-driven phosphoproteomics data")) {
      # Remove existing PTM nodes (nodes labeled "p").
      selectNodes(nodes = "p", by.col = "name", preserve.current.selection = FALSE) 
      deleteSelectedNodes()
    
    ## -----------------------------
    ## Common Steps: Node & Data Setup
    ## -----------------------------
      
    node.table.prot <- node.table %>% 
      filter(Type %in% c("Protein", "GeneProduct")) %>%
        dplyr::select(SUID, name, XrefId, Ensembl)
    
    node.table.prot.mapped <- merge(node.table.prot, biomart, by.x = "Ensembl", by.y = "ensembl_gene_id") %>%
      filter(ensembl_peptide_id != "") %>%
      filter(uniprotswissprot != "")
    
    if (mode == "Data-driven PROGENy"){
      print(mode)
    ## Get the relevant phospho sites from PROGENy that match the proteins in the node table
    matching.nodes.phospho <- cptac.progeny.pos %>%
      filter(protein %in% node.table.prot.mapped$ensembl_peptide_id) %>%
      dplyr::select(symbol, protein, site, prot_site)
    }
    
    if (mode == "Data-driven phosphoproteomics data"){
      print(mode)
      matching.nodes.phospho <- cptac.phospho.threshold %>%
        filter(protein %in% node.table.prot.mapped$ensembl_peptide_id) %>%
        dplyr::select(symbol, protein, site, prot_site)
      }
    
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
    
    ## For kinase mode, filter for only those phospho sites that have kinases on the pathway
    if (analysisMode == "kinase") {
      psp.data.human <- psp.data %>%
        filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human") %>% 
        dplyr::select(GENE, KINASE, KIN_ACC_ID, SUB_ACC_ID, SUB_MOD_RSD)
      
      ## Data mapping PSP data
      ## Add Ensembl prot id for substrate
      psp.data.human.mapped <- merge(psp.data.human, biomart, by.x = "SUB_ACC_ID", by.y = "uniprotswissprot") %>%
        mutate(substrate_ensembl_gene_id = ensembl_gene_id) %>%
        mutate(substrate_peptide_id = ensembl_peptide_id) %>%
        mutate(prot_site = paste0(substrate_peptide_id, "_", SUB_MOD_RSD)) %>%
        dplyr::select(GENE, KINASE, KIN_ACC_ID, substrate_ensembl_gene_id, SUB_ACC_ID, substrate_peptide_id, prot_site, SUB_MOD_RSD)
      
      ## Add Ensembl gene id for kinase
      psp.data.human.mapped <- merge(psp.data.human.mapped, biomart, by.x = "KIN_ACC_ID", by.y = "uniprotswissprot") %>%
        mutate(kin_ensembl_gene_id = ensembl_gene_id)
      
      ## Get filtered list of phospho sites by filtering for those with kinases on the pw.
      filtered.ptms <- psp.data.human.mapped %>%
        filter(prot_site %in% matching.nodes.phospho$prot_site) %>%
        filter(kin_ensembl_gene_id %in% node.table.prot$Ensembl) %>%
        dplyr::select(prot_site)
      
      filtered.ptms <- unique(filtered.ptms)
    
      matching.nodes.phospho <- matching.nodes.phospho %>%
        filter(prot_site %in% filtered.ptms$prot_site)
    }
    
    ## -----------------------------
    ## Visualization Branching
    ## -----------------------------
    
    ## Get protein nodes in the node table that correspond to the phospho sites. 
    ## We can then use this list of existing nodes to add new nodes
    matching.nodes.prot <- node.table.prot.mapped %>% 
      filter(ensembl_peptide_id %in% matching.nodes.phospho$protein) %>%
      dplyr::select(SUID, name, ensembl_peptide_id)
    
    message("Matching protein nodes: ", nrow(matching.nodes.prot))
    message("Matching phospho nodes: ", nrow(matching.nodes.phospho))
    
    if (nrow(matching.nodes.phospho) > 0){
      
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
        prot <- matching.nodes.prot$ensembl_peptide_id[matching.nodes.prot$SUID == p]
        message("Processing protein (Traditional): ", prot)
        phospho.nodes <- site.table %>% 
          filter(protein == prot) %>% 
          dplyr::select(prot_site)
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
      loadTableData(cptac.protein, data.key.column = "ensembl", 
                    table = "node", table.key.column = 'Ensembl')
      loadTableData(cptac.phospho, data.key.column = "prot_site", 
                    table = "node", table.key.column = 'shared name')
      RCy3::setNodeColorMapping(table.column = type.val, colors = paletteColorBrewerRdBu, 
                                style.name = style.name)
      ptms.all.data <- inner_join(ptms.all, cptac.phospho, by = join_by(name == prot_site))
      ptms.all.data.site <- ptms.all.data %>%
        dplyr::select(SUID, site) %>%
        rename(name = site)
      loadTableData(ptms.all.data.site, data.key.column = "SUID", 
                    table = "node", table.key.column = "SUID")
      setNodeLabelMapping('name', style.name = style.name)
      output$status <- renderText({
        paste("Traditional PTM visualization complete for WikiPathway", wpid)
      })
     } 
    else if (vizMode == "Pie Chart") {
      node.layout.pie <- node.layout.pie %>%
        mutate(x.ptm = (x_location + node.width / 2 + 20)) %>%
        mutate(y.ptm = y_location)
      
      cptac.phospho.full.ov <- cptac.phospho %>%
        filter(prot_site %in% matching.nodes.phospho$prot_site) %>%
        mutate(symbol_ptm = paste0(symbol, "_ptm"))
      
      write.table(cptac.phospho.full.ov,
                  paste0(cytoscapeSampleDataPath, "cptac.phospho.full.txt"),
                  sep = "\t", row.names = FALSE)

      matching.nodes.prot.pie <- matching.nodes.prot %>%
        mutate(name = paste0(name, "_ptm")) %>%
        dplyr::select(SUID, name)

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
                          paste0('file="', cytoscapeSampleDataPath, "cptac.phospho.full.txt\""),
                          'newTableName="cptac.phospho.full"',
                          'startLoadRow="2"')
      commandsRun(ovload.cmd)
      ovconnect.cmd <- paste('ov connect',
                             'mappingColNet="shared name"',
                             'mappingColTable="symbol_ptm"')
      commandsRun(ovconnect.cmd)
      ovviz.cmd <- paste('ov viz apply inner continuous',
                         paste0('attributes="', type.val, '"'),
                         'paletteName="Red-Blue"','rangeMax="1.3"', 'rangeMin="-1.3"',
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

      loadTableData(cptac.protein, data.key.column = "ensembl",
                    table = "node", table.key.column = 'Ensembl')
      RCy3::setNodeColorMapping(table.column = type.val, colors = paletteColorBrewerRdBu,
                                style.name = style.name)
      output$status <- renderText({
        paste("Pie Chart visualization complete for WikiPathway", wpid)
      })
    }
    } ##end if non-zero ptms
    else {
      output$status <- renderText({
        "There are no matching phospho sites."
      })
    }
    } ## end if data-driven
    
    else if ((mode == "Manually curated")){
      print(mode)
      cptac.phospho <- read.csv("../datasets/phospho/CPTAC_phospho_tn.txt", stringsAsFactors = F, sep = "\t")
    
      ## Add a new column, data_mapping_id, for data mapping. A new data frame is created with the new column and then read into the Cytoscape node table.
      ## The new data mapping column will contain a new id which is a composite of the parent id (Uniprot) and ptm site information, for example P27361_T202.
      ## The proteomics data has the ptm site information in a different format, i.e. S233 vs ser233.
      node.table.ptm <- node.table %>% 
        filter(ptm == 'p') %>%
        mutate(data_mapping_id=paste0(parentid, "_", position)) %>% 
        mutate(data_mapping_id = str_replace(data_mapping_id, "ser", "S"),
               data_mapping_id = str_replace(data_mapping_id, "thr", "T"),
               data_mapping_id = str_replace(data_mapping_id, "tyr", "Y")) %>%
        mutate(name=position) %>%
        mutate(name = str_replace(name, "^ser", "S"),
               name = str_replace(name, "^thr", "T"),
               name = str_replace(name, "^tyr", "Y")) %>%
        dplyr::select(SUID, data_mapping_id, name)
      
      loadTableData(node.table.ptm, data.key.column="SUID", "node", table.key.column = 'SUID') ## load the new column into Cytoscape
      
      ## Add Ensembl ids to data_mapping_id column. For the proteomics data, this will be the mapping id. We are simply moving it from the Ensembl column to the data_mapping_id
      node.table.ensembl <- node.table %>% 
        filter(Type == 'GeneProduct') %>%
        mutate(data_mapping_id=Ensembl) %>% 
        dplyr::select(SUID, data_mapping_id)
      
      loadTableData(node.table.ensembl, data.key.column="SUID", "node", table.key.column = 'SUID') ## load the new column into Cytoscape
      
      ## Remove bypasses for ptm nodes that will interfere with data visualization
      setNodeFontSizeBypass(node.table.ptm$SUID, 9)
      
      for (s in node.table.ptm$SUID) {
        clearNodePropertyBypass(s, visual.property = "NODE_FILL_COLOR")
        clearNodePropertyBypass(s, visual.property = "NODE_LABEL")
        setNodeFontSizeBypass(s, 9)
        setNodeWidthBypass(s, 35)
        setNodeHeightBypass(s, 20)
      }
    
      cptac.phospho.mapped <- merge(cptac.phospho, biomart, by.x = "protein", by.y = "ensembl_peptide_id") %>%
        mutate(prot_site=paste0(uniprotswissprot, "_", site))
      
      ## Load phosphoproteomics and proteomics data into Cytoscape, using the newly added data_mapping_id column.
      loadTableData(cptac.phospho.mapped, data.key.column="prot_site", "node", table.key.column = 'data_mapping_id') 
      loadTableData(cptac.protein, data.key.column="ensembl", "node", table.key.column = "Ensembl")
      
      ## Set visual style
      style.name = "WikiPathways"
      setNodeColorDefault('#FFFFFF', style.name = style.name)
      setNodeBorderColorDefault("#737373", style.name = style.name)
      
      ## The above for loop is a workaround due to a bug in RCy3. Once the bug is fixed, it should be updated to:
      ##clearNodePropertyBypass(node.names = ptm.nodes$SUID, visual.property = "NODE_FILL_COLOR") ## this doesnt work
      ##clearNodePropertyBypass(node.names = ptm.nodes$SUID, visual.property = "NODE_LABEL") ## this doesnt work
      
      ##RCy3::setNodeColorMapping('CCRCC.val', colors=paletteColorBrewerRdBu, style.name = style.name) 
      
      RCy3::setNodeColorMapping(table.column = type.val, colors = paletteColorBrewerRdBu, 
                                style.name = style.name)
      
    }
  })
}

## Run the Shiny App
shinyApp(ui, server)