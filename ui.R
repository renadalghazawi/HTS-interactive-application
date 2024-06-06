library(shiny)
library(reticulate)
library(openxlsx)
library(dplyr)
library(heatmaply)
library(plotly)

# Ensure the correct Conda environment is used
use_condaenv("rdkit-env", required = TRUE)

# Source the Python script
source_python("convert_smiles.py")

# Add resource path for custom directory (if not using 'www')
addResourcePath("static", "www")

ui <- fluidPage(
  titlePanel("SMILES to Molecule Converter"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose Excel File", accept = c(".xls", ".xlsx")),
      numericInput("page_size", "Molecules per Page:", value = 10, min = 1, max = 100),
      actionButton("prev_page", "Previous Page"),
      actionButton("next_page", "Next Page"),
      textOutput("page_info"),
      actionButton("convert", "Convert SMILES"),
      hr(),
      actionButton("create_heatmap", "Create Heatmap"),
      numericInput("threshold_value", "Threshold Value for Heatmap:", value = 100),
      textOutput("heatmap_info")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Molecules", uiOutput("molecule_images")),
        tabPanel("Heatmap", plotlyOutput("heatmap"), uiOutput("selected_molecule"))
      )
    )
  )
)