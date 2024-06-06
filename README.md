# HTS-interactive-application
Using R &amp; Python, this R Shiny application was developed to visualize and interact with high-throughput data

# HTS Interactive Application

## Overview
This Shiny application allows users to convert SMILES strings from an Excel file into molecule images, visualize the data in a heatmap, and interactively view detailed molecular structures. The app reads an Excel file containing SMILES strings and associated metadata, generates molecular images, and displays a heatmap with clickable cells for detailed inspection.

## Features
- Convert SMILES strings to molecule images.
- Visualize data in a heatmap.
- Interactive exploration of molecular structures and associated metadata.
- Customizable threshold for heatmap normalization.

## Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)
- Python (version 3.7 or higher)
- [rdkit](https://www.rdkit.org/) Python package
- Required R packages:
  - shiny
  - reticulate
  - openxlsx
  - dplyr
  - heatmaply
  - plotly

## Setup

### Step 1: Install R and RStudio
Download and install R from [CRAN](https://cran.r-project.org/).
Download and install RStudio from [RStudio](https://rstudio.com/products/rstudio/download/).

### Step 2: Install Python and rdkit
Ensure you have Python installed. You can download it from [python.org](https://www.python.org/).

Create a Conda environment for rdkit:

```sh
conda create -c conda-forge -n rdkit-env rdkit
conda activate rdkit-env
```

### Step 3: Install Required R Packages
Open R or RStudio and run the following commands:

```sh
install.packages(c("shiny", "reticulate", "openxlsx", "dplyr", "heatmaply", "plotly"))
```
### Step 4: Clone the Repository
Clone this repository to your local machine:

```sh
git clone https://github.com/YOUR-USERNAME/HTS-interactive-application.git
cd HTS-interactive-application
```

### Step 5: Run the Shiny App
Open ui.R in RStudio and click the "Run App" button, or run the following command in R:

```
shiny::runApp()
```

### Usage
1. Upload an Excel File: Use the file input to upload an Excel file containing SMILES strings and associated metadata.
2. Convert SMILES to Molecules: Click the "Convert SMILES" button to generate molecule images.
3. Navigate Pages: Use the "Previous Page" and "Next Page" buttons to navigate through the pages of molecules.
4. Create Heatmap: Click the "Create Heatmap" button to generate and display a heatmap.
5. Interact with Heatmap: Click on cells in the heatmap to view detailed molecular structures and metadata.
