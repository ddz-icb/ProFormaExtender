# An R installation is required!

This repository provides tools for generating ProForma-annotated output from phosphoproteomics data. It supports input files directly from Proteome Discoverer (PD) as well as MaxQuant (MQ). 
All computational logic is implemented in a single shared file called core_calculation.R. Both the PD and the MQ paths
source this file and call the same core function.

# Shiny App
To run the Shiny version, you need to open RStudio. Please install if you haven't.
The files RShiny_ui.R and core_calculation.R must be located in the same directory. Open the RShiny_ui.R file in RStudio. It is important that you have all the packages needed:


Needed packages:
  install.packages(c(
    "shiny",
    "shinyjs",
    "shinyBS",
    "readxl",
    "stringr",
    "httr",
    "jsonlite"
  ))

# File Requirements:
Ensure the following files are located in the same working directory:
RShiny_ui_v2.R   
core_calculation.R   
mq_to_pd_converter.py 

You can either start the application by using the run button in RStudio or navigating in the terminal to the directory of the files and then use the following command:

# Run app:
shiny::runApp("RShiny_ui_v2.R")

Once the app is startet, the workflow is as follows: 


# 1. Select Input Type:
- Proteome Discoverer: Directly upload a standard PD Excel export (.xlsx)
<img width="484" height="813" alt="Bildschirmfoto 2026-07-24 um 11 37 42" src="https://github.com/user-attachments/assets/da07d97b-377c-40e9-8852-cea7545d4235" />
MaxQuant: Upload your MaxQuant output files:  modificationSpecificPeptides.txt, Phospho (STY) Sites.txt, peptides.txt and click Convert MQ → Excel
# 2. Run ProForma Pipeline:
     Click Start Calculation
     The pipeline will parse the data and query the UniProt REST API to generate ProForma annotations and phosphosite mapping
# 3. Download Results:
     Preview the result table in the app and download the final annotated dataset as an Excel file (.xlsx)






Rscript proforma_cli.R --input input.xlsx --output result.tsv
