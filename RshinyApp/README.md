# An R installation is required!

This repository provides tools for generating ProForma-annotated output from phosphoproteomics data. It supports input files directly from Proteome Discoverer (PD) as well as MaxQuant (MQ). The project offers two different ways to run the same analysis pipeline: an Shiny-based user interface and an command-line interface.
All computational logic is implemented in a single shared file called core_calculation.R. Both the Shiny application and the command-line script source this file and call the same core function.

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
     Proteome Discoverer: Directly upload a standard PD Excel export (.xlsx)
     MaxQuant: Upload your MaxQuant output files:  modificationSpecificPeptides.txt (Required), Phospho (STY) Sites.txt (Required), peptides.txt (Optional, provides exact genomic           coordinates)
   Adjust filter options (e.g., remove/keep reverse hits or contaminants) and click Convert MQ → Excel
# 2. Run ProForma Pipeline:
     Click Start Calculation
     The pipeline will parse the data and query the UniProt REST API to generate ProForma annotations and phosphosite mapping
# 3. Download Results:
     Preview the result table in the app and download the final annotated dataset as an Excel file (.xlsx)





# Commandline
To use the CLI, the files proforma_cli.R and core_calculation.R must be present in the same working directory.
Please install following packages If they are not installed yet.

Needed packages:
install.packages(c("optparse", "readxl", "stringr", "httr", "jsonlite"))

The script is executed from a terminal or command prompt using Rscript. The command must be run from the directory containing the scripts or the script paths must be specified explicitly. The minimum required argument is an input Excel file, which must again be a Proteome Discoverer output file in .xlsx format. An optional output filename can be provided. Otherwise, the script automatically generates a date-based filename.


Please navigate to your directory including the files before running the following command:

Run commandline:
Rscript proforma_cli.R --input input.xlsx

OR

Rscript proforma_cli.R --input input.xlsx --output result.tsv
