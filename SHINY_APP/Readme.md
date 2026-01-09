This repository provides tools for generating ProForma-annotated output from Proteome Discoverer phospho data. The project offers two different ways to run the same analysis pipeline: a Shiny-based user interface and a command-line interface.
All computational logic is implemented in a single shared file called core_calculation.R. Both the Shiny application and the command-line script source this file and call the same core function.

# Shiny App
To run the Shiny version, the files app_commandLine.R and core_calculation.R must be located in the same directory.


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

Run app:
shiny::runApp("app_commandLine.R")

Once the app is running, users can upload a Proteome Discoverer Excel output file (.xlsx) using the file upload field in the sidebar.



# Commandline
To use the CLI, the files proforma_cli.R and core_calculation.R must be present in the same working directory.

Needed packages:
install.packages(c("optparse", "readxl", "stringr", "httr", "jsonlite"))

The script is executed from a terminal or command prompt using Rscript. The command must be run from the directory containing the scripts or the script paths must be specified explicitly. The minimum required argument is an input Excel file, which must again be a Proteome Discoverer output file in .xlsx format. An optional output filename can be provided. Otherwise, the script automatically generates a date-based filename.

Run commandline:
Rscript proforma_cli.R --input input.xlsx

OR

Rscript proforma_cli.R --input input.xlsx --output result.tsv

