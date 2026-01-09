This repository provides tools for generating ProForma-annotated output from Proteome Discoverer phospho data. The project offers two different ways to run the same analysis pipeline: a Shiny-based user interface and a command-line interface.
All computational logic is implemented in a single shared file called core_calculation.R. Both the Shiny application and the command-line script source this file and call the same core function.

# Shiny App
You need core_calculation.R and proforma_cli.R for the Shiny App version. 

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



# Commandline
Needed files: core_calculation.R and app_commandLine.R

Needed packages:
install.packages(c("optparse", "readxl", "stringr", "httr", "jsonlite"))

Run commandline:
Rscript proforma_cli.R --input input.xlsx

OR

Rscript proforma_cli.R --input input.xlsx --output result.tsv

