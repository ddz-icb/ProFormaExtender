This repository provides tools for generating ProForma-annotated output from Proteome Discoverer phospho data. The project offers two different ways to run the same analysis pipeline: a Shiny-based user interface and a command-line interface.
All computational logic is implemented in a single shared file called core_calculation.R. Both the Shiny application and the command-line script source this file and call the same core function.

# Please install R!

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

You can either start the application by using the run button in RStudio or you navigate in the terminal to the directory of the files and then use the following command:

Run app:
shiny::runApp("RShiny_ui.R")

Once the app is running, users can upload a Proteome Discoverer Excel output file (.xlsx) using the file upload field in the sidebar.
The UNIPROT API is then called and maps you iput data.



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

