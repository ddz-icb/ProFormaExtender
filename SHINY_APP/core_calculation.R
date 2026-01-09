############################################################
# core_calculation.R
# Zentrale Berechnung (inkl. UniProt API)
# Wird von Shiny und Command-Line verwendet
############################################################

library(httr)
library(readxl)
library(jsonlite)
library(stringr)

############################################################
# 1) UniProt ID Mapping (API)
############################################################

isJobReady <- function(jobId) {
  pollingInterval = 3
  nTries = 40
  
  for (i in 1:nTries) {
    url <- paste0("https://rest.uniprot.org/idmapping/status/", jobId)
    r <- GET(url, accept_json())
    status <- content(r, "parsed")
    
    if (!is.null(status[["results"]]) ||
        !is.null(status[["failedIds"]]) ||
        !is.null(status[["redirectURL"]])) {
      return(TRUE)
    }
    if (!is.null(status[["messages"]])) {
      print(status[["messages"]])
      return(FALSE)
    }
    Sys.sleep(pollingInterval)
  }
  return(FALSE)
}

getResultsURL <- function(redirectURL) {
  if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
    gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
  } else {
    gsub("/results/", "/results/stream/", redirectURL)
  }
}

uniprot_map_ids <- function(ids,
                            from = "UniProtKB_AC-ID",
                            to   = "Gene_Name") {
  
  ids <- unique(na.omit(as.character(ids)))
  if (length(ids) == 0) {
    return(data.frame(From=character(0), To=character(0)))
  }
  
  # Submit job
  r <- POST(
    url = "https://rest.uniprot.org/idmapping/run",
    body = list(
      ids  = paste(ids, collapse = ","),
      from = from,
      to   = to
    ),
    encode = "multipart",
    accept_json()
  )
  
  jobId <- content(r, "parsed")$jobId
  if (is.null(jobId))
    stop("UniProt error: no jobId returned.")
  
  # Poll job ready
  if (!isJobReady(jobId))
    stop("UniProt mapping timed out.")
  
  # Get redirect
  rdet <- GET(paste0("https://rest.uniprot.org/idmapping/details/", jobId), accept_json())
  details <- content(rdet, "parsed")
  
  url <- getResultsURL(details$redirectURL)
  url <- paste0(url, "?format=tsv")
  
  rr <- GET(url)
  df <- read.table(text = content(rr, "text"), sep = "\t", header = TRUE)
  names(df) <- c("From", "To")
  df
}

############################################################
# 2) Hauptfunktion: run_proforma_pipeline()
############################################################

run_proforma_pipeline <- function(excel_path) {
  
  dat.orig <- as.data.frame(read_excel(excel_path, na = "NA"))
  
  # Extract UniProt IDs
  master_ids_raw <- dat.orig$`Master Protein Accessions`
  all_ids <- unique(unlist(strsplit(master_ids_raw, "; ")))
  
  # Perform API mapping
  id.mapping <- uniprot_map_ids(all_ids)
  
  # Identify sequence field
  if("Annotated Sequence" %in% colnames(dat.orig)) {
    seq.orig <- dat.orig$`Annotated Sequence`
  } else if("Sequence" %in% colnames(dat.orig)) {
    seq.orig <- dat.orig$Sequence
  } else {
    stop("Sequence column not found")
  }
  
  mod.orig <- dat.orig$Modifications
  mod.master.orig <- if("Modifications in Master Proteins" %in% colnames(dat.orig))
    dat.orig$`Modifications in Master Proteins`
  else
    rep(NA, nrow(dat.orig))
  
  master.orig <- dat.orig$`Master Protein Accessions`
  
  n <- nrow(dat.orig)
  proforma.vec <- rep(NA, n)
  proforma.extend.vec <- rep(NA, n)
  phossite.id.vec <- rep(NA, n)
  phossite.id.vec2 <- rep(NA, n)
  phossite.id.vec3 <- rep(NA, n)
  
  ############################################################
  # Hauptschleife über alle Zeilen
  ############################################################
  
  for(i in seq_len(n)) {
    
    seq_i <- seq.orig[i]
    if(is.na(seq_i)) seq_i <- ""
    seq <- gsub("^\\[.+\\]\\.","", seq_i)
    seq <- gsub("\\.\\[.+\\]$","", seq)
    
    mod_i <- mod.orig[i]
    if(is.na(mod_i)) mod_i <- ""
    mod.vec <- if(nzchar(mod_i)) unlist(strsplit(mod_i, "]; ")) else character(0)
    
    if(length(mod.vec) > 1) {
      mod.vec[1:(length(mod.vec)-1)] <- paste0(mod.vec[1:(length(mod.vec)-1)], "]")
    }
    
    nterm.mod <- ""
    seqpos.mod <- matrix(ncol=2, nrow=0)
    colnames(seqpos.mod) <- c("Modification","Position")
    
    # ----- Parse modifications -----
    for(j in seq_along(mod.vec)) {
      mv <- mod.vec[j]
      
      if(grepl("N-Term", mv)) {
        nterm.mod <- gsub(" \\[N-Term\\]","", mv)
        nterm.mod <- gsub("\\d+x","", nterm.mod)
        nterm.mod <- paste0("[", nterm.mod, "]-")
      }
      else if(grepl("Carbamidomethyl", mv)) {
        cm <- gsub("^\\d+xCarbamidomethyl \\[","", mv)
        cm <- gsub("\\]$","", cm)
        for(x in unlist(strsplit(cm, "; "))) {
          pos <- gsub("C","", x)
          if(pos=="") pos <- "-Inf"
          seqpos.mod <- rbind(seqpos.mod, c("Carbamidomethyl", pos))
        }
      }
      else if(grepl("Oxidation", mv)) {
        ox <- gsub("^\\d+xOxidation \\[","", mv)
        ox <- gsub("\\]$","", ox)
        for(x in unlist(strsplit(ox, "; "))) {
          pos <- gsub("M","", x)
          if(pos=="") pos <- "-Inf"
          seqpos.mod <- rbind(seqpos.mod, c("Oxidation", pos))
        }
      }
      else if(grepl("Phospho", mv)) {
        phospho.count <- as.numeric(gsub("^(\\d+)xPhospho.+","\\1", mv))
        ph <- gsub("^\\d+xPhospho \\[","", mv)
        ph <- gsub("\\]$","", ph)
        ph <- gsub("\\(\\d+\\.*\\d*\\)","", ph)
        
        if(!is.na(phospho.count) && phospho.count > 1 && !grepl(";", ph)) {
          for(k in 1:(phospho.count-1)) {
            ph <- paste0(ph, "; ", ph)
          }
        }
        
        ph.vec <- unlist(strsplit(ph, "; "))
        for(x in ph.vec) {
          pos <- gsub("S|T|Y|\\/","", x)
          if(pos=="") pos <- "-Inf"
          seqpos.mod <- rbind(seqpos.mod, c("Phospho", pos))
        }
      }
    }
    
    # Sort descending
    if(nrow(seqpos.mod) > 0) {
      seqpos.mod <- seqpos.mod[order(as.numeric(seqpos.mod[,"Position"]), decreasing=TRUE),,drop=FALSE]
    }
    
    # ----- Build ProForma -----
    proforma <- seq
    for(j in seq_len(nrow(seqpos.mod))) {
      posj <- seqpos.mod[j,"Position"]
      modj <- seqpos.mod[j,"Modification"]
      
      if(posj == "-Inf") {
        if(nchar(nterm.mod)>0) {
          proforma <- paste0(nterm.mod, proforma)
          nterm.mod <- ""
        }
        proforma <- paste0("[",modj,"]?", proforma)
      } else {
        pnum <- as.numeric(posj)
        proforma <- paste0(
          substr(proforma, 1, pnum),
          paste0("[",modj,"]"),
          substr(proforma, pnum+1, nchar(proforma))
        )
      }
    }
    if(nchar(nterm.mod)>0) {
      proforma <- paste0(nterm.mod, proforma)
    }
    
    proforma <- gsub("pho\\]\\?\\[Pho","pho\\]\\[Pho", proforma)
    proforma.vec[i] <- proforma
    
    # ----- Phosphosite ID handling -----
    
    master.tmp <- gsub("-\\d+","", master.orig[i])
    master.tmp.vec <- unique(unlist(strsplit(master.tmp, "; ")))
    master.tmp.vec2 <- unlist(strsplit(master.orig[i], "; "))
    
    mm <- mod.master.orig[i]
    if(is.na(mm) || !nzchar(mm)) mm <- mod.orig[i]
    if(is.na(mm)) mm <- ""
    
    mod.master.tmp.vec <- if(nzchar(mm)) unlist(strsplit(mm, "]; ")) else character(0)
    
    if(length(mod.master.tmp.vec) > 1) {
      mod.master.tmp.vec[1:(length(mod.master.tmp.vec)-1)] <-
        paste0(mod.master.tmp.vec[1:(length(mod.master.tmp.vec)-1)], "]")
    }
    
    # keep only phospho
    mod.master.tmp.vec <- mod.master.tmp.vec[grep("Phospho", mod.master.tmp.vec)]
    mod.master.tmp.vec <- gsub("^.+\\[","", mod.master.tmp.vec)
    mod.master.tmp.vec <- gsub("\\].*$","", mod.master.tmp.vec)
    mod.master.tmp.vec <- gsub("\\/","", mod.master.tmp.vec)
    mod.master.tmp.vec <- gsub("\\(\\d+\\.*\\d*\\)","", mod.master.tmp.vec)
    mod.master.tmp.vec <- gsub(";",",", mod.master.tmp.vec)
    
    genes.tmp <- rep(NA, length(master.tmp.vec))
    for(j in seq_along(master.tmp.vec)) {
      mmatch <- id.mapping[id.mapping$From == master.tmp.vec[j], "To"]
      genes.tmp[j] <- if(length(mmatch)>=1) mmatch[1] else NA
    }
    
    if(length(master.tmp.vec2)>1 && length(mod.master.tmp.vec)>0 &&
       (length(mod.master.tmp.vec) %% length(master.tmp.vec2) != 0)) {
      master.tmp.vec2 <- rep(master.tmp.vec2, each=length(mod.master.tmp.vec))
    } else if(length(master.tmp.vec2)>1 && length(mod.master.tmp.vec)>0) {
      master.tmp.vec2 <- rep(master.tmp.vec2, each=length(mod.master.tmp.vec)/length(master.tmp.vec2))
    }
    
    if(length(mod.master.tmp.vec)>0) {
      phossite.id.tmp  <- paste0(master.tmp.vec2, "_", mod.master.tmp.vec)
      phossite.id.tmp2 <- unique(paste0(genes.tmp, "_", mod.master.tmp.vec))
      phossite.id.tmp3 <- paste0(master.tmp.vec2, "_", phossite.id.tmp2)
    } else {
      phossite.id.tmp <- NA
      phossite.id.tmp2 <- NA
      phossite.id.tmp3 <- NA
    }
    
    phossite.id.vec[i]  <- paste(phossite.id.tmp,  collapse="; ")
    phossite.id.vec2[i] <- paste(phossite.id.tmp2, collapse="; ")
    phossite.id.vec3[i] <- paste(phossite.id.tmp3, collapse="; ")
    
    genes.tmp2 <- paste(genes.tmp, collapse="; ")
    proforma.extend.vec[i] <- paste0(genes.tmp2, "_", proforma)
  }
  
  # Clean master desc
  if("Master Protein Descriptions" %in% colnames(dat.orig)) {
    dat.orig$`Master Protein Descriptions` <- gsub("\r\n", "; ", dat.orig$`Master Protein Descriptions`)
    dat.orig$`Master Protein Descriptions` <- gsub("\\'", "", dat.orig$`Master Protein Descriptions`)
  }
  
  out.df <- cbind(
    proforma.vec,
    proforma.extend.vec,
    phossite.id.vec,
    phossite.id.vec2,
    phossite.id.vec3,
    dat.orig
  )
  
  return(out.df)
}
