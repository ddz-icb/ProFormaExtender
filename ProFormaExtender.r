##########################################################################
#                  NECESSARY SETTINGS & DEPENDENCIES                     #
##########################################################################

options(scipen = 999)

library(readxl)
library(stringr)

##########################################################################
#                     FUNCTION DEFINITIONS                               #
##########################################################################

source(file="E:/Projekte/Phospho2.1/ddz-pd-phospho-lib.r")



##########################################################################
#                               MAIN                                     #
##########################################################################

cwd <- "E:/Projekte_small/ProFormaExtender"

dat.orig <- as.data.frame(read_excel(paste0(cwd, "/Phospho2.1_Quan_with_isoforms-(1)_phosphopeptides.xlsx")))
id.mapping <- read.table(paste0(cwd,"/20220702_Phospho2.1_UniProt_to_GeneName.tsv"),header=T,sep="\t")

mod.orig <- dat.orig$Modifications
seq.orig <- dat.orig$`Annotated Sequence`
master.orig <- dat.orig$`Master Protein Accessions`
proforma.vec <- rep(NA,nrow(dat.orig))
proforma.extend.vec <- rep(NA,nrow(dat.orig))
for(i in 1:nrow(dat.orig)){
#for(i in 11900:11977){
  print(seq.orig[i])
  seq <- gsub("^\\[.+\\]\\.","",seq.orig[i])
  seq <- gsub("\\.\\[.+\\]$","",seq)
  print(seq)
  
  
  
  mod.vec <- unlist(strsplit(mod.orig[i], "]; "))
  print(mod.vec)
  if(length(mod.vec) > 1){
    for(j in 1:length(mod.vec)-1){
      mod.vec[j] <- paste0(mod.vec[j], "]")
    }
  }
  print(mod.vec)
  
  nterm.mod <- ""
  seqpos.mod <- matrix(nrow=0, ncol=2)
  colnames(seqpos.mod) <- c("Modification","Position")
  for(j in 1:length(mod.vec)){
    if(grepl("N-Term",mod.vec[j])){
      nterm.mod <- gsub(" \\[N-Term\\]","",mod.vec[j])
      nterm.mod <- gsub("\\d+x","",nterm.mod)
      nterm.mod <- paste0("[",nterm.mod,"]-")
      print(nterm.mod)
    }else if(grepl("Carbamidomethyl",mod.vec[j])){
      carbamidomethyl.mod <- gsub("^\\d+xCarbamidomethyl \\[","",mod.vec[j])
      carbamidomethyl.mod <- gsub("\\]$","",carbamidomethyl.mod)
      print(carbamidomethyl.mod)
      carbamidomethyl.mod.vec <- unlist(strsplit(carbamidomethyl.mod, "; "))
      print(carbamidomethyl.mod.vec)
      for(k in 1:length(carbamidomethyl.mod.vec)){
        mod.pos.tmp <- gsub("C","",carbamidomethyl.mod.vec[k])
        if(mod.pos.tmp=="") mod.pos.tmp <- "-Inf"
        seqpos.mod <- rbind(seqpos.mod, c("Carbamidomethyl", mod.pos.tmp))
      }
    }else if(grepl("Oxidation",mod.vec[j])){
      oxidation.mod <- gsub("^\\d+xOxidation \\[","",mod.vec[j])
      oxidation.mod <- gsub("\\]$","",oxidation.mod)
      print(oxidation.mod)
      oxidation.mod.vec <- unlist(strsplit(oxidation.mod, "; "))
      print(oxidation.mod.vec)
      for(k in 1:length(oxidation.mod.vec)){
        mod.pos.tmp <- gsub("M","",oxidation.mod.vec[k])
        if(mod.pos.tmp=="") mod.pos.tmp <- "-Inf"
        seqpos.mod <- rbind(seqpos.mod, c("Oxidation", mod.pos.tmp))
      }
    }else if(grepl("Phospho",mod.vec[j])){
      phospho.count <- gsub("^(\\d+)xPhospho.+","\\1", mod.vec[j])
      print(paste0("phospho.count: ", phospho.count))
      phospho.mod <- gsub("^\\d+xPhospho \\[","",mod.vec[j])
      phospho.mod <- gsub("\\]$","",phospho.mod)
      phospho.mod <- gsub("\\(\\d+\\.*\\d*\\)","",phospho.mod)
      print(phospho.mod)
      if(as.numeric(phospho.count)>1 & !grepl(";", phospho.mod)){
        print("HERE!!!")
        for(k in 1:(as.numeric(phospho.count)-1)){
          phospho.mod <- paste0(phospho.mod,"; ",phospho.mod)  
        }
      }
      phospho.mod.vec <- unlist(strsplit(phospho.mod, "; "))
      print(phospho.mod.vec)
      for(k in 1:length(phospho.mod.vec)){
        mod.pos.tmp <- gsub("S|T|Y|\\/","",phospho.mod.vec[k])
        if(mod.pos.tmp=="") mod.pos.tmp <- "-Inf"
        seqpos.mod <- rbind(seqpos.mod, c("Phospho", mod.pos.tmp))
      }
    }  
  }
  print(seqpos.mod)
  seqpos.mod <- seqpos.mod[order(as.numeric(seqpos.mod[,"Position"]),decreasing=T),,drop=F]
  print(seqpos.mod)
  
  
  proforma <- seq
  for(j in 1:nrow(seqpos.mod)){
    if(seqpos.mod[j,"Position"] == "-Inf"){
      if(nterm.mod != ""){
        proforma <- paste0(nterm.mod, proforma)
        nterm.mod <- ""
      }
      proforma <- paste0("[",seqpos.mod[j,"Modification"],"]?",proforma)
    }else{
      proforma <- paste0(substr(proforma,1,as.numeric(seqpos.mod[j,"Position"])),
                         paste0("[", seqpos.mod[j,"Modification"], "]"),
                         substr(proforma,as.numeric(seqpos.mod[j,"Position"])+1,nchar(proforma)))
    }
  }
  if(nterm.mod != ""){
    proforma <- paste0(nterm.mod, proforma)
    nterm.mod <- ""  
  }
  print(proforma)
  proforma.vec[i] <- proforma
  
  
  master.tmp <- gsub("-\\d+","",master.orig[i])
  master.tmp.vec <- unique(unlist(strsplit(master.tmp, "; ")))
  genes.tmp <- rep(NA,length(master.tmp.vec))
  for(j in 1:length(master.tmp.vec)){
    genes.tmp[j] <- (id.mapping[id.mapping$From==master.tmp.vec[j],"To"])[1]  
  }
  print(genes.tmp)
  genes.tmp2 <- paste0(genes.tmp, collapse="; ")
  proforma.extend.vec[i] <- paste0(genes.tmp2, "_", proforma)
  print(proforma.extend.vec[i])
  print("---------------------------------------------------------------------")
}

dat.orig$`Master Protein Descriptions` <- gsub("\r\n", "; ", dat.orig$`Master Protein Descriptions`)
write.table(cbind(proforma.vec, proforma.extend.vec, dat.orig),
            file=paste0(cwd,"/proforma-ext.output.txt"),
            col.names=T,
            row.names=F,
            quote=F,
            sep="\t")