##########################################################################
#                  NECESSARY SETTINGS & DEPENDENCIES                     #
##########################################################################

options(scipen = 999)

library(readxl)
library(stringr)

##########################################################################
#                     FUNCTION DEFINITIONS                               #
##########################################################################

source(file="C:/Workspace/Projekte/Phospho2.1/ddz-pd-phospho-lib.r")



##########################################################################
#                               MAIN                                     #
##########################################################################

cwd <- "C:/Workspace/Projekte_small/ProFormaExtender"

#dat.orig <- as.data.frame(read_excel(paste0(cwd, "/Phospho2.1_Quan_with_isoforms-(1)_phosphopeptides.xlsx")))
#id.mapping <- read.table(paste0(cwd,"/20220702_Phospho2.1_UniProt_to_GeneName.tsv"),header=T,sep="\t")

#dat.orig <- as.data.frame(read_excel(paste0(cwd, "/humPhosRes_PhosphoPeptides_OE480_Aurora_2hDDA_nested_LFQ_Sequest_PhosPep.xlsx")))
#id.mapping <- read.table(paste0(cwd,"/20240628_humPhosRes_UniProt_to_GeneName.tsv"),header=T,sep="\t")

#dat.orig <- as.data.frame(read_excel(paste0(cwd, "/huMuPhosProt_PhosPep_n9-lean-obese-T2D_basal-clamp_PEPTIDES.xlsx")))
#id.mapping <- read.table(paste0(cwd,"/20250514_Harry_UniProt_to_GeneName_with_imputation.tsv"),header=T,sep="\t")

#dat.orig <- as.data.frame(read_excel(paste0(cwd, "/20250401_huMuPhosProt_PhosPep_n9-lean-obese-T2D_basal-clamp_noIMP_PEPTIDES.xlsx")))
#id.mapping <- read.table(paste0(cwd,"/20250514_Harry_UniProt_to_GeneName_with_imputation.tsv"),header=T,sep="\t")

dat.orig <- as.data.frame(read_excel(paste0(cwd, "/20250901_harry_9v9v9_imp-abundances_export.xlsx"), na="NA"))
id.mapping <- read.table(paste0(cwd,"/20250714_harry_noIMP_9vs9vs9_idmapping_2025_07_14.tsv"),header=T,sep="\t")



mod.orig <- dat.orig$Modifications
mod.master.orig <- dat.orig$`Modifications in Master Proteins`
if(any(grepl("Annotated Sequence", colnames(dat.orig)))){
  seq.orig <- dat.orig$`Annotated Sequence`
}else if(any(grepl("Sequence", colnames(dat.orig)))){
  seq.orig <- dat.orig$Sequence
}else{
  stop("Error! Sequence column not found!")
}
master.orig <- dat.orig$`Master Protein Accessions`
#all.uniprot.ids <- gsub("-\\d+","",master.orig)
#all.uniprot.ids <- unique(unlist(strsplit(all.uniprot.ids, "; ")))
#write.table(x=all.uniprot.ids,file=paste0(cwd,"/all.uniprot.ids.txt"),col.names=F,row.names=F,quote=F,sep="\n")
proforma.vec <- rep(NA,nrow(dat.orig))
proforma.extend.vec <- rep(NA,nrow(dat.orig))
phossite.id.vec <- rep(NA,nrow(dat.orig))
phossite.id.vec2 <- rep(NA,nrow(dat.orig))
phossite.id.vec3 <- rep(NA,nrow(dat.orig))
for(i in 1:nrow(dat.orig)){
#for(i in 1:154){
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
      oxidation.mod <- gsub("\\(\\d+\\.*\\d*\\)","",oxidation.mod)
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
  proforma <- gsub("pho\\]\\?\\[Pho","pho\\]\\[Pho",proforma)
  print(proforma)
  proforma.vec[i] <- proforma
  
  
  master.tmp <- gsub("-\\d+","",master.orig[i])
  master.tmp.vec <- unique(unlist(strsplit(master.tmp, "; ")))
  master.tmp.vec2 <- unlist(strsplit(master.orig[i], "; "))
  # NEW
  if(is.na(mod.master.orig[i])) mod.master.orig[i] <- mod.orig[i]
  if(grepl("Phospho", mod.master.orig[i]) == F) mod.master.orig[i] <- mod.orig[i]
  print(mod.master.orig[i])
  #mod.master.tmp.vec <- unique(unlist(strsplit(mod.master.orig[i], "]; ")))
  mod.master.tmp.vec <- unlist(strsplit(mod.master.orig[i], "]; "))
  print(mod.master.tmp.vec)
  mod.master.tmp.vec <- mod.master.tmp.vec[grep("Phospho", mod.master.tmp.vec)]
  if(length(mod.master.tmp.vec) > 1){
    mod.master.tmp.vec[1:(length(mod.master.tmp.vec)-1)] <- paste0(mod.master.tmp.vec[1:(length(mod.master.tmp.vec)-1)],"]")
  }
  print(mod.master.tmp.vec)
  mod.master.tmp.vec <- gsub("^.+\\[","",mod.master.tmp.vec)
  mod.master.tmp.vec <- gsub("\\].*$","",mod.master.tmp.vec)
  mod.master.tmp.vec <- gsub("\\(\\d+\\.*\\d*\\)","",mod.master.tmp.vec)
  mod.master.tmp.vec <- gsub("\\/","",mod.master.tmp.vec)
  print("mod.master.tmp.vec:")
  print(mod.master.tmp.vec)
  mod.master.tmp.vec <- gsub(";",",",mod.master.tmp.vec)
  print("mod.master.tmp.vec:")
  print(mod.master.tmp.vec)
  print("master.tmp.vec2: ")
  print(master.tmp.vec2)
  # NEW
  genes.tmp <- rep(NA,length(master.tmp.vec))
  for(j in 1:length(master.tmp.vec)){
    genes.tmp[j] <- (id.mapping[id.mapping$From==master.tmp.vec[j],"To"])[1]
  }
  print(genes.tmp)
  
  #if(length(master.tmp.vec2) > 1) master.tmp.vec2 <- rep(master.tmp.vec2,
  #                                                       each=length(mod.master.tmp.vec)/length(master.tmp.vec2))
  if(length(master.tmp.vec2) > 1 && length(mod.master.tmp.vec) %% length(master.tmp.vec2)!=0){
    master.tmp.vec2 <- rep(master.tmp.vec2,
                           each=length(mod.master.tmp.vec))  
  }else if(length(master.tmp.vec2) > 1){
    master.tmp.vec2 <- rep(master.tmp.vec2,
                           each=length(mod.master.tmp.vec)/length(master.tmp.vec2))
  }
  phossite.id.tmp <- paste0(master.tmp.vec2, "_", mod.master.tmp.vec)
  phossite.id.tmp2 <- unique(paste0(genes.tmp, "_", mod.master.tmp.vec))
  phossite.id.tmp3 <- paste0(master.tmp.vec2, "_", phossite.id.tmp2)
  
  phossite.id.vec[i] <- paste0(phossite.id.tmp, collapse="; ")
  phossite.id.vec2[i] <- paste0(phossite.id.tmp2, collapse="; ")
  phossite.id.vec3[i] <- paste0(phossite.id.tmp3, collapse="; ")
  genes.tmp2 <- paste0(genes.tmp, collapse="; ")
  proforma.extend.vec[i] <- paste0(genes.tmp2, "_", proforma)
  
  print(proforma.extend.vec[i])
  print(phossite.id.vec[i])
  print(phossite.id.vec2[i])
  print(phossite.id.vec3[i])
  print("---------------------------------------------------------------------")
}

timestamp <- gsub("^(\\d{4}-\\d{2}-\\d{2}).+", "\\1", Sys.time())
timestamp <- gsub("-", "", timestamp)
if(!is.null(dat.orig$`Master Protein Descriptions`)){
  dat.orig$`Master Protein Descriptions` <- gsub("\r\n", "; ", dat.orig$`Master Protein Descriptions`)
  dat.orig$`Master Protein Descriptions` <- gsub("\\'", "", dat.orig$`Master Protein Descriptions`)
}
write.table(cbind(proforma.vec, proforma.extend.vec, phossite.id.vec, phossite.id.vec2, phossite.id.vec3, dat.orig),
            file=paste0(cwd,"/", timestamp, "_proforma-ext.output.txt"),
            col.names=T,
            row.names=F,
            quote=F,
            sep="\t")