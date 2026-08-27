install.packages("Rcpp", type = "binary")
install.packages("installr")
install.packages("tidyverse", type = "binary")
install.packages("purrr")
install.packages("magrittr")
install.packages("ggplot2")
install.packages("readxl")
install.packages("Hmisc")
install.packages("gplots")
install.packages("dplyr")
install.packages("stringr")
install.packages("data.table")
install.packages("gtsummary")
install.packages("tidyr")
install.packages("writexl")
install.packages("openxlsx", type = "binary")
install.packages("BiocManager")
install.packages("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v110")
BiocManager::install("AnnotationHub")


library(Rcpp)
library(purrr)
library(installr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(readxl)
library(Hmisc)
library(gplots)
library(dplyr)
library(stringr)
library(data.table)
library(gtsummary)
library(tidyr)
library(writexl)
library(openxlsx)
library(BiocManager)
library(GenomicRanges)
library(Biostrings)
library(AnnotationHub)
library(ensembldb)

ah <- AnnotationHub()

#fixes path that is copied from files into the proper format for read_excel
path_fix <- function(){
  out <- clipr::read_clip()
  out <- gsub("^\" | \" $","", out)
  out <- gsub("\\\\", "/", out)
  clipr::write_clip(out)
  message("The formatted path has been copied to your clipboard")
  return(out)
}

outputdir <- choose.dir(getwd(), "choose the desired folder :)")
clipr::write_clip(outputdir)
path_fix()
outputdir <- clipr::read_clip()
outputdir <- gsub('^"|"$', '', outputdir)
fileName = "SUM159_15-7-5_middlePEP.xlsx"

#write the name you want your output file to have in the selected folder in " "
outfile = paste(outputdir, fileName, sep = "/")

#choose the CryEX file to analyse and convert it into a data frame
cryEx_table<- choose.files(caption = "select desired RNAseq sheet to analyse")
clipr::clear_clip()
cryEx_table <- clipr::write_clip(cryEx_table)
path_fix()
cryEx_table_path <- clipr::read_clip()
cryEx_table_path <- gsub('^"|"$', '', cryEx_table)
cryEx_table <- read_excel(cryEx_table_path, sheet = "middle_exon") %>% as.data.frame()



library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
edb <-ah[["AH119325"]]

tx <- transcripts(edb, columns = c("gene_id", "tx_id", "tx_support_level", "tx_is_canonical", "gene_name" ,"protein_sequence","gene_biotype", "description", "tx_biotype"))
tx_canon <- tx[tx$tx_is_canonical == TRUE, ] %>% as.data.table()



flanklen = 45 #15 amino acids  will use this to check against protein_sequence for RF determination 260203

#minpeptidelength = 14 #8 is for HLA purposes but can be changed if desired
min_aa_len <- flanklen/3-1

#fix flanking junction (right side) shift -> comes from bash extract from IGV to original excel spreadsheet [260109]
cryEx_table$flanking_jxns <- sapply(cryEx_table$flanking_jxns, function(x){
  flank_split <- strsplit(x, ":")[[1]]
  nums <- as.integer(flank_split)
  nums[2] <- nums[2] + 1
  paste(nums, collapse = ":")
})

attach(cryEx_table)

n_cry <- nrow(cryEx_table)
n_rows <- n_cry * 3

pepcolnames <- c(
  "CryEx_row","Coords","Flanking","Strand","Frame",
  "Name","True_pep","RF_check","True_RF","Peptide", "Novel_true" ,"Length"
)

peptable <- as.data.frame(
  matrix(nrow = n_rows, ncol = length(pepcolnames)),
  stringsAsFactors = FALSE
)
colnames(peptable) <- pepcolnames

peptable$CryEx_row <- rep(seq_len(n_cry), each = 3)
peptable$Frame     <- rep(c(0,1,2), times = n_cry)
peptable$Coords    <- rep(cryEx_table$coords, each = 3)
peptable$Flanking  <- rep(cryEx_table$flanking_jxns, each = 3)
peptable$Strand    <- rep(cryEx_table$strand, each = 3)
peptable$Name      <- rep(cryEx_table$GENCODE_GeneName, each = 3)

#Find true reading frame

#Identify upstream exons and give full aa sequence of upstream exon
genome_row <- match(peptable$Name, tx_canon$gene_name) 
peptable$True_pep <- tx_canon$protein_sequence[genome_row]

#construct aa sequence of flanking exons to check against true aa sequence to find true RF

RF_check <- character(n_rows)

idx <- 1
for (rown in seq_len(n_cry)) {
  
  jx <- as.integer(strsplit(cryEx_table$flanking_jxns[rown], ":")[[1]])
  strand <- cryEx_table$strand[rown]
  
  if (strand == "+") {
    start_pos <- jx[1] - flanklen + 1
    end_pos   <- jx[1]
  } else {
    start_pos <- jx[2]
    end_pos   <- jx[2] + flanklen - 1 + 1
  }
  
  interval <- GRanges(
    seqnames = cryEx_table$chrom[rown],
    ranges   = IRanges(start = start_pos, end = end_pos),
    strand   = strand
  )
  
  dna <- if (strand == "+") {
    paste(getSeq(genome, interval), collapse = "")
  } else {
    paste(rev(getSeq(genome, interval)), collapse = "")
  }
  
  rf <- c(
    strsplit(as.character(translate(DNAString(dna))), "\\*")[[1]][1],
    strsplit(as.character(translate(DNAString(substr(dna,2,nchar(dna))))), "\\*")[[1]][1],
    strsplit(as.character(translate(DNAString(substr(dna,3,nchar(dna))))), "\\*")[[1]][1]
  )
  
  RF_check[idx:(idx+2)] <- rf
  idx <- idx + 3
}

peptable$RF_check <- RF_check



#checking 7aa in upstream flanking exon
True_RF <- character(n_rows)

for (i in seq_len(n_rows)) {
  
  chk <- peptable$RF_check[i]
  tp  <- peptable$True_pep[i]
  
  if (is.na(chk) || is.na(tp) || chk == "") {
    True_RF[i] <- "no reference aa"
  } else if (nchar(chk) < min_aa_len) {
    True_RF[i] <- "FALSE"
  } else if (grepl(str_sub(chk, 8 , 14), tp, fixed = TRUE)) {
    True_RF[i] <- "TRUE"
  }  else {
    True_RF[i] <- "FALSE"
  }
}

peptable$True_RF <- True_RF

Peptide <- character(n_rows)

idx <- 1
for (rown in seq_len(n_cry)) {
  
  jx <- as.integer(strsplit(cryEx_table$flanking_jxns[rown], ":")[[1]])
  strand <- cryEx_table$strand[rown]
  
  interval <- GRanges(
    seqnames = cryEx_table$chrom[rown],
    ranges = IRanges(
      start = c(jx[1]-flanklen+1, cryEx_table$start[rown], jx[2]),
      end   = c(jx[1], cryEx_table$end[rown], jx[2]+flanklen-1+1)
    ),
    strand = strand
  )
  
  dna <- if (strand == "+") {
    paste(getSeq(genome, interval), collapse = "")
  } else {
    paste(rev(getSeq(genome, interval)), collapse = "")
  }
  
  pep <- c(
    strsplit(as.character(translate(DNAString(dna))), "\\*")[[1]][1],
    strsplit(as.character(translate(DNAString(substr(dna,2,nchar(dna))))), "\\*")[[1]][1],
    strsplit(as.character(translate(DNAString(substr(dna,3,nchar(dna))))), "\\*")[[1]][1]
  )
  
  Peptide[idx:(idx+2)] <- pep
  idx <- idx + 3
}

peptable$Peptide <- Peptide
peptable$Length  <- nchar(peptable$Peptide)


#check right flanking exon aa sequence for novel readthroughs

peptable <- peptable %>% mutate(NovelRead_check = str_sub( Peptide , -14), .before = Novel_true )

Novel_true <- character(n_rows)


for (i in seq_len(n_rows)) {
  
  chk <- peptable$NovelRead_check[i]
  tp  <- peptable$True_pep[i]
  upstm <- peptable$RF_check[i]
  trslPep <- peptable$Peptide[i]
  
  if (is.na(chk) || is.na(tp) || chk == "") {
    Novel_true[i] <- "no reference aa"
  } else if (nchar(chk) < min_aa_len || grepl(chk, upstm, fixed = TRUE) || nchar(trslPep) <= flanklen/3*2-1) {
    Novel_true[i] <- "FALSE"
  } else if (grepl(chk, tp, fixed = TRUE)) {
    Novel_true[i] <- "TRUE"
  } else {
    Novel_true[i] <- "FALSE"
  }
}

peptable$Novel_true <- Novel_true

stopifnot(nrow(peptable) == length(peptable$Peptide))
stopifnot(nrow(peptable) == length(peptable$RF_check))
stopifnot(all(table(peptable$Frame) == n_cry))



peptable <- peptable %>%
  separate(Coords, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  mutate(cry_lengthDiv3 = (end + 1 - start)/3)   #cryptic length divided by 3 (+1 added to include last nuc)

#peptable = peptable[peptable$Length >= minpeptidelength,]
peptable$UniqueID=1:dim(peptable)[1] 

wb <-createWorkbook()
addWorksheet(wb, "peptable" )
writeData(wb, "peptable", peptable)
saveWorkbook(wb, file = outfile, overwrite = TRUE)

Yellow <- createStyle(fgFill = "yellow")
Blue <- createStyle(fgFill = "lightblue")

cryLeng.list <- as.list(peptable$cry_lengthDiv3)
Length.list <- as.list(peptable$Length)
Novel.list <- as.list(peptable$Novel_true)

lapply(seq_along(cryLeng.list), function(i){
  if(cryLeng.list[[i]] == Length.list[[i]] - (flanklen/3*2)){
    addStyle(wb, sheet = "peptable" ,style = Yellow, rows = i + 1, cols = 1:14 ,gridExpand = TRUE)
  }
})

lapply(seq_along(cryLeng.list), function(i){
  if(cryLeng.list[[i]] == Length.list[[i]] - (flanklen/3*2)-1){
    addStyle(wb, sheet = "peptable" ,style = Yellow, rows = i + 1, cols = 1:14 ,gridExpand = TRUE)
  }
})

lapply(seq_along(Novel.list), function(i){
  if(Novel.list[[i]] == "TRUE"){
    addStyle(wb, sheet = "peptable", style = Blue, rows = i + 1, cols = 15, gridExpand = TRUE)
  }
})
pep_filterd <- peptable %>% dplyr::filter(peptable$True_RF == "TRUE")

filtered_cryLeng.list <- as.list(pep_filterd$cry_lengthDiv3)
filtered_Length.list <- as.list(pep_filterd$Length)
filtered_Novel.list <- as.list(pep_filterd$Novel_true)

addWorksheet(wb, "peptable_filtered")
writeData(wb, "peptable_filtered", pep_filterd)

lapply(seq_along(filtered_cryLeng.list), function(i){
  if(filtered_cryLeng.list[[i]] == filtered_Length.list[[i]] - (flanklen/3*2)){
    addStyle(wb, sheet = "peptable_filtered" ,style = Yellow, rows = i + 1, cols = 1:14 ,gridExpand = TRUE)
  }
})


lapply(seq_along(filtered_cryLeng.list), function(i){
  if(filtered_cryLeng.list[[i]] == filtered_Length.list[[i]] - (flanklen/3*2-1)){
    addStyle(wb, sheet = "peptable_filtered" ,style = Yellow, rows = i + 1, cols = 1:14 ,gridExpand = TRUE)
  }
})

lapply(seq_along(filtered_Novel.list), function(i){
  if(filtered_Novel.list[[i]] == "TRUE"){
    addStyle(wb, sheet = "peptable_filtered", style = Blue, rows = i + 1, cols = 15, gridExpand = TRUE)
  }
})


saveWorkbook(wb, file = outfile, overwrite = TRUE)

