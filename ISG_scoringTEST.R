#need to convert ENSMBL to gene symbols for ISG scoring

ISG_sig <- read.delim("~/Downloads/ISG_manual.gmt.txt", header=FALSE)

# formating into a vector
ISG_sig <- as.data.frame(t(ISG_sig))
ISG_sig <- ISG_sig$V1[-grep('ISG_sig|http:',ISG_sig$V1)]

#########
#make expr_mat matrix
########
MATR3_seq = read.xlsx("/media/bree/Expansion/RNASeqAnalysis_GB_LM2_MATR3/DESeq2_GB_LM2_MATR3/LM2_MATR3_DESeq2.xlsx", sheet = "DESeq2 output")

#remove NAs in gene symbol and TPM columns
MATR3_seq = MATR3_seq[-which(is.na(MATR3_seq$Gene)),]
MATR3_seq = MATR3_seq[-which(is.na(MATR3_seq$GB_LM2_siCtrl_Rep_1_Aligned.sortedByCoord.out.bam_)),]
# MATR3_seq = MATR3_seq[-which(is.na(MATR3_seq$GB_LM2_siCtrl_Rep_2_Aligned.sortedByCoord.out.bam_)),]
# MATR3_seq = MATR3_seq[-which(is.na(MATR3_seq$GB_LM2_siMATR3_Rep_1_Aligned.sortedByCoord.out.bam_)),]
# MATR3_seq = MATR3_seq[-which(is.na(MATR3_seq$GB_LM2_siMATR3_Rep_2_Aligned.sortedByCoord.out.bam_)),]

expr_mat = list(
  gene_symbol = MATR3_seq$Gene,
  L2f_TPM_ctrl1 = log2(MATR3_seq$GB_LM2_siCtrl_Rep_1_Aligned.sortedByCoord.out.bam_),
  L2f_TPM_ctrl2 = log2(MATR3_seq$GB_LM2_siCtrl_Rep_2_Aligned.sortedByCoord.out.bam_),
  L2f_TPM_KD1 = log2(MATR3_seq$GB_LM2_siMATR3_Rep_1_Aligned.sortedByCoord.out.bam_),
  L2f_TPM_KD2 = log2(MATR3_seq$GB_LM2_siMATR3_Rep_2_Aligned.sortedByCoord.out.bam_),
  ENSMBL = MATR3_seq$ENSEMBL
)

#remove na rowname
# sum(is.na(expr_mat$gene_symbol))
# expr_mat <- expr_mat[-which(is.na(expr_mat$gene_symbol)), ]

# Set gene symbols as rownames
expr_mat <- as.data.frame(expr_mat) #rearanged
dupes = expr_mat$gene_symbol[duplicated(expr_mat$gene_symbol)]
fix_symbol = c()
for (i in 1:length(expr_mat$gene_symbol)) {
  ENSMBL = expr_mat$ENSMBL[i]
  Symbol = expr_mat$gene_symbol[i]
  
  if(Symbol %in% dupes){
     fix_symbol[i] = paste0(Symbol, ";", ENSMBL) 
  }else{
    fix_symbol[i] = Symbol
  }
}
  
expr_mat$gene_symbol = fix_symbol
rownames(expr_mat) <- expr_mat$gene_symbol
expr_mat$gene_symbol <- NULL
expr_mat$ENSMBL <- NULL

# Convert to numeric matrix
expr_mat <- as.matrix(expr_mat)
mode(expr_mat) <- "numeric"

gene_sets <- list(ISG_sig = ISG_sig)

library(GSVA)

param <- gsvaParam(
  exprData = expr_mat,
  geneSets = gene_sets,
  kcdf = "Gaussian",   # correct for log2 TPM)
)

gsva_scores <- gsva(param)
gsva_scores <- as.data.frame(t(gsva_scores))
gsva_scores$tumor_biospecimen_id <- rownames(gsva_scores)