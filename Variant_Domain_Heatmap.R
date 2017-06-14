#import packages
library(readr)
args<-commandArgs(TRUE)

#read in Range_Finder_ex_Files and extract info
df <- read_csv("~/Desktop/Salk_Internship/Variant_Finder/ex_TFT_ALL.csv")

prot_dom = unique(df[["Protein Domain"]]) 
#print(prot_dom)
prot_dom <- prot_dom[-1]
conseqs = unique(df[["consequence"]])
print(conseqs[5])
conseqs <- conseqs[-length(conseqs)-1]
print((conseqs))

#correlate protein domains and consequences
corr <- data.frame(table(df$`Protein Domain`, df$consequence)[,])
corr <- corr[!(corr$Var2=="[]"),]

#label corr
colnames(corr)[1] <- "Protein_domain"
colnames(corr)[2] <- "Variant_type"
#make personal color for heatmap
my_palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 120)
#make matrix of corr
corr_matrix <- data.matrix(corr)
corr_matrix = matrix(corr$Freq,nrow=length(conseqs),ncol=length(prot_dom))
dimnames(corr_matrix) = list(c(conseqs),c(prot_dom))
corr_matrix <- t(corr_matrix)
corr_heatmap <- heatmap(corr_matrix, 
                        Rowv=NA, 
                        Colv=NA, 
                        col = my_palette, 
                        scale="column",
                        cexCol=1,
                        cexRow=0.5,
                        margins=c(10,10))
print(corr_matrix)
write.table(corr_matrix,
            file="heatmap_corr_matrix.csv", 
            sep=",", 
            quote="", 
            col.names=TRUE, 
            row.names=TRUE)

dev.off()

