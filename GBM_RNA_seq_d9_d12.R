library(dplyr)
library(patchwork)
library(stringr)
library(readr)
library(ggplot2)
library(readxl)

library(pheatmap)
library(grid)

#################Data Processing########################3
#Read Excel Files Week9
Ct_CD_df<- read_excel("C:/Users/willi/OneDrive/Desktop/Chandramohan Lab/GBM Data/Differential_expression_analysis_table_AB_IPA.xlsx")
colnames(Ct_CD_df)[2]<-"logFC_CD"
colnames(Ct_CD_df)[3]<-"pval_CD"
colnames(Ct_CD_df)[4]<-"padj_CD"
View(Ct_CD_df)

Ct_D2C7_df<- read_excel("C:/Users/willi/OneDrive/Desktop/Chandramohan Lab/GBM Data/Differential_expression_analysis_table_AC_IPA.xlsx")
colnames(Ct_D2C7_df)[2]<-"logFC_D2C7"
colnames(Ct_D2C7_df)[3]<-"pval_D2C7"
colnames(Ct_D2C7_df)[4]<-"padj_D2C7"
View(Ct_D2C7_df)

Ct_D_C_df<- read_excel("C:/Users/willi/OneDrive/Desktop/Chandramohan Lab/GBM Data/Differential_expression_analysis_table_AD_IPA.xlsx")
colnames(Ct_D_C_df)[2]<-"logFC_D_C"
colnames(Ct_D_C_df)[3]<-"pval_D_C"
colnames(Ct_D_C_df)[4]<-"padj_D_C"
View(Ct_D_C_df)

#Merge to include all treatments
Intermediate_df<- merge(Ct_CD_df, Ct_D2C7_df,by="ID")
Wk9_df<-merge(Intermediate_df,Ct_D_C_df,BY="ID")

View(Wk9_df)
dim(Wk9_df)

#Convert ENSEMBL to Gene name
library(org.Mm.eg.db)
annots <- select(org.Mm.eg.db, keys=as.vector(Wk9_df[,"ID"]), columns="SYMBOL", keytype="ENSEMBL")
GBM_wk9_df <- merge(Wk9_df, annots, by.x="ID", by.y="ENSEMBL")
write.csv(GBM_wk9_df,'C:/Users/willi/OneDrive/Desktop/Chandramohan Lab/GBM_DEG_wk9.csv', row.names = FALSE)

###################Load Data Table for gene expression###########################3
GBM_wk9_df<-read.csv('C:/Users/willi/OneDrive/Desktop/Chandramohan Lab/GBM_DEG_wk9.csv')
GBM_wk12_df<-read.csv('C:/Users/willi/OneDrive/Desktop/Chandramohan Lab/GBM_DEG_wk12.csv')
View(GBM_wk9_df)


#Function to generate heatmap dataframe
Generate_exp_df <- function(df,GOI_list){
  GOI_df<-df[df$SYMBOL %in% GOI_list,]
  row.names(GOI_df)<-GOI_df$SYMBOL
  GOI_df["Vehicle"]<-rep(1,length(GOI_df$SYMBOL))
  GOI_df["CD40"]<- 2**GOI_df$logFC_CD
  GOI_df["D2C7"]<- 2**GOI_df$logFC_D2C7
  GOI_df["D+C"]<- 2**GOI_df$logFC_D_C
  Exp_Table<-GOI_df[c("Vehicle","CD40","D2C7","D+C")]
  return(Exp_Table)
}

Generate_sig_df <- function(df,GOI_list){
  GOI_df<-df[df$SYMBOL %in% GOI_list,]
  row.names(GOI_df)<-GOI_df$SYMBOL
  GOI_df["Vehicle"]<-rep(1,length(GOI_df$SYMBOL))
  GOI_df["CD40"]<- 2**GOI_df$logFC_CD
  GOI_df["D2C7"]<- 2**GOI_df$logFC_D2C7
  GOI_df["D+C"]<- 2**GOI_df$logFC_D_C
  Sig_Table<-GOI_df[c("pval_CD","pval_D2C7","pval_D_C")]
  Sig_Table[is.na(Sig_Table)] <- 1
  return(Sig_Table)
}

#Function to generate heatmap with sig label
Generate_heatmap <- function(Exp_df){
  heat <- t(scale(t(as.matrix(Exp_df))))
  pheatmap(heat,cluster_rows = T,cluster_cols = F,border_color = NA, cellwidth = 80, cellheight = 16)
  
}

length(gene_list)

dim(Wk9_Adh_df)[1]

######################################Expression Table + Heatmap###############################


#Adhesion Molecules
adh_gene_list<-c("Icam1","Icam2","Vcam1","Itga4","Ackr1","Pecam1","Ninj1","Sell","Sele","Selp","Lyve1",
                 "Vegfa","Vegfb","Cd99l2","Cav1","Cxcl12","Lama4","Alcam","F11r","Jam2","Jaml","Ccl19","Itgb1","Abcb1",'Itgal')

Wk9_Adh_df<-Generate_exp_df(GBM_wk9_df,adh_gene_list)
Wk9_Adh_Sig_Table<- Generate_sig_df(GBM_wk9_df,adh_gene_list)
Generate_heatmap(Wk9_Adh_df)
grid.text("Heatmap of Adhesion Molecule Gene Exp from d9" , x = 0.5, y = 0.9)

Wk12_Adh_df<-Generate_exp_df(GBM_wk12_df,adh_gene_list)
Wk12_Adh_Sig_Table<- Generate_sig_df(GBM_wk12_df,adh_gene_list)
Generate_heatmap(Wk12_Adh_df)
grid.text("Heatmap of Adhesion Molecule Gene Exp from d12" , x = 0.5, y = 0.9)

#Arrest Chemokines
Che_gene_list<-c("Cxcl1","Cxcl2","Cxcl8","Cxcl9","Cxcl10","Cxcl12","Ccl3","Ccl5","Ccl11","Ccl19","Ccl21","Ackr1","Darc","Rap1a")

Wk9_Che_df<-Generate_exp_df(GBM_wk9_df,Che_gene_list)
Wk9_Che_Sig_Table<- Generate_sig_df(GBM_wk9_df,Che_gene_list)
Generate_heatmap(Wk9_Che_df)
grid.text("Heatmap of Chemokine Molecule Gene Exp from Week 9" , x = 0.5, y = 0.75)

Wk12_Che_df<-Generate_exp_df(GBM_wk12_df,Che_gene_list)
Wk12_Che_Sig_Table<- Generate_sig_df(GBM_wk12_df,Che_gene_list)
Generate_heatmap(Wk12_Che_df)
grid.text("Heatmap of Chemokine Molecule Gene Exp from Week 12" , x = 0.5, y = 0.75)


#Diapdesis paracellular
Diapdesis_Paracellular_gene_list<-c('Cx43', 'Cldn5', 'Esam', 'Iqgap1', 'JamA', 'JamB', 'JamC', 'Jaml', 'Lck', 'Mmp2', 'OcLN',
                                    'Pecam1', 'Tiam1', 'Timp1', 'Timp4', 'Zo1', 'Zo2','Zo3')


Wk9_Diapdesis_Paracellular_df<-Generate_exp_df(GBM_wk9_df,Diapdesis_Paracellular_gene_list)
Wk9_Diapdesis_Paracellular_Sig_Table<- Generate_sig_df(GBM_wk9_df,Diapdesis_Paracellular_gene_list)
Generate_heatmap(Wk9_Diapdesis_Paracellular_df)
grid.text("Heatmap of Diapdesis Paracellular Gene Exp from Week 9" , x = 0.5, y = 0.75)

Wk12_Diapdesis_Paracellular_df<-Generate_exp_df(GBM_wk12_df,Diapdesis_Paracellular_gene_list)
Wk12_Diapdesis_Paracellular_Sig_Table<- Generate_sig_df(GBM_wk12_df,Diapdesis_Paracellular_gene_list)
Generate_heatmap(Wk12_Diapdesis_Paracellular_df)
grid.text("Heatmap of Diapdesis Paracellular Gene Exp from Week 12" , x = 0.5, y = 0.75)



#Diapdesis transcellular
Diapdesis_Trans_gene_list<-c('Cav1', 'Cav2', 'Cavin1','Plvap', 'Snap23', 'Snap25', 'Vamp1', 'Vamp2', 'Vamp3','Vamp8')


Wk9_Diapdesis_Trans_df<-Generate_exp_df(GBM_wk9_df,Diapdesis_Trans_gene_list)
Wk9_Diapdesis_Trans_Sig_Table<- Generate_sig_df(GBM_wk9_df,Diapdesis_Trans_gene_list)
Generate_heatmap(Wk9_Diapdesis_Trans_df)
grid.text("Heatmap of Diapdesis Transcellular Gene Exp from Week 9" , x = 0.5, y = 0.75)

Wk12_Diapdesis_Trans_df<-Generate_exp_df(GBM_wk12_df,Diapdesis_Trans_gene_list)
Wk12_Diapdesis_Trans_Sig_Table<- Generate_sig_df(GBM_wk12_df,Diapdesis_Trans_gene_list)
Generate_heatmap(Wk12_Diapdesis_Trans_df)
grid.text("Heatmap of Diapdesis Transcellular Gene Exp from Week 12" , x = 0.5, y = 0.75)


#BBB Barrier


BBB_gene_list<-c('Zo1', 'Zo2', 'Cldn5', 'Marveld2', 'Ocln' ,
                 'Jam2', 'Cdh5', 'Pecam1',
                 'Gja1', 'Gja2' ,'Gja5',
                 'Plvap', 'Mfsd2a', 'Abca2')
Wk9_BBB_df<-Generate_exp_df(GBM_wk9_df,BBB_gene_list)
Wk9_BBB_Sig_Table<- Generate_sig_df(GBM_wk9_df,BBB_gene_list)
Generate_heatmap(Wk9_BBB_df)
grid.text("Heatmap of BBB Molecule Gene Exp from Week 9" , x = 0.5, y = 0.8)

Wk12_BBB_df<-Generate_exp_df(GBM_wk12_df,BBB_gene_list)
Wk12_BBB_Sig_Table<- Generate_sig_df(GBM_wk12_df,BBB_gene_list)
Generate_heatmap(Wk12_BBB_df)
grid.text("Heatmap of BBB Molecule Gene Exp from Week 12" , x = 0.5, y = 0.8)


#Clock genes
Clock_gene_list<-c('Clock', 'Arntl', 'Bmal1', 'Per1', 'Per2', 'Cry1', 'Rora', 'Rorb', 'Rorc', 'Nr1d1')

Wk9_Clock_df<-Generate_exp_df(GBM_wk9_df,Clock_gene_list)
Wk9_Clock_Sig_Table<- Generate_sig_df(GBM_wk9_df,Clock_gene_list)
Generate_heatmap(Wk9_Clock_df)
grid.text("Heatmap of Clock Gene Exp from Day 9" , x = 0.5, y = 0.9)

Wk12_Clock_df<-Generate_exp_df(GBM_wk12_df,Clock_gene_list)
Wk12_Clock_Sig_Table<- Generate_sig_df(GBM_wk12_df,Clock_gene_list)
Generate_heatmap(Wk12_Clock_df)
grid.text("Heatmap of Clock Gene Exp from Day 12" , x = 0.5, y = 0.9)


#Microglia Genes
Microglia_genes <- c('Il1b', 'Retsdr8', 'Spp1', 'Folr2', 'Nlrp3', 'Lyve1', 'Apoe')

D9_Microglia_df<-Generate_exp_df(GBM_wk9_df, Microglia_genes)
D9_Microglia_Sig_Table<- Generate_sig_df(GBM_wk9_df, Microglia_genes)
Generate_heatmap(D9_Microglia_df)
grid.text("Heatmap Day 9" , x = 0.5, y = 0.8)

D12_Microglia_df<-Generate_exp_df(GBM_wk12_df, Microglia_genes)
D12_Microglia_Sig_Table<- Generate_sig_df(GBM_wk12_df, Microglia_genes)
Generate_heatmap(D12_Microglia_df)
grid.text("Heatmap Day 12" , x = 0.5, y = 0.8)

#NLRP3 Inflammasome pathway
Inflammasome_genes <- c('Il1b', 'Nfkb1', 'Nlrp3', 'Casp1', 'Pycard', 'Tlr4', 'Myd88','Bpifd2','Casp8','Il18')

D9_Inflammasome_df<-Generate_exp_df(GBM_wk9_df, Inflammasome_genes)
D9_Inflammasome_Sig_Table<- Generate_sig_df(GBM_wk9_df, Inflammasome_genes)
Generate_heatmap(D9_Inflammasome_df)
grid.text("Heatmap Day 9" , x = 0.5, y = 0.8)

D12_Inflammasome_df<-Generate_exp_df(GBM_wk12_df, Inflammasome_genes)
D12_Inflammasome_Sig_Table<- Generate_sig_df(GBM_wk12_df, Inflammasome_genes)
Generate_heatmap(D12_Inflammasome_df)
grid.text("Heatmap Day 12" , x = 0.5, y = 0.8)



#####################################Heatmap#####################################
library(pheatmap)
library(grid)

#Z scale
heat <- t(scale(t(as.matrix(Wk9_Exp_Table))))
heat
#Heatmap
pheatmap(heat,cluster_rows = T,cluster_cols = F,border_color = NA, cellwidth = 80, cellheight = 16)
grid.text("Heatmap of Adhesion Molecule Gene Exp from Week 9" , x = 0.5, y = 0.98)

Generate_sig_heatmap(Wk9_Che_df)

Wk9_Che_Sig_Table
##################################Label Significance#########################################3
library(grid)

#Define X and Y value
x<-c(233.5,313.5,393.5)
y<-rep(1,22)
for (i in 1:22){y[i]<-420-i*16}


#For loop to put * on significant observations
for (i in 1:22){
  for  (j in 1:3){
    if(as.double(Wk9_Sig_Table[i,j])<0.01){
      grid.text("*" , x = unit(x[j], "point"), y = unit(y[i], "point"))
      print(c(i,j))
    }
  }
}






