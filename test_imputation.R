---
title: "test_imputation"
output: html_notebook
---



```{r libraries et blabla}
####objectif : imputation sur donnees manquantes et test de probabilite pour tester la qualitÃ© de l'imputation
#methode : sauvegarder le dataset de base (on garde les filtres prealables MAF et NA), remplacer NA sur des donnees connues, sauvegarde ce dataset avec NA, faire une imputation sur ces donnees, test le rÃ©sulat de l'imputation entre dataset_donnees_reelles et dataset_donnees_imputed. Calcul du pourcentage d'erreur de l'imputation : entre valeur data_real et data_imputed, meme valeur ou non ? TEst sur 1000 inferences

##libraries
library(data.table)
library(LEA)
library(dplyr)

#####on utilise 3 datasets : 
##data_raw (donnes non modifiees), 
##data_WK (on va supprimer les valeurs 0 ou 1 en retenant les positions row, col)
##data_imputed: on utilise le fichier data_WK (avec NA), on change le fichier pour un format LFMM, puis on impute les donnees (cf scipt GWAS )

##attention a avoir la meme taille de fichier (206 IP etudiees): sauvegarder les chr, pos, ref dans un autre fichier. FAire la vérif pour chaque fichier créer : data_imputed, data_raw_ et data_WK
```


```{r prepare dataset}
#setwd(dir = "C:/Users/2021lg003/Documents/")
#genotype<-"C:/Users/2021lg003/Documents/Pearl millet panel PMIGAP/chr_association_vcf/chr2_association_vcf.txt"
genotype<-"C/Users\2021lg003\Documents\Pearl millet panel PMIGAP\chr_association_vcf/chr_association_all_vcf.txt"
data_raw<-fread(genotype, header=TRUE)

#save chr, pos, ref
data_ref<-data_raw[,1:3]
##raw data (et seulement data : le dataset LFMM des donnees imputees ne contient pas de metadonnees)
data_raw<-data_raw[,-c(1:3)]
##creer le doc de travail et lire le dataset
data_raw <- data_raw[,lapply(.SD,as.numeric)]#, by=c("#CHROM","POS", "REF")]  ## all numeric except chrom & pos 

data_WK<-fread(genotype, header=TRUE)
#sapply(data_WK, class)   
#convert character to num
data_WK <- data_WK[,lapply(.SD,as.numeric)]#, by=c("#CHROM","POS", "REF")]  ## all numeric except chrom & pos 
data_WK_save<-data_WK

#output data_save
data_save<- "C:/Users/2021lg003/Documents/R/PMIGAP/test_imputation/data_save_chr2.txt"
data_save<-as.data.frame(t(c("#CHROM","POS","REF","row_nb","col_nb", "IP")))
colnames(data_save)<-data_save[1,]
data_save<-data_save[-c(1),]

##une premier position au hasard pour lancer la boucle 
row_nb<-sample(nrow(data_WK), 1)
col_nb<-sample(ncol(data_WK[4:ncol(data_WK)]), size = 1)
```


```{r boucle remplace donnees par NA}
##boucle 1000 inferences 
l <- 0  #count nb valeurs traités
while (l<1000) { 
##prendre au hasard dans data_WK une valeur connue 0 ou 1 
##selctionner row
#row_nb<-sample(nrow(data_WK), 1)
row_nb<-as.vector(row_nb)
row_WK<-data_WK[row_nb,]  

##selectionner col
#col_nb<-sample(ncol(data_WK), size = 1)
col_nb<-as.vector(col_nb)
row_WK[[col_nb]] 

data_save<-as.data.frame(c(row_WK[,1:3], row_nb, col_nb, row_WK[[col_nb]]), col.names = c("#CHROM","POS","REF","row_nb","col_nb", "IP"))

#Si  valeur ==! NA, sauvegarder la valeur dans un nouveau dataset data_save
## apprend passe a la ligne suivante
  if (is.na(data_save$IP) == FALSE) {
    data_save<-write.table(c(row_WK[,1:3], row_nb, col_nb, row_WK[[col_nb]]), file = "data_save_chr2.txt",col.names = !file.exists("data_save_chr2.txt"), row.names = F, append = TRUE) #### a modifier : le header est repetÃ© car numero IP change, essayer de dÃ©caler l'IP par colonne

  #####remplacer NA sur les memes donnees connues de data_WK
  data_WK[[row_nb,col_nb]]<-NA

  } 

print(data_WK[[row_nb,col_nb]])
row_nb<-sample(nrow(data_WK), 1)
col_nb<-sample(ncol(data_WK), size = 1)
l<-l+1
print(l)
}

####fin boucle : on a un dataset data_WK avec NA + data_raw avec les valeurs a comparer
```


```{r lfmm et imputation}
#### Keep only usefull information and replace missing data code (from "." to "9")

# vcf_haplo[,3:ncol(vcf_haplo)] <- apply(vcf_haplo[,3:ncol(vcf_haplo)],2,as.numeric)
#data_WK <- data_WK[,lapply(.SD,as.numeric),by=c("#CHROM","POS")]

#faire une imputation sur ces donnees de data_WK (meme methode que script GWAS African Rice)
####lfmm file 
genotype <- data_WK[,4:ncol(data_WK)]
genotype <- t(genotype)
write.lfmm(R = genotype,output.file = "genotype_chr2_imputation.lfmm")
# imputation (using 1000 reps, might be more), K choice to be based on PCA, snmf,...
obj.snmf = snmf(input.file="genotype_chr2_imputation.lfmm", K=1:15, entropy=T, ploidy = 1, project="new", rep=10)
plot(obj.snmf) ### K value 
best = which.min(cross.entropy(obj.snmf, K = 4))
impute(object = obj.snmf, input.file="genotype_chr2_imputation.lfmm", K = 4, run = best)
#################################"

#test le rÃ©sulat de l'imputation entre data_WK et data_raw 
##Calcul du pourcentage d'erreur de l'imputation
data_imputed <- read.lfmm("genotype_chr5_imputation.lfmm_imputed.lfmm")
data_imputed <- t(data_imputed) ##var=col=IP et row=obs=SNPs
data_imputed<-as.data.frame(data_imputed)
colnames(data_imputed)<-colnames(data_raw[4:ncol(data_raw)]) ###??
```


```{r faire sauvegardes, rajouter CHROM, POS, REF }
##change data_raw pour qu'il corresponde au format 
#data_raw_save <- data_raw[,lapply(.SD,as.numeric)]
data_WK_save<-data_WK
data_imputed_save<-data_imputed
data_raw_save<-data_raw
#rajouter CHROM, POS, REF pour chaque dataset 3 col CHROM, POS, REF, + IP
data_raw<-c(data_ref, data_raw_save)
data_imputed<-c(data_ref, data_imputed_save)
data_imputed<-as.data.frame(data_imputed)
data_raw<-as.data.frame(data_raw)

## attention bien reprendre les memes positions row et col : sauvegarder les positions dans data_save
```


```{r }
inputFile<- "C:/Users/2021lg003/Documents/R/PMIGAP/test_imputation/data_save_chr5.txt"
data_save <- file(inputFile, open = "r")

#data_save<-read.table("data_save.txt")
#colnames(data_save)<-c("#CHROM","POS","REF","row_nb","col_nb", "IP")

i<-1
while (length(oneLine <- readLines(data_save, n = 1, warn = FALSE)) > 0) {

   oneLine<-unlist(strsplit(oneLine, split=" "))
    row_nb<-oneLine[4]
    col_nb<-oneLine[5]
    col_nb<-as.integer(col_nb)
    row_nb<-as.integer(row_nb)
  oneLine
  save_count<-data_raw[[row_nb,col_nb]]== data_imputed[[row_nb,col_nb]]
  write.table(save_count, file = "count_chr5_raw_imputed.txt", col.names = F, append = TRUE)
  i<-i+1
  print(i)
}
##nb FALSE/nb total*100 : pourcentage d'erreur imputation
count_df<-read.table("count_chr5_raw_imputed.txt")
count_df<-as.character(count_df$V2)
   percent<-length(which(count_df == "TRUE"))/(length(which(count_df == "FALSE"))+length(which(count_df == "TRUE")))
  
###########
```



