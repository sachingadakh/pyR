library(reshape2) # for unique function
library(idr)

# to extract non overlapping data frame
notin  <- function(a1,a2)
{
  a1.vec <- apply(a1, 1, paste, collapse = "")
  a2.vec <- apply(a2, 1, paste, collapse = "")
  a1.without.a2.rows <- a1[!a1.vec %in% a2.vec,]
  return(a1.without.a2.rows)
}

extract_feature_ids <- function(sample_types)
{
  for(s_t in sample_types) {
    #sample_types_1 = list(x="count_C1R1", y="count_C1R3")
    #df <- genes_cpm[c(sample_types_1$x, sample_types_1$y)]
    df <- genes_cpm[c(s_t$x, s_t$y)]
    #print (colnames(df))
    cat(sprintf("Performing IDR analysis of %s, %s \n", s_t$x, s_t$y), file = stderr())
    cat(sprintf("Total number of features in samples %g \n", nrow(df)), file = stderr())
    #df = df[rowSums(df != 0) != 0,]
    #head(df)
    idr_out = est.IDR(df, mu = 1, sigma = 0.2, rho = 0.99, p = 0.7)
    selected = select.IDR(df, idr_out$idr, IDR.level = 0.1)
    filtered_count = nrow(df)-selected$n
    
    cat(sprintf("Total number of features with CPM >0 %g \n", nrow(df)), file = stderr())
    cat(sprintf("Total number of features with IDR <0.1 %g \n", selected$n), file = stderr())
    cat(sprintf("Total number of features with IDR >0.1 %g \n", filtered_count), file = stderr())
    
    p <- cor(df, method='pearson') #, use='na.or.complete'
    s <- cor(df, method='spearman') #, use='na.or.complete'
    
    cat(sprintf("Pearson Correlation of features is %f \n", p[1, 2]), file = stderr())
    cat(sprintf("Spearman Correlation of features is %f \n", s[1, 2]), file = stderr())
    
    filtered = selected$x
    #head(filtered)
    p <- cor(filtered, method='pearson') #, use='na.or.complete'
    s <- cor(filtered, method='spearman') #, use='na.or.complete'
    
    cat(sprintf("Pearson Correlation of filtered features is %f \n",   p[1, 2]), file = stderr())
    cat(sprintf("Spearman Correlation of filtered features is %f \n",   s[1, 2]), file = stderr())
    
    for (v in row.names(filtered))
      feature_ids <- c(feature_ids, v)
    cat(sprintf("Number of genes filtered %g \n", length(feature_ids), file = stderr()))
    
    feature_ids = as.list(feature_ids)
    if(length(feature_ids_main)) {
      feature_ids_main = Reduce(intersect, list(feature_ids_main, feature_ids))
    }else {
      feature_ids_main = feature_ids
    }
    cat(sprintf("Number of genes filtered totally %g \n", length(feature_ids_main), file = stderr()))
    write.table(feature_ids_main,file ="%s, %s",sep = "\t", row.names = T)
    filtered = as.data.frame(filtered)
    red_df<-notin(df,filtered)
    #head(red_df)
    if(nrow(red_df)>0) {
      red_df=as.data.frame(red_df)
      red_df= log2(red_df)
      is.na(red_df) <- sapply(red_df, is.infinite)
    }
    filtered=log2(filtered)
    is.na(filtered) <- sapply(filtered, is.infinite)
    #cat(sprintf("Number of genes filtered with IDR %g \n", nrow(filtered), file = stderr()))
    #cat(sprintf("Number of genes not filtered with IDR %g \n", nrow(red_df), file = stderr()))
    #R1_f = filtered[,sample_types_1$x]
    #R2_f = filtered[,sample_types_1$y]
    #xlabel <- paste(sample_types_1$x," CPM (log2)")
    #ylabel <- paste(sample_types_1$y," CPM (log2)")
    #head(R1_f)
    #head(filtered)
    R1_f <- filtered[,s_t$x]
    R2_f <- filtered[,s_t$y]
    xlabel <- paste(s_t$x," CPM (log2)")
    ylabel <- paste(s_t$y," CPM (log2)")
    plot(R1_f, R2_f, xlab=xlabel, ylab =ylabel, pch=".", cex=0.1, col='black', xlim =c(-5,15), ylim =c(-5,15))
    #plot(R1_f, R2_f, xlab='R1 CPM (log2)', ylab = 'R2 CPM (log2)', pch=".", cex=0.1, col='black', xlim = c(0,10), ylim=c(0,10))
    if(nrow(red_df)>0) {
      #R1_rdf = red_df[,sample_types_1$x]
      #R2_rdf = red_df[,sample_types_1$y]
      R1_rdf <- red_df[,s_t$x]
      R2_rdf <- red_df[,s_t$y]
      #plot(R1_rdf, R2_rdf, pch=".", cex=0.1, col='violet')
      points(R1_rdf, R2_rdf, pch=".", cex=0.1, col='violet')
    }
    print(plot)
    feature_ids <- vector(mode="character")
  }
  return(feature_ids_main)
}


genes_cpm <- read.table("genesCPM", sep="\t", header = T, row.names = 1)

# vector to store unique feature_ids
feature_ids <- vector(mode="character")
feature_ids_main <- list()
# to combine plots in three columns
par(mfrow=c(1,2))
#head(genes_cpm)
#SET COLUMN HEADERS AS PER FILE "genes_cpm_v1.txt" FOR PAIRWISE IDR ANALYSIS
sample_types_1 = list(x="C1R1", y="C1R2")
sample_types_2 = list(x="C2R1", y="C2R2")
sample_types_3 = list(x="C3R1", y="C3R2")
sample_types_4 = list(x="C4R1", y="C4R2")
sample_types = list(sample_types_1, sample_types_2, sample_types_3, sample_types_4)
feature_ids_1<-extract_feature_ids(sample_types)
#sample_types_1 = list(x="count_C3", y="count_C4")
#sample_types_2 = list(x="count_C2R2", y="count_C2R3")
#sample_types_3 = list(x="count_C2R1", y="count_C2R3")
#sample_types = list(sample_types_1)
#feature_ids_2<-extract_feature_ids(sample_types)
#feature_ids <- union(feature_ids_1, feature_ids_2)
#feature_ids =Reduce(intersect, list(feature_ids_1, feature_ids_2))
feature_ids = unique(feature_ids_1)
cat(sprintf("Number of genes filtered from all samples %g \n", length(feature_ids), file = stderr()))
write(unlist(feature_ids), file = "idr_filtered_numatData.txt", append = FALSE, sep = "\t", ncolumns=1)
