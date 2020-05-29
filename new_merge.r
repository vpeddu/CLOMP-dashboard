library(dplyr)

files<-list.files('/Users/vikas/Downloads/scratch/', pattern = '*.tsv', full.names = TRUE)

reports<-c()
for(i in 1:length(files)){ 
  reports[[i]]<-read.csv(files[i], sep = '\t')
}
reports


for(i in 1:length(files)){ 
  if( i == 1){ 
    final_df<-read.csv(files[i], sep = '\t', col.names = c('percent_covered', 'fragments_at_below_clade', 'fragments_at_clade', 'code', 'taxid', 'name'))
    final_df<-final_df[,c(3:6)]
    colnames(final_df)[1]<-basename(files[i])
  }
  else{ 
    temp_df<-read.csv(files[i], sep = '\t', col.names = c('percent_covered', 'fragments_at_below_clade', 'fragments_at_clade', 'code', 'taxid', 'name'))
    temp_df<-temp_df[,c(3:6)]
    colnames(temp_df)[1]<-basename(files[i])
    final_df<-merge(temp_df, final_df, by = c('taxid', 'name', 'code'), sort = F)  
      
    }
}
