# Script for identifying rows whihc overlap with a set range (cCOMP)
# Developed for using on metabolomics dataset where features min and max peak areas were compared those found in the 
# blake wash sample. The porblem was that there was a significant carry over of the products from sample to sample and
# to remove those features from samples where they shoudl be, John Chodkowski imagined to remove those if they are,
# outside their range when detected in blank wash sample.


install.packages('DescTools')
library(DescTools)

x=as.numeric(as.character(c(1,2,9,1,3,7)))
y=c(10,8,15,3,4,20)
id=c(1,1,1,2,2,2)
g=c('a','e','b','a','e','b')

df=data.frame(cbind(id,x,y,g))
df$x=as.numeric(as.character(df$x))
df$y=as.numeric(as.character(df$y))

df_result=NULL
output=NULL

for(f in unique(id)){
  
  df_sub=df[df$id == f,]
  e_val<- df_sub[df_sub$g == 'e',]
  first=e_val$x
  sec=e_val$y
  cCOMP=c(first, sec)
  
  for(i in c('a','e','b')){
    
    df_sub_plus=df_sub[df_sub$g==i,]
    c1=c(df_sub_plus$x,df_sub_plus$y)
    output=(c1 %overlaps% cCOMP)
    output=cbind(output,i,f)
    df_result=rbind(df_result,output)
    
  }
}

