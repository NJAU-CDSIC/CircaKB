rm(list=ls())
# 1.Init parameters-------------------------------------------
library(jsonlite)
library(dplyr)
args=commandArgs(trailingOnly = T)
# Read parameters from the JSON file
params <- jsonlite::fromJSON(args[1])

dup=as.numeric(params$dup)# Duplicated sampling counts
sTime=as.numeric(params$startTime) # Starting sampling time
interval=as.numeric(params$interval)# Sampling interval
samplingCounts=as.numeric(params$samplingCount)# Sampling counts
idenCounts=as.numeric(params$idenCounts)# Identifier column number
nRow=as.numeric(params$geneNum)# gene number

dataMode=params$dataMode#  data mode

uGroup=params$singleGroup# Single groups
uGroup=uGroup[!duplicated(uGroup)]
uComGroup=params$compareGroup# Compare groups
uComGroup=uComGroup[!duplicated(uComGroup)]

tsMethods=as.numeric(params$singleMethods)# Circadian oscillation methods id
tsMethods=tsMethods[!duplicated(tsMethods)]
tdMethods=as.numeric(params$diffMethods)# Differential rhythmicity methods id
tdMethods=tdMethods[!duplicated(tdMethods)]

sMethods=data.frame(id = 1:4,methods = c("Meta2d","Cosinor", "Rain", "GeneCycle"))
dMethods=data.frame(id = 1:4,methods = c("diffCircadian", "CircaCompare", "LimoRhyde", "DODR"))

u_sMethods=sMethods[sMethods$id %in% tsMethods,]$methods
u_dMethods=dMethods[dMethods$id %in% tdMethods,]$methods


saveAddr=params$saveAddr
fileDir=params$fileDir



library(stringr)
preGroup=c()
if(length(uComGroup)==0){
  preGroup=uGroup
}else{
  sPre=str_split(uComGroup,"Vs.",simplify = T)
  
  for(i in 1:nrow(sPre)){
    preGroup=c(preGroup,sPre[i,])
  }
  preGroup=union(preGroup,uGroup)
}

cat("Resampling Count:",dup,"\n")
cat("Start Time:",sTime,"\n")
cat("Sampling Interval:",interval,"\n")
cat("Sampling Counts:",samplingCounts,"\n")
cat("idenCounts:",idenCounts,"\n")
cat("Data Mode:",dataMode,"\n")
cat("Gene Number:",nRow,"\n")

cat("Single Groups:",uGroup,"\n")
cat("Compare Groups:",uComGroup,"\n")
cat("Circadian Oscillation Models:",u_sMethods,"\n")
cat("Differential Rhythmicity Models:",u_dMethods,"\n")


cat("File Path:",fileDir,"\n")
cat("Save Path:",saveAddr,"\n")

# 2.Get data
dat=as.data.frame(read.csv(file=fileDir,header = T),stringAsFactors=F) 
# 2.1 Set data rows
if(nrow(dat)<nRow){
  nRow=nrow(dat)
}
dat=dat[c(1:nRow),]
# 2.2 Transfer data to numeric
dat=cbind(dat[,1:idenCounts],as.data.frame(lapply(dat[,(idenCounts+1):ncol(dat)],as.numeric)))

if(idenCounts>1){# Remove rows where all identifier columns are empty
  dat=dat[!(rowSums( dat[,1:idenCounts]=="" |is.na(dat[,1:idenCounts]) )==idenCounts),]
}else{
  dat=dat[!(dat[,1:idenCounts]=="" |is.na(dat[,1:idenCounts])),]
}
# 2.3 Identify the unique identifier column
ndat=dat
if(idenCounts==1){
  colnames(ndat)[1]="mark"
  
  if(length(names(table(duplicated(ndat[,1]))))==1 && names(table(duplicated(ndat[,1])))=="FALSE"){
  }else{# Merge duplicate entries
    ndat=aggregate(.~mark,mean,data=ndat)# Calculate the average value for duplicate genes
    ndat=ndat[order(ndat$mark),]
  }
}else{
  for(i in c(1:idenCounts)){
    m=table(duplicated(ndat[,c(1:idenCounts)]))
    if(length(names(m)=="TRUE")==1){
      oldIdenty=c(i,colnames(ndat)[i])
      colnames(ndat)[i]="mark"
      break
    }
  }
  if(!("mark" %in% colnames(ndat))){# Merge duplicate entries in the mark column
    colnames(ndat)[1]="mark"
    temp=ndat[,c(1,(idenCounts+1):ncol(ndat))]
    expr_mean=aggregate(.~mark,mean,data=temp)
    expr_mean=expr_mean[order(expr_mean$mark),]
    
    new_dat=ndat[!duplicated(ndat$mark),]
    new_dat=new_dat[order(new_dat$mark),]
    
    expr_mean=cbind(expr_mean[,1],new_dat[,2:idenCounts],expr_mean[,c(2:ncol(expr_mean))])
    
    colnames(expr_mean)[1:idenCounts]=colnames(new_dat)[1:idenCounts]
    ndat=expr_mean
  }
  
}
# 2.3 Split data according to the sampling mode
switch(dataMode,
       "mode1"={# No duplicate sampling
         start=idenCounts+1
         for(i in seq_along(preGroup)){
           temp=cbind(ndat[,1:idenCounts],ndat[,start:(start+samplingCounts-1)])
           assign(preGroup[i],temp)
           start=start+samplingCounts
         }
       },
       "mode2"={# CT00_dup1,CT04_dup1,...CT00_dup2
         start=idenCounts+1
         for(i in seq_along(preGroup)){
           temp=cbind(ndat[,1:idenCounts],ndat[,start:(start+samplingCounts*dup-1)])
           assign(preGroup[i],temp)
           start=start+samplingCounts*dup
         }
       },
       "mode3"={# CT00_dup1,CT00_dup2,CT04_dup1,...
         start=idenCounts+1
         for(i in seq_along(preGroup)){
           gtemp=ndat[,start:(start+samplingCounts*dup-1)]
           for (j in c(1:dup)) {
             s=seq(j,ncol(gtemp)-dup+j,dup)
             if(j==1){
               x=gtemp[,s]
             }else{
               x=cbind(x,gtemp[,s])
             }
           }
           temp=cbind(ndat[,1:idenCounts],x)
           assign(preGroup[i],temp)
           start=start+samplingCounts*dup
         }
       })

# 3.Data processing---------------------------------------------------------
for (i in seq_along(preGroup)) {
  temp=get(preGroup[i])
  temp=temp[order(temp$mark),]
  assign(preGroup[i],temp)
}
Times=as.numeric(seq(sTime,sTime+interval*(samplingCounts-1),by=interval))
rowName=get(preGroup[1])$mark
# geneList
if(idenCounts>1){
  geneList=as.data.frame(get(preGroup[1])[,1:idenCounts],stringAsFactors=F)
  colnames(geneList)[as.numeric(oldIdenty[1]) ]=oldIdenty[2]
}else{# idenCounts=1
  
  geneList=as_tibble(get(preGroup[1])[,1:idenCounts])
  colnames(geneList)=colnames(get(preGroup[1]))[1]
  colnames(geneList)[as.numeric(oldIdenty[1])]=oldIdenty[2]
}

for (i in seq_along(uGroup)) {
  temp=get(uGroup[i])
  temp=temp[order(temp$mark),]
  rownames(temp)=temp$mark
  temp=temp[,-1:-idenCounts]
  assign(uGroup[i],temp)
}

source("batchProcessing.R")
# Execute the Circadian Oscillation model

if(length(u_sMethods)>0 && length(uGroup)>0){
  sResult=rhyOscillation()

}
# Execute the differential rhythmicity model

if(length(u_dMethods)>0 && length(uComGroup)>0){
  dResult=diffRhythm()
}


