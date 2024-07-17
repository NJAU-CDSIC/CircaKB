library(dplyr)
library(tidyr)
library(fdrtool)
library(stringr)
#' @title Create Directory
#' @description This function checks if a directory exists at the specified address, and if it does not exist, it creates the directory, including any necessary but nonexistent parent directories.
#' @param addr A character string representing the path to the directory to be created.
#' @return The function returns the address of the created directory.
#' @examples
#' # Example usage of createDir
#' dir_path <- createDir("path/to/directory")
#' @export
createDir <-function(addr){
  
  if(! dir.exists(addr)){
    
    dir.create(addr,recursive = TRUE)
  }
  return (addr)
}

#' @title Project Acrophase Function
#' @description Adjusts the given acrophase angle to the range [-π, π]. This function ensures the acrophase angle is within the standard range for circular data analysis. If the input angle exceeds π or is less than -π, it adjusts it by adding or subtracting 2π accordingly.
#' @param acr Numeric value representing the acrophase angle in radians.
#' @return Numeric value representing the adjusted acrophase angle within the range [-π, π].
#' @examples
#' # Example usage of project_acr
#' acrophase <- project_acr(3*pi/2)
#' @export
project_acr<-function(acr){
  acr =acr %% (2*pi)
  if (acr > pi){
    acr =acr- 2*pi
  } 
  else if (acr < -pi){
    acr =acr+ 2*pi  
  }
  return (acr)
}

#' @title Convert Acrophase to Hours
#' @description Converts the acrophase angle (in radians) to hours based on the given period.This function calculates the corresponding hour value based on the acrophase angle and the period of the circadian rhythm.
#' @param acrophase Numeric value representing the acrophase angle in radians.
#' @param period Numeric value representing the period of the circadian rhythm in hours (default: 24).
#' @return Numeric value representing the hour value corresponding to the acrophase.
#' @examples
#' # Example usage of acrophase_to_hours
#' acrophase_hours <- acrophase_to_hours(1.5, period = 24)
#' @export
acrophase_to_hours<-function(acrophase,period=24){
  acrophase = project_acr(acrophase)
  hours = -period * acrophase/(2*pi)
  if (hours < 0){
    hours =hours+ period
  }
  return (hours)
}

#' @title Perform Cosinor Analysis
#'
#' @description Applies the Cosinor analysis to time-series data in each group provided. The function calculates
#' periodicity, amplitude, acrophase, and statistical significance using the cosinor2 package.
#'
#' @param group Character vector specifying the names of data frames containing time-series data.
#' @param dup Numeric value indicating duplication factor for time points.
#' @param Times Numeric vector representing time points or initial time.
#' @param interval Numeric value representing time interval between measurements.
#' @param saveAddr Character string specifying the directory path to save results.
#' @param geneList Character vector of gene identifiers.
#' @param nRow Numeric value indicating the number of rows in each data frame.
#' @param rowName Character vector specifying row names for result matrices.
#' @return Save results to CSV file.
#' @export
Cosinor <- function(group, dup, Times, interval, saveAddr, geneList, nRow, rowName){
  print("Cosinor")
  library(cosinor2)
  if(dup==0||dup==1){
    c_time=Times
  }else{
    c_time=seq(Times[1],(dup*length(Times)-1)*interval+Times[1],by=interval)
  }
  
  colName=c("p","q","amplitude","acrophase","mesor","amp.p","amp.q","acr.p","acr.q","acrophase_h","period")
  # Res=list()
  for (i in c(1:length(group))) {
    print(paste("Cosinor",sep=":",i))
    res=data.frame(matrix(ncol = 11, nrow = nRow,dimnames = list(rowName,colName)))
    Gtemp=get(group[i])
    aPer=c()
    for (j in c(1:nRow)) {
      
      temp_table=Gtemp[j,] %>%  gather("time","value")
      temp_table$time=c_time
      temp_table$value=as.numeric(temp_table$value)
      
      n0 <- apply(Gtemp[j,] == 0, 1, sum)# Count the number of zeros in each row
      n0=as.numeric(n0)
      if(n0>=ncol(Gtemp[j,])*0.3){
        besttime=0
        res[j,]=c(1,1,0,0,0,1,1,1,1,NA)#p,q,amplitude,acrophase,mesor,amp.p,amp.q,acr.p,acr.q,acrophase_h
      }else{
        print(j)
        peri=cosinor2::periodogram(data = Gtemp[j,], time = c_time,periods = 20:28,alpha=0.05)
        l=length(peri[["plot_env"]][["besttime"]])
        per=as.numeric(peri[["plot_env"]][["besttime"]][l])
        res[j,]$period=per
        aPer=c(aPer,per)
        fit.cosinor <- cosinor.lm(value ~ time(time), period = per,data =temp_table)
        
        res[j,1] <- cosinor.detect(fit.cosinor)[4]#p
        if(is.na(fit.cosinor$coefficients[2]) || is.na(fit.cosinor$coefficients[1])){
          res[j,3]=0
          res[j,4]=0
          res[j,10]=NA
        }else{
          res[j,3] <- fit.cosinor$coefficients[2]#amplitude
          res[j,4] <- 2*pi+correct.acrophase(fit.cosinor)#acrophase
          res[j,10]=acrophase_to_hours(res[j,4],per)#acrophase_h
        }
        
        s <- summary(fit.cosinor)$transformed.table
        res[j,5]=fit.cosinor$coefficients[1]#mesor
        if(is.na(s$p.value[2]) || is.na(s$p.value[3])){
          res[j,6]=1
          res[j,8]=1
        }else{
          res[j,6]<- s$p.value[2]#amp.p
          res[j,8] <- s$p.value[3]#acr.p
        }
        
      }
      
    }
    res$p[res$p>1]=0.9999
    res$p[res$p<0]=0.9999
    
    res$q=fdrtool(res$p, statistic="pvalue",plot = F,color.figure = F,verbose = F)$qval
    res$amp.q = fdrtool(res$amp.p, statistic="pvalue",plot = F,color.figure = F,verbose = F)$qval
    res$acr.q = fdrtool(res$acr.p, statistic="pvalue",plot = F,color.figure = F,verbose = F)$qval
    temLag=as.numeric(res$acrophase_h)
    modifyAcrophase_h=ifelse(temLag-Times[1]<0,temLag+aPer-Times[1],temLag-Times[1]) 
      
    RRres=data.frame(period=res$period,amp=res$amplitude,lag=modifyAcrophase_h,mesor=res$mesor,p=res$p,q=res$q,acrophase=res$acrophase,acrophase_h=res$acrophase_h)
    
    RRres=as.data.frame(lapply(RRres,as.character),stringsAsFactors=F)
    rownames(RRres)=rownames(res)
    RRres[is.na(RRres)]="NA"
    
    RRres=cbind(geneList,RRres)
    addr=createDir(paste(saveAddr,sep = "/",group[i]))
    write.csv(RRres,file = paste0(addr,"/","Cosinor",".csv"),row.names = F)
    
    # Res[[group[i]]]=RRres
  }
  print("The Cosinor algorithm has been completed!")
  return(T)
}

#' @title Perform Rain Analysis
#'
#' @description  Applies the Rain analysis to time-series data in each group provided. The function calculates
#' periodicity, lag, and statistical significance using the rain package.
#'
#' @param group Character vector specifying the names of data frames containing time-series data.
#' @param dup Numeric value indicating duplication factor for time points.
#' @param Times Numeric vector representing time points or initial time.
#' @param interval Numeric value representing time interval between measurements.
#' @param saveAddr Character string specifying the directory path to save results.
#' @param geneList Character vector of gene identifiers.
#' @param nRow Numeric value indicating the number of rows in each data frame.
#' @param rowName Character vector specifying row names for result matrices.
#' @return Save results to CSV file.
#' @export
Rain<-function(group, dup, Times, interval, saveAddr, geneList, nRow, rowName){
  print("Rain")
  library(rain)
  
  if(dup==0||dup==1){
    nr.series=1
  }else{
    nr.series=dup
    colName=c()
    for (i in c(1:dup)) {
      colName=c(colName,paste("ZT",Times,sep = "_",i))
    }
    for(i in seq_along(group)){
      temp=get(group[i])
      colnames(temp)=colName
      assign(group[i],temp)
    }
  }
  # Res=list()
  for (i in seq_along(group)) {
    Rtemp=rain(t(get(group[i])[c(1:nRow),]),deltat=interval,
               period = 24,
               method="independent",#independent: Multiple cycles are interpreted as repetitions of a single cycle
               period.delta=4,
               verbose = T,
               peak.border=c(0,1),
               nr.series=nr.series,
               adjp.method="ABH")
    colnames(Rtemp)[2]="Lag"
    RRtemp=data.frame(period=Rtemp$period,amp=NA,lag=Rtemp$Lag,mesor=NA,p=Rtemp$pVal,q=NA)
    RRtemp[is.na(RRtemp)]="NA"
    RRtemp=as.data.frame(lapply(RRtemp,as.character),stringsAsFactors=F)
    rownames(RRtemp)=rownames(Rtemp)
    
    RRtemp=cbind(geneList,RRtemp)

    addr=createDir(paste(saveAddr,sep = "/",group[i]))
    write.csv(RRtemp,file = paste0(addr,"/","Rain",".csv"),row.names = F)
    # Res[[group[i]]]=RRtemp
  }
  print("The Rain algorithm has been completed!")
  return(T)
  
}

#' @title Perform Meta2d Analysis
#'
#' @description Applies Meta2d analysis to time-series data in each group provided. The function calculates
#' periodicity, amplitude, lag, and statistical significance using MetaCycle package algorithms.
#'
#' @param group Character vector specifying the names of data frames containing time-series data.
#' @param dup Numeric value indicating duplication factor for time points.
#' @param Times Numeric vector representing time points or initial time.
#' @param interval Numeric value representing time interval between measurements.
#' @param saveAddr Character string specifying the directory path to save results.
#' @param geneList Character vector of gene identifiers.
#' @param nRow Numeric value indicating the number of rows in each data frame.
#' @param rowName Character vector specifying row names for result matrices.
#' @return Save results to CSV file.
#' @export
Meta2d <- function(group, dup, Times, interval, saveAddr, geneList, nRow, rowName){
  print("Meta2d")
  library(MetaCycle)
  if(dup==0||dup==1){
    c_time=Times
    dup=1
  }else{
    c_time=seq(Times[1],(dup*length(Times)-1)*interval+Times[1],by=interval)
  }
  for(i in seq_along(group)){
    temp=get(group[i])
    temp=cbind(rownames(temp),temp)
    colnames(temp)=c("mark",c_time)
    assign(group[i],temp)
  }
  # Res=list()
  # cores <- detectCores(logical=F)
  for(i in seq_along(group)){
    temp_res=meta2d(inDF=as.data.frame(get(group[i])),timepoints=as.numeric(c_time),
                    infile = "NULL",
                    minper = 20,
                    maxper = 28,
                    ARSdefaultPer = 24,
                    parallelize=F,
                    # nCores=cores-1,
                    weightedPerPha=T,
                    adjustPhase = "predictedPer",
                    filestyle="csv",outputFile = F)
    
    colnames(temp_res[["ARS"]])[7]="Lag"
    temp_res[["meta"]]=NULL
    
    if(is.null(temp_res[["JTK"]])){
      JTK=as.data.frame(matrix(data="NA",nrow = nrow(genelist),ncol = 5,dimnames = list(NULL,c("period","amp","lag","p","q"))),stringsAsFactors = F)
    }else{
      JTK=data.frame(period=temp_res[["JTK"]]$PER,amp=temp_res[["JTK"]]$AMP,lag=temp_res[["JTK"]]$LAG,mesor=NA,p=temp_res[["JTK"]]$ADJ.P,q=temp_res[["JTK"]]$BH.Q)
      JTK=as.data.frame(lapply(JTK,as.character),stringsAsFactors=F)
      JTK[is.na(JTK)]="NA"
      rownames(JTK)=temp_res[["JTK"]]$CycID
    }
    
    JTK=cbind(geneList,JTK)
    
    addr=createDir(paste(saveAddr,sep = "/",group[i]))
    write.csv(JTK,file = paste0(addr,"/","JTK_CYCLE",".csv"),row.names = F)
    # temp_res[["JTK"]]=JTK
    
    if(is.null(temp_res[["ARS"]])){
      ARS=as.data.frame(matrix(data="NA",nrow = nrow(genelist),ncol = 5,dimnames = list(NULL,c("period","amp","lag","p","q"))),stringsAsFactors = F)
    }else{
      lg=as.numeric(temp_res[["ARS"]]$Lag)
      lp=as.numeric(temp_res[["ARS"]]$period)
      tempLag=ifelse(lg-Times[1]<0,
                     lg+lp-Times[1],
                     lg-Times[1])
      ARS=data.frame(period=temp_res[["ARS"]]$period,amp=temp_res[["ARS"]]$amplitude,lag=temp_res[["ARS"]]$Lag,mesor=temp_res[["ARS"]]$mean,p=temp_res[["ARS"]]$pvalue,q=temp_res[["ARS"]]$fdr_BH)
      ARS=as.data.frame(lapply(ARS,as.character),stringsAsFactors=F)
      ARS[is.na(ARS)]="NA"
      rownames(ARS)=temp_res[["ARS"]]$CycID
    }
    ARS=cbind(geneList,ARS)
    
    addr=createDir(paste(saveAddr,sep = "/",group[i]))
    write.csv(ARS,file = paste0(addr,"/","ARSER",".csv"),row.names = F)

    # temp_res[["ARS"]]=ARS
    
    if(is.null(temp_res[["LS"]])){
      LS=as.data.frame(matrix(data="NA",nrow = nrow(genelist),ncol = 5,dimnames = list(NULL,c("period","amp","lag","p","q"))),stringsAsFactors = F)
    }else{
      lg=as.numeric(temp_res[["LS"]]$PhaseShift)
      lg=lg/dup
      lp=as.numeric(temp_res[["LS"]]$Period)
      lg=lg-Times[1]
      tempLag=ifelse(lg<0 | is.na(lg) |lg=="",NA,
                     ifelse(lg-lp>0,lg-lp,lg))
      LS=data.frame(period=temp_res[["LS"]]$Period,amp=NA,lag=tempLag,mesor=NA,p=temp_res[["LS"]]$p,q=temp_res[["LS"]]$BH.Q)
      LS=as.data.frame(lapply(LS,as.character),stringsAsFactors=F)
      LS[is.na(LS)]="NA"
      rownames(LS)=temp_res[["LS"]]$CycID
    }
    LS=cbind(geneList,LS)
    
    addr=createDir(paste(saveAddr,sep = "/",group[i]))
    write.csv(LS,file = paste0(addr,"/","Lomb-Scargle",".csv"),row.names = F)

    # temp_res[["LS"]]=LS
    
    # Res[[group[i]]]=temp_res
  }
  print("The JTK_CYCLE algorithm has been completed!")
  print("The ARSER algorithm has been completed!")
  print("The Lomb-scargle algorithm has been completed!")
  return(T)
}

#' @title Perform differential expression analysis using Fisher's G-test and Robust G-test methods.
#'
#' @description Applies GeneCycle analysis to time-series data in each group provided. The function calculates
#' periodicity, and statistical significance using Fisher's G-test and Robust G-test.
#'
#' @param group Character vector specifying the names of data frames containing time-series data.
#' @param dup Numeric value indicating duplication factor for time points.
#' @param Times Numeric vector representing time points or initial time.
#' @param interval Numeric value representing time interval between measurements.
#' @param saveAddr Character string specifying the directory path to save results.
#' @param geneList Character vector of gene identifiers.
#' @param nRow Numeric value indicating the number of rows in each data frame.
#' @param rowName Character vector specifying row names for result matrices.
#' @return Save results to CSV file.
#' @export
GeneCycle <- function(group, dup, Times, interval, saveAddr, geneList, nRow, rowName){
  print("GeneCycle")
  library(GeneCycle)
  colName=c("p","q")
  if(dup==0||dup==1){
    c_time=Times
  }else{
    c_time=seq(Times[1],(dup*length(Times)-1)*interval+Times[1],by=interval)
  }
  for(i in seq_along(group)){
    temp=get(group[i])[c(1:nRow),]
    colnames(temp)=c_time
    assign(group[i],as.data.frame(t(temp)))
  }
  
  # Res=list()
  #fisher.g.test
  for(i in seq_along(group)){
    res=data.frame(matrix(ncol = 2, nrow = nRow,dimnames = list(rowName,colName)))
    temp_res=fisher.g.test(get(group[i]))
    fdr.out = fdrtool(temp_res, statistic="pvalue",plot = F,color.figure = F,verbose = F)
    res[,1]=fdr.out$pval
    res[,2]=fdr.out$qval
    
    RRtemp=data.frame(period=NA,amp=NA,lag=NA,mesor=NA,p=res$p,q=res$q)
    RRtemp=as.data.frame(lapply(RRtemp,as.character),stringsAsFactors=F)
    rownames(RRtemp)=rownames(res)
    RRtemp[is.na(RRtemp)]="NA"
    fisher=list(RRtemp)
    fisher <- setNames(fisher, "fisher")
    
    fisher=cbind(geneList,fisher)
    
    addr=createDir(paste(saveAddr,sep = "/",group[i]))
    write.csv(fisher,file = paste0(addr,"/","Fisher's G-test",".csv"),row.names = F)

    # Res[[group[i]]]=fisher
    
  }
  print("The Fisher's G-test algorithm has been completed!")
  #robust.g.test
  for (i in seq_along(group)) {
    res=data.frame(matrix(ncol = 2, nrow = nRow,dimnames = list(rowName,colName)))
    spe5 = robust.spectrum(get(group[i]))
    pval = robust.g.test(spe5)  
    tempfile=paste("g_pop_length_",nrow(spe5),sep = "",".txt")
    unlink(tempfile)
    
    num.na=is.na(pval)
    num.overrange=pval>1
    pval[num.na]=0
    pval[num.overrange]=1
    
    fdr.out = fdrtool(pval, statistic="pvalue",plot = F,color.figure = F,verbose = F)
    fdr.out$pval[num.na]=NaN
    fdr.out$qval[num.na]=NaN
    res[,1]=fdr.out$pval
    res[,2]=fdr.out$qval
    
    RRtemp=data.frame(period=NA,amp=NA,lag=NA,mesor=NA,p=res$p,q=res$q)
    RRtemp=as.data.frame(lapply(RRtemp,as.character),stringsAsFactors=F)
    rownames(RRtemp)=rownames(res)
    RRtemp[is.na(RRtemp)]="NA"
    robust=list(RRtemp)
    robust <- setNames(robust, "robust")
    
    robust=cbind(geneList,robust)
    
    addr=createDir(paste(saveAddr,sep = "/",group[i]))
    write.csv(robust,file = paste0(addr,"/","Robust G-test",".csv"),row.names = F)
    
    # Res[[group[i]]]=append(Res[[group[i]]],robust)
  }
  print("The Robust G-test algorithm has been completed!")
  return(T)
}

#' @title Perform differential expression analysis using HANOVA and robust DODR methods.
#'
#' @description This function conducts differential expression analysis on pre-processed data using the HANOVA and robust DODR methods. It transposes dataframes to match specified time points, performs statistical tests, adjusts p-values using the Benjamini-Hochberg method, and saves results to CSV files.
#'
#' @param PreUcomgroup List of dataframes containing expression data for each group.
#' @param dup Integer, number of replicates.
#' @param Times Numeric vector, time points for measurements.
#' @param interval Numeric, interval between time points.
#' @param saveAddr Character, directory path to save results.
#' @param geneList Character dataframe, list of gene names or identifiers.
#' @param nRow Integer, number of rows (genes) to process.
#' @param rowName Character vector, row names for the output dataframe.
#' @param cmb Dataframe specifying contrasts between groups for analysis.
#'
#' @return TRUE if the analysis completes successfully.
#'
#' @export
DODR <-function(PreUcomgroup, dup, Times, interval, saveAddr, geneList, nRow, rowName,cmb){
  print("DODR")
  library(DODR)
  colName=c("p","q")
  if(dup==0||dup==1){
    c_time=Times
  }else{
    c_time=seq(Times[1],(dup*length(Times)-1)*interval+Times[1],by=interval)
  }
  
  for(i in seq_along(PreUcomgroup)){
    temp=get(PreUcomgroup[i])
    colnames(temp)=c_time
    assign(PreUcomgroup[i],as.data.frame(t(temp)))
  }
  
  
  #HANOVA & robustDODR
  # Res=list()
  for (i in c(1:nrow(cmb))) {
    res_HANOVA=data.frame(matrix(ncol = 2, nrow = nRow,dimnames = list(rowName,colName)))
    res_robustDODR=data.frame(matrix(ncol = 2, nrow = nRow,dimnames = list(rowName,colName)))
    
    temp=dodr(get(cmb[i,1]),get(cmb[i,2]),c_time,c_time,24,method = c("HANOVA","robustDODR"))
    
    HANOVA.p.value=temp[["p.value.table"]]$HANOVA
    HANOVA.q.value=p.adjust(HANOVA.p.value, method = "BH")
    res_HANOVA[,1]=HANOVA.p.value
    res_HANOVA[,2]=HANOVA.q.value
    res_HANOVA=as.data.frame(lapply(res_HANOVA,as.character),stringsAsFactors=F)
    
    res_HANOVA=cbind(geneList,res_HANOVA)
    addr=createDir(paste(saveAddr,sep = "/",cmb$contrast[i]))
    write.csv(res_HANOVA,file = paste0(addr,"/","HANOVA",".csv"),row.names = F)
    
    # Res[["HANOVA"]][[cmb[i,]$contrast]]=res_HANOVA
    
    robustDODR.p.value=temp[["p.value.table"]]$robustDODR
    robustDODR.q.value=p.adjust(robustDODR.p.value, method = "BH")
    res_robustDODR[,1]=robustDODR.p.value
    res_robustDODR[,2]=robustDODR.q.value
    res_robustDODR=as.data.frame(lapply(res_robustDODR,as.character),stringsAsFactors=F)
    
    res_robustDODR=cbind(geneList,res_robustDODR)
    addr=createDir(paste(saveAddr,sep = "/",cmb$contrast[i]))
    write.csv(res_robustDODR,file = paste0(addr,"/","robust DODR",".csv"),row.names = F)
    # Res[["robustDODR"]][[cmb[i,]$contrast]]=res_robustDODR
  }
  print("The HANOVA algorithm has been completed!")
  print("The robust DODR algorithm has been completed!")
  return(T)
  
}

#' @title Perform comparative analysis using CircaCompare method.
#'
#' @description This function conducts comparative analysis between two groups using the CircaCompare method.
#' It calculates differential measures such as amplitude, phase, and mesor differences,
#' adjusts p-values, and saves results to CSV files.
#'
#' @param PreUcomgroup List of dataframes containing expression data for each group.
#' @param dup Integer, number of replicates.
#' @param Times Numeric vector, time points for measurements.
#' @param interval Numeric, interval between time points.
#' @param saveAddr Character, directory path to save results.
#' @param geneList Character dataframe, list of gene names or identifiers.
#' @param nRow Integer, number of rows (genes) to process.
#' @param rowName Character vector, row names for the output dataframe.
#' @param cmb Dataframe specifying contrasts between groups for analysis.
#'
#' @return TRUE if the analysis completes successfully.
CircaCompare <-function(PreUcomgroup, dup, Times, interval, saveAddr, geneList, nRow, rowName,cmb){
  library(circacompare)
  library(tidyr)
  print("Circacompare")
  if(dup==0||dup==1){
    c_time=Times
  }else{
    c_time=seq(Times[1],(dup*length(Times)-1)*interval+Times[1],by=interval)
  }
  
  for(i in seq_along(PreUcomgroup)){
    temp=get(PreUcomgroup[i])
    colnames(temp)=c_time
    assign(PreUcomgroup[i],as.data.frame(t(temp)))
  }
  result_list <- list()
  
  for (i in 1:nrow(cmb)) {
    dat1 <- get(cmb[i, 1])
    dat2 <- get(cmb[i, 2])
    
    # Initialize a temporary result matrix
    temp_res <- as.data.frame(matrix(NA, nrow = nRow, ncol = 15,
                                     dimnames = list(rowName,
                                                     c(cmb[i, ]$V2, cmb[i, ]$V1,
                                                       paste0(cmb[i, ]$V2, "_mesor"), paste0(cmb[i, ]$V1, "_mesor"), "logdiff_mesor", "diff_mesor.p",
                                                       paste0(cmb[i, ]$V2, "_amp"), paste0(cmb[i, ]$V1, "_amp"), "logdiff_amp", "diff_amp.p",
                                                       paste0(cmb[i, ]$V2, "_phase"), paste0(cmb[i, ]$V1, "_phase"), "logdiff_phase", "diff_phase.p", "Reserve"))),stringAsFactors=F)
    
    temp_res$diff_mesor.p <- 1
    temp_res$diff_amp.p <- 1
    temp_res$diff_phase.p <- 1
    
    # Loop through each row index j
    for (j in 1:nRow) {
      df1 <- data.frame(time = c_time, measure = dat1[, j], group = cmb[i, 1])
      df2 <- data.frame(time = c_time, measure = dat2[, j], group = cmb[i, 2])
      df <- rbind(df1, df2)
      
      out <- try(circacompare(x = df, col_time = "time", col_group = "group", col_outcome = "measure", alpha_threshold = 0.1), silent = TRUE)
      
      if (class(out) == "try-error") {
        print("out is N/A!")
      } else {
        temp <- out[[2]]
        temp_res[j, ]$Reserve <- -temp[13, ]$value
        
        temp[5, ]$value <- log2(temp[4, ]$value / temp[3, ]$value)
        temp[9, ]$value <- log2(temp[8, ]$value / temp[7, ]$value)
        
        # The difference between phases cannot exceed half a cycle
        if (abs(temp[12, ]$value - temp[11, ]$value) > 12) {
          if (temp[12, ]$value > temp[11, ]$value) {
            temp[12, ]$value <- abs(24 - temp[12, ]$value)
          } else {
            temp[11, ]$value <- abs(24 - temp[11, ]$value)
          }
        }
        
        temp[13, ]$value <- log2(temp[12, ]$value / temp[11, ]$value) * (temp_res[j, ]$Reserve / abs(temp_res[j, ]$Reserve))
        
        temp_res[j, 1:14] <- temp[c(1:14), 2]
      }
    }
    
    rownames(temp_res) <- rowName
    result_list[[i]] <- temp_res
  }
  
  Res=list()
  for(i in c(1:length(result_list))){
    temp=result_list[[i]]
    ntemp=as.data.frame(lapply(temp,as.character),stringsAsFactors=F)
    rownames(ntemp)=rownames(temp)
    RRtemp=data.frame(Amp.p=ntemp$diff_amp.p,Phase.p=ntemp$diff_phase.p,Mesor.p=ntemp$diff_mesor.p,
                      Amp.diff=ntemp$logdiff_amp,Phase.diff=ntemp$Reserve,Mesor.diff=ntemp$logdiff_mesor)
    rownames(RRtemp)=rownames(ntemp)
    RRtemp[is.na(RRtemp)]="NA"
    
    RRtemp=cbind(geneList,RRtemp)
    addr=createDir(paste(saveAddr,sep = "/",cmb$contrast[i]))
    write.csv(RRtemp,file = paste0(addr,"/","CircaCompare",".csv"),row.names = F)
    # 
    # Res[[cmb$contrast[i]]]=RRtemp
  }
  print("The CircaCompare algorithm has been completed!")
  return(T)
}

#' @title Perform differential circadian analysis using the diffCircadian method.
#'
#' @description  This function conducts differential circadian analysis between two groups using the diffCircadian method.
#' It calculates differential measures such as amplitude, phase, mesor, and fit differences,
#' adjusts p-values, and saves results to CSV files.
#'
#' @param PreUcomgroup List of dataframes containing expression data for each group.
#' @param dup Integer, number of replicates.
#' @param Times Numeric vector, time points for measurements.
#' @param interval Numeric, interval between time points.
#' @param saveAddr Character, directory path to save results.
#' @param geneList Character dataframe, list of gene names or identifiers.
#' @param nRow Integer, number of rows (genes) to process.
#' @param rowName Character vector, row names for the output dataframe.
#' @param cmb Dataframe specifying contrasts between groups for analysis.
#'
#' @return TRUE if the analysis completes successfully.
diffCircadian<-function(PreUcomgroup, dup, Times, interval, saveAddr, geneList, nRow, rowName,cmb){
  library(tidyr)
  library(nloptr)
  library(diffCircadian)
  if(dup==0||dup==1){
    c_time=Times
  }else{
    c_time=seq(Times[1],(dup*length(Times)-1)*interval+Times[1],by=interval)
  }
  
  for(i in seq_along(PreUcomgroup)){
    temp=get(PreUcomgroup[i])
    colnames(temp)=c_time
    assign(PreUcomgroup[i],as.data.frame(temp))
  }
  result_list <- list()

  for (i in 1:nrow(cmb)) {
    dat1 <- as.matrix(get(cmb[i, 1]))
    dat2 <- as.matrix(get(cmb[i, 2]))
    
    temp_res <- as.data.frame(matrix(NA, nrow = nRow, ncol = 16, 
                                     dimnames = list(rowName, 
                                                     c(paste0(cmb[i, ]$V1, "_mesor"), paste0(cmb[i, ]$V2, "_mesor"), "diff_mesor", "diff_mesor.p",
                                                       paste0(cmb[i, ]$V1, "_amp"), paste0(cmb[i, ]$V2, "_amp"), "diff_amp", "diff_amp.p",
                                                       paste0(cmb[i, ]$V1, "_phase"), paste0(cmb[i, ]$V2, "_phase"), "diff_phase", "diff_phase.p",
                                                       paste0(cmb[i, ]$V1, "_fit"), paste0(cmb[i, ]$V2, "_fit"), "diff_fit", "diff_fit.p"))),stringAsFactors=F) 
    temp_res$diff_mesor.p <- 1
    temp_res$diff_amp.p <- 1
    temp_res$diff_phase.p <- 1
    temp_res$diff_fit.p <- 1
    
    for (j in 1:nRow) {
      adiff <- try(LR_diff(c_time, dat1[j, ], c_time, dat2[j, ], type = "all", period = 24, method = "LR", FN = TRUE), silent = TRUE)
      if (class(adiff) == "try-error") {
        print("out is N/A!")
      } else {
        tres <- as.data.frame(unlist(adiff, recursive = FALSE))
        # tres$phase_1 <- ifelse(tres$phase_1 > 0, tres$phase_1, tres$phase_1 + 24)
        # tres$phase_2 <- ifelse(tres$phase_2 > 0, tres$phase_2, tres$phase_2 + 24)
        tres$logdiff_amp <- log2(tres$amp_2 / tres$amp_1)
        tres$diff_phase <- tres$phase_2 - tres$phase_1
        tres$logdiff_mesor <- log2(tres$offset_2 / tres$offset_1)
        tres$logdiff_fit <- log2(tres$sigma2_2 / tres$sigma2_1)
        
        temp_res[j, ] <- c(tres$offset_1, tres$offset_2, tres$logdiff_mesor, tres$pvalue.2,
                           tres$amp_1, tres$amp_2, tres$logdiff_amp, tres$pvalue,
                           tres$phase_1, tres$phase_2, tres$diff_phase, tres$pvalue.1,
                           tres$sigma2_1, tres$sigma2_2, tres$logdiff_fit, tres$pvalue.3)
      }
    }
    
    
    temp_res=as.data.frame(lapply(temp_res, as.character),stringsAsFactors = F)
    rownames(temp_res) <- rowName
    result_list[[i]] <- temp_res
  }
  
  # Res=list()
  for(i in c(1:length(result_list))){
    ntemp=result_list[[i]]
    
    RRtemp=data.frame(Amp.p=ntemp$diff_amp.p,Phase.p=ntemp$diff_phase.p,Mesor.p=ntemp$diff_mesor.p,Fit.p=ntemp$diff_fit.p,
                      Amp.diff=ntemp$diff_amp,Phase.diff=ntemp$diff_phase,Mesor.diff=ntemp$diff_mesor,Fit.diff=ntemp$diff_fit)
    rownames(RRtemp)=rownames(ntemp)
    
    RRtemp=cbind(geneList,RRtemp)
    addr=createDir(paste(saveAddr,sep = "/",cmb$contrast[i]))
    write.csv(RRtemp,file = paste0(addr,"/","diffCircadian",".csv"),row.names = F)
    
    # Res[[cmb$contrast[i]]]=RRtemp
  }
  print("The diffCircadian algorithm has been completed!")
  return(T)
}

#' @title Perform differential expression analysis using LimoRhyde method.
#'
#' @description This function conducts differential expression analysis between two conditions using the LimoRhyde method.
#' It calculates statistical measures such as p-values and adjusted p-values (q-values),
#' adjusts for periodicity in time series data using limorhyde, and saves results to CSV files.
#' @param PreUcomgroup List of dataframes containing expression data for each group.
#' @param dup Integer, number of replicates.
#' @param Times Numeric vector, time points for measurements.
#' @param interval Numeric, interval between time points.
#' @param saveAddr Character, directory path to save results.
#' @param geneList Character dataframe, list of gene names or identifiers.
#' @param nRow Integer, number of rows (genes) to process.
#' @param rowName Character vector, row names for the output dataframe.
#' @param cmb Dataframe specifying contrasts between groups for analysis.
#'
#' @return TRUE if the analysis completes successfully.
LimoRhyde<-function(PreUcomgroup, dup, Times, interval, saveAddr, geneList, nRow, rowName,cmb){
  library(limorhyde)
  library(limma)
  if(dup==0||dup==1){
    c_time=Times
  }else{
    c_time=rep(Times,dup)
  }
  
  
  # Res=list()
  for (i in c(1:nrow(cmb))) {
    dat=cbind(get(cmb[i,]$V1),get(cmb[i,]$V2))
    
    m_time1=data.frame(cond=cmb[i,]$V1,time=as.numeric(c_time))
    m_time2=data.frame(cond=cmb[i,]$V2,time=as.numeric(c_time))
    m_time=rbind(m_time1,m_time2)
    m_time=cbind(m_time,limorhyde(m_time$time, 'time_'))
    
    
    design = model.matrix(~ cond * (time_cos + time_sin), data = m_time)
    
    fit = lmFit(dat, design)
    fit = eBayes(fit, trend = TRUE)
    drLimma = topTable(fit, coef = 5:6, number = Inf)
    drLimma=cbind(rownames(drLimma),drLimma)
    colnames(drLimma)[1]="gene_id"
    rownames(drLimma)=c(1:nrow(drLimma))
    
    res_limorhyde=data.frame(p=drLimma$P.Value,q=drLimma$adj.P.Val)
    rownames(res_limorhyde)=drLimma$gene_id
    res_limorhyde[is.na(res_limorhyde)]="NA"
    res_limorhyde=as.data.frame(lapply(res_limorhyde,as.character),stringsAsFactors=F)
    
    res_limorhyde=cbind(geneList,res_limorhyde)
    addr=createDir(paste(saveAddr,sep = "/",cmb$contrast[i]))
    write.csv(res_limorhyde,file = paste0(addr,"/","LimoRhyde",".csv"),row.names = F)
    # Res[[cmb$contrast[i]]]=res_limorhyde
  }
  print("The LimoRhyde algorithm has been completed!")
  return(T)
  
}






