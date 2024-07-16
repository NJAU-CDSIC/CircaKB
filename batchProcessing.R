source("Methods.R")

#' @title Perform single-group circadian rhythm analysis using specified methods.
#'
#' @description This function executes circadian rhythm analysis using multiple methods specified in u_sMethods. It handles errors gracefully and returns a list of results from each method.
#'
#' @return NULL, or an error message if any method fails.
rhyOscillation <-function(){
  cat("Executing Circadian Oscillation Analysis---------------------------\n")
  # tGroup=paste0("t_",uGroup)
  result_list <- mclapply(u_sMethods, function(name) {
    tryCatch(
      do.call(name,list(group=uGroup, dup = dup, Times=Times, interval=interval, saveAddr=saveAddr, geneList=geneList, nRow=nRow, rowName=rowName)),
      error = function(e) {
        message("An error occurred with function:", name, "Error message:", e$message)
        e$message
      }
    )
  })
  # result_list=setNames(result_list, u_sMethods)
  return(result_list)
}

#' @title Perform multi-group differential rhythmicity analysis using specified methods.
#'
#' @description This function executes differential rhythmicity analysis using multiple methods specified in u_dMethods. It handles errors gracefully and returns a list of results from each method.
#'
#' @return NULL, or an error message if any method fails.
diffRhythm <- function(){
    cat("Executing Differential Rhythmicity Analysis---------------------------\n")
    PreUcomgroup=c()
    sp=str_split(uComGroup,"Vs.",simplify = T)
    for(i in c(1:nrow(sp))){
      PreUcomgroup = union(PreUcomgroup,sp[i,])
    }
    
    cmb = as.data.frame(matrix(nrow=nrow(sp),ncol=3,dimnames = list(NULL, c("V1","V2","contrast"))))
    for(i in c(1:nrow(cmb))){
      cmb[i,1]=sp[i,1]
      cmb[i,2]=sp[i,2]
      cmb[i,3]=paste(sp[i,1],sep = "Vs.",sp[i,2])
    }
    cmb = cmb

    Com_list <- mclapply(u_dMethods, function(name) {
      tryCatch(
        do.call(name, list(PreUcomgroup=PreUcomgroup,dup = dup, Times=Times, interval=interval, 
                           saveAddr=saveAddr, geneList=geneList, nRow=nRow, rowName=rowName, cmb=cmb)),
        error = function(e) {
          message("An error occurred with function:", name, "Error message:", e$message)
          e$message
        }
      )
    })
    
    Com_list <<- setNames(Com_list, u_dMethods)
    
}


