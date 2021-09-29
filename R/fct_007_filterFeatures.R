
#' filterFeatures
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param ... A sequence of arguments comprising of the filters to applied, each followed by the respective threshold value/values.
#'
#' @return
#'
#' @export
#'
#' @importMethodsFrom patRoon as.data.frame 
#' 
#' @examples
#'

filterFeatures <- function(obj, 
                           filterList = list(filterMinInt = 2000, filterBlank = 3)) {
  
  if (length(filterList) == 0) {
    warning("No filters selected. Please select at least one filter.")
    return(obj)
  }
  
  listOfViableFilters <- c("filterMinInt", "filterBlank")
  
   filters <- names(filterList)
   
    for(i in 1:length(filters)) {
      
      switch(names(filterList)[i],
             filterMinInt = (obj <- filterMinInt(obj, unlist(filterList[i]))),
             filterBlank = (obj <- filterBlank(obj, unlist(filterList[i]))) )
      
    }
 
  #dots <- as.data.frame(list(...))
  #return(dots)
   return(obj)
   
}




tdots <- filterFeatures("test", "test2", list(5,7), "test3", 1000, "test4", 500)

ttdots <- filterFeatures("test", list("test2", 7), list("test3", 1000), list("test4", 500))


tdots
tdots[[1]][1]
class(tdots[1])
unlist(tdots[[1]])


# feats2_f3 <- feats2_f3[order(feats2_f3$mz),]