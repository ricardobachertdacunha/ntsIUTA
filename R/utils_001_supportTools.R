

#' @title getScreeningListTemplate
#'
#' @param projPath The project folder location. Default is \code{setup$projPath}.
#'
#' @return Pastes a template .csv file of the screeningList into the projPath.
#' 
#' @export
#'
#' @examples
#' 
#' 
#' 
getScreeningListTemplate <- function(projPath = setup$projPath) {
  base::file.copy(from = base::paste0(base::system.file(package = "ntsIUTA", dir = "extdata"),"/screeningList_template.csv"),
                  to = setup$projPath,
                  overwrite = FALSE)
}



#' @title getColors
#'
#' @param x An \linkS4class{OnDiskMSnExp} object with one or more files or the number of colors to be produced.
#' @param which Possible entries are \code{groups} and \code{samples} for getting group or samples colors, respectively.
#'
#' @return A vector of colors names according to the given \linkS4class{OnDiskMSnExp} object.
#' 
#' @export
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr count
#'
#' @examples
#' 
#' 
#' 
getColors <- function(x, which = c("groups","samples")) {
  
  #Examples
  # getColors(ntsIUTA::rawDataExample, "groups")
  
  
  colors <- c(RColorBrewer::brewer.pal(8, "Greys")[6],
              RColorBrewer::brewer.pal(8, "Greens")[6],
              RColorBrewer::brewer.pal(8, "Blues")[6],
              RColorBrewer::brewer.pal(8, "Oranges")[6],
              RColorBrewer::brewer.pal(8, "Purples")[6],
              RColorBrewer::brewer.pal(8, "PuRd")[6],
              RColorBrewer::brewer.pal(8, "YlOrRd")[6],
              RColorBrewer::brewer.pal(8, "PuBuGn")[6],
              RColorBrewer::brewer.pal(8, "GnBu")[6],
              RColorBrewer::brewer.pal(8, "BuPu")[6],
              RColorBrewer::brewer.pal(8, "Dark2"))
  
  if (!base::class(x) == "numeric" & !base::class(x) == "integer")
  {
    numberOfGroups <- base::length(base::unique(x$sample_group))
    
    if (numberOfGroups > 18)
    {
      require(grDevices)
      colors <- grDevices::colorRampPalette(colors)(numberOfGroups)
    }
    group_colors <- colors[1:numberOfGroups]
    base::names(group_colors) <- base::unique(x$sample_group)
    count <- dplyr::count(x@phenoData@data, sample_group)
    sample_group_colors <- base::rep(group_colors, times = count[,"n"])
    base::names(sample_group_colors) <- base::unique(x$sample_name)
    
    if (which == "groups") return(group_colors)
    if (which == "samples") return(sample_group_colors)
    
  } else {
    numberOfGroups <- x
    
    if (numberOfGroups > 18)
    {
      require(grDevices)
      colors <- grDevices::colorRampPalette(colors)(numberOfGroups)
    }
    colors <- colors[1:numberOfGroups]
    
    return(colors)
    
  }
  
}

