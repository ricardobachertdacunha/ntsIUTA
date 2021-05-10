

#' @title getInstParam
#' @description Tranforms the list of parameters for peak picking (PP), alignment (\emph{i.e.} retention time adjustment across samples)
#' and grouping into the respetive class objects.
#' As input a list as follows should be given where for each data processing step the input should be given as described in
#' the documentation of the aimed function of the \pkg{xcms} package.
#' \itemize{
#'   \item \code{paramList <- list(PP = list("..."), preGrouping = list("..."), alignment = list("..."), grouping = list("..."))}
#' }
#' \code{preGrouping} is used for \code{\link[xcms]{adjustRtime}} function
#' when using the method \code{PeakGroups} as it requires peaks groups (\emph{i.e.} features) for alignment.
#' Note that currently for PP only \code{MassifquantParam} is possible, for alignment only \code{PeakGroupsParam} is possible
#' and for grouping only \code{PeakDensityParam}. Other methods will be integrated in the future.
#' 
#' @param paramList A list of parameters for PP, alignment and grouping as described above.
#' See \code{\link{paramListExample}} object for a structure example.
#'
#' @return A \code{list} with parameters for PP, aligment and grouping to be passed to \code{\link{peakPicking}} and \code{\link{makeFeatures}}. 
#' 
#'  @references
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#' 
#' @export
#'
#' @examples
#' paramList <- ntsIUTA::paramListExample
#' getInstParam(paramList)
#' 
getInstParam <- function(paramList = paramList) {
  
  instParam <- list()
  
  instParam$PP <- NULL
  
  # Make instParam for PP
  if(!is.null(paramList$PP)) {
    
    if(paramList$PP$funcName == class(xcms::MassifquantParam())) {
      
      instParam$PP <- xcms::MassifquantParam(ppm = paramList$PP$ppm,
                                             peakwidth = paramList$PP$peakwidth,
                                             snthresh = paramList$PP$snthresh,
                                             prefilter = paramList$PP$prefilter,
                                             mzCenterFun = paramList$PP$mzCenterFun,
                                             integrate = paramList$PP$integrate,
                                             mzdiff = paramList$PP$mzdiff,
                                             fitgauss = paramList$PP$fitgauss,
                                             noise = paramList$PP$noise,
                                             criticalValue = paramList$PP$criticalValue,
                                             consecMissedLimit = paramList$PP$consecMissedLimit,
                                             unions = paramList$PP$unions,
                                             checkBack = paramList$PP$checkBack,
                                             withWave = paramList$PP$verboseColumns,
                                             verboseColumns = paramList$PP$verboseColumns)
      
    }
    
    #TODO Add parameters for other methods
    
  }
  
  instParam$preGrouping <- NULL
  instParam$alignment <- NULL
  
  # Make instParam for alignment, including preGrouping for PeakGroupsParam method
  if(!is.null(paramList$alignment)) {
    
    if(paramList$alignment$funcName == class(xcms::PeakGroupsParam())) {
      
      instParam$preGrouping <-  xcms::PeakDensityParam(sampleGroups = "placeHolder",
                                                       bw = paramList$preGrouping$bw,
                                                       minFraction = paramList$preGrouping$minFraction,
                                                       minSamples = paramList$preGrouping$minSamples,
                                                       binSize = paramList$preGrouping$binSize,
                                                       maxFeatures = paramList$preGrouping$maxFeatures)
      
      instParam$alignment <- xcms::PeakGroupsParam(minFraction = paramList$alignment$minFraction,
                                                   extraPeaks = paramList$alignment$extraPeaks,
                                                   smooth = paramList$alignment$smooth,
                                                   span = paramList$alignment$span,
                                                   family = paramList$alignment$family,
                                                   peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
                                                   subset = integer(),
                                                   subsetAdjust = "average")
      
    }
    
    #Add parameters for other methods
    
  }
  
  instParam$grouping <- NULL
  
  #Make instParam for grouping.
  if(!is.null(paramList$grouping)) {
    
    if(paramList$grouping$funcName == class(xcms::PeakDensityParam(sampleGroups = "placeHolder"))) {
      
      instParam$grouping <- PeakDensityParam(sampleGroups = "placeHolder",
                                             bw = paramList$grouping$bw,
                                             minFraction = paramList$grouping$minFraction,
                                             minSamples = paramList$grouping$minSamples,
                                             binSize = paramList$grouping$binSize,
                                             maxFeatures = paramList$grouping$maxFeatures)
      
    }
    
    #Add parameters for other methods
    
  }
  
  return(instParam)
  
}



### paramList -----

# paramListExample <- list(
#   
#   instName = "Generic",
#   
#   PP = list(
#     funcName = "MassifquantParam",
#     ppm = 20,
#     peakwidth = c(8, 60),
#     snthresh = 3,
#     prefilter = c(6, 500),
#     mzCenterFun = "wMean",
#     integrate = 1,
#     mzdiff = -0.001,
#     fitgauss = TRUE,
#     noise = 200,
#     criticalValue = 1.5,
#     consecMissedLimit = 2,
#     unions = 1,
#     checkBack = 1,
#     withWave = TRUE,
#     verboseColumns = TRUE
#   ),
#   
#   preGrouping = list(
#     funcName = "PeakDensityParam",
#     bw = 5,
#     minFraction = 0.5,
#     minSamples = 1,
#     binSize = 0.003,
#     maxFeatures = 100
#   ),
#   alignment = list(
#     funcName = "PeakGroupsParam",
#     minFraction = 1,
#     extraPeaks = 0,
#     smooth = "loess",
#     span = 0.2,
#     family = "gaussian"
#   ),
#   grouping = list(
#     funcName = "PeakDensityParam",
#     bw = 3,
#     minFraction = 0.5,
#     minSamples = 1,
#     binSize = 0.003,
#     maxFeatures = 100
#   ),
#   componentization = list(
#     
#   )
# )


#additional function for ploting adepted XIC from MSnbase 
.vertical_sub_layout_ex <- function(x, sub_plot = 2) {
  sqrt_x <- base::sqrt(base::length(x))
  ncol <- base::ceiling(sqrt_x)
  nrow <- base::round(sqrt_x)
  rws <- base::split(1:(ncol * nrow * sub_plot), f = base::rep(1:nrow, each = sub_plot * ncol))
  base::do.call(rbind, base::lapply(rws, matrix, ncol = ncol))
}



#additional function for ploting adepted XIC from MSnbase
.plotXIC_ex <- function(x, main = "", col = "white", colramp = colorRampPalette(c("#383E47", "#5E8CAA", "#16B9E5", "#16E5C9", "#16E54C")), #topo.colors colorRampPalette(rev(RColorBrewer::brewer.pal(7,"RdGy")))
                        grid.color = "lightgrey", pch = 21,
                        layout = base::matrix(1:2, ncol = 1), plotTargetMark = plotTargetMark,...) {
  # start edit
  plot_strip <- function(..., v, h) base::plot(...)
  dots <- base::list(...)
  #print(list(...))
  # end edit
  
  if (!plotTargetMark)
  {
    dots$h <- NULL
    dots$v <- NULL
  }
  
  targetColor <- "#0644E9"
  
  if (base::is.matrix(layout))
    graphics::layout(layout)
  
  ## Chromatogram.
  bpi <- base::unlist(base::lapply(base::split(x$i, x$rt), max, na.rm = TRUE))
  brks <- lattice::do.breaks(base::range(x$i), nint = 256)
  graphics::par(mar = c(0, 4, 2, 1))
  
  # start edit
  plot_strip(base::as.numeric(base::names(bpi)), bpi, xaxt = "n", col = col, main = main,
             bg = lattice::level.colors(bpi, at = brks, col.regions = colramp), xlab = "",
             pch = pch, cex = 1.5, ylab = "", las = 2, ...)
  # end edit
  
  if(!is.null(dots$v)) graphics::abline(v=dots$v,col=targetColor,lty=3)
  
  graphics::mtext(side = 4, line = 0, "Intensity", cex = graphics::par("cex.lab"))
  graphics::grid(col = grid.color)
  graphics::par(mar = c(3.5, 4, 0, 1))
  
  # start edit
  plot_strip(x$rt, x$mz, main = "", pch = pch, col = col, xlab = "", ylab = "", cex = 1.5,
             yaxt = "n", bg = lattice::level.colors(x$i, at = brks, col.regions = colramp),
             ...)
  
  if(!is.null(dots$h)) graphics::abline(h=dots$h,col=targetColor,lty=3)
  if(!is.null(dots$v)) graphics::abline(v=dots$v,col=targetColor,lty=3)
  
  if(!is.null(dots$v) & !is.null(dots$h))
  {
    graphics::rect(dots$v-10, dots$h-((5/1E6)*dots$h), dots$v+10, dots$h+((5/1E6)*dots$h),
                   col = NA, lty = 2, border = targetColor)
  }
  # end edit
  
  graphics::axis(side = 2, las = 2)
  graphics::grid(col = grid.color)
  graphics::mtext(side = 1, line = 2.5, "Retention time (sec.)", cex = par("cex.lab"))
  graphics::mtext(side = 4, line = 0, "m/z", cex = par("cex.lab"))
  
}


#' @title plotRawChrom_old
#' @description Plot BPC or TIC chromatograms of an \linkS4class{OnDiskMSnExp} object.
#'
#' @param x An \linkS4class{OnDiskMSnExp} object with one or more files.
#' @param fileIndex The index of the file/s to extract the centroids or profile data.
#' @param mz Optional target \emph{m/z} to obtain an extracted ion chromatogram (EIC).
#' @param ppm The mass deviation to extract the data for the EIC in \code{ppm}.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below.
#' @param rtWindow The time deviation to collect centroids or profile data. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{sec}.
#' @param msLevel The MS level to extract the data. POssible values are 1 or 2.
#' @param type The type of chromatogram. Possible entries are "bpc" of base peak chromatogram or "tic" for total ion chromatogram.
#' The default is "tic".
#'
#' @return A plot.
#'
#' @examples
#' 
#' plotRawChrom_old(ntsIUTA::rawDataExample)
#' 
#' 
plotRawChrom_old <- function(x = rawData, fileIndex = NULL,
                             mz = NULL, ppm = 20,
                             rt = NULL, rtWindow = NULL,
                             rtUnit = "min",
                             msLevel = 1, type = "tic") {
  
  require(magrittr)
  
  if (rtUnit == "min") if (!is.null(rt)) rt <- rt*60
  if (rtUnit == "min") if (!is.null(rtWindow)) rtWindow <- rtWindow*60
  
  if (!base::is.null(fileIndex)) {x <- MSnbase::filterFile(x, fileIndex)}
  
  rtr <- c(base::min(MSnbase::rtime(x)), base::max(MSnbase::rtime(x)))
  mzr <- NULL
  if (!is.null(mz)) {
    
    if (length(mz) == 1) { mzr <- c(mz - ((ppm/1E6)*mz), mz + ((ppm/1E6)*mz)) }
    if (length(mz) == 2) { mzr <- c(mz[1], mz[2]) }
    
  }
  
  if (!is.null(rt)) { rtr <- c((rt) - ifelse(!is.null(rtWindow), rtWindow, 1), (rt) + ifelse(!is.null(rtWindow), rtWindow, 1)) }
  if (is.null(rt)) if (unique(!is.null(rtWindow))) if (length(rtWindow) == 2) { rtr <- c(rtWindow[1], rtWindow[2]) }
  
  colors <- getColors(x, "samples")
  
  if (is.null(mzr)) {
    chrom <- MSnbase::chromatogram(x, aggregationFun = base::ifelse(type == "bpc", "max", base::ifelse(type == "tic", "sum", stop("Unkown type, please use bpc or tic."))), rt = rtr, msLevel = msLevel)
  } else {
    chrom <- MSnbase::chromatogram(x, aggregationFun = base::ifelse(type == "bpc", "max", base::ifelse(type == "tic", "sum", stop("Unkown type, please use bpc or tic."))), rt = rtr, mz = mzr, msLevel = msLevel)
  }
  
  if (!base::is.null(mzr))
  {
    main <- base::paste0("EIC of ", round(mz, digits = 4)," +/- ", round(ppm, digits = 0)," ppm")
  } else {
    main <- base::toupper(type)
  }
  
  # graphics::layout(base::matrix(1:1))
  # graphics::par(mar=c(5,4,4,2)+0.1)
  
  plot <- MSnbase::plot(chrom, col = colors, lwd = 0.5, main = main, xlab = "Retention time (sec.)", ylab = "Intensity",)
  
  plot <- plot %>% graphics::legend("topleft" , legend = x$sample_name, col = colors, lty = 1, lwd = 2, cex = 0.6, bty = "n", xjust = 0)
  
  
  return(plot)
  
}




#' @title makeSetup
#' @description Project setup by screening a given folder for MS files or
#' by creating a new folder where the MS files can be added using \code{addFiles}.
#' Setting \code{createRproject} to \code{TRUE} will create and open an R project in the selected or given folder.
#' When \code{convertFiles} is \code{TRUE} the function uses \code{\link{mzMLconverter}} to automatically convert specified raw MS files to mzML.
#' 
#'
#' @param projPath The directory for the project. The default is the \code{utils::choose.dir()} to select or create a folder.
#' Note that the function will look into subfolders for mzML files or specified raw files when \code{convertFiles} is set to \code{TRUE}.
#' @param date The project date. The default is the system date generated by \code{\link[base]{Sys.time}}.
#' @param groups A vector with the identifier/name of each sample replicate group. If \code{NULL},
#' a grouping tentative will be made using the name of the MS files found in the given \code{projPath}.
#' @param blanks The name of the sample replicate group to be used for blank subtraction.
#' Different blank groups can be assigned to replicate sample groups by editing the \code{sampleInfo} in the output \code{setup} list.
#' @param polarity A vector specifying the polarity mode used in each MS file.
#' Possible values are \code{positive} or \code{negative} for positive and negative mode, respectively.
#' @param convertFiles Logical, set to \code{TRUE} for non mzML MS files being converted to centroided mzML.
#' The default is \code{FALSE} for no conversion even if files are present.
#' @param convertFrom The name of the file format or vendor format to be found and converted.
#' Possible entries are: \emph{thermo}, \emph{bruker}, \emph{agilent}, \emph{ab} (from AB Sciex) and \emph{waters}.
#' @param convertToCentroid Logical, set to \code{TRUE} to convert profile data to centroided. The default and recommended is \code{TRUE}.
#' @param save Logical, set to \code{TRUE} to save the object setup in the \strong{rData} folder.
#' If \strong{rData} folder does not exist in given \code{projPath}, the folder will be created.
#' @param makeNewProject Logical, set to TRUE to create an R project in the given \code{projPath} and open a new R session.
#'
#' @return The output is a simple list, containing: (1) \code{projPath}, the chosen project path; (2) the \code{date}, the given project date; and 
#' (3) the \code{sampleInfo}, a \code{\link[base]{data.frame}} with the following columns:
#' (1) \code{filePath}, the path of each mzML file, including any converted MS files;
#' (2) \code{sample}, the retrieved file name to be used as sample identifier;
#' (3) \code{group}, the name of each replicate sample group
#' (e.g. samples with the same name but with different numbering, such as in the case of replicate samples);
#' (4) \code{blank}, the specified blank group names or the blank group detected by file name (\emph{i.e.} blank, Blank and/or B as file name);
#' and (5) \code{polarity}, the polarity mode used for each sample, possible entries are \code{positive} and \code{negative};
#' 
#' 
#' @details The data.frame \code{sampleInfo} in the output list can always be edited for correction of
#' sample replicate group names, blank replicate groups and assigned polarity, which is set to \code{positive} by default when not given as argument.
#' 
#' @export
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom dplyr mutate
#' @importFrom rstudioapi initializeProject openProject
#'
#' @examples
#' 
#' 
#' 
setupProject <- function(projPath = utils::choose.dir(base::getwd(), "Select or create a project folder"),
                         date = base::Sys.Date(),
                         groups = NULL,
                         blanks = NULL,
                         polarity = "positive",
                         convertFiles = FALSE,
                         convertFrom = NULL,
                         convertToCentroid = TRUE,
                         save = TRUE,
                         makeNewProject = FALSE) {
  
  #Examples
  # projPath <-  system.file(package = "ntsIUTA", dir = "extdata")
  # setup <- setupProject(projPath = projPath, save = FALSE)
  # setup
  
  
  setup <- base::list()
  setup$projPath  <- projPath
  setup$date <- date
  
  #Create holder for list of samples
  sampleInfo <- base::data.frame(filePath = base::as.character(), sample = base::as.character(),
                                 group = base::as.character(), blank = base::as.character())
  
  
  if (convertFiles)
  {
    if (!base::is.null(convertFrom))
      ntsIUTA::mzMLconverter(path = projPath, convertFrom = convertFrom, centroidData = centroidData)
    else { stop("Vendors should be specified for file recognition. Possible entries are: thermo, bruker, agilent,
                 ab (from AB Sciex) and waters.") }
  }
  
  
  #Screen for MS files in project folder and add info to sampleInfo
  msFiles <- base::list.files(path = projPath, pattern = ".mzML|.mzXML", recursive = TRUE, full.names = TRUE, no.. = TRUE)
  
  if (base::length(msFiles) == 0)
  {
    warning("No mzML or mzXML files were found in selected project path. Use addFiles() to add files to project.")
    
    #Adds files info to sample data frame
  } else {
    sampleInfo[1:base::length(msFiles),] <- NA
    sampleInfo$filePath <- msFiles #base::dirname(msFiles)
    sampleInfo$sample <- tools::file_path_sans_ext(base::basename(msFiles))
    sampleInfo$blank <- blanks
    
    if (base::length(groups) < 2 & base::is.null(groups))
    {
      sampleInfo$group <- tools::file_path_sans_ext(base::basename(msFiles))
      #tentative to group samples and add blank group
      sampleInfo <- dplyr::mutate(sampleInfo, group = base::ifelse(base::grepl("A-r|qc|QC", sampleInfo$sample), "QC", group))
      sampleInfo <- dplyr::mutate(sampleInfo, group = base::ifelse(base::grepl("Blank|blank|B-r", sampleInfo$sample), "Blank", group))
      
    } else sampleInfo$group <- groups
    
    if (base::length(blanks) < 2 & base::is.null(blanks)) sampleInfo$blank <- base::ifelse("Blank" %in% sampleInfo$group, "Blank", NA)
    
    sampleInfo$polarity <- polarity
    
    setup$sampleInfo <- sampleInfo
  }
  
  
  if (save) #TODO Add possibility to move the R project file 
  {
    rData <- base::paste0(projPath,"\\rData")
    
    if (base::dir.exists(rData))
    {
      base::saveRDS(setup, file = base::paste0(rData,"\\setup.rds"))
    } else {
      base::dir.create(rData)
      base::saveRDS(setup, file = base::paste0(rData,"\\setup.rds"))
    }
  }
  
  
  if (makeNewProject)
  {
    #TODO add template for main script file
    sp <- base::file.path(projPath,"mainScript.R")
    base::cat(
      
      "#Copy template from template.R using x and/or use ?ntsIUTA for a tutorial.\n
#Run the following code to load the project setup:\n
library(ntsIUTA) #/n
setup <- readRDS('rData/setup.rds')",
      
      file = sp, sep = "")
    rstudioapi::initializeProject(projPath)
    rstudioapi::openProject(projPath, newSession = TRUE)
    base::print("Run  rstudioapi::navigateToFile('mainScript.R')  in the new project to open the mainScript.R file.")
  }
  
  return(setup)
  
}

