
#logo header
sourcepath <- system.file(package = "ntsIUTA", dir = "extdata")
logo <- file.path(sourcepath, "logo.png")
img <- htmltools::img(
  src = knitr::image_uri(logo),
  alt = "logo",
  style =
    "position:absolute;
     top:3%;
     margin-left:55.5%;
     width:350px;
     height:157px;"
)

htmlhead <- paste0('<script> document.write(\'<div class = "logo">', img, '</div>\') </script> ')
packagepath <- "C:\\Users\\Ricardo\\Documents\\CodeProjects\\ntsIUTA\\inst\\extdata"
readr::write_lines(htmlhead, file = paste0(packagepath, "\\logo.html"))
readr::write_lines(htmlhead, file = paste0(sourcepath, "\\logo.html"))

# .logo { 
#   width: 350px;
#   height: 157px;
#   text-align: right;
# }
# .logo img {
#   display: inline-block;
# }
# <div class="logo">
#   <img src=`r file.path(system.file(package = "ntsIUTA", dir = "extdata"), "logo.png")` alt = "example" class="logo"/>
# </div>
 