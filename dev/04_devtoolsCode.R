
#updates documentation and installs the pkg
devtools::document()
devtools::document()
devtools::install(upgrade = "never")
library(ntsIUTA)

#check the pkg
devtools::check()

#load
devtools::document()
devtools::document()
devtools::load_all()

#show loaded packages
search()

#see todo list
todor::todor()

#show s4 methods
showMethods("a method")


#update the wd to for the package
setwd("C:/Users/Ricardo/Documents/CodeProjects/ntsIUTA")

#ntsIUTA system path
system.file(package = "ntsIUTA", dir = "extdata")
