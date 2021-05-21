
#Updates documentation and installs the pkg
devtools::document()
devtools::document()
devtools::install(upgrade = "never")
library(ntsIUTA)

#Check the pkg
devtools::check()
