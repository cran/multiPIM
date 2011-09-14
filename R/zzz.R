################################################################################
######################################### zzz.R ################################
################################################################################

## Created: 2009.05.26
## Author: Stephan Ritter

.onAttach <- function(libname, pkgname) {

  packageStartupMessage('This is multiPIM version ',
                        packageDescription('multiPIM')$Version, "\n")

}
