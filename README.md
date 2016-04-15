# LowMACAApp
Shiny app around the Bioconductor package LowMACA

LowMACA app is a GUI extension of LowMACA package with new interactive plots and analysis types.
For our Bioconductor package, visit http://www.bioconductor.org/packages/release/bioc/html/LowMACA.html

## HOW TO INSTALL
LowMACA was tested on Windows, Linux and MacOS but there are some things you have to do before running the app according to your operating system.
If you are a Windows user, you should be able to run the app without any installation required. Just skip to "how to run the App"

If you are a MacOS or Linux user, you are required to install Clustal Omega (our trusted aligner) and Ghostscript. The latter is often already present in the major unix distribution but check if the version is > 9.1. 
* For Clustal Omega, we wrote a simple installer for Unix OS, just run "source Download_and_Install_ClustlO.sh" from a terminal.
* For Ghostscript, we found a useful pkg here http://pages.uoregon.edu/koch/Ghostscript-9.16.pkg for macOS
The Linux binaries can be found here:
	http://downloads.ghostscript.com/public/binaries/ghostscript-9.16-linux-x86_64.tgz (64bit)  
	http://downloads.ghostscript.com/public/binaries/ghostscript-9.16-linux-x86.tgz (32 bit)
After the installation, check if the command gs is in the PATH by running which gs from command line terminal.

## HOW TO RUN THE APP
LowMACA relies on some external libraries. Let's start by installing them.
```{r}
myLibraries <- c("devtools","googleVis","reshape2","plyr","d3Network","ggplot2"
                ,"dplyr","shinythemes","shinydashboard","DT","gmp")
sapply(myLibraries , function(x) if(!x %in% installed.packages()[,1]) install.packages(x))

source("https://bioconductor.org/biocLite.R")
biocLite("LowMACAAnotation")
biocLite("LowMACA")
```
If everything is ok, run the app in this way:
```{r}
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("transformPhenotype", "gmelloni" , ref="production")
```

## KNOWN ISSUES
* Package "DT" could be not compatible with R version 3.2.2 . In case runApp.R stops at downloading DT from CRAN, try the development version:
```{r}
library(devtools)
devtools::install_github('rstudio/DT')
```
* Iceweasel Debian browser could have some problems with google charts plots. The LowMACA app was tested successfully on Firefox, Chrome and Safari
* Beware that R 3.2.0 is not correct. Biostrings and MotIV packages have known issues with this R version
