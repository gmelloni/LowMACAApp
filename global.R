# detach("package:LowMACA" , unload=TRUE)
# detach("package:LowMACAAnnotation" , unload=TRUE)

library(shiny)
library(gmp)
# library(shinyjs)
# shinyjs::useShinyjs()
# library(devtools)

# While waiting that the LowMACA package passes to version 1.1, we have to call it as an external package
library(LowMACA)
# load_all(file.path("data" , "v0.99.5_bis" , "LowMACAAnnotation"))
# load_all(file.path("data" , "v0.99.5_bis" , "LowMACA"))
load_all(file.path("data" , "cooccur"))
# detach("package:devtools", unload=TRUE)
source(file.path("data" , "custom_functions" , "allPfamAnalysis_special.R"))
source(file.path("data" , "custom_functions" , "conditionalDisabledPanel.R"))

# Under windows, both cluslomega and ghostscript are inside the LowMACA app
# Under Unix, ghostscript is launched as gs from the PATH library
# ClustalO must be installed first
if(Sys.info()['sysname']=="Windows") {
    clustalo_cmd <- file.path("data" , "clustal-omega-1.2.0-win32" , "clustalo.exe")
    if(grepl("64" , Sys.info()['release'])){
        Sys.setenv(R_GSCMD = file.path(getwd() , "data" , "Ghostscript" , "bin" , "gswin64c.exe"))
    }else{
        Sys.setenv(R_GSCMD = file.path(getwd() , "data" , "Ghostscript" , "bin" , "gswin32c.exe"))
    }
} else {
    clustalCommand <- Sys.which("clustalo")
    if(clustalCommand=="") {
        if(dir.exists(file.path("data" , "ClustalForUnix"))) {
            clustalo_cmd <- file.path("data" , "ClustalForUnix" , "bin" , "clustalo")
        } else {
            stop("ClustalOmega MUST be installed under Unix systems before running the app!")
        }
    } else {
        clustalVersion <- system(paste(clustalCommand ,  "--version") , intern=TRUE)
        if(!grepl("^1.2" , clustalVersion)) {
            if(dir.exists(file.path("data" , "ClustalForUnix"))) {
                clustalo_cmd <- file.path("data" , "ClustalForUnix" , "bin" , "clustalo")
            } else {
                warning("Clustal Omega version could be not compatible. LowMACA was tested on 1.2.x. Use our installer to download the correct version" , immediate.=TRUE)
            }
        } else {
            message("It looks like you already installed Clustal Omega, good boy!")
            clustalo_cmd <- "clustalo"
        }
    }
    gs_cmd <- Sys.which("gs")
    if(gs_cmd==""){
        warning("Ghostscript is not in your PATH. You will not be able to visualize the Logo Plot" , immediate.=TRUE)
    } else {    
        gsversion <- system("gs --version" , intern=TRUE)
        if(! as.numeric(gsversion) >= 9.1) {
            warning("Ghostscript version is < 9.0 and could not be compatible. The Logo plot could not be shown" , immediate.=TRUE)
        }
    }
}

myPfam <- readRDS(file.path("data" ,"custom_data", "myPfam.RData"))
myPfam_red <- readRDS(file.path("data" ,"custom_data", "myPfam_red.RData"))
myUni <- getMyUni()
tumor_type <- readRDS(file.path("data" , "custom_data", "tumor_type.RData"))
repos <- readRDS(file.path("data" , "custom_data", "cBioAllMutations05052015.RData"))

start <- unique(sort(myPfam$Pfam_ID))[1]
pfam_list <- unique(sort(myPfam$Pfam_ID))
miss_type <- c("Missense_Mutation"
                    ,"In_Frame_Del"
                    ,"In_Frame_Ins"
                    )
trunc_type <- c("Frame_Shift_Del"
                ,"Nonsense_Mutation"
                ,"Translation_Start_Site"
                ,"Frame_Shift_Ins"
                ,"Nonstop_Mutation"
                ,"Splice_Site"
                ,"Indel"
                )
