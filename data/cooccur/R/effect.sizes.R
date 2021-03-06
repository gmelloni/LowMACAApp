effect.sizes <-
function(mod, standardized=TRUE, matrix=FALSE){
  ptab <- mod$results
  if ("sp1_name" %in% colnames(ptab)){
      sp1 <- as.character(ptab$sp1_name)
      sp2 <- as.character(ptab$sp2_name)
  }else{
      sp1 <- ptab$sp1
      sp2 <- ptab$sp2
  }
  
  if (standardized==T){
      effs <- data.frame(sp1,sp2,effect=as.numeric(((ptab$obs_cooccur - ptab$exp_cooccur)/mod$sites)))
  }else{
      effs <- data.frame(sp1,sp2,effect=as.numeric((ptab$obs_cooccur - ptab$exp_cooccur)))
  }
  
  if (matrix==F){
    effs
  }else{
    effs_1 <- effs
    effs_2 <- effs[,c(2,1,3)]
    colnames(effs_2) <- c("sp1","sp2","effect")
    effs <- rbind(effs_1,effs_2)
    effs <- recast(data=effs,formula=sp1~sp2,measure.var="effect",id.var=c("sp1","sp2"))
    #m <- as.matrix(effs[[1]])
    
    #rns <- effs[[2]][[1]]
    #cns <- effs[[2]][[2]]
    
    #row.names(m) <- as.character(rns$sp1)
    #colnames(m) <- as.character(cns$sp2)
    
    #m <- effs[mod$spp.names,mod$spp.names]
    #as.dist(as.matrix(m))
        row.names(effs) <- effs$sp1
    effs$sp1 <- NULL
        effs <- as.matrix(effs)
    effs <- effs[mod$spp.names,mod$spp.names]
    as.dist(as.matrix(effs))
    
  }
  
}
