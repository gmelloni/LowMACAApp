library(googleVis)
library(reshape2)
library(plyr)
library(d3Network)
library(ggplot2)
library(dplyr)
library(DT)
# tmp <- file.path("data" , "repos.sqlite3")
# repos_sqlite <- src_sqlite(tmp, create = FALSE)

TableOR_corr <- function(x , y){
  xtab <- table(x , y)+.5
  n00 <- xtab[1,1]
  n01 <- xtab[1,2]
  n10 <- xtab[2,1]
  n11 <- xtab[2,2]
  OR <- (n00 * n11)/(n01 * n10)
  return(OR)
}


shinyServer(function(input, output , session) {
  #Put The Logo
  output$myImage <- renderImage({
    width  <- session$clientData$output_myImage_width
    height <- as.character(as.numeric(session$clientData$output_myImage_height)/2)
    # For high-res displays, this will be greater than 1
    pixelratio <- session$clientData$pixelratio
    filename <- normalizePath(file.path('data', 'images' , 'lowmaca_logo.png'))
    list(src = filename, alt="Coolest Logo Ever!", width=width , height=height)
  }, deleteFile = FALSE)
  # Pur LowMACA cheatsheet
  output$cheatsheet <- renderImage({
    width  <- session$clientData$output_myImage_width*2
    height <- session$clientData$output_myImage_height*2
    # For high-res displays, this will be greater than 1
    pixelratio <- session$clientData$pixelratio
    filename <- normalizePath(file.path('data', 'images' , 'LowMACA_cheatsheet.png'))
    list(src = filename, alt="Coolest Logo Ever!", width=width , height=height)
  }, deleteFile = FALSE)
  #Embed a youtube video
  output$myVideo <- renderUI({
      h3("Video Tutorial" 
      , br() 
      , tags$iframe(width="560" 
                  , height="315"
                  , src="https://www.youtube.com/embed/FJaxzsFFBTg"
                  , frameborder="0" 
                  , allowfullscreen=NA)
      )
  })
  #Change gene selector according to pfam
  observe({
    #This updater will select only the genes of a particular Pfam
    #In case No Pfam is selected, all the genes become available
    updateSelectInput(session , "genes" , "Choose Your Genes:"
                  , choices= if(input$pfam=="No Pfam") { 
                                c("All Genes" , sort(unique(myUni[ , "Gene_Symbol"]))) 
                              } else { 
                                c("All Genes" 
                                  , sort(unique(myPfam[ myPfam$Pfam_ID==strsplit(input$pfam , "_")[[1]][1], "Gene_Symbol"])) 
                                )
                              }
                  , selected="All Genes")
  })
  #Change method if No Pfam is selected
  observe({
    #If No Pfam is selected, the method is forced to full protein
    updateSelectInput(session , "fullProtein" , "Running Full Protein:"
                      , choices=if(input$pfam=="No Pfam") {
                                  c("Yes") 
                                } else { 
                                  c("No", "Yes")
                                }
                      , selected=if(input$pfam=="No Pfam") {
                                  "Yes" 
                                } else {
                                  "No"
                                })
  })
  #Create LowMACA obj
  lowMACAObj <- eventReactive(input$goButton, {
  withProgress(message = 'LowMACA Analysis', value = 0, {
    n <- 6
    incProgress(1/n, detail = "Creating a LowMACA Object...")
    bw <- input$bw
    tumor_type <- input$tumor_type
    genes <- if("All Genes" %in% input$genes) sort(unique(myPfam[ myPfam$Pfam_ID==strsplit(input$pfam , "_")[[1]][1], "Gene_Symbol"])) else input$genes
    #print(genes)
    if(input$fullProtein=="Yes")
      lowMACAObj <- newLowMACA(genes=genes)
    else
      lowMACAObj <- newLowMACA(pfam=strsplit(input$pfam , "_")[[1]][1] , genes=genes)
    lmParams(lowMACAObj)$density_bw <- bw
    if(!"All Tumors" %in% tumor_type)
      lmParams(lowMACAObj)$tumor_type <- tumor_type
    lmParams(lowMACAObj)$clustal_cmd <- clustalo_cmd
    mut_type <- switch(input$mut_type , 
                        "All Mutations" = "all"
                        ,"Missense Type" = "missense"
                        , "Truncating Type" = "truncating")
    lmParams(lowMACAObj)$mutation_type <- mut_type
    incProgress(1/n, detail = "Align Sequences...")
    lowMACAObj <- alignSequences(lowMACAObj)
    incProgress(1/n, detail = "Get Mutations...")
    # myIn <- paste0("('" , paste(genes , collapse="','") , "')")
    # repos <- as.data.frame(tbl(repos_sqlite, sql(paste("SELECT * FROM repos_tbl WHERE Gene_Symbol IN" , myIn))))
    repos <- repos[ repos$Gene_Symbol %in% genes , ]
    lowMACAObj <- getMutations(lowMACAObj , repos=repos)
    incProgress(1/n, detail = "Map Mutations...")
    lowMACAObj <- mapMutations(lowMACAObj)
    incProgress(1/n, detail = "Mutation Entropy...")
    lowMACAObj <- entropy(lowMACAObj)
  })
    return(lowMACAObj)
  })
  #LowMACA plot
  output$plot <- renderPlot({
      if(as.logical(input$goButton)){
        withProgress(message = 'LowMACA Plot', value = 0, {
          n <- 3
          incProgress(1/n, detail = "Calculating Statistics...")
          bw <- input$bw
          lowMACAObj <- lowMACAObj()
          lmParams(lowMACAObj)$density_bw <- bw
          lowMACAObj <- entropy(lowMACAObj)
          incProgress(1/n, detail = "Rendering...")
          lmPlot(lowMACAObj , conservation=input$cons)
          })
      }else{
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("Select A Pfam Domain, for example, start typing ras\n",
                              "Select the genes of your interest, for example KRAS, NRAS and HRAS\n",
                              "Choose the tumor types you want to analyze\n",
                              "To find single mutation hotspots, leave the bandwidth at 0\n",
                              "Select the minmum accepted level of conservation\n",
                              "Do you want to align the entire proteins or just the domains?\n",
                              "######## LOWMACAIZE! ########"), 
                              cex = 1.6, col = "black")
      }
  })
  #LowMACA Barplot
  output$bpAll <- renderGvis({
    if(as.logical(input$goButton)){
      withProgress(message = 'LowMACA Google Barplot', value = 0, {
          n <- 3
          incProgress(1/n, detail = "Munging Data...")
          df <- as.data.frame(t(lowMACAObj()@mutations$aligned))
          df$Multiple_Aln_pos <- as.character(seq(1,nrow(df)))
          df <- df[ , c("Multiple_Aln_pos" , colnames(df)[colnames(df)!="Multiple_Aln_pos"])]
          mut <- lowMACAObj()@mutations$data
          alignment <- lowMACAObj()@alignment$ALIGNMENT
          mut_aligned <- merge( mut , alignment[!is.na(alignment$Ref) , ]
              , by.x=c("Entrez" , "Gene_Symbol" , "Amino_Acid_Position") 
              , by.y=c("Entrez" , "Gene_Symbol" , "Ref") 
              , all.x=TRUE)
          if(lowMACAObj()@arguments$mode=="genes")
            mut_aligned$domainID <- mut_aligned$Gene_Symbol
          mut_aligned <- mut_aligned[!is.na(mut_aligned$Align), ]
          mut_aligned_agg <- aggregate(Amino_Acid_Change ~ domainID + Align 
                , mut_aligned , FUN=function(x) paste(sort(unique(x)) , collapse=","))
          for(i in colnames(mut_aligned_agg))
            mut_aligned_agg[ , i] <- as.character(mut_aligned_agg[,i])
          mut_aligned_agg_cast <- dcast(data= mut_aligned_agg 
            , formula = Align ~ domainID , value.var="Amino_Acid_Change")
          colnames(mut_aligned_agg_cast)[colnames(mut_aligned_agg_cast)!="Align"] <- 
              paste0(colnames(mut_aligned_agg_cast)[colnames(mut_aligned_agg_cast)!="Align"] , ".html.tooltip")
          df2 <- merge(df , mut_aligned_agg_cast , by.x="Multiple_Aln_pos" , by.y="Align" , all.x=TRUE)
          df2 <- df2[ order(as.numeric(df2$Multiple_Aln_pos)) , ]
          for(i in grep("tooltip" , colnames(df2) , value=TRUE)) {
            df2[is.na(df2[ , i]) , i] <- ""
            df2[ , i] <- paste(sub(".html.tooltip" , "" , i)
                            , paste("Position:" , df2$"Multiple_Aln_pos")
                            , paste("Number Of Mutations:" , df2[ , sub(".html.tooltip" , "" , i)])
                            , paste("Mutation Types:" , df2[ , i] , sep="\\n")
                            , sep="\\n\\n"
                            )
          }
          incProgress(1/n, detail = "Rendering Plot...")
          return(gvisColumnChart(df2 
            , xvar="Multiple_Aln_pos" 
            , yvar=sort(grep("Multiple_Aln_pos" , colnames(df2) , invert=TRUE , value=TRUE)) 
            , options=list(
                    isStacked=TRUE
                    , explorer="{actions: ['dragToZoom','rightClickToReset'], maxZoomIn:0.05 , keepInBounds: true , axis: 'both'}"
                    #, crosshair="{trigger:'both'}"
                    , title="LowMACA tooltip Barchart: drag an area to zoom, right click to go back, tooltip for more information"
                    , height=800
                    , width=3000
                    , vAxis.title="Number Of Mutations"
                    , hAxis.title="Positions Of The Consensus Sequence"
                    #, hAxis.viewWindowMode="pretty"
                    #, vAxis.viewWindowMode="pretty"
                    , chartArea= "{width: '90%', height: '90%'}"
                    , legend= "{position: 'in' , textStyle: {color: 'black', fontSize: 10} }"
                    , titlePosition= 'in', axisTitlesPosition= 'in'
                    , hAxis= "{textPosition: 'out'}", vAxis= "{textPosition: 'out'}"
                    , bar.groupWidth= '75%'
                    #, legend.position="bottom"
                    #, tooltip="{isHtml:'true'}"
                    , tooltip.trigger="both"
                    #, chartArea='{left:100,top:500 ,right:50}'
                    ,enableScrollWheel=TRUE
            )
          ))
        })
    } else {
      return(gvisBarChart(data.frame(myGene="no data" , myMutations=0)
        , options=list(isStacked=TRUE 
        #, height=1000 
        #, width=1500
        #, chartArea='{left:0,top:0,width:"400%",height:"400%"}'
        )))
    }
  })
  networkObj <- reactive({
          align <- lowMACAObj()@mutations$aligned
          alignmelt <- melt(align)
          alignmelt2 <- alignmelt
          colnames(alignmelt) <- c("GeneA" , "Position" , "MutationsA")
          colnames(alignmelt2) <- c("GeneB" , "Position" , "MutationsB")
          if(length(unique(alignmelt$GeneA))<3)
            return(NULL)
          alignmerge <- merge(alignmelt , alignmelt2)
          alignmerge <- with(alignmerge , alignmerge[ GeneA!=GeneB , ])
          alignmerge <- with(alignmerge , alignmerge[ MutationsA!=0 & MutationsB!=0 , ])
          alignmerge$Mutations <- as.numeric(apply(alignmerge , 1 , function(x) min(as.numeric(x['MutationsA']) , as.numeric(x['MutationsB']))))
          #alignmerge$Mutations <- as.numeric(apply(alignmerge , 1 , function(x) sqrt(as.numeric(x['MutationsA']) * as.numeric(x['MutationsB']))))       
          gm_mean = function(a){prod(a)^(1/length(a))}
          alignagg <- aggregate(Mutations ~ GeneA + GeneB , alignmerge , FUN=length )
          alignagg$GeneA <- as.character(alignagg$GeneA)
          alignagg$GeneB <- as.character(alignagg$GeneB)
          alignagg$dup <- apply(alignagg , 1 , function(x) paste(sort(x[c('GeneA' , 'GeneB')]) , collapse=""))
          MisLinks <- alignagg[ !duplicated(alignagg$dup) , ]
          MisLinks$dup <- NULL
          MisNodes1 <- data.frame(Sequence=unique(alignagg[ , "GeneA"]) 
            , Gene=sapply( strsplit( unique(alignagg[ , "GeneA"]) , "\\|" ) , function(x) x[1]))
          MisNodes2 <- data.frame(Sequence=unique(alignagg[ , "GeneB"]) 
            , Gene=sapply( strsplit( unique(alignagg[ , "GeneB"]) , "\\|" ) , function(x) x[1]))
          MisNodes <- rbind(MisNodes1 , MisNodes2)
          MisNodes$Gene <- as.integer(as.numeric(factor(MisNodes$Gene , levels=unique(MisNodes$Gene))))
          MisNodes <- unique(MisNodes)
          MisNodes <- MisNodes[ order(MisNodes$Sequence) , ]
          rownames(MisNodes) <- as.character(1:nrow(MisNodes))
          MisLinks <- MisLinks[ order(MisLinks$GeneA) , ]
          MisLinks$GeneA <- suppressMessages(suppressWarnings(
                          as.integer(mapvalues(MisLinks$GeneA , from=MisNodes$Sequence , to=as.integer(rownames(MisNodes))-1))
                          ))
          MisLinks$GeneB <- suppressMessages(suppressWarnings(
                          as.integer(mapvalues(MisLinks$GeneB , from=MisNodes$Sequence , to=as.integer(rownames(MisNodes))-1))
                          ))
          MisLinks$Mutations <- as.integer(MisLinks$Mutations)
          return(list(MisLinks=MisLinks , MisNodes=MisNodes))
    })
  #Network
  output$networkPlot <- renderPrint({
      if(as.logical(input$goButton)){
        if(is.null(networkObj()))
          return(d3SimpleNetwork(data.frame(GeneA=c("Sorry" , "No Connections") 
                                        , GeneB=c("No plot with no connections" , "No plot with no connections")) 
                                      , standAlone = FALSE 
                                      , parentElement='#networkPlot'
                                      , width=1000
                                      , height=1000
                                      #, zoom=TRUE
                                      , charge=-300
                                      , fontsize=30
                                      ,linkDistance=400))
        withProgress(message = 'LowMACA Network', value = 0, {
          n <- 3
          incProgress(1/n, detail = "Munging Data...")
          MisNodes <- networkObj()[["MisNodes"]]
          MisLinks <- networkObj()[["MisLinks"]][ networkObj()[["MisLinks"]]$Mutations>=input$CommonPositions , ]
          incProgress(1/n, detail = "Rendering Network...")
          if(nrow(MisLinks)<2)
            return(d3SimpleNetwork(data.frame(GeneA=c("Sorry" , "No Connections") 
                                        , GeneB=c("Move the selector to a lower value" , "Move the selector to a lower value")) 
                                      , standAlone = FALSE 
                                      , parentElement='#networkPlot'
                                      , width=1000
                                      , height=1000
                                      #, zoom=TRUE
                                      , charge=-300
                                      , fontsize=30
                                      ,linkDistance=400))
          return(d3ForceNetwork(Links = MisLinks, Nodes = MisNodes,
                         Source = "GeneA", Target = "GeneB",
                         Value = "Mutations", NodeID = "Sequence",
                         Group = "Gene", width = 1000, height = 800
                         , fontsize=20
                         , linkWidth="function(d) { return Math.log(d.value); }"
                         , opacity = 0.9 
                         , zoom=TRUE 
                         , standAlone = FALSE
                         , charge=-100
                         , parentElement = '#networkPlot' , linkDistance = 400 ))
          })
      } else {
        return(d3SimpleNetwork(data.frame(GeneA=c("No Data Yet" , "Welcome") 
                                        , GeneB=c("Run a new analysis" , "Run a new analysis")) 
                                      , standAlone = FALSE 
                                      , parentElement='#networkPlot'
                                      , width=1000
                                      , height=1000
                                      #, zoom=TRUE
                                      , charge=-300
                                      , fontsize=30
                                      ,linkDistance=400))
      }
  })
  #Protter plot
  output$protter <- renderImage({
    if(as.logical(input$goButton)){
      withProgress(message = 'LowMACA Protter', value = 0, {
            n <- 3
            incProgress(1/n, detail = "Munging Data...")
            tmp <- tempfile()
            protter(lowMACAObj() , tmp , conservation=input$cons)
            incProgress(1/n, detail = "Render Image...")
            list(src = tmp, alt="Protter Plot", width="800" , height="600")
            })
    } else {
      tmp <- tempfile(fileext='.png')
      png(tmp)
      plot.new()
      dev.off()
      list(src = tmp, contentType = 'image/png' , alt="Protter Plot", width="800" , height="600")
    }
  }, deleteFile = TRUE)
  #Add text on top of plot
  output$text1 <- renderText({
    if(!as.logical(input$goButton)){
        "Ready to Plot!"
    } else {
      pfam_name <- paste(unique(myPfam[ myPfam$Pfam_ID==strsplit(input$pfam , "_")[[1]][1] , "Pfam_Name"]) , "Domain")
      if("All Genes" %in% input$genes)
        paste("LowMACA Analysis of:" , pfam_name)
      else
        paste("LowMACA Analysis of:" , pfam_name , "choosing the following genes:" , paste(paste(input$genes , collapse=",")) )
    }
    })
  #Pfam Helper table
  output$PfamTable <- DT::renderDataTable({
      return(myPfam_red[ order(myPfam_red$Gene_Symbol) , ])
  } , rownames = FALSE  ,options = list(searchHighlight = TRUE , 
    scrollX = TRUE,
    scrollCollapse = FALSE), filter = 'top')
  #HGNC Helper Table
  output$geneTable <- DT::renderDataTable({
      return(myUni[ order(myUni$Gene_Symbol) , colnames(myUni)!="AMINO_SEQ"])
  } , rownames = FALSE  ,options = list(searchHighlight = TRUE , 
    scrollX = TRUE,
    scrollCollapse = FALSE), filter = 'top')
  lfmObj <- reactive({
    #lfm(lowMACAObj() , conservation=input$cons , metric="pvalue")
    lfm(lowMACAObj() , conservation=input$cons)
    })
  #Significant mutation table
  output$table <- DT::renderDataTable({
        if(as.logical(input$goButton)){
          out <- tryCatch( unique( lfmObj()[ , c("Gene_Symbol" , "Amino_Acid_Change" , "Tumor_Type" , "Sample" , "Multiple_Aln_pos" , "metric")]) 
                , error=function(e) data.frame( Gene_Symbol="No Data" , Amino_Acid_Change="No Data" , Tumor_Type="No Data" 
                                                    , Multiple_Aln_pos="No Data" , metric="No Data"))
          if("Sample" %in% colnames(out)){
            out$Tumor_Type <- with(out , paste(Sample , Tumor_Type) )
            out$Sample <- NULL
          }
        } else {
          out <- data.frame( Gene_Symbol="No Data" , Amino_Acid_Change="No Data" , Tumor_Type="No Data" , Multiple_Aln_pos="No Data" , metric="No Data")
        }
        return(out)
  } , rownames = FALSE  ,options = list(searchHighlight = TRUE ,
    	scrollX = TRUE,
    	scrollCollapse = FALSE), filter = 'top')

  #Cooccurrence analysis position based
  cooccurObjPos <- reactive({
    df <- lfmObj()
    if(nrow(df)==0)
      return(NULL)
    df$gs_mut <- paste(df$Gene_Symbol, df$Multiple_Aln_pos, sep='_')
    mat <- acast(df, Sample ~ gs_mut, fun.aggregate=length)
    mat[mat>=2] <- 1
    cooccur_myData <- cooccur(mat=mat , type="site_spp", thresh=FALSE, spp_names=TRUE , prob="comb")
    return(list(cooccur_myData=cooccur_myData , mat=mat))
    })
  output$coocPos <- renderPlot({
    if(!as.logical(input$goButton))
      return(plot.new())
    if(is.null(cooccurObjPos()))
      return(plot.new())
    plot(cooccurObjPos()[["cooccur_myData"]] ,plotrand=FALSE ) #+ ggtitle("Cooccurrence\nMutual Exclusivity\nPlot")

  })
  output$coocDataPos <- renderTable({
    if(!as.logical(input$goButton))
      return(data.frame(GeneA.PosX="No Data Yet" , GeneB.PosY="No Data Yet" 
                        #, Mutations.in.A="No Data Yet" , Mutations.in.B="No Data Yet" 
                        , pVal.MutEx="No Data Yet" 
                        , pVal.Cooc="No Data Yet" , Effect.Size="No Data Yet"))
    df <- lfmObj()
    if(nrow(df)==0)
      return(data.frame(GeneA.PosX="No Sign Mutations" , GeneB.PosY="No Sign Mutations" 
                        #, Mutations.in.A="No Sign Mutations" , Mutations.in.B="No Sign Mutations" 
                        , pVal.MutEx="No Sign Mutations" , pVal.Cooc="No Sign Mutations" , Log10OddsRatio="No Sign Mutations"))
    results <- cooccurObjPos()[["cooccur_myData"]]$results
    if(nrow(results)==0)
      return(data.frame(GeneA.PosX="No Sign Mutations" , GeneB.PosY="No Sign Mutations" 
                        #, Mutations.in.A="No Sign Mutations" , Mutations.in.B="No Sign Mutations" 
                        , pVal.MutEx="No Sign Mutations" , pVal.Cooc="No Sign Mutations" , Log10OddsRatio="No Sign Mutations"))
    #effectSize <- effect.sizes(cooccurObjPos()[["cooccur_myData"]])
    #results <- merge(results , effectSize , by.x=c("sp1_name" , "sp2_name") , by.y=c("sp1" , "sp2") , all.x=TRUE)
    results <- results[ results$p_lt<=0.05 | results$p_gt<=0.05 , ]
    results <- results[ , c("sp1_name" ,  "sp2_name" , "sp1_inc" , "sp2_inc" 
                          , "obs_cooccur" , "prob_cooccur" , "exp_cooccur"
                          , "p_lt" , "p_gt" 
                          #, "effect" 
                          )]
    colnames(results) <- c("GeneA.PosX" , "GeneB.PosY" , "Mutations.in.A" , "Mutations.in.B" 
                              , "Observed.Cooc" , "Cooc.Prob" , "Expected.Cooc" , "pVal.MutEx" 
                              , "pVal.Cooc" 
                              #, "Effect.Size"
                              )
    results$Mutations.in.A <- as.integer(results$Mutations.in.A)
    results$Mutations.in.B <- as.integer(results$Mutations.in.B)  
    mat <- cooccurObjPos()[["mat"]]
    results$Log10OddsRatio <- log10(sapply(1:nrow(results) , function(i) {
                                x <- factor(as.vector(mat[ , as.character(results[i,"GeneA.PosX"]) ] ) , levels=c(0,1))
                                y <- factor(as.vector(mat[ , as.character(results[i,"GeneB.PosY"]) ] ) , levels=c(0,1))
                                if(length(x)==0 || length(y)==0)
                                  return(NA)
                                else
                                  return(TableOR_corr(x=x , y=y))
                                }))
    return(results[ , c("GeneA.PosX" , "GeneB.PosY" 
                        #, "Mutations.in.A" , "Mutations.in.B" 
                        , "pVal.MutEx" , "pVal.Cooc" , "Log10OddsRatio")])
  }
  , include.rownames=FALSE , digits=4)

  cooccurObjGene <- reactive({
    df <- lfmObj()
    if(nrow(df)==0)
      return(NULL)
    mat <- acast(df, Sample ~ Gene_Symbol, fun.aggregate=length)
    mat[mat>=2] <- 1
    cooccur_myData <- cooccur(mat=mat , type="site_spp", thresh=FALSE, spp_names=TRUE , prob="comb")
    return(list(cooccur_myData=cooccur_myData , mat=mat))
  })
  output$coocGene <- renderPlot({
    if(!as.logical(input$goButton))
      return(plot.new())
    if( is.null(cooccurObjGene()) )
      return(plot.new())
    plot(cooccurObjGene()[["cooccur_myData"]] ,plotrand=FALSE ) #+ ggtitle("Cooccurrence\nMutual Exclusivity\nPlot")
  })
  output$coocDataGene <- renderTable({
    if(!as.logical(input$goButton))
      return(data.frame(GeneA="No Data Yet" , GeneB="No Data Yet" 
                      #, Mutations.in.A="No Data Yet" , Mutations.in.B="No Data Yet" 
                      , pVal.MutEx="No Data Yet" 
                      , pVal.Cooc="No Data Yet" , Effect.Size="No Data Yet"))
    df <- lfmObj()
    if(nrow(df)==0)
      return(data.frame(GeneA="No Sign Mutations" , GeneB="No Sign Mutations" 
                        #, Mutations.in.A="No Sign Mutations" , Mutations.in.B="No Sign Mutations" 
                        , pVal.MutEx="No Sign Mutations" 
                        , pVal.Cooc="No Sign Mutations" , Log10OddsRatio="No Sign Mutations"))
    results <- cooccurObjGene()[["cooccur_myData"]]$results
    if(nrow(results)==0)
      return(data.frame(GeneA.PosX="No Sign Mutations" , GeneB.PosY="No Sign Mutations" 
                        #, Mutations.in.A="No Sign Mutations" , Mutations.in.B="No Sign Mutations" 
                        , pVal.MutEx="No Sign Mutations" , pVal.Cooc="No Sign Mutations" , Log10OddsRatio="No Sign Mutations"))
    #effectSize <- effect.sizes(cooccurObjGene()[["cooccur_myData"]] , standardized=FALSE)
    #results <- merge(results , effectSize , by.x=c("sp1_name" , "sp2_name") , by.y=c("sp1" , "sp2") , all.x=TRUE)
    results <- results[ results$p_lt<=0.05 | results$p_gt<=0.05 , ]
    results <- results[ , c("sp1_name" ,  "sp2_name" , "sp1_inc" , "sp2_inc" 
                          , "obs_cooccur" , "prob_cooccur" , "exp_cooccur"
                          , "p_lt" , "p_gt" 
                          #, "effect" 
                          )]
    colnames(results) <- c("GeneA" , "GeneB" , "Mutations.in.A" , "Mutations.in.B" 
                            , "Observed.Cooc" , "Cooc.Prob" , "Expected.Cooc" , "pVal.MutEx" 
                            , "pVal.Cooc" 
                            #, "Effect.Size"
                            )
    results$Mutations.in.A <- as.integer(results$Mutations.in.A)
    results$Mutations.in.B <- as.integer(results$Mutations.in.B)
    mat <- cooccurObjGene()[["mat"]]
    results$Log10OddsRatio <- log10(sapply(1:nrow(results) , function(i) {
                                x <- factor(as.vector(mat[ , as.character(results[i,"GeneA"]) ] ) , levels=c(0,1))
                                y <- factor(as.vector(mat[ , as.character(results[i,"GeneB"]) ] ) , levels=c(0,1))
                                if(length(x)==0 || length(y)==0)
                                  return(NA)
                                else
                                  return(TableOR_corr(x=x , y=y))
                                }))
    return(results[ , c("GeneA" , "GeneB"
                        #, "Mutations.in.A" , "Mutations.in.B" 
                        ,"pVal.MutEx" , "pVal.Cooc" , "Log10OddsRatio")])
  }
  , include.rownames=FALSE , digits=3)
  
  ###################################
  ##### For datadriven analysis #####
  ###################################

  # Show the uploaded table
  output$ddcontents <- DT::renderDataTable({
    if(!as.logical(input$example))
      inFile <- input$myfile
    else
      inFile <- list(name="custom_analysis_example.txt" 
        , datapath=file.path("data" , "custom_data" , "custom_analysis_example.txt"))
    if(is.null(inFile))
      return(NULL)
    read.table(inFile$datapath, header=TRUE, sep="\t")
  	} , rownames=FALSE  ,options = list(searchHighlight = TRUE , 
    	scrollX = TRUE,
    	scrollCollapse = FALSE), filter = 'top')

  # output$controlAction <- reactive({
  #   inFile <- input$myfile
  #   if(!is.null(inFile))
  #     return("enable")
  #   else
  #     return("disable")
  # })
  
  # Enable action button
  observe({
    if(!as.logical(input$example))
      inFile <- input$myfile
    else
      inFile <- list(name="custom_analysis_example.txt" , datapath=file.path("data" , "custom_data" , "custom_analysis_example.txt"))
    if(!is.null(inFile))
      session$sendCustomMessage("go_with_the_analysis" , list(mex='Go with the analysis'))
  })

	# Change text under the input data table
  output$changingtext <- renderText({
    	if(!as.logical(input$goDDButton))
      		return("Waiting for the analysis...")
    	else
      		return("")
  })
  
  # Perform the allPfamAnalysis
  ddAnalysis <- eventReactive(input$goDDButton, {
    progress <- shiny::Progress$new()
    progress$set(message = "Custom Analysis...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    n <- 5
    updateProgress <- function(detail = NULL) {
      progress$inc(amount = 1/n, detail = detail)
    }
    if(!as.logical(input$example))
      inFile <- input$myfile
    else
      inFile <- list(name="custom_analysis_example.txt" , datapath=file.path("data" , "custom_data" , "custom_analysis_example.txt"))
    if (is.null(inFile))
      return(NULL)
    myfile <- read.table(inFile$datapath, header=TRUE, sep="\t" , as.is=TRUE)
    ddmut_type <- switch(input$ddmuttype , 
                      "All Mutations" = "all"
                      ,"Missense Type" = "missense"
                      , "Truncating Type" = "truncating")
    return(allPfamAnalysis_special(repos=myfile , allLowMACAObjects=NULL , parallel=FALSE
        , mutation_type=ddmut_type
        , NoSilent=TRUE , mail=NULL
        , perlCommand="perl",verbose=TRUE , conservation=input$ddcons
        , use_hmm=FALSE, datum=FALSE , updateProgress=updateProgress
        , clustal_cmd=clustalo_cmd)
    )
  })
  # Print the results as two separated tables
  output$ddalignedresult <- DT::renderDataTable({
    if( is.null(ddAnalysis()) )
      return(NULL)
    ddAnalysis()$AlignedSequence
  } , rownames=FALSE  ,options = list(searchHighlight = TRUE , 
    scrollX = TRUE,
    scrollCollapse = FALSE), filter = 'top')
  output$ddsingleseqresult <- DT::renderDataTable({
    if (is.null(ddAnalysis()))
      return(NULL)
    ddAnalysis()$SingleSequence
  } , rownames=FALSE  ,options = list(searchHighlight = TRUE,
    scrollX = TRUE,
    scrollCollapse = FALSE), filter = 'top')
  # observe({
  #   #Take a dependency on calculation
  #   ddAnalysis()
  #   # enable the download button
  #   shinyjs::enable("downloadDataAlign")
  #   sinyjs::enable("downloadDataSingle")
  # })
  # This command keep the download button unavailable until the computation is running
  # shinyjs::disable("data_file")

  observe({
    ddAnalysis()
    # notify the browser that the data is ready to download
    session$sendCustomMessage("download_ready" , list(mex=""))
  })


  #Download alignment results as tab delimited text
  output$downloadDataAlign <- downloadHandler(
    filename=function() {
      if(!as.logical(input$example))
        inFile <- input$myfile
      else
        inFile <- list(name="custom_analysis_example.txt" , datapath=file.path("data" , "custom_data" , "custom_analysis_example.txt"))
      tobepaste <- sub("\\....$","",inFile$name)
      paste(tobepaste , "AlignmentResults.txt" , sep="_")
    }
    , content=function(file) {
        write.table(ddAnalysis()$AlignedSequence
                , file=file
                , sep = "\t"
                , row.names = FALSE
                , quote=FALSE)
    }
  )
  #Download single sequence results as tab delimited text
  output$downloadDataSingle <- downloadHandler(
    filename=function() {
      if(!as.logical(input$example))
        inFile <- input$myfile
      else
        inFile <- list(name="custom_analysis_example.txt" , datapath=file.path("data" , "custom_data" , "custom_analysis_example.txt"))
      tobepaste <- sub("\\....$","",inFile$name)
      paste(tobepaste , "SingleSeqResults.txt" , sep="_")
    }
    , content=function(file) {
        write.table(ddAnalysis()$SingleSequence
                , file=file
                , sep = "\t"
                , row.names = FALSE
                , quote=FALSE)
    }
  )

  #Close the R process if the browser is closed
  session$onSessionEnded(function() {
        stopApp()
   })
})


