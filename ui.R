library(shinythemes)
library(shinydashboard)

# This function is a button disabler (triggered by some event)
tagfordisable <- singleton(tags$head(HTML(
'
  <script type="text/javascript">
    $(document).ready(function() {
      // disable download at startup. downloadDataSingle is the id of the downloadButton
      $("#downloadDataSingle").attr("disabled", "true").attr("onclick", "return false;");
      $("#downloadDataAlign").attr("disabled", "true").attr("onclick", "return false;");

      Shiny.addCustomMessageHandler("download_ready", function(message) {
        $("#downloadDataSingle").removeAttr("disabled").removeAttr("onclick").html(
          "<i class=\\"fa fa-download\\"></i> Download of Single Sequence Results is ready" + message.mex);
        $("#downloadDataAlign").removeAttr("disabled").removeAttr("onclick").html(
          "<i class=\\"fa fa-download\\"></i> Download of Aligned Results is ready" + message.mex);
      });
    })
  </script>
'
)))

# Define UI
shinyUI(
  dashboardPage(
    dashboardHeader(title = "LowMACA App")
    ,dashboardSidebar(
      #Logo
      imageOutput("myImage")
      ,sidebarMenu(
        menuItem("App Panel" , tabName="apppanel" , icon=icon("area-chart"))
        ,menuItem("Custom Analysis" , tabName="customanalysis" , icon=icon("area-chart"))
        ,menuItem("Pfam Reference" , tabName="pfamHelper" , icon=icon("question-sign", lib = "glyphicon"))
        ,menuItem("Gene/Protein Reference" , tabName="geneHelper" , icon=icon("question-sign", lib = "glyphicon"))
        ,menuItem("Help" , tabName="help" , icon=icon("info-sign" , lib="glyphicon"))
        ,menuItem("Contacts" , tabName="contact" , icon=icon("list"))
      )
    )
    ,dashboardBody(
      tabItems(
        tabItem("apppanel" , 
          fluidPage(
            theme = shinytheme("spacelab")
            ,tags$head(tags$script(src = 'http://d3js.org/d3.v3.min.js'))
            ,sidebarLayout(
              sidebarPanel(
                  #Pfam Selector
                  selectInput("pfam" 
                            , "Select A Pfam Domain:"
                            , selected=paste( start
                                      , unique( myPfam[ myPfam$Pfam_ID==start , "Pfam_Name"] ) , sep="_")
                            , choices = c("No Pfam" 
                                      , paste(pfam_list , unique( myPfam[ myPfam$Pfam_ID %in% pfam_list , "Pfam_Name"]) , sep="_")
                                      )
                            )
                  #Gene Selector
                  ,selectInput("genes", "Choose Your Genes:"
                            ,choices=c("All Genes" 
                                  , sort(unique(myPfam[ myPfam$Pfam_ID==start , "Gene_Symbol"])))
                            ,multiple=TRUE
                            ,selected="All Genes"
                            )
                  #Tumor Selector
                  ,selectInput("tumor_type" , "Choose A Tumor Type:"
                            ,choices=c("All Tumors" , tumor_type)
                            ,multiple=TRUE
                            ,selected="All Tumors"
                            )
                  #Mutation Type Selector
                  ,selectInput("mut_type" , "Choose A Mutation Type Category:"
                            ,choices=c("All Mutations" , "Missense Type" , "Truncating Type")
                            ,multiple=FALSE
                            ,selected="Missense Type"
                            )
                  #Bandwidth Controller
                  ,sliderInput("bw", 
                            "Bandwidth:", 
                             value = 0,
                             min = 0, 
                             max = 5,
                             step=.1)
                  #Conservation Controller
                  ,sliderInput("cons", 
                            "Trident Conservation Score Threshold:", 
                             value = 0.1,
                             min = 0.05, 
                             max = 1,
                             step=0.05)
                  # ,radioButtons("fullProtein", "Running Full Protein:",
                  #Choose if the you want to run in full protein or just pfam
                  ,selectInput("fullProtein", "Running Full Protein:",
                            ,choices=c("No" , "Yes")
                            ,multiple=FALSE
                            ,selected="No"
                            )
                  #Action Button
                  ,actionButton("goButton", "LowMACAize!")
              )
              ,mainPanel(
                textOutput("text1")
                ,tabsetPanel(type = "tabs"
                  ,tabPanel("LowMACA Plot" , plotOutput("plot" , height='1000px'))
                  ,tabPanel("LowMACA Barplot" , htmlOutput("bpAll") , style = "overflow-x:scroll")
                  ,tabPanel("LowMACA Network" ,  
                    sliderInput("CommonPositions", 
                            "Number of common position for Network Plot:", 
                             value = 1,
                             min = 1, 
                             max = 50,
                             step=1
                    )
                    ,htmlOutput("networkPlot") , style = "overflow-x:scroll")
                  ,tabPanel("Protter Plot" , imageOutput("protter"))
                  ,tabPanel("Significant Mutations Table", DT::dataTableOutput("table"))
                  ,tabPanel("Cooc - MutEx Analysis" ,
                        h2("Cooccurence and Mutual Exclusivity analysis of significant hotspots"),
                        fluidRow(
                          column(6,
                              h3("Position based Analysis"),
                              plotOutput("coocPos" , height="500px"),
                              hr(),
                              tableOutput("coocDataPos")
                              ),
                          column(6,
                              h3("Gene based Analysis"),
                              plotOutput("coocGene" , height="500px"),
                              hr(),
                              tableOutput("coocDataGene")
                              )
                              )
                  )
                )
              )
            )
          )
        )
        ,tabItem("customanalysis" , 
          fluidPage(
            tagfordisable 
            , conditionalDisabledPanel(element="goDDButton" , tagmessage="go_with_the_analysis") 
            , sidebarLayout(
              sidebarPanel(
              fileInput('myfile', 'Choose a tab delimited file',
                accept=c('text/comma-separated-values' , 'text/plain'))
              ,actionButton("example" , "Load Example Dataset")
              ,helpText("LowMACA accepts only tab delimited files with a specific format.")
              ,helpText("Check how the example was created inside data/custom_data/custom_analysis_example.txt")
              ,tags$hr()
              ,sliderInput("ddcons", 
                "Trident Conservation Score Threshold:", 
                value = 0.1,
                min = 0.05, 
                max = 1,
                step=0.05)
              ,selectInput("ddmuttype" , "Choose A Mutation Type Category:"
                ,choices=c("All Mutations" , "Missense Type" , "Truncating Type")
                ,multiple=FALSE
                , selected="Missense Type"
              )
              ,tags$hr()
              # ,conditionalPanel( condition = " input.myfile == null " , actionButton("goDDButton", "Go with the analysis"))
              ,actionButton("goDDButton", "Go with the analysis")
              ,tags$hr()
              ,tags$hr()
              # ,downloadButton('downloadDataAlign', 'Download Align Base Results')
              # ,downloadButton('downloadDataSingle', 'Download Single Sequence Results')
              ,downloadButton('downloadDataAlign')
              ,downloadButton('downloadDataSingle')
              ,helpText("Download will be available once the processing is completed.")              
            )
            ,mainPanel(
                tabsetPanel(
                  type = "tabs"
                  ,tabPanel("Uploaded Data" , fluidPage(DT::dataTableOutput("ddcontents") , textOutput("changingtext")))
                  ,tabPanel("Alignment Results" , DT::dataTableOutput("ddalignedresult"))
                  ,tabPanel("Single Sequence Results" , DT::dataTableOutput("ddsingleseqresult"))
                )
            )
          ))
        )
        ,tabItem("pfamHelper" , fluidPage(h3("Search your Genes, UNIPROT Proteins or Pfam Domain") , DT::dataTableOutput("PfamTable")))
        ,tabItem("geneHelper" , fluidPage(h3("Search your Genes, even if they don't have any domain") , DT::dataTableOutput("geneTable")))
        # ,tabItem("help" , fluidPage(htmlOutput("myVideo")))
        # ,tabItem("help" , fluidPage(imageOutput("cheatsheet") , htmlOutput("myVideo")))
        ,tabItem("help" , fluidRow(
                          column(4,
                              # h3("What's this app?"),
                              imageOutput("cheatsheet")
                              ),
                          column(8,
                              htmlOutput("myVideo")
                              )
                              )
                )
        ,tabItem("contact" 
          ,fluidPage(sidebarLayout(
            sidebarPanel(h2("LowMACA Contacts..."))
            ,mainPanel(
              p(h3("The LowMACA App was created by Giorgio Melloni at the Center for Genomic Science of IIT@SEMM, Milan, Italy")),
              p(h3("If you want to cite us, please use the cgsb website at" 
                  , span(a("CGSB website" , href="https://cgsb.genomics.iit.it/wiki/projects/LowMACA" , target="_blank") , style = "color:blue"))),
              p(h3("For any question or suggestion regarding the app, write to" 
                  , span("giorgio.melloni [a] iit.it" , style = "color:blue"))),
              p(h3("For our Bioconductor package, please refer to" 
                  , span(a("LowMACA Package" , href="http://www.bioconductor.org/packages/release/bioc/html/LowMACA.html" , target="_blank") , style = "color:blue"))),
              br(),
              p(h3(strong("Thank you for using our app :)")))
            ))))
      )
    )
  )
)