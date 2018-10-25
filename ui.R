#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
rm(list = ls())
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)
options(shiny.maxRequestSize = 6000*1024^2)
#.rs.restartR() # Restart R session
# check if pkgs are installed already, if not, install automatically:
source("installPkgsR.R")
source("SwitchButton.R")
source("helpers.R") # Load all the code needed to show feedback on a button click
selectSteps <- list("Quality control and cell filtering", "Gene variability across single cells", "Linear dimensional reduction (PCA)", "Non-linear dimensional reduction (UMAP/tSNE)", "Differentially expressed genes", "Discriminating marker genes")

header <- dashboardHeader(title = "iS-CellR",
                          tags$li(class = "dropdown",
                                  tags$a(href="http://www.immunocore.com", target="_blank", 
                                  tags$img(height = "20px", src="logo.png"),style = "text-align: center")
                          #tags$li(a(href = 'http://www.immunocore.com',
                          #            img(src = 'logo.png',
                          #                height = "50px"),
                          #            style = "text-align: center"),
                          #          class = "dropdown"))
                          ))

sidebar <- dashboardSidebar( uiOutput("userName"),
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             ),
                            sidebarMenu(
                                menuItem("Home", href = NULL, tabName = "home", icon = icon("home"), selected = T),
                                menuItem("Get data", icon = icon("folder-open"), tabName = "getdata"),
                                menuItem("Tools", tabName = "tools", icon = icon("cubes"), 
                                        menuSubItem("iS-CellR", tabName = "iS-CellR", icon = icon("dot-circle-o")))
                            ),
                            helpText("Developed by ", 
                                     a("Mitul Patel", href = ""), br(), a("Immunocore Limited", href = "https://immunocore.com/", target='_blank'),
                                     style = "padding-left:1em; padding-right:1em;position:absolute; bottom:1em; ")
                            )
                            
body <- dashboardBody(tags$head(includeCSS(rel = "stylesheet", type = "text/css", path = "style.css")),
                      tabItems(
                        tabItem(
                          tabName = "home",
                          fluidRow(
                            box(
                              title = h2("iS-CellR - Interactive platform for Single-cell RNAseq"), solidHeader = FALSE,
                              collapsible = TRUE, width = 12, status="info", 
                              includeMarkdown("about.md"),
                              #bsAlert("alertCITE"),
                              #bsButton("citeBton", label = "Mitulkumar V Patel (2018). iS-CellR: a user-friendly tool for analyzing and visualizing single-cell RNA sequencing data, Bioinformatics, 1-2. doi: 10.1093/bioinformatics/bty517. https://doi.org/10.1093/bioinformatics/bty517", block = TRUE, class = "citebutton", style = "success", icon = icon("bullhorn"), size = "default") 
                              bsButton("citeBton", label = HTML(paste("<strong>","Mitulkumar V Patel","</strong>", " (2018). iS-CellR: a user-friendly tool for analyzing and visualizing single-cell RNA sequencing data. ", paste("<em>", "Bioinformatics,", "</em>", sep="")," 1-2. doi: 10.1093/bioinformatics/bty517. <a href=https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty517/5048937, target='blank'>[Article]</a>", sep="")), block = TRUE, class = "citebutton", icon = icon("bullhorn"), type = "toggle", value = FALSE, size = "default")                             
                            ),
                            infoBoxOutput("pipeline"),
                            tags$div(id="BOX1",
                            bsModal(id="workflow", title="iS-CellR workflow", trigger="btn", tags$img(src='iS-CellR_workflow.png', align="center")),
                            tags$head(tags$style("#BOX1 .modal-body {padding: 10px}
                              #BOX1 .modal-content  {-webkit-border-radius: 6px !important;-moz-border-radius: 6px !important;border-radius: 6px !important;}
                              #BOX1 .modal-dialog { width: 65%; display: inline-block; text-align: left; vertical-align: top;}
                              #BOX1 .modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important;}
                              #BOX1 .modal { text-align: center; padding-right:10px; padding-top: 24px;}
                              #BOX1 .close { font-size: 18px}"
                            ))),
                            infoBoxOutput("tools"),
                            tags$head(tags$style(type = "text/css", "#ToolList th {display:none;}")),
                            tags$div(id="BOX2",
                            bsModal(id="wktools", title="iS-CellR workflow packages", trigger="btn", formattableOutput("ToolList")),
                            tags$head(tags$style("#BOX2 .modal-body {padding: 10px}
                              #BOX2 .modal-content  {-webkit-border-radius: 6px !important;-moz-border-radius: 6px !important;border-radius: 6px !important;}
                              #BOX2 .modal-dialog { width: 65%; display: inline-block; text-align: left; vertical-align: top;}
                              #BOX2 .modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important;}
                              #BOX2 .modal { text-align: center; padding-right:10px; padding-top: 24px;}
                              #BOX2 .close { font-size: 18px}"
                            ))),
                            infoBoxOutput("example"),
                            tags$div(id="BOX3",
                            bsModal(id="wkDemo", title="Run iS-CellR demo", trigger="btn", htmlOutput("DemoInfo")),
                            tags$head(tags$style("#BOX3 .modal-body {padding: 10px}
                              #BOX3 .modal-content  {-webkit-border-radius: 6px !important;-moz-border-radius: 6px !important;border-radius: 6px !important;}
                              #BOX3 .modal-dialog { width: 65%; display: inline-block; text-align: left; vertical-align: top;}
                              #BOX3 .modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important;}
                              #BOX3 .modal { text-align: center; padding-right:10px; padding-top: 24px;}
                              #BOX3 .close { font-size: 18px}"
                            )))
                              )
                            ),
                        #-------------------------------------------------------------#
                        tabItem(
                          tabName = "dashboard",
                          # infoBoxes with fill=FALSE
                          fluidRow(
                            box(
                              title = "Histogram", status = "primary", solidHeader = TRUE,
                              collapsible = TRUE,
                              plotOutput("plot3", height = 250)
                            ),
                            
                            box(
                              title = "Inputs", status = "warning", solidHeader = TRUE,
                              "Box content here", br(), "More box content",
                              sliderInput("slider", "Slider input:", 1, 100, 50),
                              textInput("text", "Text input:")
                            )
                          )
                        ),
                        #-------------------------------------------------------------#
                        tabItem(
                          tabName = "widgets",
                          fluidRow(
                            box( id ="greetbox",
                                 width  = 12, 
                                 height = "100%",
                                 solidHeader = TRUE, 
                                 status = "info",
                                 div(id="greeting", "Greeting here")
                          ))
                        ),
                        #-------------------------------------------------------------#
                          tabItem(
                            tabName = "getdata",
                            fluidRow(
                              box( 
                                useShinyjs(),
                                tags$style(appCSS),
                                title = "Choose file to upload",
                                collapsed = TRUE, collapsible = FALSE, width=6,
                                prettyRadioButtons(inputId = 'sep', label = 'Separator', 
                                             choices = c(Comma=',',
                                               Semicolon=';',
                                               Tab='\t'),
                                             thick = FALSE, shape = "round", fill = TRUE,
                                             animation = "pulse", status = "warning",
                                             selected = ',', inline=T),
                                prettyRadioButtons(inputId = 'quote', label = 'Quote',
                                             choices = c(None='',
                                               'Double Quote'='"',
                                               'Single Quote'="'"),
                                             thick = FALSE, shape = "round", fill = TRUE, 
                                             animation = "pulse", status = "warning",
                                             selected = '"', inline=T),
                                #tags$head(tags$style(type = "text/css", "#useheader {display:inline-block;}")),
                                fluidRow(column(3,
                                  prettyCheckbox(inputId = "header", label = "Header", icon = icon("check-square"), bigger = TRUE, outline = TRUE, plain = TRUE, animation = "pulse", status = "warning", value=TRUE)),
                                #checkboxInput('header', 'Header', TRUE),
                                column(3, prettyCheckbox(inputId = 'load10x', label = '10X genomics', bigger = TRUE, outline = TRUE, plain = TRUE, icon = icon("share-square"), status = "success", value=FALSE))),
                                switchButton(inputId = "SwitchUpload", label = "Local file", value = TRUE, col = "GB", type = "YN"),
                                useShinyjs(),
                                hidden(fileInput('file1', '',
                                          accept = c(
                                            'text/csv',
                                            'text/comma-separated-values',
                                            'text/tab-separated-values',
                                            'text/plain',
                                            '.csv',
                                            '.tsv',
                                            '.Rds',
                                            '.mtx'),
                                            multiple = TRUE
                                )),
                                hidden(shinyFilesButton('sfile', label = 'Browse...', title = 'Load file from sever', multiple = TRUE, icon=icon("upload")))
                                ),
                              box(
                                title = "Project summary",
                                collapsed = FALSE, collapsible = TRUE, width=6,
                                tags$head(tags$style(type = "text/css", "#DataSummary th {display:none;}")),
                                formattableOutput("DataSummary"),
                                tags$hr(),
                                useShinyjs(),  # Set up shinyjs
                                bsAlert("alert")
                              )
                              ),
                            fluidRow(
                            box(title = "Data", status = "primary", solidHeader = FALSE,
                                collapsed = FALSE, collapsible = TRUE, width=12,
                                #tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; -moz-box-shadow: none;box-shadow: none;}'))), 
                                uiOutput("DTshow"),
                                div(DT::dataTableOutput('dt') %>% withSpinner(type = 6, color="#bf00ff", size = 1))
                            )
                            )
                            ),
                            #-------------------------------------------------------------#
                            tabItem(
                              tabName = "iS-CellR",
                                fluidRow(tags$head(tags$style(".checkbox-inline {margin: 0 !important;} div.checkbox {margin-top: 10px;}")),
                                  box(title = "iS-CellR", status = "primary", solidHeader = TRUE,
                                      collapsed = FALSE, collapsible = TRUE, width=12,
                                      useShinyjs(),  # Include shinyjs
                                      box(id = "checkList", title="iS-CellR analysis steps", status = "warning", collapsible = FALSE, solidHeader = TRUE, width=6, 
                                      fluidRow(column(8,
                                      prettyCheckboxGroup(inputId = 'seuratSteps', label = "Select iS-CellR step to perform:", thick = TRUE, choices = selectSteps, animation = "pulse", status = "success"),
                                      actionButton("selectall", "Select All", class = "taskDFbutton")),
                                      column(4,
                                        useShinyjs(),
                                        prettyRadioButtons(inputId = 'clustLabels', label = 'Lables for clusters:', 
                                              choices = c("Labels from header"='useheader',
                                               "Custom labels"='customLabels',
                                               "Auto"='defLabels'), thick = FALSE, fill = TRUE, shape = "round",
                                             animation = "pulse", status = "success",
                                              selected = 'defLabels',inline = F))),
                                      useShinyjs(),  # Set up shinyjs
                                      tags$div(id="BOX4",
                                        bsModal(id="ClusterLabels", title="Define cluster labels", trigger = "customLabels", withSpinner(uiOutput("defineLable", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1), 
                                          tags$head(tags$style("#ClusterLabels .modal-footer{ display:none}")),
                                          footer = tagList(# modalButton("Cancel"),
                                                    div(style="display:inline-block", withBusyIndicatorUI(actionButton("changeLabels", "Apply", icon("user-circle"), disabled = TRUE,
                                                       style="float:left; color: #fff; background-color: #009933; border-color: #2e6da4"))))
                                          ),
                                        #wellPanel(actionButton("changeLables", "Apply"))),
                                        tags$head(tags$style("#BOX4 .modal-body {padding: 10px}
                                        #BOX4 .modal-content  {-webkit-border-radius: 6px !important;-moz-border-radius: 6px !important;border-radius: 6px !important;}
                                        #BOX4 .modal-dialog { width: 65%; display: inline-block; text-align: left; vertical-align: top;}
                                        #BOX4 .modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important;}
                                        #BOX4 .modal { text-align: center; padding-right:10px; padding-top: 24px;}
                                        #BOX4 .close { font-size: 18px}"
                                      )))),
                                      box(title="Running selected iS-CellR steps", status = "success", collapsible = FALSE, solidHeader = TRUE, width=6, 
                                      "Selected steps:", br(),
                                      textOutput("chsteps"), tags$hr(),
                                      #wellPanel(
                                      actionButton("seuratRUN", "Run iS-CellR", width = "165px", icon("play-circle"), class = "runbutton", disabled = TRUE),
                                      withBusyIndicatorUI(downloadButton("iSCellRreport", "Generate report", class = "actbutton", disabled = TRUE)),
                                      tags$style(type='text/css', "#seuratRUN { margin-bottom:10px;margin-right: -7px;}"),
                                      tags$style(type='text/css', "#iSCellRreport { margin-top:0px;margin-right: -7px;}"), br(),
                                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                      tags$div("Running iS-CellR...",id="loadmessage"),
                                                      img(src="ajax-loader-bar.gif")), run = "completed"
                                      )
                                    )
                                  ),
                              fluidRow(
                                uiOutput("seuratQC"),
                                uiOutput("seuratVargenes"),
                                uiOutput("seuratPCA"),
                                uiOutput("seurattSNE"),
                                uiOutput("seuratDiffExpr"),
                                uiOutput("seuratMakers"))
                      )#tabitem seurat
                            
    ))

ui <- dashboardPage(skin = "purple", header, sidebar, body)

