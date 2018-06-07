# 
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
rm(list = ls())
#.rs.restartR() # Restart R session
# check if pkgs are installed already, if not, install automatically:
source("installPkgsR.R")
options(shiny.maxRequestSize = 6000*1024^2)
source("helpers.R") # Load all the code needed to show feedback on a button click
selectSteps <- list("Quality control and cell filtering", "Gene variability across single cells", "Linear dimensional reduction (PCA)", "Non-linear dimensional reduction (tSNE)", "Differentially expressed genes", "Discriminating marker genes")

generateReport <- function() {
  out.markdown <- ''
  withProgress(message = 'Generating Report',
               value = 0, {
                 out.markdown <- rmarkdown::render(
                   input = "Tirosh2016_seurat.Rmd",
                   output_format = "html_document")
                 
                 setProgress(1)
               })
  
  read_file(out.markdown)
}

server <- function(input, output, session) {

  # Define variables for each user session
  source("setglobs.R")
  #source("sessionCommon.R", local=TRUE)


  dataModal <- function(failed = FALSE) {
    tags$div(id="loginBOX",
    modalDialog(
      title="Information",
      textInput("username", "User*"),
      textInput("projectid", "Project id"),
      textInput("projectname", "Project name"),
      textInput("sampleinfo", "Sample info"),
      #passwordInput("password", "Password:"),
      size = "s",
      fade = TRUE,
      footer = tagList(
        # modalButton("Cancel"),
        div(style="display:inline-block", actionButton("ok", "OK", icon("user-circle"),
                                                       style="color: #fff; background-color: #009933; border-color: #2e6da4"))
      )
     ),
    tags$head(tags$style("#loginBOX .modal-body {padding: 10px}
                        #loginBOX .modal-content {-webkit-border-radius: 8px !important;-moz-border-radius: 8px !important;}
                        #loginBOX .modal-dialog { display: inline-block; text-align: left; vertical-align: top;}
                        #loginBOX .modal-header {background-color: #337ab7; color: white; font-weight: bold !important;border-top-left-radius: 8px;border-top-right-radius: 8px;}
                        #loginBOX .modal { text-align: center; padding-left:10px; font-weight: bold !important; padding-right:10px; padding-top: 24px;}"
                       ))
    )
  }
  
  # Show modal when button is clicked.  
  # This `observe` is suspended only whith right user credential
  
  obs1 <- observe({
    showModal(dataModal())
  })
  
  # When OK button is pressed, attempt to authenticate. If successful,
  # remove the modal. 
  
  obs2 <- observeEvent(input$ok,{
    isolate({
      Username <- input$username
      values$ProjectID <- input$projectid
      values$ProjectName <- input$projectname
      values$SampleInfo <- input$sampleinfo 
      #Password <- input$password
    })
    if (input$username != "") {
      Logged <<- TRUE
      values$authentication <- Username
      obs1$suspend()
      removeModal()
    } #else {
      #Logged <<- TRUE
      #values$authentication <- Username
      #obs1$suspend()
      #removeModal()
    #}     
  })
  
  output[["userName"]] <- renderUI({
    sidebarUserPanel(values$authentication, subtitle = a(href = "#", icon("user", class = "text-success"), "Current user"))
  })
  
  isRequiredSTEP <- reactive({ # reactive part = thos code is repeated when user input changes
   
    validate( # define error messages if user doesn't choose anything
      need(input$seuratSteps != "",
           "Please select at least one step to run the pipeline"),
      need(input$file1 != "", "Please upload a count file")
    )
    
  })
 
  output$showbusy <-
    renderText({
      c(
        '<img align="middle", 
        src="',
        "ajax-loader-bar.gif",
        '">'
      )
    })

  #observeEvent(c(input$file1,input$sep,input$quote,input$header),{
   observe({
    if (!is.null(input$file1)) {
      createAlert(session, "alert", "info", title = "Reading file contents .....",
                  content = "Please wait...", style="warning", dismiss = FALSE, append = FALSE)
    }
  })

  # read in the data
  count_data <- reactive({
    #validate(need(input$file1 != "", "Please upload a count file"))
    inFile <- input$file1
    if (is.null(inFile)){ 
      #return(NULL)
      createAlert(session, "alert", "warning", title = "File Not Loaded Yet!",
                  content = "Please upload file...", style="danger", dismiss = FALSE, append = FALSE)
    } else {
    countMatrix <- read.delim(inFile$datapath, header = input$header,
                           sep = input$sep, quote = input$quote, stringsAsFactors=FALSE, row.names = 1)
    countFile$val <- countMatrix
    geneID <- rownames(countMatrix)
    GeneNames$val <- geneID
    #rownames(countMatrix) <- sapply(rownames(countMatrix), function(geneID) 
    #       toString(tags$a(href=paste0(link2Ens(geneID)), geneID, target="_blank")))

    if (!is.null(countMatrix)){ 
    createAlert(session, "alert", "success", title = "File Loading Complete",
                content = "You can now proceed...", style="success", dismiss = FALSE, append = FALSE)    
    countMatrix
    }
   }
  })
  
  output$dt <- DT::renderDataTable({
       #DT::datatable(count_data(), escape=FALSE, selection = 'none', options = list(scrollX = TRUE)) %>% formatStyle(0, cursor = 'pointer')
    DT::datatable(count_data(), escape=FALSE, selection = 'none', options = list(searchHighlight = TRUE, scrollX = TRUE))
  })
  
  observeEvent(c(input$file1,input$sep,input$quote,input$header),{
  output$DataSummary <- renderFormattable({
    if(length(input$file1) != 0) { 
      countFile <- data.matrix(countFile$val)
      #genes50 <- rowSums(countFile>0) # genes found in > 50 cells
      expr0 <- sum(colSums(countFile == 0))
      expr <- sum(colSums(countFile > 0))
      percentZero <- expr0*100/(expr0 + expr)
      ZeroExpr <- countFile[apply(countFile[, -1], MARGIN = 1, function(x) all(x == 0)), ] # Select all genes that has 0expr in all cells
      SummaryInfo$val$nCells <- ncol(countFile)
      SummaryInfo$val$nGenes <- nrow(countFile)
      SummaryInfo$val$nZero <- nrow(ZeroExpr)

      summaryProjectID <- data.table("Project ID :", values$ProjectID)
      summaryProjectName <- data.table("Project name :", values$ProjectName)
      summarySampleInfo <- data.table("Sample info :", values$SampleInfo)
      summaryCells <- data.table("Total cells :", ncol(countFile))
      summaryGenes <- data.table("Total genes :", nrow(countFile))
      #summarygenes50 <- data.table("Genes found in >50 cells", genes50)
      summaryZero <- data.table("Genes 0Expr :", nrow(ZeroExpr))
      
      summaryTab1 <- rbind(summaryProjectID, summaryProjectName, summarySampleInfo)
      summaryTab2 <- rbind(summaryCells, summaryGenes, summaryZero)
      summary <- cbind(summaryTab1, summaryTab2)
      colnames(summary) <- c("col1", "col2", "col3", "col4")
      align_column=c("l","l","l","r")

      formattable(summary, align=align_column, 
        list(col1 = formatter("span",
                          style = x ~ ifelse(x == c("Project ID :", "Project name :", "Sample info :"), style(color = "Gray", font.size = "14px", font.weight = "bold"), NA)),
        col2 = formatter(.tag = "span", style = function(x) style(font.size = "14px")),       
        col3 = formatter("span",
                          style = x ~ ifelse(x == c("Total cells :", "Total genes :", "Genes 0Expr :"), style(color = "Gray", font.size = "14px", font.weight = "bold"), NA)),
        #col4 = color_bar("orange")
        col4 = formatter("span", style = function(x) style(display = "bar", direction = "rtl", font.size = "15px", `border-radius` = "5px", `padding-left` = "4px", `padding-right` = "4px", width = paste(proportion(x),"px",sep=""), `background-color` = csscolor("#ffa31a")))
        ))
      #colnames(summary) <- c("", "", "", "")
      #summary
    }
  })
  })

 observeEvent(input$seuratSteps,{
  output$chsteps <- renderText({
    isRequiredSTEP()
    runsteps <- paste(input$seuratSteps, collapse = ", ")
    runsteps
    #Selected$steps <<- runsteps
  })
})
 ####### Enable/Disable UI elements #########
 
  observe({
   if(!is.null(input$seuratSteps) && !is.null(input$file1)){
     shinyjs::enable("seuratRUN")
   }
   else{
     shinyjs::disable("seuratRUN")
   }
 })

 observe({
   if((!is.null(input$seuratSteps) || !is.null(input$file1)) && (!is.null(DownloadPlot$val$PCAplot2D) || !is.null(DownloadPlot$val$tSNEplot2D) || !is.null(DownloadPlot$val$Joyplot) || !is.null(DownloadPlot$val$Vlnplot) || !is.null(DownloadPlot$val$Featureplot) || !is.null(DownloadPlot$val$Dotplot) || !is.null(DownloadPlot$val$Heatmap) || !is.null(DownloadPlot$val$CoExprplot) || !is.null(DownloadPlot$val$InterHet) || !is.null(DownloadPlot$val$IntraHetplotly))){
     shinyjs::enable("clustLabels")
   }
   else{
     shinyjs::disable("clustLabels")
   }
 })
 
 observe({
   if(!is.null(DownloadPlot$val$Summaryplot)){
     shinyjs::enable("summaryplotfull")
     shinyjs::enable("downsummaryplot")
     shinyjs::enable("ApplyFilter")
   }
   else{
     shinyjs::disable("summaryplotfull")
     shinyjs::disable("downsummaryplot")
     shinyjs::disable("ApplyFilter")
   }
 })
 
  observeEvent(input$ApplyFilter,{
     shinyjs::enable("RemoveFilter")
 })
 
 observe({
   if(!is.null(DownloadPlot$val$Summaryplot) || !is.null(DownloadPlot$val$PCAplot2D) || !is.null(DownloadPlot$val$Vizplot) || !is.null(DownloadPlot$val$tSNEplot2D) || !is.null(DownloadPlot$val$Joyplot) || !is.null(DownloadPlot$val$Vlnplot) || !is.null(DownloadPlot$val$Featureplot) || !is.null(DownloadPlot$val$Dotplot) || !is.null(DownloadPlot$val$Heatmap) || !is.null(DownloadPlot$val$CoExprplot) || !is.null(DownloadPlot$val$InterHet) || !is.null(DownloadPlot$val$IntraHetplot)){
     shinyjs::enable("iSCellRreport")
   }
   else{
     shinyjs::disable("iSCellRreport")
   }
 })
 
 observe({
   if(!is.null(DownloadPlot$val$PCAplot)){
     shinyjs::enable("PCA2Dfull")
     shinyjs::enable("downPCA2D")
     shinyjs::show("PrintLabelPCA")
     }
   else{
     shinyjs::disable("PCA2Dfull")
     shinyjs::disable("downPCA2D")
     shinyjs::hide("PrintLabelPCA")
   }
 })
 
 observe({
   if(!is.null(DownloadPlot$val$Vizplot)){
     shinyjs::enable("Vizplotfull")
     shinyjs::enable("downVizplot")
   }
   else{
     shinyjs::disable("Vizplotfull")
     shinyjs::disable("downVizplot")
   }
 })
 
 observe({
   if(!is.null(DownloadPlot$val$PCAHeatmap)){
     shinyjs::enable("PCAHeatmapfull")
     shinyjs::enable("downPCAHeatmap")
   }
   else{
     shinyjs::disable("PCAHeatmapfull")
     shinyjs::disable("downPCAHeatmap")
   }
 })
 
 observe({
   if(!is.null(DownloadPlot$val$tSNEplot)){
     shinyjs::enable("tSNE2Dfull")
     shinyjs::enable("downtSNE2D")
     shinyjs::show("PrintLabeltSNE")
   }
   else{
     shinyjs::disable("tSNE2Dfull")
     shinyjs::disable("downtSNE2D")
     shinyjs::hide("PrintLabeltSNE")
   }
 })
 
 observe({
   if(length(input$brush1_rows_all) > 0){
     shinyjs::enable("tSNEbrushdata")
   }
   else{
     shinyjs::disable("tSNEbrushdata")
   }
 })
 
 observe({
   if(length(input$JoyplotGenes) == 0){
     shinyjs::disable("updateJoyplot")
   }
   else{
     shinyjs::enable("updateJoyplot")
   }
 })
 observe({
   if(!is.null(DownloadPlot$val$Joyplot)){
     shinyjs::enable("Joyplotfull")
     shinyjs::enable("downJoyplot")
   }
   else{
     shinyjs::disable("Joyplotfull")
     shinyjs::disable("downJoyplot")
   }
 })
 observe({
   if(length(input$VlnplotGenes) == 0){
     shinyjs::disable("updateVlnplot")
   }
   else{
     shinyjs::enable("updateVlnplot")
   }
 })
 observe({
   if(!is.null(DownloadPlot$val$Vlnplot)){
     shinyjs::enable("Vlnplotfull")
     shinyjs::enable("downVlnplot")
   }
   else{
     shinyjs::disable("Vlnplotfull")
     shinyjs::disable("downVlnplot")
   }
 })
 
 observe({
   if(length(input$DotplotGenes) == 0){
     shinyjs::disable("updateDotplot")
   }
   else{
     shinyjs::enable("updateDotplot")
   }
 })
 observe({
   if(!is.null(DownloadPlot$val$Dotplot)){
     shinyjs::enable("Dotplotfull")
     shinyjs::enable("downDotplot")
   }
   else{
     shinyjs::disable("Dotplotfull")
     shinyjs::disable("downDotplot")
   }
 })
 
 observe({
   if(length(input$FeatureGenes) == 0){
     shinyjs::disable("updateFeatplot")
   }
   else{
     shinyjs::enable("updateFeatplot")
   }
 })
 observe({
   if(!is.null(DownloadPlot$val$Featureplot)){
     shinyjs::enable("Featureplotfull")
     shinyjs::enable("downFeatureplot")
     shinyjs::show("PrintLabel_Featplot")
   }
   else{
     shinyjs::disable("Featureplotfull")
     shinyjs::disable("downFeatureplot")
     shinyjs::hide("PrintLabel_Featplot")
   }
 })
 
 observe({
   if(length(input$HeatmapGenes) == 0){
     shinyjs::disable("updateHeatmap")
   }
   else{
     shinyjs::enable("updateHeatmap")
   }
 })
 observe({
   if(!is.null(DownloadPlot$val$Heatmap)){
     shinyjs::enable("Heatmapfull")
     shinyjs::enable("downHeatmap")
   }
   else{
     shinyjs::disable("Heatmapfull")
     shinyjs::disable("downHeatmap")
   }
 })
 
 observe({
   if(length(input$CoExprGenes) <= 1){
     shinyjs::disable("updateCoExprplot")
   }
   else{
     shinyjs::enable("updateCoExprplot")
   }
 })
 
 observe({
   if(!is.null(DownloadPlot$val$CoExprplot)){
     shinyjs::enable("CoExprplotfull")
     shinyjs::enable("downCoExprplot")
     shinyjs::show("PrintLabel")
     shinyjs::show("GeneNotFound_CoExpr")
   }
   else{
     shinyjs::disable("CoExprplotfull")
     shinyjs::disable("downCoExprplot")
     shinyjs::hide("PrintLabel")
     shinyjs::show("GeneNotFound_CoExpr")
   }
 })
 
 observe({
   if(!is.null(input$geneList)){
     shinyjs::enable("updateInterHetplot")
   }
   else{
     shinyjs::disable("updateInterHetplot")
   }
 })
 observe({
   if(!is.null(DownloadPlot$val$InterHet)){
     shinyjs::enable("downInterHet")
   }
   else{
     shinyjs::disable("downInterHet")
   }
 })
 
 observe({
   if(length(input$IntraHetGenes) == 0){
     shinyjs::disable("updateIntraHetplot")
   }
   else{
     shinyjs::enable("updateIntraHetplot")
   }
 })
 
 observe({
   if(!is.null(DownloadPlot$val$IntraHetplot)){
     shinyjs::enable("IntraHetplotfull")
     shinyjs::enable("downIntraHetplot")
   }
   else{
     shinyjs::disable("IntraHetplotfull")
     shinyjs::disable("downIntraHetplot")
   }
 })


 observe({
 if(input$selectall == 0) return(NULL)
   else if(input$selectall %% 2 == 0) 
    {
       #myChoices = list("PCA analysis" = "PCA", "tSNE analysis" = "tSNE", "Marker genes" = "Genes", "Features" = "Feat")
       updateCheckboxGroupInput(session,"seuratSteps","Select iS-CellR step to perform:", choices=selectSteps)
    } 
    else #if(input$selectall > 0)
    {
       #myChoices = list("PCA analysis" = "PCA", "tSNE analysis" = "tSNE", "Marker genes" = "Genes", "Features" = "Feat")
       updateCheckboxGroupInput(session,"seuratSteps","Select iS-CellR step to perform:", choices=selectSteps, selected=unlist(selectSteps))
    }
})  


observeEvent(input$clustLabels,{
    if(input$clustLabels == "customLabels") {
       toggleModal(session, modalId = "ClusterLabels", toggle = "open")
    }
  })

observeEvent(input$changeLabels,{
       toggleModal(session, modalId = "ClusterLabels", toggle = "close")
  })

  output$defCluster <- DT::renderDataTable({
       #DT::datatable(count_data(), escape=FALSE, selection = 'none', options = list(scrollX = TRUE)) %>% formatStyle(0, cursor = 'pointer')
    DT::datatable(ClusterLabInfo$val, escape=FALSE, selection = 'none', 
      editable = TRUE, extensions = 'KeyTable', options = list(searchHighlight = TRUE, scrollX = TRUE, keys = TRUE))
  })

  proxy2 = dataTableProxy('defCluster')
  
  #observeEvent(input$reset_barcodeDT, {
  #  withBusyIndicatorServer("reset_barcodeDT", {
  #  proxy2 %>% selectRows(NULL)
  #    })
  #})
  
  observeEvent(input$defCluster_cell_edit, {
    info = input$defCluster_cell_edit
      str(info)
      i = info$row
      j = info$col
      v = info$value
      df <- ClusterLabInfo$val
      df[i, 2] <- DT::coerceValue(v, df[i, 2])
      replaceData(proxy2, df, resetPaging = FALSE)  # important
      ClusterLabInfo$val <<- df
  })
  
observeEvent(input$defCluster_cell_edit,{
    df <- data.frame(ClusterLabInfo$val)
    is.na(df$new.cluster.ids) <- df$new.cluster.ids == ''

    if(any(is.na(df$new.cluster.ids))){
      shinyjs::disable("changeLabels")
    } else {
      shinyjs::enable("changeLabels")
    }
})

  observe({
    if(input$clustLabels == "customLabels") {
       output[["defineLable"]] <- renderUI({ 
        fluidRow(
          box(title = "Current clusters", status = "primary", solidHeader = FALSE,
                                      collapsed = FALSE, collapsible = FALSE, width=12,
            DT::dataTableOutput('defCluster')                          
            )
          )
        })
    }
 })

 #####################################################################################################
 ###### info boxes
 
  output$pipeline <- renderInfoBox({
    box1<-infoBox(
      "Pipeline",  
      HTML(paste("<b>Pipeline steps</b>",
                 "Find out more...",
                 sep = "<br/>")),
      icon = icon("gears"),
      href = "#",
      color = "aqua", fill = TRUE)
    box1$children[[1]]$attribs$class<-"action-button"
    box1$children[[1]]$attribs$id<-"button_box_01"
    return(box1)
  })

  observeEvent(input$button_box_01, {
    toggleModal(session,"workflow","open")
  })

  output$tools <- renderInfoBox({
    box2<-infoBox(
      "Tools",
      HTML(paste("Tools in pipeline",
                 "Find out more...",
                 sep = "<br/>")),
      icon = icon("cubes"),
      href = "#",
      color = "green", fill = TRUE)
    box2$children[[1]]$attribs$class<-"action-button"
    box2$children[[1]]$attribs$id<-"button_box_02"
    return(box2) 
  })

  observeEvent(input$button_box_02, {
    toggleModal(session,"wktools","open")
  })

  output$example <- renderInfoBox({
    box3<-infoBox(
      "Run demo",
      HTML(paste("Run demo",
                 "Find out more...",
                 sep = "<br/>")),
      icon = icon("edit"),
      href = "#",
      color = "yellow", fill = TRUE)
    box3$children[[1]]$attribs$class<-"action-button"
    box3$children[[1]]$attribs$id<-"button_box_03"
    return(box3) 
  })

  observeEvent(input$button_box_03, {
    toggleModal(session,"wkDemo","open")
  })

output$ToolList <- renderFormattable({
      pkgLIST <- data.matrix(sort(pkgs, decreasing=FALSE))
      p1 <- data.frame()
      p2 <- data.frame()
      
      for (i in seq(1,length(pkgLIST), 3)) {
        p1 <- cbind(pkgLIST[i], pkgLIST[i+1], pkgLIST[i+2])
        p2 <- rbind(p2, p1)
      }
      align_column=c("l","l","l")
      colnames(p2) <- c("col1","col2","col3")

      formattable(p2, align=align_column, 
        list(
        col1 = formatter("span", style = function(x) style(font.size = "14px")),
        col2 = formatter("span", style = function(x) style(font.size = "14px")),       
        col3 = formatter("span", style = function(x) style(font.size = "14px"))
        ))
  })

  output$DemoInfo <- renderText({
    paste("<b> Run demo with the data from Tirosh et al., 2016.</b><br/><br/>How to Get Started and Load Data<br/>
      It is important that the input file should follow the same format as descibed below.<br/><br/>
      To get started, please load in a CSV/TSV separated values file. The file should contain the count data where:<br/>
      1. Columns are Cells<br/>
      2. Rows are Genes<br/>
      3. Cell names or Cell type in Column header followed by Cell id. eg. Tcell_TC1_fcount_1.<br/><br/>
      <strong>NOTE:</strong> An example CSV file can be accessed <a href='data/Maligant50.csv', target='blank', download = 'Maligant50.csv'><strong>HERE</strong>.</a> (right click save as).<br/><br/>
      Notes on Plots and Processed Data<br/>
      All plots can be downloaded at a high resolution in PDF. This functionality works best in <strong>Chrome</strong> and <strong>Firefox</strong>. Within the Shiny app, user  can control the labelling of clusters.  User can also provide list of genes to compare their expression levels.<br/>
")
    })
  ############################## RUN Seurat and render plots ################################
  
  observeEvent(input$seuratRUN,{
    source("scripts/SeuratRUN.R", local=TRUE)
  })
  
  observeEvent(c(input$seuratRUN,input$ApplyFilter,input$RemoveFilter),{
  output$seuratresTab1 <- renderPlot({
    req(c(input$seuratRUN,input$ApplyFilter,input$RemoveFilter))
    isolate({
    source("scripts/variablegenes.R", local = TRUE)$value
    })
  })
  })
  observeEvent(c(input$seuratRUN,input$ApplyFilter,input$RemoveFilter),{
  output$summaryplot <- renderPlot({
    req(c(input$seuratRUN,input$ApplyFilter,input$RemoveFilter))
    isolate({
    source("scripts/SummaryPlot.R", local = TRUE)$value
    })
  })
  })
  
  observeEvent(input$ApplyFilter, {
    # When the button is clicked, wrap the code in a call to `withBusyIndicatorServer()`
    withBusyIndicatorServer("ApplyFilter", {
      Sys.sleep(1)
          source("scripts/QCrun.R", local = TRUE)
          source("scripts/SummaryPlot.R", local = TRUE)$value
     })
  })
  
  observeEvent(input$RemoveFilter,{
    # When the button is clicked, wrap the code in a call to `withBusyIndicatorServer()`
    withBusyIndicatorServer("RemoveFilter", {
      Sys.sleep(1)
      # if (input$select == "error") {
      #isolate({
        shinyjs::reset("QCoptions")
        source("scripts/SeuratRUN.R", local=TRUE)
        source("scripts/SummaryPlot.R", local = TRUE)$value
    })
  })
  
  output$summaryplotZoom <- renderPlot({
    req(input$summaryplotfull)
    isolate({
      print(DownloadPlot$val$Summaryplot)
      })
  })
  
  output$seuratVizPCA <- renderPlot({
    req(input$updateVizplot)
    isolate({
    source("scripts/PCAViZplot.R", local = TRUE)$value
    })
  })
  output$seuratVizPCAZoom <- renderPlot({
    req(input$Vizplotfull)
    isolate({
      print(DownloadPlot$val$Vizplot)
    })
  })
  
  output$seuratPCAHeatmap <- renderPlot({
    req(input$updatePCAHeatmap)
    #isolate({
      source("scripts/PCAHeatmap.R", local = TRUE)$value
    #})
  })
  output$seuratPCAHeatmapZoom <- renderPlot({
    req(input$PCAHeatmapfull)
    #isolate({
      print(DownloadPlot$val$PCAHeatmap)
    #})
  })

observeEvent(input$changeLabels,{
  output$seuratPCAplot <- renderPlotly({
    req(c(input$changeLabels,input$PrintLabelPCA))
      isolate({
          source("scripts/PCAplot.R", local = TRUE)$value
      })
  })
})

toListen <- eventReactive(input$clustLabels, {
#  toListen <- reactive({
  if(input$clustLabels == "defLabels" || input$clustLabels == "useheader"){
    list(input$seuratRUN,input$clustLabels,input$ApplyFilter,input$RemoveFilter) 
  } 
  })

  observeEvent(toListen(),{
  output$seuratPCAplot <- renderPlotly({
    obsList <- list()
    obsList[[length(toListen())+1]] <- input$PrintLabelPCA
      req(obsList)
      isolate({
          source("scripts/PCAplot.R", local = TRUE)$value
      })
  })
})
  
output$PCA2DplotZoom <- renderPlotly({
    req(input$PCA2Dfull)
    #isolate({
      print(DownloadPlot$val$PCAplot)
    #})
  })

  observeEvent(toListen(),{
  output$seuratPCAplot3D <- renderPlotly({
    obsList <- list()
    obsList[[length(toListen())+1]] <- input$changeLabels
    req(obsList)

    isolate({
      source("scripts/PCAplot3D.R", local = TRUE)$value
    })
  })
  })

  output$seuratPCAplot3DZoom <- renderPlotly({
    req(input$PCAfull)
    isolate({
      print(DownloadPlot$val$PCA3D)
    })
  })
  
  observeEvent(input$changeLabels,{
  output$seurattSNEplot <- renderPlotly({
    req(c(input$changeLabels,input$PrintLabeltSNE))
      isolate({
          source("scripts/tSNEplot.R", local = TRUE)$value
      })
  })
  })

  observeEvent(toListen(),{
  output$seurattSNEplot <- renderPlotly({
    obsList <- list()
    obsList[[length(toListen())+1]] <- input$PrintLabeltSNE
    req(obsList)

    isolate({
        source("scripts/tSNEplot.R", local = TRUE)$value
    })
  })
  })

  output$tSNE2DplotZoom <- renderPlotly({
    req(input$tSNE2Dfull)
    #isolate({
      print(DownloadPlot$val$tSNEplot)
    #})
  })

  observeEvent(input$changeLabels,{
  output$seurattSNEplot2 <- renderPlot({
    req(input$changeLabels)
    requireTSNE()
    #isolate({
    source("scripts/tSNEplot2.R", local = TRUE)$value
      })
    #})
  })

  observeEvent(toListen(),{
  output$seurattSNEplot2 <- renderPlot({
    req(toListen())
    requireTSNE()
    #isolate({
    source("scripts/tSNEplot2.R", local = TRUE)$value
      })
    #})
  })

  observeEvent(input$changeLabels,{
  output$seurattSNEplot3 <- renderPlot({
    req(input$changeLabels)
    requireTSNE()
    #isolate({
    source("scripts/tSNEplot3.R", local = TRUE)$value
    #  })
    })
  })

  observeEvent(toListen(),{
  output$seurattSNEplot3 <- renderPlot({
    req(toListen())
    requireTSNE()
    #isolate({
    source("scripts/tSNEplot3.R", local = TRUE)$value
    #  })
    })
  })

toListenDiff <- reactive({
    list(input$seuratRUN,input$clustLabels,input$changeLabels,input$ApplyFilter,input$RemoveFilter) 
  })

  observeEvent(toListen(),{
  output$seurattSNEplot3D <- renderPlotly({
    req(toListen())
    isolate({
      source("scripts/tSNEplot3D.R", local = TRUE)$value
    })
  })
  })
  observeEvent(input$changeLabels,{
  output$seurattSNEplot3D <- renderPlotly({
    req(input$changeLabels)
    isolate({
      source("scripts/tSNEplot3D.R", local = TRUE)$value
    })
  })
  })
  output$seurattSNEplot3DZoom <- renderPlotly({
    req(input$tSNEfull)
    isolate({
      print(DownloadPlot$val$tSNE3D)
    })
  })


observeEvent(toListenDiff(),{ 
withBusyIndicatorServer("changeLabels", {
Sys.sleep(1) 
    output$seuratJoyplot <- renderPlot({
      req(input$updateJoyplot)
        isolate({
          source("scripts/Joyplot_Feature.R", local = TRUE)$value
        }) 
    })
    output$seuratJoyplotZoom <- renderPlot({
      req(input$Joyplotfull)
        isolate({
          print(DownloadPlot$val$Joyplot)
        })
    })
  
    output$seuratVlnplot <- renderPlot({
      req(input$updateVlnplot)
        isolate({
          source("scripts/Vlnplot_Feature.R", local = TRUE)$value
        })
    })
    output$seuratVlnplotZoom <- renderPlot({
      req(input$Vlnplotfull)
        isolate({
          print(DownloadPlot$val$Vlnplot)
        })
    })

    output$seuratDotplot <- renderPlot({
      req(input$updateDotplot)
        isolate({
          source("scripts/Dotplot_Feature.R", local = TRUE)$value
        })
    })
    output$seuratDotplotZoom <- renderPlot({
      req(input$Dotplotfull)
        isolate({
          print(DownloadPlot$val$Dotplot)
        })
    })
  
    output$seuratFeatureplot <- renderPlotly({
      req(input$updateFeatplot)
        isolate({
          source("scripts/Featureplot_Feature.R", local = TRUE)$value
        })
    })
    output$seuratFeatureplotZoom <- renderPlotly({
      req(input$Featureplotfull)
        isolate({
          print(DownloadPlot$val$Featureplot)
        })
    })
  
    output$seuratHeatmap <- renderPlot({
      req(input$updateHeatmap)
        isolate({
          source("scripts/Heatmap_Feature.R", local = TRUE)$value
        })
    })
    output$seuratHeatmapZoom <- renderPlot({
      req(input$Heatmapfull)
        isolate({
          print(DownloadPlot$val$Heatmap)
        })
    })
#})
 # 
#observeEvent(toListenDiff(),{  
  CoExprPLOT_Err1 <- eventReactive(input$updateCoExprplot,{
    validate(
      need(DownloadPlot$val$CoExprplot, paste(unlist(GenesAbsent$val), "Not found. Please check the gene names.", sep = ": "))
    )
    })

  CoExprPLOT_Err2 <- eventReactive(input$updateCoExprplot,{
    validate(
      need(DownloadPlot$val$CoExprplot, "Expression threshold higher then maximum. Please check the expression threshold.")
    )
  })
  
  output$seuratCoExprplot <- renderPlotly({
    req(input$updateCoExprplot)
      isolate({
        source("scripts/CoExpr_Feature.R", local = TRUE)$value
      })
    
    if(length(GenesAbsent$val) > 0){
      CoExprPLOT_Err1()
    } else {
      CoExprPLOT_Err2()
      print(DownloadPlot$val$CoExprplot)
    }
  })

  output$seuratCoExprplotZoom <- renderPlotly({
    req(input$CoExprplotfull)
    isolate({
      print(DownloadPlot$val$CoExprplot)
    })
  })
  output$seuratInterHetplot <- renderPlot({
    req(input$updateInterHetplot)
    GeneList()
    isolate({
    source("scripts/InterHet_plot.R", local = TRUE)$value
    })
  })
  output$seuratIntraHetplot <- renderPlotly({
    requireGeneList()
    req(input$updateIntraHetplot)
    isolate({
    source("scripts/IntraHet_plot.R", local = TRUE)$value
    })
  })
  output$seuratIntraHetplotZoom <- renderPlotly({
    req(input$IntraHetplotfull)
    isolate({
      print(DownloadPlot$val$IntraHetplotly)
    })
  })
  })
})

observeEvent(c(input$changeLabels,input$PrintLabel_Featplot),{
output$seuratFeatureplot <- renderPlotly({
      req(c(input$changeLabels,input$PrintLabel_Featplot))
        isolate({
          source("scripts/Featureplot_Feature.R", local = TRUE)$value
        })
    })
})

CoExprPLOT_Err1 <- eventReactive(input$updateCoExprplot,{
    validate(
      need(DownloadPlot$val$CoExprplot, paste(unlist(GenesAbsent$val), "Not found. Please check the gene names.", sep = ": "))
    )
    })

  CoExprPLOT_Err2 <- eventReactive(input$updateCoExprplot,{
    validate(
      need(DownloadPlot$val$CoExprplot, "Expression threshold higher then maximum. Please check the expression threshold.")
    )
  })

observeEvent(c(input$changeLabels,input$PrintLabel),{
output$seuratCoExprplot <- renderPlotly({
    req(c(input$changeLabels,input$PrintLabel))
      isolate({
        source("scripts/CoExpr_Feature.R", local = TRUE)$value
      })
      if(length(GenesAbsent$val) > 0){
      CoExprPLOT_Err1()
    } else {
      CoExprPLOT_Err2()
      print(DownloadPlot$val$CoExprplot)
    }
  })
})


##################################################################################################
  
  observeEvent(c(input$nGenes1,input$Expr1),{
  output$numGenes <- renderText({
    countFiletmp <- as.data.frame(countFileQC$val)
    freq <- as.numeric((rowSums(countFiletmp > input$Expr1)*100)/ncol(countFiletmp)) # Calculate frequency of cells in which each gene expressed
    mean <- rowMeans(countFiletmp) # Calculate mean expression across all cells for each gene
    df.freqmean <- cbind(freq, mean)
    newdata <- df.freqmean[order(-freq),]
    topnGenes <- head(newdata, input$nGenes1) # select number of first n genes (df sorted descending)
    #d25 <- colSums(topnGenes >= 25) # Count genes which expressed in >25% of cells
    minpercent <- min(topnGenes[,1]) # Minumum percentage of cells in which input$nGenes genes expressed 
    ngenes <- paste0(input$nGenes1," genes are expressed in at least ", "<b>", format(round(minpercent, 2), nsmall = 2), "%</b>"," of cells with expression > ", input$Expr1)
    HTML(paste(ngenes))
    })
  })
  
  observeEvent(c(input$nGenes2,input$nCells,input$Expr2),{
    output$QCinfo <- renderText({
      # Genes must be present in n Cells
      countFiletmp <- as.data.frame(countFileQC$val)
      g.sub <- countFiletmp[apply(countFiletmp[, -1], MARGIN = 1, function(x) any(x > input$Expr2)), ]
      g.freq <- as.numeric(rowSums(g.sub > input$Expr2))
      g.mean <- rowMeans(g.sub)
      g.df <- cbind(g.freq,g.mean)
      g.df.count <- colSums(g.df >= input$nCells)
      
      # Cells must have n genes
      #c.sub <- countFiletmp[apply(countFiletmp[, -1], MARGIN = 1, function(x) any(x > input$Expr2)), ]
      c.freq <- as.numeric(colSums(g.sub > input$Expr2))
      c.mean <- colMeans(g.sub)
      c.df <- cbind(c.freq,c.mean)
      c.df.count <- colSums(c.df >= input$nGenes2)
      QCfilter$val$nGenes2 <- input$nGenes2
      QCfilter$val$nCells <- input$nCells
      QCfilter$val$Expr2 <- input$Expr2
    
      qcMSG <- paste("<b>",g.df.count[[1]],"</b>", " genes are expressed in at least ", input$nCells, " cells with expression >", input$Expr2,"<br/>", "<b>", c.df.count[[1]], "</b>", "cells have at least ", input$nGenes2, "genes with expression >", input$Expr2)

      QCfilter$val$qcMSG <- qcMSG
      
      HTML(paste(qcMSG))
    })
  })
  
    # read in the data GeneList1
  GeneList <- reactive({
    validate(need(input$geneList, "Please upload a Gene list file"))
    inFile <- input$geneList
    if (is.null(inFile)) return(NULL)
    GeneList <- read.delim(inFile$datapath, header = input$geneheader,
                           sep = input$genesep, stringsAsFactors=FALSE)
    GeneListglob$val <- GeneList
    GeneList
  })
  
  requireGeneList <- reactive({ # reactive part = thos code is repeated when user input changes
    validate( # define error messages if user doesn't choose anything
      need(input$geneList, "Please upload a Gene list file")
    )
  })
  
  requireTSNE <- reactive({ # reactive part = thos code is repeated when user input changes
    validate( # define error messages if user doesn't choose anything
      need(tSNEmatrix$val, "Please wait until tSNE finished...")
    )
  })
  
  tSNEplotDT <- reactive({
    brush1Data <- as.data.frame(tSNEmatrix$val)
    brush1Data$Celltype <- brush1Data$Cell
    brush1Data$Celltype <- gsub("\\..*|_.*|-.*", "", brush1Data$Celltype)
    brush1Data
  })
  
  output$brush1 <- DT::renderDataTable(brushedPoints(tSNEplotDT(), input$tSNEbrush, xvar = "tSNE_1", yvar = "tSNE_2"), options = list(scrollX = TRUE), escape = FALSE)

  output$tSNEplot3 <- DT::renderDataTable(tSNEplotDT(), options = list(scrollX = TRUE), server = FALSE, escape = FALSE, selection = 'none')

  output$tSNE_select <- renderText({
      #req(DownloadPlot$val$tSNEplot2)
        HTML(paste0("The data for each point on the plot can be viewed in the table next via selecting points on the plot. Just drag the cursor across points to make selection and corresponding data should be displayed in the table.</br></br>"))
     })
  
  output$tSNE_highlight <- renderText({
      #req(DownloadPlot$val$tSNEplot3)
      HTML(paste0("The points on the plot are highlighted in BLACK that corresponds to the data currently listed in table. The change in current listing datatable will automatically highlight corresponding points. The table can also be searched using the serachbox to highligh specific points.</br></br>"))
    })
    
  # download the brushed data
  observeEvent(input$tSNEbrush,{
  output$tSNEbrushdata <- downloadHandler(
    "tSNE_selected.txt",
    content = function(file) {
    rows = input$brush1_rows_all #download rows on all pages after being filtered
    #rows = input$brush1_rows_current # download rows on current page
    write.table(brushedPoints(tSNEplotDT(), input$tSNEbrush, xvar = "tSNE_1", yvar = "tSNE_2")[rows,  , drop = FALSE], file, sep="\t", row.names = F)
  })    
})

observe({
  output$header <- renderText({
    if(length(CoExprValue$val$Gene1) > 0 && length(CoExprValue$val$Gene2) > 0){
      HTML(paste0("<h4>","Gene expression values","</h4>"))
    }
  })
  
  output$Gene12ExprTable <- renderTable({
    if(length(CoExprValue$val$Gene1) > 0 && length(CoExprValue$val$Gene2) > 0){ 
      outputGene1 <- data.table(CoExprValue$val$Gene1, CoExprValue$val$Gene1Min, CoExprValue$val$Gene1Max, CoExprValue$val$Gene1Mean)
      outputGene2 <- data.table(CoExprValue$val$Gene2, CoExprValue$val$Gene2Min, CoExprValue$val$Gene2Max, CoExprValue$val$Gene2Mean)
      outputGene12 <- rbind(outputGene1, outputGene2)
      colnames(outputGene12) <- c("", "Min", "Max", "Mean")
      outputGene12
    }
  })
})
  
  observeEvent(input$updateJoyplot,{
     selectedGenes <- input$JoyplotGenes
    output$GeneNotFound_Joy <- renderValueBox({
      if(is.null(GenesAbsent$val)) {
          valueBox(
            unlist(length(selectedGenes)), HTML(paste("<strong>","Success .!!  ","</strong>"," All Genes found!")), icon = icon("thumbs-up", "fa-6x")
          ,color="green")
          } else if(length(GenesAbsent$val) == length(selectedGenes)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Failed .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("ban", "fa-6x")
          ,color="red")  
          } else if(!is.null(GenesAbsent$val) && !is.null(DownloadPlot$val$Joyplot)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Warning .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("exclamation-triangle", "fa-6x")
          ,color="yellow")
          }
       })
  })
        
  observeEvent(input$updateVlnplot,{
    selectedGenes <- input$VlnplotGenes
   output$GeneNotFound_Vln <- renderValueBox({
      if(is.null(GenesAbsent$val)) {
          valueBox(
            unlist(length(selectedGenes)), HTML(paste("<strong>","Success .!!  ","</strong>"," All Genes found!")), icon = icon("thumbs-up")
          ,color="green")
          } else if(length(GenesAbsent$val) == length(selectedGenes)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Failed .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("ban", "fa-6x")
          ,color="red")  
          } else if(!is.null(GenesAbsent$val) && !is.null(DownloadPlot$val$Vlnplot)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Warning .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("exclamation-triangle", "fa-6x")
          ,color="yellow")
          }
       })
  })
  
  observeEvent(input$updateFeatplot,{
    selectedGenes <- input$FeatureGenes
    output$GeneNotFound_Feat <- renderValueBox({
      if(is.null(GenesAbsent$val)) {
          valueBox(
            unlist(length(selectedGenes)), HTML(paste("<strong>","Success .!!  ","</strong>"," All Genes found!")), icon = icon("thumbs-up")
          ,color="green")
          } else if(length(GenesAbsent$val) == length(selectedGenes)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Failed .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("ban", "fa-6x")
          ,color="red")  
          } else if(!is.null(GenesAbsent$val) && !is.null(DownloadPlot$val$Featureplot)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Warning .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("exclamation-triangle", "fa-6x")
          ,color="yellow")
          }
       })
  })
  
  observeEvent(input$updateDotplot,{
    selectedGenes <- input$DotplotGenes
   output$GeneNotFound_Dot <- renderValueBox({
      if(is.null(GenesAbsent$val)) {
          valueBox(
            unlist(length(selectedGenes)), HTML(paste("<strong>","Success .!!  ","</strong>"," All Genes found!")), icon = icon("thumbs-up")
          ,color="green")
          } else if(length(GenesAbsent$val) == length(selectedGenes)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Failed .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("ban", "fa-6x")
          ,color="red")  
          } else if(!is.null(GenesAbsent$val) && !is.null(DownloadPlot$val$Dotplot)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Warning .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("exclamation-triangle", "fa-6x")
          ,color="yellow")
          }
       })
  })
  
  observeEvent(input$updateHeatmap,{
    selectedGenes <- input$HeatmapGenes
    output$GeneNotFound_Heat <- renderValueBox({
      if(is.null(GenesAbsent$val)) {
          valueBox(
            unlist(length(selectedGenes)), HTML(paste("<strong>","Success .!!  ","</strong>"," All Genes found!")), icon = icon("thumbs-up")
          ,color="green")
          } else if(length(GenesAbsent$val) == length(selectedGenes)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Failed .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("ban", "fa-6x")
          ,color="red")  
          } else if(!is.null(GenesAbsent$val) && !is.null(DownloadPlot$val$Heatmap)) {
          valueBox(
            unlist(length(GenesAbsent$val)), HTML(paste("<strong>","Warning .!!  ","</strong>", paste("<em>", unlist(GenesAbsent$val), "</em>", collapse = ", "), " NOT FOUND")), icon = icon("exclamation-triangle", "fa-6x")
          ,color="yellow")
          }
       })
  })

  observeEvent(input$updateCoExprplot,{
    selectedGenes <- input$CoExprGenes
    gene1 <- CoExprValue$val$Gene1
    gene2 <- CoExprValue$val$Gene2
    #Expr <- input$ExprCutoff
   output$GeneNotFound_CoExpr <- renderValueBox({
      if(!is.null(GenesAbsent$val)) {
          valueBox(
            unlist(GenesAbsent$val), HTML(paste("<strong>","Failed .!!  ","</strong>", "Genes not found")), icon = icon("ban", "fa-6x")
          ,color="red")  
          } else if(CoExprValue$val$Gene1Max < input$ExprCutoff || !is.null(CoExprValue$val$Gene1Ab)){
          valueBox(
            CoExprValue$val$Gene1, HTML(paste("<strong>","Warning .!!  ","</strong>", " Expression threshold too high")), icon = icon("ban", "fa-6x")
          ,color="yellow")  
          } else if(CoExprValue$val$Gene2Max < input$ExprCutoff || !is.null(CoExprValue$val$Gene2Ab)){
          valueBox(
            CoExprValue$val$Gene2, HTML(paste("<strong>","Warning .!!  ","</strong>", " Expression threshold too high")), icon = icon("ban", "fa-6x")
         ,color="yellow")  
          } else if(!is.null(CoExprValue$val$Genes)){
          valueBox(
            paste(CoExprValue$val$Genes[1],",",CoExprValue$val$Genes[2]), HTML(paste("<strong>","Warning .!!  ","</strong>", " Co-expression threshold too high")), icon = icon("ban", "fa-6x")
         ,color="yellow")  
          }
       })
  })
  
  observeEvent(input$updateIntraHetplot,{
  output$header_celltypes <- renderText({
    req(input$updateIntraHetplot)
    if(!is.null(GenesAbsent$val)) {
      HTML(paste0("<h4>","Sample id not found","</h4>"))
    } else { 
      HTML(paste0("<h4>","All sample id found","</h4>"))
    }
  })
  output$CelltypeNotFound <- renderText({
    ngenes <- paste(unlist(GenesAbsent$val), collapse = ", ")
    ngenes
  })
  })
  

  #############################################################################################
  ###### Download plts as PDF
  
  ## call the plot function when downloading the image
  output$downsummaryplot <- downloadHandler(
    filename =  function() {
      paste('SummaryPlot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$Summaryplot)
      dev.off()
    } 
  )

  output$downPCA2D <- downloadHandler(
    filename =  function() {
      paste('PCAplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$PCAplot2D)
      dev.off()
    } 
  )
  
  output$downtSNE2D <- downloadHandler(
    filename =  function() {
      paste('tSNEplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$tSNEplot2D)
      dev.off()
    } 
  )
  
   output$downInterHet <- downloadHandler(
    filename =  function() {
      paste('InterHetplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$InterHet)
      dev.off()
    } 
  )
  output$downCoExprplot <- downloadHandler(
    filename =  function() {
      paste('CoExprplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$CoExprplot1)
      dev.off()
    } 
  )
  
  output$downHeatmap <- downloadHandler(
    filename =  function() {
      paste('Features_Heatmap', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$Heatmap)
      dev.off()
    } 
  )
  output$downJoyplot <- downloadHandler(
    filename =  function() {
      paste('Features_Joyplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$Joyplot)
      dev.off()
    } 
  )
  output$downDotplot <- downloadHandler(
    filename =  function() {
      paste('Features_Dotplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$Dotplot)
      dev.off()
    } 
  )
  output$downFeatureplot <- downloadHandler(
    filename =  function() {
      paste('Featuresplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=16, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$Featureplot2D)
      dev.off()
    } 
  )
  output$downVlnplot <- downloadHandler(
    filename =  function() {
      paste('Features_Vlnplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$Vlnplot)
      dev.off()
    } 
  )
  output$downIntraHetplot <- downloadHandler(
    filename =  function() {
      paste('IntraHetplot', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=10, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$IntraHetplot)
      dev.off()
    } 
  )
  output$downVizplot <- downloadHandler(
    filename =  function() {
      paste('VizplotPCA', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=15, height=8, paper='special', pointsize=16, family="Helvetica")
      print(DownloadPlot$val$Vizplot)
      dev.off()
    } 
  )
  
  output$downPCAheatmap <- downloadHandler(
    filename =  function() {
      paste('PCAHeatmap', 'pdf', sep='.')
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      pdf(file, width=15, height=8, paper='special', pointsize=16, family="Helvetica")
      plot.new()
      print(DownloadPlot$val$PCAHeatmap)
      dev.off()
    } 
  )
  
  ###########################################################################################
  
  output$iSCellRreport <- downloadHandler(
    filename = function() {
      paste('iS-CellR-report', sep = '.', 'html')
      #switch(
      #        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      #))
    },
    content = function(file) {
      if (file.exists('iS-CellR-report.html')) file.remove('iS-CellR-report.html')
      src <- normalizePath('iS-CellR-report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'iS-CellR-report.Rmd', overwrite = TRUE)
      withBusyIndicatorServer("iSCellRreport", {
        Sys.sleep(1)
      out <- rmarkdown::render('iS-CellR-report.Rmd',
                               #params = list(text = input$text),
                               envir = new.env(parent = globalenv()
                               #switch(input$format,
                              #        PDF = pdf_document(), 
                               #       HTML = html_document(), 
                              #        Word = word_document()
                               ))
      file.rename(out, file)
      })
    }
  )

  ###########################################################################################
  
  hideBoxes <- reactive({
    perm.vector <- as.vector(input$seuratSteps)
    perm.vector
  }) 

################# Dropdown Gene list ###################### 

output$choose_JoyplotGenes <- reactiveUI(function() {

  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 16px;} .selectize-dropdown { font-size: 12px; line-height: 12px; }")
        selectizeInput(inputId = "JoyplotGenes", label = "Genes of interest", multiple = TRUE, choices = NULL,
        options = list(placeholder = 'Please select gene', maxOptions = 1000))
})

output$choose_VlnplotGenes<- reactiveUI(function() {

  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 16px;} .selectize-dropdown { font-size: 12px; line-height: 12px; }")
        selectizeInput(inputId = "VlnplotGenes", label = "Genes of interest", multiple = TRUE, choices = NULL,
        options = list(placeholder = 'Please select gene', maxOptions = 1000))                          
})

output$choose_DotplotGenes <- reactiveUI(function() {

  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 16px;} .selectize-dropdown { font-size: 12px; line-height: 12px; }")
        selectizeInput(inputId = "DotplotGenes", label = "Genes of interest", multiple = TRUE, choices = NULL,
        options = list(placeholder = 'Please select gene', maxOptions = 1000))                        
})

output$choose_FeatureGenes <- reactiveUI(function() {

  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 16px;} .selectize-dropdown { font-size: 12px; line-height: 12px; }")
        selectizeInput(inputId = "FeatureGenes", label = "Genes of interest", multiple = TRUE, choices = NULL,
          options = list(placeholder = 'Please select gene', maxOptions = 1000))
})

output$choose_HeatmapGenes <- reactiveUI(function() {

  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 16px;} .selectize-dropdown { font-size: 12px; line-height: 12px; }")
        selectizeInput(inputId = "HeatmapGenes", label = "Genes of interest", multiple = TRUE, choices = NULL,
          options = list(placeholder = 'Please select gene', maxOptions = 1000))
})

output$choose_CoExprGenes <- reactiveUI(function() {

  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 16px;} .selectize-dropdown { font-size: 12px; line-height: 12px; }")
        selectizeInput(inputId = "CoExprGenes", label = "Genes of interest", multiple = TRUE, choices = NULL,
          options = list(placeholder = 'Please select two genes', maxItems = 2))                          
})

observe({
 onclick('JoyplotGenes', function(){ 
  if(is.null(input$JoyplotGenes)){
    updateSelectizeInput(session, 'JoyplotGenes', choices = c(unlist(GeneNames$val)), server=TRUE)
    }
  })
 
 onclick('VlnplotGenes', function(){ 
  if(is.null(input$VlnplotGenes)){
    updateSelectizeInput(session, 'VlnplotGenes', choices = c(unlist(GeneNames$val)), server=TRUE)
    }
  })

 onclick('DotplotGenes', function(){ 
  if(is.null(input$DotplotGenes)){
    updateSelectizeInput(session, 'DotplotGenes', choices = c(unlist(GeneNames$val)), server=TRUE)
    }
  })

 onclick('FeatureGenes', function(){ 
  if(is.null(input$FeatureGenes)){
    updateSelectizeInput(session, 'FeatureGenes', choices = c(unlist(GeneNames$val)), server=TRUE)
    }
  })

 onclick('HeatmapGenes', function(){ 
  if(is.null(input$HeatmapGenes)){
    updateSelectizeInput(session, 'HeatmapGenes', choices = c(unlist(GeneNames$val)), server=TRUE)
    }
  })

 onclick('CoExprGenes', function(){ 
  if(is.null(input$CoExprGenes)){
    updateSelectizeInput(session, 'CoExprGenes', choices = c(unlist(GeneNames$val)), server=TRUE)
    }
  })
})

################# Dropdown Gene list ###################### 


observeEvent(input$seuratRUN,{
    listSteps <- hideBoxes()
    if("Features" %in% listSteps){
      #hideElement(id = "QCbox1", anim = TRUE)
      #shinyjs::disable("box1")
    }
    else
    {
      #showElement(id = "QCbox1")
    }
  })
  
  observe({
    if(input$seuratRUN == 0) return() 
    isolate({
      output[["seuratQC"]] <- renderUI({
        listSteps <- hideBoxes()
        if("Quality control and cell filtering" %in% isolate(listSteps)){
            box(id = "QCbox1", title = "iS-CellR results: quality control and cell filtering", status = "success", solidHeader = FALSE,
                collapsed = TRUE, collapsible = TRUE, width=12,
                tabBox(
                  # The id lets us use input$tabset1 on the server to find the current tab
                  id = "seuratQC", width=12,
                  tabPanel("Summary plot", status = "success",
                  fluidRow(
                    box(
                      title = "Summary from plot", status = "info", solidHeader = FALSE,
                      collapsed = FALSE, collapsible = TRUE, width=5,
                      fluidRow(
                      column(6, numericInput("nGenes1", "Number of genes", 
                                               min = 1, max = length(GeneNames$val), value = 100, step = 500, width = "120px")),
                      column(6, numericInput("Expr1", "Expression value", 
                                               min = 0, value = 0, step = 0.5, width = "120px"))),
                      fluidRow(
                        column(12, htmlOutput("numGenes")))),
                    box(
                      title = "QC and filtering", status = "warning", solidHeader = FALSE,
                      collapsed = FALSE, collapsible = TRUE, width=7,
                      fluidRow(useShinyjs(), id = "QCoptions",
                      column(4, numericInput("nCells", "Minimum cells", width = "120px", 
                                             min = 1, max = length(GeneNames$val), value = 3, step = 25)),
                      column(4, numericInput("nGenes2", "Minimum genes", 
                                             min = 1, max = length(GeneNames$val), value = 3, step = 500, width = "120px")),
                      column(4, numericInput("Expr2", "Expression value", 
                                             min = 0, value = 0, step = 0.5, width = "120px"))),
                      fluidRow(useShinyjs(),
                        column(4, 
                      withBusyIndicatorUI(actionButton("ApplyFilter", "Apply filter", class = "taskQC1button", width = "113px", disabled = TRUE)),
                      withBusyIndicatorUI(actionButton("RemoveFilter", "Remove filter", class = "taskQC2button", disabled = TRUE))),
                      column(8, htmlOutput("QCinfo"))),
                      tags$style(type='text/css', "#ApplyFilter { margin-bottom: 10px; margin-right: -7px;}"),
                      tags$style(type='text/css', "#RemoveFilter { margin-right: -7px;}"))), 
                  useShinyjs(),
                  actionButton("summaryplotfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE), 
                  downloadButton("downsummaryplot", label = "Save plot", class = "taskDFbutton", disabled = TRUE), tags$hr(),
                  withSpinner(plotOutput("summaryplot", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1),
                  tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                  tags$head(tags$style(".modal{height:95%;color:purple;}")),
                  tags$head(tags$style(".modal-dialog{width:95%}")),
                  tags$head(tags$style(".modal-body{min-height:95%}")),
                  tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                  bsModal("summaryplotZOOM", "Summary plot", "summaryplotfull", withSpinner(plotOutput("summaryplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))
                  )
              
                  #tabPanel("Geneplot", status = "success", plotOutput("seuratGeneplot", width = "100%", height = "500px"))
                )
            )
        }
      })
    #}
  #else if(input$seuratRUN > 0){
  #observeEvent(input$seuratRUN,{
      output[["seuratVargenes"]] <- renderUI({
        listSteps <- hideBoxes()
        if("Gene variability across single cells" %in% listSteps){
          box(title = "iS-CellR results: variable genes across the single cells", status = "danger", solidHeader = FALSE,
              collapsed = TRUE, collapsible = TRUE, width=12,
              tabBox(
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "seuratVargenes", width=12,
                tabPanel("Variability in gene expresion", status = "danger", 
                         withSpinner(plotOutput("seuratresTab1", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1))
                #downloadButton("seuratTab1", "Download plot")),
              ))
        }
      })
  #}
   #  else if(input$seuratRUN > 0){
  #observeEvent(input$seuratRUN,{
      output[["seuratPCA"]] <- renderUI({
        listSteps <- hideBoxes()
        if("Linear dimensional reduction (PCA)" %in% listSteps){
          box(title = "iS-CellR results: linear dimensional reduction (PCA)", status = "info", solidHeader = FALSE,
              collapsed = TRUE, collapsible = TRUE, width=12,
              tabBox(
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "seuratPCA", width=12, 
                #conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                #                tags$div("Running Seurat...",id="loadmessage"),
                #                img(src="spinner.gif")),
                tabPanel("PCA 2D", status = "info",
                         actionButton("PCA2Dfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE), 
                         downloadButton("downPCA2D", label = "Save plot", class = "taskDFbutton", disabled = TRUE), 
                         hidden(checkboxInput('PrintLabelPCA', 'Label clusters', value = FALSE)), tags$hr(),
                         withSpinner(plotlyOutput("seuratPCAplot", width = "100%", height = "600px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("PCA2DplotZOOM", "PCA 2D", "PCA2Dfull", withSpinner(plotlyOutput("PCA2DplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))),
                tabPanel("PCA 3D", status = "warning", 
                         #actionButton("PCAplot3D", "Plot 3D", icon("hand-o-right")),
                         conditionalPanel(condition = "output.seuratPCAplot3D",
                                          actionButton("PCAfull", "Full screen", icon("desktop"), class = "taskDFbutton")),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("PCAplot3DZOOM", "PCA plot", "PCAfull", withSpinner(plotlyOutput("seuratPCAplot3DZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1)),
                         withSpinner(plotlyOutput("seuratPCAplot3D", width = "100%", height = "600px"), type = 6, color="#bf00ff", size = 1)),
                tabPanel("PCA Vizplot", status = "info", #htmlOutput("showbusyVizPCA"), 
                         fluidRow(column(2,
                                         numericInput("VizPCAuse1", "PCA to plot:",
                                                      min = 1, max = 12, value = 1, step = 1, width = "120px"),
                                         numericInput("VizPCAuse2", "to",
                                                      min = 1, max = 12, value = 3, step = 1, width = "120px")),
                                  column(2,
                                         numericInput("PCAvizGenes", "Number of Genes", width = "150px", min = 1, value = 30, step = 5),
                                         numericInput("PCAvizfont", "Label size", width = "150px", min = 0, value = 10, step = 1))),
                                         useShinyjs(),
                                         actionButton("updateVizplot", "Plot it", icon("hand-o-right"),
                                                      class = "taskRUNbutton"), 
                                         actionButton("Vizplotfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE), 
                                         downloadButton("downVizplot", label = "Save plot", class = "taskDFbutton", disabled = TRUE), tags$hr(),
                         withSpinner(plotOutput("seuratVizPCA", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("VizplotZOOM", "PCA Vizplot", "Vizplotfull", withSpinner(plotOutput("seuratVizPCAZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))),
                tabPanel("PCA heatmap", status = "info", 
                         fluidRow(column(2,
                                         numericInput("PCAHeatuse1", "PCA to plot:",
                                                      min = 1, max = 12, value = 1, step = 1, width = "120px")),
                                  column(2,
                                         numericInput("PCAHeatCells", "Number of Cells", width = "150px", min = 1, value = 500, step = 10)),
                                  column(2,
                                         numericInput("PCAHeatfont", "Label size", width = "150px", min = 0, value = 10, step = 1))),
                                         useShinyjs(),
                                         actionButton("updatePCAHeatmap", "Plot it", icon("hand-o-right"),
                                                      class = "taskRUNbutton"), 
                                         actionButton("PCAHeatmapfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE), 
                                         downloadButton("downPCAHeatmap", label = "Save plot", class = "taskDFbutton", disabled = TRUE), tags$hr(),
                         withSpinner(plotOutput("seuratPCAHeatmap", width = "100%", height = "650px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("PCAHeatmapZOOM", "PCA Heatmap", "PCAHeatmapfull", withSpinner(plotOutput("seuratPCAHeatmapZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1)))
                )
              )
        }
      })
    # }
    #else if(input$seuratRUN > 0){
 # observeEvent(input$seuratRUN,{
      output[["seurattSNE"]] <- renderUI({
        listSteps <- hideBoxes()
        if("Non-linear dimensional reduction (tSNE)" %in% listSteps){
          box(title = "iS-CellR results: non-linear dimensional reduction (tSNE)", status = "warning", solidHeader = FALSE,
              collapsed = TRUE, collapsible = TRUE, width=12,
              tabBox(
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "seurattSNE", width=12, 
                tabPanel("tSNE 2D", status = "warning", 
                         actionButton("tSNE2Dfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE), 
                         downloadButton("downtSNE2D", label = "Save plot", class = "taskDFbutton", disabled = TRUE), 
                         hidden(checkboxInput('PrintLabeltSNE', 'Label clusters', value = FALSE)), tags$hr(),
                         withSpinner(plotlyOutput("seurattSNEplot", width = "100%", height = "600px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("tSNE2DplotZOOM", "tSNE 2D", "tSNE2Dfull", withSpinner(plotlyOutput("tSNE2DplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))),
                tabPanel("tSNE 3D", status = "warning", 
                         #actionButton("tSNEplot3D", "Plot 3D", icon("hand-o-right")),
                         conditionalPanel(condition = "output.seurattSNEplot3D",
                                          actionButton("tSNEfull", "Full screen", icon("desktop"), class = "taskDFbutton")),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("tSNEplot3DZOOM", "tSNE plot", "tSNEfull", withSpinner(plotlyOutput("seurattSNEplot3DZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1)),
                         withSpinner(plotlyOutput("seurattSNEplot3D", width = "100%", height = "600px"), type = 6, color="#bf00ff", size = 1)),
                tabPanel("tSNE select", status = "warning",
                         fluidRow(column(12, htmlOutput("tSNE_select"))),
                         fluidRow(column(6, #hidden(checkboxInput('PrintLabeltSNE2', 'Label clusters', value = FALSE)),
                                withSpinner(plotOutput("seurattSNEplot2", width = "100%", height = "500px", brush = brushOpts(id="tSNEbrush")), type = 6, color="#bf00ff", size = 1)), 
                         column(6, div(DT::dataTableOutput("brush1"), downloadButton("tSNEbrushdata", "Download data", class = "taskDFbutton"))))),
                tabPanel("tSNE highlight", status = "warning",
                         fluidRow(column(12, htmlOutput("tSNE_highlight"))),
                         fluidRow(column(6, DT::dataTableOutput("tSNEplot3")), 
                         column(6, #hidden(checkboxInput('PrintLabeltSNE3', 'Label clusters', value = FALSE)),
                                withSpinner(plotOutput("seurattSNEplot3", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1))))
                #downloadButton("seuratTab1", "Download plot")),
              )
                )
        }
      })
    #}
    #else if(input$seuratRUN > 0){
    #observeEvent(input$seuratRUN,{
      output[["seuratDiffExpr"]] <- renderUI({
        listSteps <- hideBoxes()
        if("Differentially expressed genes" %in% listSteps){
          box(title = "iS-CellR results: differentially expressed genes", status = "primary", solidHeader = FALSE,
              collapsed = TRUE, collapsible = TRUE, width=12,
              tabBox(
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "seuratDiffExpr", width=12, 
                tabPanel("Joy plot", status = "primary", 
                         fluidRow(column(6, #style='height:180px',
                         uiOutput("choose_JoyplotGenes"),
                            #textAreaInput("JoyplotGenes", paste0("Genes of interest e.g:", paste(sample(GeneNames,3), collapse = ',')), width = "500px"),
                         actionButton("updateJoyplot", "Plot it", icon("hand-o-right"), class = "taskRUNbutton", disabled = TRUE), 
                                     # style="color: #fff; background-color: #337ab7; border-color: #2e6da4", disabled = TRUE), 
                         actionButton("Joyplotfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE), 
                         downloadButton("downJoyplot", label = "Save plot", class = "taskDFbutton", disabled = TRUE)),
                         column(6, #div(style='height:10px;',
                              valueBoxOutput("GeneNotFound_Joy", width = NULL))
                         ), tags$hr(),
                         withSpinner(plotOutput("seuratJoyplot", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("JoyplotZOOM", "Joy plot", "Joyplotfull", withSpinner(plotOutput("seuratJoyplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))),
                tabPanel("Violin plot", status = "primary", 
                         fluidRow(column(6, #style='height:180px',
                            #textAreaInput("VlnplotGenes", paste0("Genes of interest e.g:", paste(sample(GeneNames,3), collapse = ',')), width = "500px"), 
                         uiOutput("choose_VlnplotGenes"),
                         actionButton("updateVlnplot", "Plot it", icon("hand-o-right"),                        
                                      class = "taskRUNbutton", disabled = TRUE), 
                         actionButton("Vlnplotfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE),
                         downloadButton("downVlnplot", label = "Save plot", class = "taskDFbutton", disabled = TRUE)),
                         column(6, #div(style='height:10px;'),
                                valueBoxOutput("GeneNotFound_Vln", width = NULL))
                         ), tags$hr(),
                         withSpinner(plotOutput("seuratVlnplot", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("VlnplotZOOM", "Violin plot", "Vlnplotfull", withSpinner(plotOutput("seuratVlnplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))),
                #downloadButton("seuratTab1", "Download plot")),
                tabPanel("Dot plot", status = "primary", 
                         fluidRow(column(6, #style='height:180px',
                            #textAreaInput("DotplotGenes", paste0("Genes of interest e.g:", paste(sample(GeneNames,3), collapse = ',')), width = "500px"), 
                         uiOutput("choose_DotplotGenes"),
                         actionButton("updateDotplot", "Plot it", icon("hand-o-right"),
                                      class = "taskRUNbutton", disabled = TRUE), 
                         actionButton("Dotplotfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE), 
                         downloadButton("downDotplot", label = "Save plot", class = "taskDFbutton", disabled = TRUE)),
                         column(6, #style='height:180px',
                                valueBoxOutput("GeneNotFound_Dot", width = NULL))
                         ), tags$hr(),
                         withSpinner(plotOutput("seuratDotplot", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("DotplotZOOM", "Dot plot", "Dotplotfull", withSpinner(plotOutput("seuratDotplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))),
                tabPanel("Scatter plot", status = "primary", 
                         fluidRow(column(6, #style='height:180px',
                            #textAreaInput("FeatureGenes", paste0("Genes of interest e.g:", paste(sample(GeneNames,3), collapse = ',')), width = "500px"), 
                         uiOutput("choose_FeatureGenes"),
                            useShinyjs(),
                         actionButton("updateFeatplot", "Plot it", icon("hand-o-right"),
                                      class = "taskRUNbutton", disabled = TRUE),  
                         actionButton("Featureplotfull", "Full screen", class = "taskDFbutton", disabled = TRUE),
                         downloadButton("downFeatureplot", label = "Save plot", class = "taskDFbutton", disabled = TRUE),
                         hidden(checkboxInput('PrintLabel_Featplot', 'Label clusters', value = FALSE))), br(),
                         column(6, #style='height:180px',
                                valueBoxOutput("GeneNotFound_Feat", width = NULL))
                         ), tags$hr(),
                         withSpinner(plotlyOutput("seuratFeatureplot", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("FeatureplotZOOM", "Scatter plot", "Featureplotfull", withSpinner(plotlyOutput("seuratFeatureplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))),
                tabPanel("Heatmap", status = "primary", 
                         fluidRow(column(6, #style='height:180px',
                            #textAreaInput("HeatmapGenes", paste0("Genes of interest e.g:", paste(sample(GeneNames,3), collapse = ',')), width = "500px"), 
                         uiOutput("choose_HeatmapGenes"),
                            useShinyjs(),
                         actionButton("updateHeatmap", "Plot it", icon("hand-o-right"),
                                      class = "taskRUNbutton", disabled = TRUE),
                         actionButton("Heatmapfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE),
                         downloadButton("downHeatmap", label = "Save plot", class = "taskDFbutton", disabled = TRUE)),
                         column(6, #style='height:180px',
                                valueBoxOutput("GeneNotFound_Heat", width = NULL))
                         ), tags$hr(), 
                         withSpinner(plotOutput("seuratHeatmap", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("HeatmapZOOM", "Heatmap", "Heatmapfull", withSpinner(plotOutput("seuratHeatmapZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1)))
              ))
        }
      })
    #}
    #else if(input$seuratRUN > 0){
  #observeEvent(input$seuratRUN,{
      output[["seuratMakers"]] <- renderUI({
        listSteps <- hideBoxes()
        if("Discriminating marker genes" %in% listSteps){
          box(title = "iS-CellR results: discriminating markers", color="purple", solidHeader = FALSE,
              collapsed = TRUE, collapsible = TRUE, width=12,
              tabBox(
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "seuratMakers", width=12,
                tabPanel("CoExpress",  color="purple", 
                         fluidRow(column(5, #style='height:180px',
                            #textAreaInput("CoExprGenes", paste0("Genes of interest (name two genes) e.g:", paste(sample(GeneNames,2), collapse = ',')), width = "250px"),
                         uiOutput("choose_CoExprGenes"),
                         useShinyjs(),
                         actionButton("updateCoExprplot", "Plot it", icon("hand-o-right"),
                                      class = "taskRUNbutton", disabled = TRUE),
                         actionButton("CoExprplotfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE),
                         downloadButton("downCoExprplot", label = "Save plot", class = "taskDFbutton", disabled = TRUE),
                         hidden(checkboxInput('PrintLabel', 'Label clusters', value = FALSE))), #br(), br(),
                         #hidden(valueBoxOutput("GeneNotFound_CoExpr", width = NULL))),
                         #column(1,
                         # htmlOutput("spacer")
                         # ),
                         column(6, #style='height:180px',
                           numericInput("ExprCutoff", "Expression threshold:",
                                      min = 0.1, max = 20, value = 0.5, step = 0.5, width = "30%"),
                         htmlOutput("header"), tableOutput("Gene12ExprTable"))), 
                         fluidRow(column(12,
                          hidden(valueBoxOutput("GeneNotFound_CoExpr", width = NULL)))), tags$hr(),
                         withSpinner(plotlyOutput("seuratCoExprplot", width = "100%", height = "600px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("CoExprplotZOOM", "CoExpress genes plot", "CoExprplotfull", withSpinner(plotlyOutput("seuratCoExprplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))),
                tabPanel("Inter-sample heterogeneity",  color="fuchsia",
                         checkboxInput('geneheader', 'Header', TRUE),
                         radioButtons('genesep', 'Separator', 
                                      c(Comma=',',
                                        Semicolon=';',
                                        Tab='\t'),
                                      '\t',inline=T),
                         fileInput('geneList', '',
                                   accept = c(
                                     'text/csv',
                                     'text/comma-separated-values',
                                     'text/tab-separated-values',
                                     'text/plain',
                                     '.csv',
                                     '.tsv',
                                     multiple = FALSE
                                   )
                         ),
                         useShinyjs(),
                         actionButton("updateInterHetplot", "Plot it", icon("hand-o-right"),
                                      class = "taskRUNbutton", disabled = TRUE), 
                         downloadButton("downInterHet", label = "Save plot", class = "taskDFbutton", disabled = TRUE), tags$hr(),
                         #div(plotOutput("seuratInterHetplot") %>% withSpinner(type = 6, color="#0dc5c1", size = 1)),
                         withSpinner(plotOutput("seuratInterHetplot", width = "100%", height = "500px"), type = 6, color="#bf00ff", size = 1)
                         #plotOutput("seuratInterHetplot"),
                         ),
                tabPanel("Intra-sample heterogeneity",  color="fuchsia",
                         fluidRow(column(6,
                         conditionalPanel(condition = "output.seuratInterHetplot",
                                          textAreaInput("IntraHetGenes", "Sample of interest (from Inter-sample heterogeneity plot)", width = "300px")),
                         useShinyjs(),
                         actionButton("updateIntraHetplot", "Plot it", icon("hand-o-right"),
                                      class = "taskRUNbutton", disabled = FALSE), 
                         actionButton("IntraHetplotfull", "Full screen", icon("desktop"), class = "taskDFbutton", disabled = TRUE),
                         downloadButton("downIntraHetplot", label = "Save plot", class = "taskDFbutton", disabled = TRUE)),       
                         column(6, #style='height:180px',
                                htmlOutput("header_celltypes"), tags$style(type='text/css', '#CelltypeNotFound {background-color: white; color: red; width = "auto";}'),
                                textOutput("CelltypeNotFound"))), tags$hr(),
                         withSpinner(plotlyOutput("seuratIntraHetplot", width = "100%", height = "600px"), type = 6, color="#bf00ff", size = 1),
                         tags$head(tags$style(".modal-lg{width: 95%;height: 95%;}")),
                         tags$head(tags$style(".modal{height:95%;color:purple;}")),
                         tags$head(tags$style(".modal-dialog{width:95%}")),
                         tags$head(tags$style(".modal-body{min-height:95%}")),
                         tags$head(tags$style(".modal-header {background-color: #0099ff; color: white; font-weight: bold !important; font-size: 38px !important; }")),
                        bsModal("IntraHetplotZOOM", "Intra-sample heterogeneity plot", "IntraHetplotfull", withSpinner(plotlyOutput("seuratIntraHetplotZoom", width = "100%", height = "700px"), type = 6, color="#bf00ff", size = 1))
                #downloadButton("seuratTab1", "Download plot")),
              )))
        }
      })
      #}
 })
})
      #file.read(out.markdown)
}

  #session$onSessionEnded(stopApp)
  

