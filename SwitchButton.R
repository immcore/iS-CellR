#' A function to change the Original checkbox of rshiny
#' into a nice true/false or on/off switch button
#' No javascript involved. Only CSS code.
#' 
#' To be used with CSS script 'button.css' stored in a 'www' folder in your Shiny app folder
#' 
#' @param inputId The input slot that will be used to access the value.
#' @param label Display label for the control, or NULL for no label.
#' @param value Initial value (TRUE or FALSE).
#' @param col Color set of the switch button. Choose between "GB" (Grey-Blue) and "RG" (Red-Green)
#' @param type Text type of the button. Choose between "TF" (TRUE - FALSE), "OO" (ON - OFF) or leave empty for no text.

switchButton <- function(inputId, label, value=FALSE, col = "GB", type="TF") {
  
  # color class
  if (col != "RG" & col != "GB") {
    stop("Please choose a color between \"RG\" (Red-Green) 
      and \"GB\" (Grey-Blue).")
  }
  if (!type %in% c("OO", "TF", "YN")){
   warning("No known text type (\"OO\", \"TF\" or \"YN\") have been specified, 
     button will be empty of text") 
  }
  if(col == "RG"){colclass <- "RedGreen"}
  if(col == "GB"){colclass <- "GreyBlue"}
  if(type == "OO"){colclass <- paste(colclass,"OnOff")}
  if(type == "TF"){colclass <- paste(colclass,"TrueFalse")}
  if(type == "YN"){colclass <- paste(colclass,"YesNo")}
  
  # No javascript button - total CSS3
  # As there is no javascript, the "checked" value implies to
  # duplicate code for giving the possibility to choose default value
  
  if(value){
    tagList(
      tags$div(class = "form-group shiny-input-container",
        tags$div(class = colclass,
          tags$label(label, class = "control-label"),
          tags$div(class = "onoffswitch",
            tags$input(type = "checkbox", name = "onoffswitch", class = "onoffswitch-checkbox",
              id = inputId, checked = ""
            ),
            tags$label(class = "onoffswitch-label", `for` = inputId,
              tags$span(class = "onoffswitch-inner"),
              tags$span(class = "onoffswitch-switch")
            )
          )
        )
      )
    )
  } else {
    tagList(
      tags$div(class = "form-group shiny-input-container",
        tags$div(class = colclass,
          tags$label(label, class = "control-label"),
          tags$div(class = "onoffswitch",
            tags$input(type = "checkbox", name = "onoffswitch", class = "onoffswitch-checkbox",
              id = inputId
            ),
            tags$label(class = "onoffswitch-label", `for` = inputId,
              tags$span(class = "onoffswitch-inner"),
              tags$span(class = "onoffswitch-switch")
            )
          )
        )
      )
    ) 
  }
}
