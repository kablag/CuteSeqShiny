library(shiny)
# library(rhandsontable)

shinyUI(fluidPage(
  tags$head(
    includeScript("www/js/allowTab.js"),
    includeScript("www/js/google-analytics.js")
  ),
  titlePanel("CuteSeq"),
  wellPanel(
    tags$h3("Input"),
    tabsetPanel(
      id = "inputTypeTabs",
      tabPanel(
        "Plain",
        fluidRow(
          column(
            6,
            HTML(paste0("<div class='form-group shiny-input-container' >",
                        "<label for='plainSequenceInput'>Sequence</label>",
                        "<textarea id='plainSequenceInput' class='form-control' style='height:200px;'>",
                        "",
                        "</textarea>",
                        "</div>")),
            fluidRow(
              column(2,
                     actionButton("plainInputExample", "Example")
              )
            )
          ),
          # uiOutput("limitSeqSliderUI"),
          column(6,
                 HTML(paste0("<div class='form-group shiny-input-container'>",
                             "<label for='plainFeaturesInput'>Features ([ID]; Name; Sequence; [Type]; [Max Mismatch])</label>",
                             "<textarea id='plainFeaturesInput' class='form-control' style='height:200px;'",
                             " onkeydown='insertTab(this, event);'></textarea>",
                             "</div>")),
                 fluidRow(
                   column(2,
                          selectInput("plainFeaturesInputSep",
                                      "Separator",
                                      c("Auto" = "auto",
                                        ";" = ";",
                                        "Tab" = "\t",
                                        "Space" = " "))
                   ),
                   column(3,
                          checkboxInput("allowInDels",
                                        HTML(
                                          paste(
                                            "Allow Insertions/Deletions",
                                            "(max number of indels = Max Mismatch option)",
                                            sep = "<br>")),
                                        FALSE)
                   )
                 )
          )

        )),
      tabPanel("GenBank",
               fileInput("gbFile", "Upload GenBank File")
      )
    )),
  uiOutput("seqNameSelectUI"),
  wellPanel(
    tags$h3("Features"),
    fluidRow(
      column(3,
             uiOutput("colorByUI")),
      column(3,
             uiOutput("labelByUI")),
      column(3,
             checkboxInput("considerStrand", "Consider Strand", TRUE)),
      column(3,
             checkboxInput("showUnmatched", "Show Unmatched", FALSE))
    ),
    DT::dataTableOutput("featuresTbl")),
  wellPanel(
    tags$h3("Result"),
    # fluidRow(
    #   column(2,
    #          numericInput("linesWidth", "Lines Width", 0, min = 0, step = 1)),
    #   column(2,
    #          selectInput("spacingEveryNth", "Spacing Every Nth",
    #                      list("None" = "0",
    #                           "3" = "3",
    #                           "5" = "5",
    #                           "10" = "10"),
    #                      selected = "10")),
    #   column(2,
    #          checkboxInput("includeLegend", "Include Legend", TRUE)),
    #   column(2,
    #          colourpicker::colourInput("mismatchColor",
    #                                    "Mismatch Color",
    #                                    "chartreuse",
    #                                    palette = "limited"))
    # ),
    fluidRow(
      column(2,

                      numericInput("linesWidth", "Lines Width", 0, min = 0, step = 1),

                      selectInput("spacingEveryNth", "Spacing Every Nth",
                                  list("None" = "0",
                                       "3" = "3",
                                       "5" = "5",
                                       "10" = "10"),
                                  selected = "10"),

                      checkboxInput("includeLegend", "Include Legend", TRUE),

                      colourpicker::colourInput("mismatchColor",
                                                "Mismatch Color",
                                                "chartreuse",
                                                palette = "limited")
             ),
      column(6,
    htmlOutput("cuteSeqHtml")))
  )
))
