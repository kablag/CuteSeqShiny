library(shiny)
# library(rhandsontable)
cuteSeqVersion <- 1.0

shinyUI(fluidPage(
  tags$head(
    includeScript("www/js/allowTab.js"),
    includeScript("www/js/google-analytics.js")
  ),
  titlePanel(
    HTML(
      paste0(
        sprintf("CuteSeq %.01f&emsp;", cuteSeqVersion),
        tags$a(href = "www.evrogen.com",
               tags$img(alt = "Evrogen",
                        src = "http://evrogen.com/img/evrogen.png",
                        border = 0))
      ))),
  wellPanel(
    tags$h3("Input"),
    tabsetPanel(
      id = "inputTypeTabs",
      tabPanel(
        "Plain",
        div(style = "margin-top:12px;",
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
        )
        )),
      tabPanel("GenBank",
               div(style = "margin-top:12px;",
                   fluidRow(
                     column(2,
               fileInput("gbFile", "Upload GenBank File")))
               )
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
           wellPanel(
             numericInput("linesWidth", "Lines Width", 0, min = 0, step = 1),

             selectInput("spacingEveryNth", "Spacing Every Nth",
                         list("None" = "0",
                              "3" = "3",
                              "5" = "5",
                              "10" = "10"),
                         selected = "10"),

             checkboxInput("includeLegend", "Include Legend", TRUE)#,
             # colourpicker::colourInput("mismatchColor",
             #                           "Mismatch Color",
             #                           "chartreuse",
             #                           palette = "limited")
           )),
    column(2,
           wellPanel(
             fluidRow(
               column(6,
                      checkboxInput("lockPalette",
                                    "Lock Palette",
                                    FALSE)),
               column(6,
                      h5("Load Palette"))),
             fluidRow(
               column(6,
                      downloadButton("savePalette",
                                     "Save Palette")),
               column(6,
                      fileInput("loadPalette",
                                # "Load Palette"))
                                NULL))
             ),
             uiOutput("changePaletteUI")
             # DT::dataTableOutput("paletteTbl")
           )),
    column(6,
           htmlOutput("cuteSeqHtml"))
  ),
  div(
    style = "margin-bottom:5px;",
    HTML(paste0(
      "Konstantin Blagodatskikh (",
      tags$a(href = "mailto:k.blag@yandex.ru?Subject=CuteSeq",
             target = "_top",
             "k.blag@yandex.ru"),
      ")"))
  )
))

