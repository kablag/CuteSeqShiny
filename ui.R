
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
# library(rhandsontable)

shinyUI(fluidPage(
  tags$head(includeScript("www/js/allowTab.js")),

  titlePanel("CuteSeq"),
  tabsetPanel(
    tabPanel("Plain",
             # textAreaInput("plainSequenceInput",
             #               "Sequence",
             #               "TAAGAGGTCCTTCACCAGCCTCCTCTCCCGGCATTATCCCATCTACCCCTCCACATTCAAGTTTTTGGAAAGATTCTACACTCCCAGTCTCTACTTCCTCACTTCTTCCTTGCTGCCCACGCCATAAACTAGCTGCTGCCTCCAGCATTGCCCTGACACCTAGTGGCTGGTGTCACCAAGACGCTAGACCCAATGGTTATTTATTTATTTATTTACTTATTTTGAGACGGAGTCTCACTCTGTCGCCCAGGCTGGAGTGCAGCGGTGCCATCTCGGCTCACTGCAACTTCCGCCTCCAGGGTTCAAGTGGTTCTCGTGCCTCAGCCTCCCAAGTAGTTTGGACTACAGGTGCCTGCCACCATGTCTGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTTGTCTTAAACTCCTGACCTCAAGTGATCCACCCACCTCGGCCTCCCAAAATGCTAGGATTATAGGCGTGAGCCACCGCACCNGGCCAATGGTTGTTTTTCAGGTCTTCTCTTGCTTGACTTCCCAGAGGGATCCCTTACTGTTGCACCTACCCTTCTGGGAACTCTCTTCCTCTGGCGTCTGTGATATTTCCCTCTCCTGCTGGCTCCTCCCTCTCCAGATGCTGTTTCTCACATCTACTCTCTTCTAGAGAGTGTGGTAGACAGAATAATGGTCACCAAAGATGTCCCTGCATGAATCCCTGGAACTTGTGAATATGATAGGTTAAATGGCCAAAAGGGAATTAAGGTTGCAGATGGAATTAAGCTGACCAATCTCCTGATTTTATTTTATTTTATTTTGTTTTTGAGGTGGAGTTTCGCTCTTGTTGCCCAACTGGAGTGCAATGGTGTGATCTCGGCTCACTGCAACCTCCGCCTGCCAGGTTCGAGAGATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGATTACAGGCACCCGCCATCATGCCTGGCTAATTTTTTAAATTTTTAGTAGAGACAGGG"),
             wellPanel(
               HTML(paste0("<div class='form-group shiny-input-container' >",
                           "<label for='plainSequenceInput'>Sequence</label>",
                           "<textarea id='plainSequenceInput' class='form-control' style='width:800px;height:220px;''>",
                           "TAAGAGGTCCTTCACCAGCCTCCTCTCCCGGCATTATCCCATCTACCCCTCCACATTCAAGTTTTTGGAAAGATTCTACACTCCCAGTCTCTACTTCCTCACTTCTTCCTTGCTGCCCACGCCATAAACTAGCTGCTGCCTCCAGCATTGCCCTGACACCTAGTGGCTGGTGTCACCAAGACGCTAGACCCAATGGTTATTTATTTATTTATTTACTTATTTTGAGACGGAGTCTCACTCTGTCGCCCAGGCTGGAGTGCAGCGGTGCCATCTCGGCTCACTGCAACTTCCGCCTCCAGGGTTCAAGTGGTTCTCGTGCCTCAGCCTCCCAAGTAGTTTGGACTACAGGTGCCTGCCACCATGTCTGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTTGTCTTAAACTCCTGACCTCAAGTGATCCACCCACCTCGGCCTCCCAAAATGCTAGGATTATAGGCGTGAGCCACCGCACCNGGCCAATGGTTGTTTTTCAGGTCTTCTCTTGCTTGACTTCCCAGAGGGATCCCTTACTGTTGCACCTACCCTTCTGGGAACTCTCTTCCTCTGGCGTCTGTGATATTTCCCTCTCCTGCTGGCTCCTCCCTCTCCAGATGCTGTTTCTCACATCTACTCTCTTCTAGAGAGTGTGGTAGACAGAATAATGGTCACCAAAGATGTCCCTGCATGAATCCCTGGAACTTGTGAATATGATAGGTTAAATGGCCAAAAGGGAATTAAGGTTGCAGATGGAATTAAGCTGACCAATCTCCTGATTTTATTTTATTTTATTTTGTTTTTGAGGTGGAGTTTCGCTCTTGTTGCCCAACTGGAGTGCAATGGTGTGATCTCGGCTCACTGCAACCTCCGCCTGCCAGGTTCGAGAGATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGATTACAGGCACCCGCCATCATGCCTGGCTAATTTTTTAAATTTTTAGTAGAGACAGGG",
                           "</textarea>",
                           "</div>")),
               # uiOutput("limitSeqSliderUI"),
               HTML(paste0("<div class='form-group shiny-input-container' style='width:800px'>",
                           "<label for='plainFeaturesInput'>Features ([ID]; Name; Sequence; [Type]; [Max Mismatch])</label>",
                           "<textarea id='plainFeaturesInput' class='form-control' style='height:200px;'",
                           " onkeydown='insertTab(this, event);'>",
                           "1; VKORC1_2-f3; tcaccaagacgctagacc; primer
2; VKORC1_2-m; attggccaggtgcgg; probe
3; VKORC1_2-wt; ccattggccgggtgc; probe
4; VKORC1_2-r2; tctgggaagtcaagcaagaga; primer
5; VKORC1_2-f2; ggcctcccaaaatgctagga; primer
</textarea>",
                           "</div>")),
               fluidRow(
                 column(2,
                        selectInput("plainFeaturesInputSep",
                                    "Separator",
                                    c("Auto" = "auto",
                                      ";" = ";",
                                      "Tab" = "\t",
                                      "Space" = " "))
                 )
               )
             )),
    tabPanel("GenBank",
             fileInput("gbFile", "Upload GenBank File")
    )
  ),
  uiOutput("seqNameSelectUI"),
  wellPanel(
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
    fluidRow(
      column(2,
             numericInput("linesWidth", "Lines Width", 0, min = 0, step = 1)),
      column(2,
             checkboxInput("includeLegend", "Include Legend", TRUE)),
      column(2,
             colourpicker::colourInput("mismatchColor",
                                       "Mismatch Color",
                                       "indianred1",
                                       palette = "limited"))
    ),
    htmlOutput("cuteSeqHtml"))
))
