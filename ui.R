
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(rhandsontable)

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
             uiOutput("limitSeqSliderUI"),
             HTML(paste0("<div class='form-group shiny-input-container' style='width:800px'>",
                         "<label for='plainFeaturesInput'>Features ([ID]; Name; Sequence; Type; [Max Mismatch])</label>",
                         "<textarea id='plainFeaturesInput' class='form-control' style='height:200px;'",
                         " onkeydown='insertTab(this, event);'>",
                         "1; VKORC1_2-f3; tcaccaagacgctagacc; primer
2; VKORC1_2-m; attggccaggtgcgg; probe
3; VKORC1_2-wt; ccattggccgggtgc; probe
4; VKORC1_2-r2; tctgggaagtcaagcaagaga; primer
5; VKORC1_2-f2; ggcctcccaaaatgctagga; primer
</textarea>",
                         "</div>")),
             selectInput("plainFeaturesInputSep",
                         "Separator",
                         c(";" = ";",
                           "Whitespace" = " ",
                           "Tab" = "\t")))),
    tabPanel("GenBank",
             fileInput("gbFile", "Upload GenBank File")
    )
  ),
  uiOutput("seqNameSelectUI"),
  wellPanel(
  fluidRow(
    column(4,
           uiOutput("colorByUI")),
    column(4,
           uiOutput("labelByUI"))
  ),
  fluidRow(
    column(4,
           checkboxInput("considerStrand", "Consider Strand", TRUE)),
    column(4,
           checkboxInput("showUnmatched", "Show Unmatched", FALSE))),
  DT::dataTableOutput("featuresTbl")),
  wellPanel(
  fluidRow(
    column(4,
           numericInput("linesWidth", "Lines Width", 0, min = 0, step = 1)),
    column(4,
           checkboxInput("includeLegend", "Include Legend", TRUE))
  ),
  htmlOutput("cuteSeqHtml"))
))
