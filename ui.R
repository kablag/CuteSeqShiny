library(shiny)
# library(rhandsontable)

shinyUI(fluidPage(
  tags$head(includeScript("www/js/allowTab.js")),

  titlePanel("CuteSeq"),
  tabsetPanel(id = "inputTypeTabs",
    tabPanel("Plain",
             # textAreaInput("plainSequenceInput",
             #               "Sequence",
             #               "TAAGAGGTCCTTCACCAGCCTCCTCTCCCGGCATTATCCCATCTACCCCTCCACATTCAAGTTTTTGGAAAGATTCTACACTCCCAGTCTCTACTTCCTCACTTCTTCCTTGCTGCCCACGCCATAAACTAGCTGCTGCCTCCAGCATTGCCCTGACACCTAGTGGCTGGTGTCACCAAGACGCTAGACCCAATGGTTATTTATTTATTTATTTACTTATTTTGAGACGGAGTCTCACTCTGTCGCCCAGGCTGGAGTGCAGCGGTGCCATCTCGGCTCACTGCAACTTCCGCCTCCAGGGTTCAAGTGGTTCTCGTGCCTCAGCCTCCCAAGTAGTTTGGACTACAGGTGCCTGCCACCATGTCTGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTTGTCTTAAACTCCTGACCTCAAGTGATCCACCCACCTCGGCCTCCCAAAATGCTAGGATTATAGGCGTGAGCCACCGCACCNGGCCAATGGTTGTTTTTCAGGTCTTCTCTTGCTTGACTTCCCAGAGGGATCCCTTACTGTTGCACCTACCCTTCTGGGAACTCTCTTCCTCTGGCGTCTGTGATATTTCCCTCTCCTGCTGGCTCCTCCCTCTCCAGATGCTGTTTCTCACATCTACTCTCTTCTAGAGAGTGTGGTAGACAGAATAATGGTCACCAAAGATGTCCCTGCATGAATCCCTGGAACTTGTGAATATGATAGGTTAAATGGCCAAAAGGGAATTAAGGTTGCAGATGGAATTAAGCTGACCAATCTCCTGATTTTATTTTATTTTATTTTGTTTTTGAGGTGGAGTTTCGCTCTTGTTGCCCAACTGGAGTGCAATGGTGTGATCTCGGCTCACTGCAACCTCCGCCTGCCAGGTTCGAGAGATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGATTACAGGCACCCGCCATCATGCCTGGCTAATTTTTTAAATTTTTAGTAGAGACAGGG"),
             wellPanel(
               HTML(paste0("<div class='form-group shiny-input-container' >",
                           "<label for='plainSequenceInput'>Sequence</label>",
                           "<textarea id='plainSequenceInput' class='form-control' style='width:800px;height:220px;''>",
                           "ATTTTTGCCACATTGAAGGAAAATTATTTCCACCAAGATTTCCCTACAGCCAAACGATCTACCAACTACAAAAATGGAAAAAATAATTTAGGACATGTAAAGTTCAAATGTTTTGCCTCCCACGTTTCNGTTTCAAGAAGCTATTCGAGATAAATCGCTCCGTGGTCACAGGACTTAGAAAGGTGGAGGTAAACACACACAAGCATTATAAGATAAGAAGTAACAGATGAATTAGTTGAAAGGGACTGATTTCGGGGGAA",
                           "</textarea>",
                           "</div>")),
               # uiOutput("limitSeqSliderUI"),
               HTML(paste0("<div class='form-group shiny-input-container' style='width:800px'>",
                           "<label for='plainFeaturesInput'>Features ([ID]; Name; Sequence; [Type]; [Max Mismatch])</label>",
                           "<textarea id='plainFeaturesInput' class='form-control' style='height:200px;'",
                           " onkeydown='insertTab(this, event);'>",
                           "1; ABO-f; TACCAACTACAAAAATGGAA; primer
2; ABO-wt; TCCCACGTTTCGGTTTC; probe
3; ABO-m; TTTCTGTTTCAAGAAGCTATT; probe
4; ABO-r; AGTCCTGTGACCACGGAG; primer
5; ABO-r2; TGCTTGTGTGTGTTTACCGCCA; primer-r2; 1
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
                                       "chartreuse",
                                       palette = "limited"))
    ),
    htmlOutput("cuteSeqHtml"))
))
