library(shiny)
library(genbankr)
library(dplyr)
library(pipeR)
library(data.table)
library(RColorBrewer)
# library(Biostrings)

source("generic.R")

clearSeq <- function(dnaseq) {
  gsub("[^ACGTMRWSYKVHDBN]", "", dnaseq, ignore.case = TRUE)
}

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })

  values <- reactiveValues()

  gbFile <- reactive({
    req(input$gbFile)
    readGenBank(input$gbFile$datapath)
  })

  output$seqNameSelectUI <- renderUI({
    req(gbFile())
    selectInput("seqNameSelect",
                "Select Sequence",
                names(gbFile()@sequence))
  })

  observe({
    req(input$seqNameSelect)
    isolate({
      values$features <- gbFile()@other_features %>>%
        as.data.frame() %>>%
        filter(seqnames == input$seqNameSelect) %>>%
        as.data.table() %>>%
        getGeneiousTypes()
      values$sequence <- gbFile()@sequence[[input$seqNameSelect]]
    })
    # ,
    #          (start < first(flatMap$seqI) &&
    #             end < first(flatMap$seqI)) ||
    #            (start < first(flatMap$seqI) &&
    #               end < first(flatMap$seqI)))
  })

  observe({
    req(input$plainSequenceInput)
    isolate({
      dnaseq <<- input$plainSequenceInput %>>%
        clearSeq() %>>%
        Biostrings::DNAString()
      values$sequence <- dnaseq
    })
  })

  observe({
    req(input$plainFeaturesInput, values$sequence)
    getSep <- function(sep) {
      switch(sep,
             "auto" = {
               if (grepl(";" , input$plainFeaturesInput))
                 return(";")
               if (grepl("\t" , input$plainFeaturesInput))
                 return("\t")
               " "
             },
             sep)
    }
    tblSeparator <- getSep(input$plainFeaturesInputSep)
    dt <- tryCatch(
      fread(input$plainFeaturesInput,
            sep = tblSeparator,
            header = FALSE,
            fill = TRUE),
      error = function(e) NULL)
    if (is.null(dt) || !nrow(dt))
      return(NULL)
    if (typeof(dt[[1]]) == "integer")
      dt <- dt[, -1]
    if (ncol(dt) >= 2) {
      switch(as.character(ncol(dt)),
             "2" = {
               dt[, type := "primer"]
               dt[, maxMismatch := 0]
             },
             "3" = {
               if (typeof(dt[[3]]) != "character")
                 return(NULL)
               dt[, maxMismatch := 0]
             },
             "4" = {
               if (typeof(dt[[4]]) != "integer") {
                 dt[[4]] <- NULL
                 dt[, maxMismatch := 0]
               }
             },
             {
               dt <- dt[, 1:4]
               if (typeof(dt[[4]]) != "integer") {
                 dt[[4]] <- NULL
                 dt[, maxMismatch := 0]
               }
             }
      )
      # if (ncol(dt) > 4)
      #   dt <- dt[, 1:4]
      setnames(dt, 1:4, c("seqnames", "seq", "type", "maxMismatch"))
      dt[type == "", type := "primer"]
      dt[is.na(maxMismatch), maxMismatch := 0]
      dt <- dt[seq != "" & seqnames != "" & type != ""]
      #
      # dt[, seq := DNAString(seq)]
      dt <- rbindlist(list(dt, dt),
                      use.names = TRUE, fill = FALSE, idcol = NULL)
      dt[, strand := rep(c("+","-"), each = nrow(dt) / 2)]
      # dt[, ID := .I]
      dt[, seq := sapply(seq, function(x) clearSeq(x))]
      dt2 <- data.table()
      # dt[, c("start", "end", "mismatchPositions") := {
      #   pattern <- switch(strand,
      #                     "+" = seq,
      #                     "-" = Biostrings::reverseComplement(Biostrings::DNAString(seq)))
      #   res <-
      #     Biostrings::matchPattern(
      #       pattern,
      #       values$sequence,
      #       max.mismatch = maxMismatch)
      #   mismatches <- Biostrings::mismatch(pattern, res)
      #
      #   # dt2 <- rbindlist(list(dt2, data.table(1, list(c(1,2)))))
      #   # print(dt2)
      #   if (length(res)) list(res[1]@ranges@start,
      #                         res[1]@ranges@start + res[1]@ranges@width,
      #                         mismatches[1])
      #   else list(-1L, -1L, list(c(0L, 0L)))},
      #   by = ID]
      for (i in 1:nrow(dt)) {
        pattern <- switch(dt[i, strand],
                          "+" = dt[i, seq],
                          "-" = Biostrings::reverseComplement(Biostrings::DNAString(dt[i, seq])))
        res <-
          Biostrings::matchPattern(
            pattern,
            values$sequence,
            max.mismatch = dt[i, maxMismatch])

        if (length(res)) {
          mismatches <- Biostrings::mismatch(pattern, res)
          sname <- {
            if (length(mismatches) > 1)
                          function(name, i) sprintf("%s~%i", name, i)
            else
              function(name, i) name
          }
          for (mmatchIndex in 1:length(mismatches)) {
            dt2 <- rbindlist(
              list(dt2,
                   data.table(seqnames = sname(dt[i, seqnames], mmatchIndex),
                              seq = dt[i, seq],
                              type = dt[i, type],
                              strand = dt[i, strand],
                              maxMismatch = dt[i, maxMismatch],
                              start = res[mmatchIndex]@ranges@start,
                              end = res[mmatchIndex]@ranges@start + res[mmatchIndex]@ranges@width,
                              mismatches = mismatches[mmatchIndex])))
            # "seqnames", "seq", "type", "maxMismatch"
          }
        }
        else {
          dt2 <- rbindlist(
            list(dt2,
                 data.table(seqnames = sname(dt[i, seqnames], mmatchIndex),
                            seq = dt[i, seq],
                            type = dt[i, type],
                            strand = dt[i, strand],
                            maxMismatch = dt[i, maxMismatch],
                            start = -1L,
                            end = -1L,
                            mismatches = list(integer()))))

        }
      }
      # print(dt2)
      isolate({
        values$features <- dt2
      })
    }

  })

  output$colorByUI <- renderUI({
    req(values$features)
    selectInput("colorBy",
                "Color by",
                colnames(values$features),
                if ("type" %in% colnames(values$features))
                  "type")
  })

  output$labelByUI <- renderUI({
    req(values$features)
    selectInput("labelBy",
                "Label by",
                colnames(values$features),
                if ("label" %in% colnames(values$features))
                  "label")
  })
  output$limitSeqSliderUI <- renderUI({
    req(values$sequence)
    sliderInput("limitSeqSlider",
                "Limit Sequence",
                min = 1 + values$sequence@offset,
                max = values$sequence@length + values$sequence@offset,
                value = c(1 + values$sequence@offset,
                          values$sequence@length + values$sequence@offset),
                step = 1)
  })

  output$featuresTbl <- DT::renderDataTable({
    req(values$features)
    dt <- {if (input$showUnmatched)
      values$features
      else values$features[start != -1]}
    DT::datatable(dt,
                  selection = "multiple"
    )
  })

  output$cuteSeqHtml <- renderUI({
    # input$processCuteSeq
    # isolate({
    req(values$sequence, values$features, input$colorBy, input$labelBy)

    ft <- values$features[start != -1]
    if (!nrow(ft))
      return(NULL)
    if (length(input$featuresTbl_rows_selected)) {
      ft <- ft[input$featuresTbl_rows_selected, ]
    }
    HTML(
      paste0(
        "<div style='font-family:monospace;overflow:hidden;width:800px;word-break:break-all;white-space:normal;'>",
        cuteSeq(values$sequence, ft,
                colorBy = input$colorBy,
                labelBy = input$labelBy,
                mismatchColor = input$mismatchColor,
                considerStrand = input$considerStrand,
                includeLegend = input$includeLegend,
                linesWidth = input$linesWidth),
        "</div>"
      ))
    # })
  })

})
