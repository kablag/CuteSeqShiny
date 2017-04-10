library(shiny)
library(genbankr)
library(dplyr)
library(pipeR)
library(data.table)
library(RColorBrewer)
# library(Biostrings)

source("generic.R")

options(digits.secs=2)
writeLogs <- TRUE

createLogEntry <- function(text) {
  if (writeLogs)
    cat(sprintf("%s: %s\n", Sys.time(), text))
}

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })

  values <- reactiveValues(
    features = list(
      GenBank = NULL,
      Plain = NULL
    ),
    sequence = list(
      GenBank = NULL,
      Plain = NULL
    ),
    paletteUpdatesCounter = 0
  )

  gbFile <- reactive({
    req(input$gbFile)
    createLogEntry("Loading GenBank file")
    readGenBank(input$gbFile$datapath)
  })

  output$seqNameSelectUI <- renderUI({
    req(gbFile())
    createLogEntry("Generating seqNameSelectUI")
    selectInput("seqNameSelect",
                "Select Sequence",
                names(gbFile()@sequence))
  })

  observe({
    req(input$seqNameSelect)
    createLogEntry("Selecting sequence from GenBank file")
    isolate({
      values$features[["GenBank"]] <- gbFile()@other_features %>>%
        as.data.frame() %>>%
        filter(seqnames == input$seqNameSelect) %>>%
        as.data.table() %>>%
        getGeneiousTypes()
      values$sequence[["GenBank"]] <- gbFile()@sequence[[input$seqNameSelect]]
    })
    # ,
    #          (start < first(flatMap$seqI) &&
    #             end < first(flatMap$seqI)) ||
    #            (start < first(flatMap$seqI) &&
    #               end < first(flatMap$seqI)))
  })

  observeEvent(input$plainInputExample, {
    createLogEntry("Loading plain input example")
    updateTextAreaInput(
      session, "plainSequenceInput",
      value = paste0("RATTTTTGCCACATTRGAAGGAAAATTATTTCCACCAAGATTTCCCT",
                     "ACAGCCAAACGATCTACCAACTACAAAAATGGAAAAAATAATTTAGGACATGTA",
                     "AAGTTCAAATGTTTTGCCDTCCCACGTTTCNGTTTCAAGAAGCTATTYCGAGATAA",
                     "ATCGCTCCGTGGTCACAGGACTTAGAAAGGTGGAGGTAAACACACACAAGCATT",
                     "ATAAGATAAGAAGTAACAGATGAATTAGTTGAAAGGGACTGATTTCGGGGGAAN"))
    updateTextAreaInput(
      session, "plainFeaturesInput",
      value = paste0(
        "1; ABO-f; TACCAACTACAAAAATGGAA\n",
        "2; ABO-wt; (FAM)-TCCCACGTTTCGGTTTC-(BHQ1); FAM\n",
        "3; ABO-m; (hex)-tttctgtttcaagaagc(t-lna)att-(bhq1); HEX\n",
        "4; ABO-r; AGTCCTGTGACCACGGAG; ROX\n",
        "5; ABO-r2; TGCTTGTGTGTGTTTACCGCCA; primer; 1"))
  })

  observe({
    req(input$plainSequenceInput)
    createLogEntry("Parsing plain sequence")
    isolate({
      # seqs <- str_match_all(input$plainSequenceInput,
      #                       ">(.*)\\n([^>]*)")[[1]]
      dnaseq <- input$plainSequenceInput %>>%
        clearSeq() %>>%
        Biostrings::DNAString()
      values$sequence[["Plain"]] <- dnaseq
    })
    # str_match_all(tstr, ">(.*)\\n([^>]*)")
  })

  observe({
    req(input$plainFeaturesInput, values$sequence[["Plain"]] )
    createLogEntry("Parsing plain features")
    isolate({
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
        setnames(dt, 1:4, c("seqnames", "seq", "type", "maxMismatch"))
        dt[type == "", type := sapply(seq, function(x) getDyeFromSeq(x))]
        dt[type == "", type := "primer"]
        dt[is.na(maxMismatch), maxMismatch := 0]
        dt <- dt[seq != "" & seqnames != "" & type != ""]
        dt <- rbindlist(list(dt, dt),
                        use.names = TRUE, fill = FALSE, idcol = NULL)
        dt[, strand := rep(c("+","-"), each = nrow(dt) / 2)]
        dt[, seq := sapply(seq, function(x) clearSeq(x))]
        dt2 <- data.table()
        for (i in 1:nrow(dt)) {
          pattern <- switch(dt[i, strand],
                            "+" = dt[i, seq],
                            "-" = Biostrings::reverseComplement(Biostrings::DNAString(dt[i, seq])))
          res <-
            Biostrings::matchPattern(
              pattern,
              values$sequence[["Plain"]],
              max.mismatch = dt[i, maxMismatch],
              with.indels = input$allowInDels,
              fixed = FALSE)

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
            }
          }
          else {
            dt2 <- rbindlist(
              list(dt2,
                   data.table(seqnames = dt[i, seqnames],
                              seq = dt[i, seq],
                              type = dt[i, type],
                              strand = dt[i, strand],
                              maxMismatch = dt[i, maxMismatch],
                              start = -1L,
                              end = -1L,
                              mismatches = list(integer()))))

          }
        }
        values$features[["Plain"]] <- dt2
      }
    })
  })

  observe({
    input$markAmbiguity
    input$considerStrand
    req(input$inputTypeTabs)
    createLogEntry("Select working seq and features by input tab")
    values$workingFeatures <- values$features[[input$inputTypeTabs]]
    values$workingSequence <- values$sequence[[input$inputTypeTabs]]
    isolate({
      values$invalidatePalette <- TRUE
    })

  })

  output$colorByUI <- renderUI({
    req(values$workingFeatures)
    createLogEntry("Creating colorBy UI")
    selectInput("colorBy",
                "Color by",
                colnames(values$workingFeatures),
                if ("type" %in% colnames(values$workingFeatures))
                  "type")
  })

  output$labelByUI <- renderUI({
    req(values$workingFeatures)
    createLogEntry("Creating labelBy UI")
    selectInput("labelBy",
                "Label by",
                colnames(values$workingFeatures),
                if ("label" %in% colnames(values$workingFeatures))
                  "label")
  })
  output$limitSeqSliderUI <- renderUI({
    req(values$workingSequence)
    createLogEntry("Creating limitSeqSlider UI")
    sliderInput("limitSeqSlider",
                "Limit Sequence",
                min = 1 + values$workingSequence@offset,
                max = values$workingSequence@length + values$workingSequence@offset,
                value = c(1 + values$workingSequence@offset,
                          values$workingSequence@length + values$workingSequence@offset),
                step = 1)
  })

  output$featuresTbl <- DT::renderDataTable({
    req(values$workingFeatures)
    createLogEntry("Updating Features table")
    dt <- {
      if (input$inputTypeTabs == "GenBank") {
        values$workingFeatures
      } else {
        dt <- {
          if (input$showUnmatched)
            values$workingFeatures
          else values$workingFeatures[start != -1]
        }
        dt[,
           mismatchesN := {
             if (length(mismatches[[1]]))
               length(mismatches)
             else
               integer(0)
           },
           by = seqnames]
        dt[, -"mismatches", with = FALSE]
      }
    }
    DT::datatable(dt,
                  selection = "multiple"
    )
  })

  flatMap <- reactive({
    # input$processCuteSeq
    # isolate({
    req(input$labelBy,
        values$workingSequence,
        values$workingFeatures[[input$labelBy]])
    createLogEntry("Generating flat map")
    ft <- {
      if (input$inputTypeTabs == "GenBank") {
        values$workingFeatures
      } else {
        values$workingFeatures[start != -1]
      }
    }
    if (!nrow(ft))
      return(NULL)
    if (length(input$featuresTbl_rows_selected)) {
      ft <- ft[input$featuresTbl_rows_selected, ]
    }

    genFlatMap(
      values$workingSequence,
      ft,
      gbSequenceStart = 1,
      colorBy = input$colorBy,
      labelBy = input$labelBy,
      considerStrand = input$considerStrand,
      markAmbiguity = input$markAmbiguity
    )
  })


  autoPalette <- reactive({
    flatMap()
    req(input$colorBy,
        input$paletteType,
        input$inputTypeTabs, values$workingFeatures)
    createLogEntry("Generating auto palette")
    isolate({
      generatePalette(
        {
          if (input$inputTypeTabs == "GenBank") {
            values$workingFeatures
          } else {
            values$workingFeatures[start != -1]
          }
        },
        input$colorBy,
        input$considerStrand,
        # input$mismatchColor
        "#7FFF00",
        {
          if (input$paletteType == "standard") {
            if (input$lockPalette)
              values$workingPalette
            else
              fread("www/palettes/standard.txt")[
                , idParam := genParamID(param)]
          } else {
            if (input$lockPalette)
              values$workingPalette
            else
              NULL
          }
        }
      )
    })
  })

  output$changePaletteUI <- renderUI({
    req(autoPalette())
    createLogEntry("Creating palette inputs")
    isolate({
      # ids <- str_match(names(input), "Color_(.*)") %>>%
      #   na.omit()
      # sapply(ids[, 1], function(x) rminput[[x]] <- NULL)
      ui <- apply(autoPalette(), 1,
                  function(el){
                    colourpicker::colourInput(sprintf("Color_%s", el["idParam"]),
                                              sprintf("%s", el["param"]),
                                              el["color"])})
      values$invalidatePalette <- FALSE
      ui
    })
  })

  # workingPalette <- reactive({
  observe({
    ids <- str_match(names(input), "Color_(.*)") %>>%
      na.omit()
    req(!values$invalidatePalette,
         nrow(ids) != 0)
    sapply(ids[, 1], function(x) input[[x]])
    isolate({
      # print(values$paletteUpdatesCounter)
      values$paletteUpdatesCounter <- values$paletteUpdatesCounter - 1
      req(all(autoPalette()[,idParam] %in% ids[,2]),
          values$paletteUpdatesCounter <= 0)
      ids <- data.table(ids)
      # req(!values$invalidatePalette)
      createLogEntry("Updating working palette")

      setnames(ids, c("V1", "V2"), c("uiId", "idParam"))
      # sapply(ids[!(idParam %in% autoPalette()[, idParam]),
      #            uiId],
      #        function(id) {
      #          # print(id)
      #          # removeUI(
      #          #   selector = sprintf("div:has(#%s)", id),
      #          #   immediate = FALSE)
      #          # rlist::list.remove(input, id)
      #          # input[[id]] <- NULL
      #        })
      ids <- ids[,#idParam %in% autoPalette()[, idParam],
                 color := sapply(ids[, uiId], function(x) input[[x]])]
      palette <- sapply(ids[, uiId], function(x) input[[x]])
      # values$workingPalette <-
      #   autoPalette()[idParam %in% ids[, idParam],
      #                 .(param, idParam, color = ids[, color])]
      values$workingPalette <- merge(autoPalette(), ids,
            all.x = TRUE, by = c("idParam"))[
              !is.na(color.y), color.x := color.y][
                ,.(param, idParam, color = color.x)]
      # print(values$workingPalette)

    })
  })


  cuteSeqResult <- reactive({
    req(nrow(values$workingPalette) != 0)#, values$workingSequence)#, values$redraw)
    createLogEntry("Generating cuteSeqResult")
    input$includeLegend
    input$linesWidth
    input$spacingEveryNth
    input$featuresTbl_rows_selected
    isolate({
      values$paletteUpdatesCounter <- 0L
      HTML(
        paste0(
          "<div style='font-family:monospace;overflow:hidden;",
          "word-break:break-all;white-space:normal;'>",
          cuteSeq(
            flatMap = flatMap(),
            seqPalette = values$workingPalette,
            # seqPalette = autoPalette(),
            includeLegend = input$includeLegend,
            linesWidth = input$linesWidth,
            spacingEveryNth = as.integer(input$spacingEveryNth)),
          "</div>"
        ))
    })
  })

  output$cuteSeqHtml <- renderUI({
    # req(input$calcCuteSeq)
    # isolate({
      req(cuteSeqResult())
      createLogEntry("Generating cuteSeqHtml UI")
      list(
        cuteSeqResult()
      )
    # })
  })


  output$savePalette <- downloadHandler(
    filename = "palette.txt",
    content = function(file) {
      write.csv(values$workingPalette[, .(param, color)],
                file, row.names = FALSE)
    },
    contentType = "text/csv"
  )

  observe({
    req(input$loadPalette$datapath)
    createLogEntry("Loading Palette")
    isolate({
      values$paletteUpdatesCounter <- 0L
      tbl <- fread(input$loadPalette$datapath)[
        ,
        {
          uiElementName <- sprintf("Color_%s", genParamID(param))
          if (!is.null(input[[uiElementName]]) &&
              input[[uiElementName]] != color) {
            colourpicker::updateColourInput(
              session,
              uiElementName,
              value = color)
            values$paletteUpdatesCounter <- values$paletteUpdatesCounter + 1L
            # anyUpdated <- TRUE
            # invalidateLater(1000)
          }
        }, by = param]
      # values$paletteLoaded <- runif(1)
    })
  })
})
