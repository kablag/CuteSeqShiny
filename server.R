library(shiny)
library(genbankr)
library(dplyr)
library(pipeR)
library(data.table)
library(RColorBrewer)
# library(Biostrings)

source("generic.R")

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
    )
  )

  gbFile <- reactive({
    req(input$gbFile)
    cat("Loading GenBank file\n")
    readGenBank(input$gbFile$datapath)
  })

  output$seqNameSelectUI <- renderUI({
    req(gbFile())
    cat("Generating seqNameSelectUI\n")
    selectInput("seqNameSelect",
                "Select Sequence",
                names(gbFile()@sequence))
  })

  observe({
    req(input$seqNameSelect)
    cat("Selecting sequence from GenBank file\n")
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
    cat("Loading plain input example\n")
    updateTextAreaInput(
      session, "plainSequenceInput",
      value = paste0("ATTTTTGCCACATTGAAGGAAAATTATTTCCACCAAGATTTCCCT",
                     "ACAGCCAAACGATCTACCAACTACAAAAATGGAAAAAATAATTTAGGACATGTA",
                     "AAGTTCAAATGTTTTGCCTCCCACGTTTCNGTTTCAAGAAGCTATTCGAGATAA",
                     "ATCGCTCCGTGGTCACAGGACTTAGAAAGGTGGAGGTAAACACACACAAGCATT",
                     "ATAAGATAAGAAGTAACAGATGAATTAGTTGAAAGGGACTGATTTCGGGGGAA"))
    updateTextAreaInput(
      session, "plainFeaturesInput",
      value = paste0(
        "1; ABO-f; TACCAACTACAAAAATGGAA; primer\n",
        "2; ABO-wt; (FAM)-TCCCACGTTTCGGTTTC-(BHQ1); probe_wt\n",
        "3; ABO-m; (hex)-tttctgtttcaagaagc(t-lna)att-(bhq1); probe_m\n",
        "4; ABO-r; AGTCCTGTGACCACGGAG; primer\n",
        "5; ABO-r2; TGCTTGTGTGTGTTTACCGCCA; primer-r2; 1"))
  })

  observe({
    req(input$plainSequenceInput)
    cat("Parsing plain sequence\n")
    isolate({
      dnaseq <- input$plainSequenceInput %>>%
        clearSeq() %>>%
        Biostrings::DNAString()
      values$sequence[["Plain"]] <- dnaseq
    })
  })

  observe({
    req(input$plainFeaturesInput, values$sequence)
    cat("Parsing plain features\n")
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
              # "seqnames", "seq", "type", "maxMismatch"
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
    req(input$inputTypeTabs)
    cat("Select working seq and features by input tab\n")
    values$workingFeatures <- values$features[[input$inputTypeTabs]]
    values$workingSequence <- values$sequence[[input$inputTypeTabs]]
    isolate({
      values$invalidatePalette <- TRUE
    })

  })

  output$colorByUI <- renderUI({
    req(values$workingFeatures)
    cat("Creating colorBy UI\n")
    selectInput("colorBy",
                "Color by",
                colnames(values$workingFeatures),
                if ("type" %in% colnames(values$workingFeatures))
                  "type")
  })

  output$labelByUI <- renderUI({
    req(values$workingFeatures)
    cat("Creating labelBy UI\n")
    selectInput("labelBy",
                "Label by",
                colnames(values$workingFeatures),
                if ("label" %in% colnames(values$workingFeatures))
                  "label")
  })
  output$limitSeqSliderUI <- renderUI({
    req(values$workingSequence)
    cat("Creating limitSeqSlider UI\n")
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
    cat("Updating Features table\n")
    dt <- {
      if (input$inputTypeTabs == "GenBank") {
        values$workingFeatures
      } else {
        if (input$showUnmatched)
          values$workingFeatures
        else values$workingFeatures[start != -1]
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
    cat("Generating flat map\n")
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
      considerStrand = input$considerStrand
    )
  })


  autoPalette <- reactive({
    req(input$colorBy, input$considerStrand,
        input$inputTypeTabs, values$workingFeatures)
    cat("Generating auto palette\n")
    isolate({
      palette <- generatePalette(
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
          if (input$lockPalette)
            values$workingPalette
          else
            NULL
        }
      )
      # removeUI(
      #   selector = "div:has(> #color_uis)",
      #   immediate = TRUE
      # )
      palette
    })
  })
  # observe({
  #   req(input$colorBy,
  #       values$workingFeatures,
  #       values$workingFeatures[[input$colorBy]],
  #       !input$lockPalette)
  #   cat("Generating auto palette\n")
  #   values$autoPalette <-
  #       generatePalette(
  #         {
  #           if (input$inputTypeTabs == "GenBank") {
  #             values$workingFeatures
  #           } else {
  #             values$workingFeatures[start != -1]
  #           }
  #         },
  #         input$colorBy,
  #         input$considerStrand,
  #         # input$mismatchColor
  #         "#7FFF00"
  #       )
  #
  # })

  output$changePaletteUI <- renderUI({
    req(autoPalette())
    cat("Creating palette inputs\n")
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
    values$updateWorkingPalette
    # req(input$Color_mismatchColor)

    # isolate({
    # print(values$invalidatePalette)
    # print("==")
    # print(autoPalette())
    # print("==")
    # print(ids)
    # print("================")
    # req(ids)
    # input[[ids[1, 1]]]
    # isolate({
    req(autoPalette(), !values$invalidatePalette)
    ids <- data.table(ids)
    # req(!values$invalidatePalette)
    cat("Updating working palette\n")

    setnames(ids, c("V1", "V2"), c("uiId", "idParam"))
    sapply(ids[!(idParam %in% autoPalette()[, idParam]),
               uiId],
           function(id) removeUI(
             selector = sprintf("div:has(> %s)", id),
             immediate = FALSE
           ))
    ids <- ids[idParam %in% autoPalette()[, idParam]]
    palette <- sapply(ids[, uiId], function(x) input[[x]])
    values$workingPalette <-
      autoPalette()[idParam %in% ids[, idParam], .(param, idParam, color = palette)]
    # })
  })


  cuteSeqResult <- reactive({
    req(flatMap(), nrow(values$workingPalette) != 0)#, values$redraw)
    print(values$workingPalette)
    # req(autoPalette())
    cat("Generating cuteSeqResult\n")
    isolate({
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
    req(cuteSeqResult())
    cat("Generating cuteSeqHtml UI\n")
    list(
      cuteSeqResult()
    )
  })


  output$savePalette <- downloadHandler(
    filename = "palette.txt",
    content = function(file) {
      write.csv(values$workingPalette[, .(param, color)] %>>% (? .),
                file, row.names = FALSE)
    },
    contentType = "text/csv"
  )

  observe({
    req(input$loadPalette)
    cat("Loading Palette\n")
    isolate({
      tbl <- fread(input$loadPalette$datapath)[
        ,
        {
          uiElementName <- sprintf("Color_%s", genParamID(param))
          if (!is.null(input[[uiElementName]])) {
            print(paste(uiElementName, "ok"))
            colourpicker::updateColourInput(
              session,
              uiElementName,
              value = color)
            values$updateWorkingPalette <- runif(1)
          }
        }, by = param]
    })
  })
})
