# addFeatureToLayers <- function(layers, featureId, from, to, compact = TRUE) {
#   to.layer <- NA
#   if (from < first(layers$seqI))
#     from <- first(layers$seqI)
#   if (to > last(layers$seqI))
#     to <- last(layers$seqI)
#   if (compact) {
#     for (layer.i in 3:ncol(layers)) {
#       if (all(is.na(layers[seqI >= from & seqI <= to][[layer.i]]))) {
#         to.layer <- layer.i - 2
#         break()
#       }
#     }
#   } else {
#     if (all(is.na(layers[seqI >= from & seqI <= to][[ncol(layers)]]))) {
#       to.layer <- ncol(layers) - 2
#     }
#   }
#   if (is.na(to.layer))
#     to.layer <- ncol(layers) - 1
#   layers[seqI >= from & seqI <= to,
#          c(sprintf("layer_%s", to.layer)) := featureId]
#   NULL
# }
library(stringr)

genParamID <- function(param) {
  str_replace_all(param,
                  c("-" = "B", "\\+" = "A", " " = "-"))
}

clearSeq <- function(dnaseq) {
  dnaseq <- str_replace_all(dnaseq, "\\(([a-zA-Z])-[^\\)]*\\)", "\\1")
  dnaseq <- str_replace_all(dnaseq, "-?\\([^\\)]*\\)-?", "")
  gsub("[^ACGTMRWSYKVHDBN]", "", dnaseq, ignore.case = TRUE)
}

getColors <- function(ncolors) {
  if (ncolors > 12)
    return(
      apply(
        col2rgb(
          sample(colors(distinct = TRUE), ncolors)
        ),
        2,
        function(color) rgb(color[1], color[2], color[3], maxColorValue = 255))
    )
  switch(as.character(ncolors),
         "0" = "",
         "1" = "#8DD3C7",
         "2" = c("#8DD3C7", "#BEBADA"),
         brewer.pal(ncolors, "Set3"))
}

addFlat <- function(flatMap, ID, start, end, mismatches = NULL) {
  if (end < data.table::first(flatMap$seqI) ||
      start > data.table::last(flatMap$seqI))
    return(NULL)
  if (start < data.table::first(flatMap$seqI))
    start <- data.table::first(flatMap$seqI)
  if (end > (data.table::last(flatMap$seqI) + 1))
    end <- data.table::last(flatMap$seqI)
  subMap <- flatMap[J(start:(end - 1))]
  hasMapI <- subMap[!is.na(map), seqI]
  flatMap[J(hasMapI), intersectionHere := TRUE]
  flatMap[J(start:(end - 1)), map := ID]
  if (length(mismatches))
    flatMap[J(unlist(mismatches) + start - 1), mismatchHere := TRUE]
  NULL
}

getGeneiousTypes <- function(featuresTable) {
  featuresTable[!is.na(note),
                type := sub("Geneious type: (.*)", "\\1", note)]
}

getDyeFromSeq <- function(dnaseq, dyes = c("FAM", "HEX", "VIC", "ROX", "Cy5")) {
  res <- str_extract(dnaseq,
                     regex(paste(dyes, collapse = "|"),
                           ignore_case = TRUE))
  if (is.na(res))
    ""
  else
    toupper(res)
}

generatePalette <- function(gbFeatures,
                            colorBy,
                            considerStrand,
                            mismatchColor,
                            currentPalette) {
  if (considerStrand) {
    gbFeatures[, c(colorBy) := paste(get(colorBy), strand)]
  }
  uniqueColorByParams <- unique(gbFeatures[, get(colorBy)])
  newPalette <- data.table(param = c(uniqueColorByParams, "Mismatch Color"),
                           idParam = c(genParamID(uniqueColorByParams),
                                       genParamID("Mismatch Color")),
                           color = c(getColors(length(uniqueColorByParams)),
                                     mismatchColor)) %>>%
    setkey(param)
  if (!is.null(currentPalette)) {
    merge(newPalette, currentPalette,
          all.x = TRUE, by = c("param"))[
            !is.na(color.y), color.x := color.y][
              ,.(param, idParam=idParam.x, color = color.x)]
  } else {
    newPalette
  }
}

genFlatMap <- function(gbSequence,
                       gbFeatures,
                       gbSequenceStart = 1,
                       colorBy,
                       labelBy,
                       considerStrand = TRUE,
                       markAmbiguity = TRUE) {
  features <- as.data.table(
    gbFeatures)[, ID := .I]
  if (considerStrand) {
    features[, c(colorBy) := paste(get(colorBy), strand)]
  }
  flatMap <- data.table(
    seqI = (gbSequenceStart + gbSequence@offset):
      (gbSequence@length + gbSequence@offset),
    gbSeq = strsplit(as.character(gbSequence), NULL)[[1]],
    map = as.integer(NA),
    intersectionHere = FALSE,
    mismatchHere = FALSE,
    ambiguityHere = FALSE
  )

  setkey(flatMap, seqI)
  if (is.null(features[["mismatches"]]))
    features[, mismatches := list(integer())]
  features[, addFlat(flatMap, ID, start, end, mismatches), by = ID]
  # features[, addFeatureToLayers(allLayers[[strand]], ID, start, end), by = ID]

  # flatMap[!is.na(map),
  #         # color := ifelse(is.na(map),
  #         #                   "",
  #                           seqPalette[map, color],
  #         # ),
  #         by = seqI]
  if (markAmbiguity) {
    flatMap[, ambiguityHere := !(gbSeq %in% c("A", "T", "G", "C"))]
  }
  flatMap[ambiguityHere == TRUE, gbSeq :=
            sapply(gbSeq, function(ambSeq)
              sprintf("[%s]",
                      paste(strsplit(Biostrings::IUPAC_CODE_MAP[ambSeq], NULL)[[1]],
                            collapse = "/")))
          ]
  flatMap[!is.na(map), label := features[map, get(labelBy)]]
  flatMap[!is.na(map), typeID := features[map, get(colorBy)]]
  flatMap[is.na(map), map := 0]
  flatMap[, dif := {
    map != data.table::shift(map, fill = FALSE) |
      intersectionHere != data.table::shift(intersectionHere, fill = FALSE) |
      mismatchHere != data.table::shift(mismatchHere, fill = FALSE) |
      ambiguityHere != data.table::shift(ambiguityHere, fill = FALSE)
  }]
  flatMap
}

cuteSeq <- function(flatMap,
                    seqPalette,
                    includeLegend = TRUE,
                    linesWidth = 60,
                    spacingEveryNth = 10) {
  setkey(seqPalette, param)
  flatMapCopy <- copy(flatMap)
  if (spacingEveryNth)
    flatMapCopy[seq(spacingEveryNth, nrow(flatMapCopy), spacingEveryNth),
                gbSeq := sprintf("<span style='letter-spacing:0.5em;'>%s</span>", gbSeq)]
  flatMapCopy[,
              coloredSeq :=
                ifelse(dif,
                       ifelse(map == 0,
                              ifelse(ambiguityHere,
                                     {
                                       sprintf("%s<span style='font-weight:bold;'>%s",
                                               ifelse(!is.na(data.table::shift(typeID, fill = FALSE)),
                                                      "</span>",
                                                      ""),
                                               gbSeq)
                                     },
                                     sprintf("</span>%s", gbSeq)),
                              sprintf("</span><span style='background-color: %s;%s%s' title='%s'>%s",
                                      ifelse(mismatchHere,
                                             seqPalette["Mismatch Color", color],
                                             seqPalette[typeID, color]),
                                      ifelse(intersectionHere,
                                             "text-decoration:underline;",
                                             ""),
                                      ifelse(ambiguityHere,
                                             "font-weight:bold;",
                                             ""),
                                      label,
                                      gbSeq)
                              #)
                       ),
                       gbSeq),
              by = seqI]
  if (linesWidth != 0) {
    flatMapCopy[seq(1, nrow(flatMapCopy), by = linesWidth),
                coloredSeq := sprintf("<br>%s", coloredSeq)]
  }

  legendTbl <- ""
  if (includeLegend) {
    legendTbl <-
      sprintf("%s%s%s%s",
              shiny::renderTable(seqPalette[!"Mismatch Color",
                                            Color := sprintf("<span style='background-color: %s;'>ATGC</span>", color)]
                                 [!"Mismatch Color", Features := paste(unique(flatMapCopy[typeID == param, label]), collapse = ", "),
                                   by = param]
                                 [!"Mismatch Color", Strand :=  str_extract(param,"[-\\+]$")]
                                 [!"Mismatch Color", Param :=  str_replace(param,"[-\\+]$", "")]
                                 # [,c("param", "strand") = list(str_match(param, "(.*) ?([+-])$?"))]
                                 [!"Mismatch Color", .(Color, Param, Strand, Features)],
                                 colnames = FALSE,
                                 sanitize.text.function = function(x) x)()
              ,
              ifelse(any(flatMapCopy[, mismatchHere] == TRUE),
                     sprintf("<span style='background-color: %s'>&emsp;&emsp;</span>&emsp;Mismatch",
                             seqPalette["Mismatch Color", color]),
                     ""),
              ifelse(any(flatMapCopy[, intersectionHere] == TRUE),
                     "<br><span style='text-decoration:underline;'>&emsp;&emsp;</span>&emsp;Intersection",
                     ""),
              ifelse(any(flatMapCopy[, ambiguityHere] == TRUE),
                     {
                       # ambNuc <- unique(flatMapCopy[ambiguityHere == TRUE, gbSeq])
                       # sprintf("<br><b>MRWS</b>&emsp;Ambiguity (%s)",
                       #         paste(paste(ambNuc, Biostrings::IUPAC_CODE_MAP[ambNuc], sep = " = "),
                       #               collapse = ", "))
                       "<br><b>[A/C]</b>&emsp;Ambiguity"
                     },
                     ""))
  }
  # styles <-
  #   paste0("<style>",
  #          "span.intersection{text-decoration:underline;}",
  #          sprintf("span.mismatch{background-color: %s;}", mismatchColor),
  #          paste0(seqPalette[, sprintf("span.param-%s{background-color:%s;}",
  #                                      .I,
  #                                      color)],
  #                 collapse = ""),
  #          "</style>\n")
  paste0(
    "<p>",
    # styles,
    paste0(c(flatMapCopy$coloredSeq, "</span>"), collapse = ""),
    "<br><br>",
    legendTbl,
    "</p>")
}

