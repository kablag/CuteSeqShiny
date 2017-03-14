library(RColorBrewer)
library(data.table)
library(pipeR)

addFeatureToLayers <- function(layers, featureId, from, to, compact = TRUE) {
  to.layer <- NA
  if (from < first(layers$seqI))
    from <- first(layers$seqI)
  if (to > last(layers$seqI))
    to <- last(layers$seqI)
  if (compact) {
    for (layer.i in 3:ncol(layers)) {
      if (all(is.na(layers[seqI >= from & seqI <= to][[layer.i]]))) {
        to.layer <- layer.i - 2
        break()
      }
    }
  } else {
    if (all(is.na(layers[seqI >= from & seqI <= to][[ncol(layers)]]))) {
      to.layer <- ncol(layers) - 2
    }
  }
  if (is.na(to.layer))
    to.layer <- ncol(layers) - 1
  layers[seqI >= from & seqI <= to,
         c(sprintf("layer_%s", to.layer)) := featureId]
  NULL
}

addFlat <- function(flatMap, ID, start, end) {
  if ((start < data.table::first(flatMap$seqI) &&
       end < data.table::first(flatMap$seqI)) ||
      (start < data.table::first(flatMap$seqI) &&
       end < data.table::first(flatMap$seqI)))
    return(NULL)
  if (start < data.table::first(flatMap$seqI))
    start <- data.table::first(flatMap$seqI)
  if (end > data.table::last(flatMap$seqI))
    end <- data.table::last(flatMap$seqI)
  flatMap[J(start:end), map := ID]
  NULL
}

getGeneiousTypes <- function(featuresTable) {
  featuresTable[!is.na(note),
                type := sub("Geneious type: (.*)", "\\1", note)]
}


cuteSeq <- function(gbSequence,
                    gbFeatures,
                    gbSequenceStart = 1,
                    colorBy,
                    labelBy,
                    considerStrand = TRUE,
                    includeLegend = TRUE,
                    linesWidth = 60) {
  # gbS <<- gbSequence
  features <- as.data.table(
    gbFeatures)[, ID := .I]
  if (considerStrand) {
    features[, c(colorBy) := paste(get(colorBy), strand)]
  }
  uniqueColorByParams <- unique(features[, get(colorBy)])

  getColors <- function(ncolors) {
    if (ncolors > 12)
      return(sample(colors(distinct = TRUE), ncolors))
    switch(as.character(ncolors),
           "1" = "#8DD3C7",
           "2" = c("#8DD3C7", "#BEBADA"),
           brewer.pal(ncolors, "Set3"))
  }

  colorMap <-
    data.table(gbtype = uniqueColorByParams,
               color = getColors(length(uniqueColorByParams)))
  setkey(colorMap, gbtype)
  flatMap <- data.table(
    seqI = (gbSequenceStart + gbSequence@offset):
      (gbSequence@length + gbSequence@offset),
    gbSeq = strsplit(as.character(gbSequence), NULL)[[1]],
    map = as.numeric(NA)
  )
  setkey(flatMap, seqI)
  features[, addFlat(flatMap, ID, start, end), by = ID]
  # features[, addFeatureToLayers(allLayers[[strand]], ID, start, end), by = ID]

  flatMap[, color := ifelse(is.na(map),
                            "",
                            colorMap[features[map, get(colorBy)], color]),
          by = seqI]
  flatMap[is.na(map), map := 0]
  flatMap[, dif := map != data.table::shift(map, fill = FALSE)]
  flatMap[,
          coloredSeq :=
            ifelse(dif,
                   ifelse(map == 0,
                          sprintf("</span>%s", gbSeq),
                          sprintf("</span><span style='background-color: %s' title='%s'>%s",
                                  color,
                                  features[map, get(labelBy)],
                                  gbSeq)),
                   gbSeq),
          by = seqI]
  if (linesWidth != 0) {
    flatMap[seq(1, nrow(flatMap), by = linesWidth),
            coloredSeq := sprintf("<br>%s", coloredSeq)]
  }
  legendTbl <- ""
  if (includeLegend) {
    legendTbl <-
      colorMap[,
               toprint := paste0(sprintf("<span style='background-color: %s'>%s</span>: ", color, gbtype),
                                 paste(features[get(colorBy) == gbtype, get(labelBy)], collapse = ", ")),
               by = gbtype][, toprint] %>>%
      paste0(collapse = "<br>")
  }
  paste0(paste0(c(flatMap$coloredSeq, "</span>"), collapse = ""),
         "<br><br>",
         legendTbl)
}

