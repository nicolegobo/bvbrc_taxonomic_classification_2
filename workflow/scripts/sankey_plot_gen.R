#!/usr/bin/env Rscript
library(argparse)
library(htmlwidgets)
library(pavian)
library(plyr)

parser <- ArgumentParser(description= 'This script creates a Sankey plot HTML from one or more Kraken2/Bracken reports')
parser$add_argument('--input', '-i', help= 'one or more paths to input report file(s), space-separated')
parser$add_argument('--output', '-o', help= 'path to the output file')
xargs<- parser$parse_args()
### Input is one or more kraken/bracken reports, space-separated
input_reports <- unlist(strsplit(xargs$input, " "))

# get sample name from path
# from microbiome/16S pipelines (bracken) or pathogen pipeline (kraken)
get_sample_name <- function(input_report) {
  if (grepl("bracken_output", input_report)) {
    file_name <- sub(".*/bracken_output/", "", input_report)
    sub("*_bracken_report.txt", "", file_name)
  } else if (grepl("kraken_output", input_report)) {
    file_name <- sub(".*/kraken_output/", "", input_report)
    sub("*_k2_report.txt", "", file_name)
  } else {
    sub("*_report.txt", "", input_report)
  }
}
sample_names <- vapply(input_reports, get_sample_name, character(1), USE.NAMES = FALSE)



### defines the build_sankey_networkfunction from Pavian ###
### this is where the analysis begins
build_sankey_network <- function(my_report, taxRanks =  c("D","K","P","C","O","F","G","S"), maxn=10,
                                 zoom = F, title = NULL,
                                 ...) {
  stopifnot("taxRank" %in% colnames(my_report))
  if (!any(taxRanks %in% my_report$taxRank)) {
    warning("report does not contain any of the taxRanks - skipping it")
    return()
  }
  my_report <- subset(my_report, taxRank %in% taxRanks)
  my_report <- plyr::ddply(my_report, "taxRank", function(x) x[utils::tail(order(x$cladeReads,-x$depth), n=maxn), , drop = FALSE])

  my_report <- my_report[, c("name","taxLineage","taxonReads", "cladeReads","depth", "taxRank")]

  my_report <- my_report[!my_report$name %in% c('-_root'), ]
  #my_report$name <- sub("^-_root.", "", my_report$name)

  splits <- strsplit(my_report$taxLineage, "\\|")

  ## for the root nodes, we'll have to add an 'other' link to account for all cladeReads
  root_nodes <- sapply(splits[sapply(splits, length) ==2], function(x) x[2])

  sel <- sapply(splits, length) >= 3
  splits <- splits[sel]

  links <- data.frame(do.call(rbind,
                              lapply(splits, function(x) utils::tail(x[x %in% my_report$name], n=2))), stringsAsFactors = FALSE)
  colnames(links) <- c("source","target")
  links$value <- my_report[sel,"cladeReads"]

  my_taxRanks <- taxRanks[taxRanks %in% my_report$taxRank]
  taxRank_to_depth <- stats::setNames(seq_along(my_taxRanks)-1, my_taxRanks)


  nodes <- data.frame(name=my_report$name,
                      depth=taxRank_to_depth[my_report$taxRank],
                      value=my_report$cladeReads,
                      stringsAsFactors=FALSE)

  for (node_name in root_nodes) {
    diff_sum_vs_all <- my_report[my_report$name == node_name, "cladeReads"] - sum(links$value[links$source == node_name])
    if (diff_sum_vs_all > 0) {
      nname <- paste("other", sub("^._","",node_name))
      #nname <- node_name
      #links <- rbind(links, data.frame(source=node_name, target=nname, value=diff_sum_vs_all, stringsAsFactors = FALSE))
      #nodes <- rbind(nodes, nname)
    }
  }

  names_id = stats::setNames(seq_len(nrow(nodes)) - 1, nodes[,1])
  links$source <- names_id[links$source]
  links$target <- names_id[links$target]
  links <- links[links$source != links$target, ]

  nodes$name <- sub("^._","", nodes$name)
  links$source_name <- nodes$name[links$source + 1]

  if (!is.null(links))
    sankeyD3::sankeyNetwork(
      Links = links,
      Nodes = nodes,
      doubleclickTogglesChildren = TRUE,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      NodeGroup = "name",
      NodePosX = "depth",
      NodeValue = "value",
      dragY = TRUE,
      xAxisDomain = my_taxRanks,
      numberFormat = "pavian",
      title = title,
      nodeWidth = 15,
      linkGradient = FALSE,
      nodeShadow = TRUE,
      nodeCornerRadius = 5,
      units = "cladeReads",
      fontSize = 12,
      iterations = maxn * 100,
      align = "none",
      highlightChildLinks = TRUE,
      orderByPath = TRUE,
      scaleNodeBreadthsByString = TRUE,
      zoom = zoom,
      ...
    )
}
## end function

reports <- pavian::read_reports(input_reports)

# Build one sankey widget per sample, saved as a standalone self-contained
# fragment in a scratch dir, then stitch the fragments into a single page
# with a dropdown to switch between samples.
tmp_dir <- tempfile("sankey_widgets_")
dir.create(tmp_dir)

fragment_paths <- character(length(reports))
for (i in seq_along(reports)) {
  sankey_widget <- build_sankey_network(reports[[i]], taxRanks = c("D","K","P","C","O","F","G","S"), maxn = 10,
                                          zoom = F)
  frag_path <- file.path(tmp_dir, paste0("sample_", i, ".html"))
  htmlwidgets::saveWidget(sankey_widget, file = frag_path, title = sample_names[i], selfcontained = TRUE)
  fragment_paths[i] <- frag_path
}

# Each fragment is a full standalone HTML document, so it's embedded via
# <iframe srcdoc="..."> rather than concatenated into the page body.
option_tags <- character(length(sample_names))
iframe_tags <- character(length(sample_names))
for (i in seq_along(sample_names)) {
  frag_html <- paste(readLines(fragment_paths[i], warn = FALSE), collapse = "\n")
  frag_html <- gsub("&", "&amp;", frag_html, fixed = TRUE)
  frag_html <- gsub("\"", "&quot;", frag_html, fixed = TRUE)

  display <- if (i == 1) "block" else "none"
  option_tags[i] <- sprintf('<option value="sample-%d">%s</option>', i, sample_names[i])
  iframe_tags[i] <- sprintf(
    '<iframe id="sample-%d" class="sankey-frame" style="display:%s;" srcdoc="%s"></iframe>',
    i, display, frag_html
  )
}

page <- paste0(
  '<!DOCTYPE html>\n<html>\n<head>\n<meta charset="utf-8">\n<title>Multisample Sankey Diagram</title>\n',
  '<style>\n',
  '  body { font-family: sans-serif; margin: 1em; }\n',
  '  .sankey-frame { width: 100%; height: 90vh; border: none; }\n',
  '  #sample-select { font-size: 1em; margin-bottom: 1em; padding: 0.3em; }\n',
  '</style>\n</head>\n<body>\n',
  '<label for="sample-select"><strong>Sample:</strong></label>\n',
  '<select id="sample-select" onchange="showSample(this.value)">\n',
  paste(option_tags, collapse = "\n"), '\n',
  '</select>\n',
  paste(iframe_tags, collapse = "\n"), '\n',
  '<script>\n',
  'function showSample(id) {\n',
  '  var frames = document.getElementsByClassName("sankey-frame");\n',
  '  for (var i = 0; i < frames.length; i++) {\n',
  '    frames[i].style.display = (frames[i].id === id) ? "block" : "none";\n',
  '  }\n',
  '}\n',
  '</script>\n</body>\n</html>\n'
)

writeLines(page, xargs$output)
unlink(tmp_dir, recursive = TRUE)
