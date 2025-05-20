
generate_taxa_tree <- function(treedata,
                               size = 0.2,
                               layout = "circular") {
  ggtree::ggtree(treedata, size = size, layout = layout)
}


create_ps_from_mm <- function(mm, only_marker = TRUE) {
  ot <- otu_table(mm)
  tt <- tax_table(mm)
  st <- sample_data(mm)
  mt <- marker_table(mm)
  sig_features <- mt$feature
  
  # extract all nodes correlated with the significant features
  # First, all parent nodes of marker
  down_nodes <- strsplit(sig_features, "|", fixed = TRUE) %>%
    purrr::map(~ purrr::map_chr(
      seq_along(.x), function(y) paste(.x[1:y], collapse = "|")))
  down_nodes <- unique(unlist(down_nodes))
  # Two, all children nodes of marker
  all_features <- tt@.Data[, 1]
  up_nodes <- purrr::map(sig_features,
                         ~ all_features[grepl(.x, all_features, fixed = TRUE)])
  up_nodes <- unique(unlist(up_nodes))
  idx <- match(unique(c(down_nodes, up_nodes)), all_features)
  
  if (only_marker) {
    ot <- ot[idx, ]
    tt <- tt[idx, ]
  }
  ps <- phyloseq(ot, tt, st)
  
  ps
}


get_treedata_phyloseq <- function(ps, sep = "|") {
  if (!taxa_are_rows(ps)) {
    stop("Requires taxa in rows of phyloseq")
  }
  
  taxa <- tax_table(ps)
  otu <- otu_table(ps)
  row.names(otu) <- taxa@.Data[, 1]
  taxa_nms <- row.names(otu)
  
  tree_table <- data.frame(
    taxa = taxa_nms,
    abd = rowMeans(otu),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      taxa = paste("r__Root", .data$taxa, sep = "|"),
      abd = .data$abd / max(.data$abd) * 100
    )
  
  taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
  nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)
  # add root node
  nodes <- c("r__Root", nodes)
  
  ## data may not contain all the seven ranks of the taxa, such as
  ## enterotypes_arumugam only contains Phylum and Genus ranks
  taxa_deepest <- taxa_split[[which.max(lengths(taxa_split))]]
  prefix <- vector("character", length(taxa_deepest))
  for (i in seq_along(taxa_deepest)) {
    if (!grepl("__$", taxa_deepest[i])) {
      prefix[i] <- gsub("(.*)__.*", "\\1", taxa_deepest[i])
    } else {
      pos <- nchar(taxa_deepest[i]) - 2
      prefix[i] <- substr(taxa_deepest[i], pos, pos)
    }
  }
  
  levels <- purrr::map_chr(nodes, ~ gsub("__.*$", "", .x)) %>%
    factor(levels = rev(prefix))
  
  nodes_parent <- purrr::map_chr(
    taxa_split,
    ~ .x[length(.x) - 1]
  )
  # root must be a parent node
  nodes_parent <- c("root", nodes_parent)
  
  ## tips comes first ?
  is_tip <- !nodes %in% nodes_parent
  index <- vector("integer", length(is_tip))
  index[is_tip] <- seq_len(sum(is_tip))
  index[!is_tip] <- (sum(is_tip) + 1):length(is_tip)
  
  edges <- cbind(
    parent = index[match(nodes_parent, nodes)],
    child = index
  )
  edges <- edges[!is.na(edges[, 1]), ]
  
  # not label the tips
  node_label <- nodes[!is_tip]
  
  phylo <- structure(
    list(
      edge = edges,
      node.label = node_label,
      tip.label = nodes[is_tip],
      edge.length = rep(1, nrow(edges)),
      Nnode = length(node_label)
    ),
    class = "phylo"
  )
  
  mapping <- data.frame(
    node = index,
    abd = c(100, tree_table$abd),
    node_label = nodes,
    stringsAsFactors = FALSE
  )
  mapping$node_class <- levels
  
  tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
}

generate_cladogram_annotation <- function(marker,
                                          color,
                                          sep = "|") {
  enrich_group <- marker$enrich_group
  if (length(color) != length(unique(enrich_group))) {
    stop("the number of colors must be equal to ",
         "the number of enriched groups.")
  }
  
  feature <- marker$feature
  label <- strsplit(feature, split = sep, fixed = TRUE) %>%
    purrr::map_chr(utils::tail, n = 1)
  label_level <- lengths(strsplit(feature, sep, fixed = TRUE))
  
  # may be no marker are identified enriched in some groups
  # drop the levels of this groups if the enrich_group is a factor
  if (inherits(enrich_group, "factor")) {
    enrich_group <- droplevels(enrich_group)
  }
  
  # named colors: set the colors based on the matched names to groups
  if (is.vector(color) && !is.null(names(color))) {
    if (!all(names(color) %in% enrich_group)) {
      stop("names of `color` muste be contained in enriched groups")
    }
    color <- color[match(enrich_group, names(color))]
  } else {
    # colors will be matched in order (usually alphabetical) with the groups
    names(color) <- sort(unique(enrich_group))
    color <- color[match(enrich_group, names(color))]
  }
  
  annotation <- data.frame(
    node = label,
    color = color,
    enrich_group = enrich_group,
    stringsAsFactors = FALSE
  )
  
  annotation
}



  
get_offset <- function(x) {
  (x * 0.2 + 0.2)^2
}


get_angle <- function(tree, node) {
  if (length(node) != 1) {
    stop("The length of `node` must be 1")
  }
  tree_data <- tree$data
  sp <- tidytree::offspring(tree_data, node)$node
  sp2 <- c(sp, node)
  sp.df <- tree_data[match(sp2, tree_data$node), ]
  mean(range(sp.df$angle))
}


get_short_label_id <- function(clade_label, clade_label_level) {
  ind <- clade_label$level < clade_label_level
  unique_id <- get_unique_id(sum(ind))
  short_label <- unique_id[seq_len(sum(ind))]
  
  short_label
}


get_unique_id <- function(n, depth = 1) {
  args <- lapply(seq_len(depth), FUN = function(x) letters)
  x <- do.call(expand.grid, args = list(args, stringsAsFactors = FALSE))
  x <- x[, rev(names(x)), drop = FALSE]
  x <- do.call(paste0, x)
  if (n <= length(x)) {
    return(x[seq_len(n)])
  }
  
  return(c(x, get_unique_id(n - length(x), depth = depth + 1)))
}



set_marker_annotation <- function(p,
                                  color,
                                  label,
                                  size = 5,
                                  shape = 22,
                                  ...) {
  dat <- data.frame(
    color = color,
    label = label,
    stringsAsFactors = FALSE
  )
  
  # suppress warning: The shape palette can deal with a maximum of 6 discrete
  # values because more than 6 becomes difficult to discriminate; you have 18.
  # Consider specifying shapes manually if you must have them.
  # using scale_shape_manual
  p <- p +
    geom_point(
      data = dat, inherit.aes = FALSE,
      aes_(x = 0, y = 0, shape = ~label),
      size = 0, stroke = 0,
    ) +
    scale_shape_manual(values = rep(shape, nrow(dat)), limits = dat$label) +
    guides(
      shape = guide_legend(
        override.aes = list(
          size = size,
          shape = shape,
          fill = dat$color
        ),
        order = 2,
        ...
      )
    )
  
  p
}



plot_cladogram_rw <- function (mm, color, only_marker = FALSE, branch_size = 0.2, 
          alpha = 0.2, node_size_scale = 1, node_size_offset = 1, clade_label_level = 4, 
          clade_label_font_size = 4, annotation_shape = 22, annotation_shape_size = 5, 
          group_legend_param = list(), marker_legend_param = list()) 
{

  ps <- create_ps_from_mm(mm, only_marker = only_marker)
  # ps <- create_ps_from_mm(mm_lefse_1, only_marker = F)
  tree <- get_treedata_phyloseq(ps) 
  #tree@data <- tree@data$[which(tree@data$node_label=="")]
  tree <- generate_taxa_tree(tree, size = branch_size)
  annotation <- generate_cladogram_annotation(mm@marker_table, 
                                              color = color)
  annotation_info <- dplyr::left_join(annotation, tree$data, 
                                      by = c(node = "label")) %>% mutate(label = .data$node, 
                                                                         id = .data$node.y, level = as.numeric(.data$node_class))
  hilight_para <- dplyr::transmute(annotation_info, node = .data$id, 
                                   fill = .data$color, alpha = alpha, extend = get_offset(.data$level))
  hilights_g <- purrr::pmap(hilight_para, geom_hilight)
  tree <- purrr::reduce(hilights_g, `+`, .init = tree)
  hilights_df <- dplyr::distinct(annotation_info, .data$enrich_group, 
                                 .data$color) %>% arrange(.data$enrich_group)
  hilights_df$x <- 0
  hilights_df$y <- 1
  # group_legend_param <- c(group_legend_param), list(title = NULL, 
  #                                                  order = 1, override.aes = list(fill = hilights_df$color)))
  #group_lgd <- do.call(guide_legend, group_legend_param)
  tree <- tree + geom_rect(aes_(xmin = ~x, xmax = ~x, ymax = ~y, 
                                ymin = ~y, fill = ~enrich_group), data = hilights_df, 
                           inherit.aes = FALSE) + guides(fill = FALSE)
  nodes_colors <- rep("white", nrow(tree$data))
  nodes_colors[annotation_info$id] <- annotation_info$color
  node_size <- node_size_scale * log(tree$data$abd) + node_size_offset
  tree$data$node_size <- node_size
  tree <- tree + geom_point2(aes(size = I(node_size)), fill = nodes_colors, 
                             shape = 21)
  clade_label <- dplyr::transmute(annotation_info, node = .data$id, 
                                  offset = get_offset(.data$level) - 0.4, angle = purrr::map_dbl(.data$id, 
                                                                                                 get_angle, tree = tree) + 90, label = .data$label, 
                                  fontsize = clade_label_font_size, barsize = 0, hjust = 0.5, 
                                  level = .data$level) %>% dplyr::arrange(desc(.data$level))
  ind <- clade_label$level < clade_label_level
  short_label <- get_short_label_id(clade_label, clade_label_level)
  clade_label_para <- mutate(clade_label, label = c(.data$label[!ind], 
                                                    short_label), level = NULL)
  clade_label_g <- purrr::pmap(clade_label_para, geom_cladelabel)
  tree <- purrr::reduce(clade_label_g, `+`, .init = tree)
  guide_label <- clade_label[ind, ] %>% mutate(label2 = paste0(short_label, 
                                                               ": ", .data$label), color = annotation_info$color[match(.data$label, 
                                                                                                                       annotation_info$label)])
  marker_legend_param <- c(marker_legend_param, list(p = tree, 
                                                     color = guide_label$color, label = guide_label$label2, 
                                                     shape = annotation_shape, size = annotation_shape_size))
  p <- do.call(set_marker_annotation, marker_legend_param) + 
    theme(legend.position = "right", legend.title = element_blank())
  p
}






# # mm, color, only_marker = FALSE, branch_size = 0.2,
# # alpha = 0.2, node_size_scale = 1, node_size_offset = 1, clade_label_level = 4,
# # clade_label_font_size = 4, annotation_shape = 22, annotation_shape_size = 5,
# # group_legend_param = list(), marker_legend_param = list()
# 
# ps <- create_ps_from_mm(mm_lefse_1, only_marker = FALSE)
# pretree <- get_treedata_phyloseq(ps)
# tree <- as.phylo(pretree)
# tree <- generate_taxa_tree(pretree, size = branch_size, layout = "circular")
# annotation <- generate_cladogram_annotation(mm_lefse_1@marker_table,
#                                             color = color)
# annotation_info <- dplyr::left_join(annotation, tree$data,
#                                     by = c(node = "label")) %>% mutate(label = .data$node,
#                                                                        id = .data$node.y, level = as.numeric(.data$node_class))
# hilight_para <- dplyr::transmute(annotation_info, node = .data$id,
#                                  fill = .data$color, alpha = alpha, extend = get_offset(.data$level))
# hilights_g <- purrr::pmap(hilight_para, geom_hilight)
# tree <- purrr::reduce(hilights_g, `+`, .init = tree)
# hilights_df <- dplyr::distinct(annotation_info, .data$enrich_group,
#                                .data$color) %>% arrange(.data$enrich_group)
# hilights_df$x <- 0
# hilights_df$y <- 1
# # group_legend_param <- c(group_legend_param), list(title = NULL,
# #                                                  order = 1, override.aes = list(fill = hilights_df$color)))
# #group_lgd <- do.call(guide_legend, group_legend_param)
# tree <- tree + geom_rect(aes_(xmin = ~x, xmax = ~x, ymax = ~y,
#                               ymin = ~y, fill = ~enrich_group), data = hilights_df,
#                          inherit.aes = FALSE) + guides(fill = FALSE)
# nodes_colors <- rep("white", nrow(tree$data))
# nodes_colors[annotation_info$id] <- annotation_info$color
# node_size <- node_size_scale * log(tree$data$abd) + node_size_offset
# tree$data$node_size <- node_size
# tree <- tree + geom_point2(aes(size = I(node_size)), fill = nodes_colors,
#                            shape = 21)
# clade_label <- dplyr::transmute(annotation_info, node = .data$id,
#                                 offset = get_offset(.data$level) - 0.4, angle = purrr::map_dbl(.data$id,
#                                                                                                get_angle, tree = tree) + 90, label = .data$label,
#                                 fontsize = clade_label_font_size, barsize = 0, hjust = 0.5,
#                                 level = .data$level) %>% dplyr::arrange(desc(.data$level))
# ind <- clade_label$level < clade_label_level
# short_label <- get_short_label_id(clade_label, clade_label_level)
# clade_label_para <- mutate(clade_label, label = c(.data$label[!ind],
#                                                   short_label), level = NULL)
# clade_label_g <- purrr::pmap(clade_label_para, geom_cladelabel)
# tree <- purrr::reduce(clade_label_g, `+`, .init = tree)
# guide_label <- clade_label[ind, ] %>% mutate(label2 = paste0(short_label,
#                                                              ": ", .data$label), color = annotation_info$color[match(.data$label,
#                                                                                                                      annotation_info$label)])
# marker_legend_param <- c(marker_legend_param, list(p = tree,
#                                                    color = guide_label$color, label = guide_label$label2,
#                                                    shape = annotation_shape, size = annotation_shape_size))
# p <- do.call(set_marker_annotation, marker_legend_param) +
#   theme(legend.position = "right", legend.title = element_blank())
# p









check_rank_names <- function(ps) {
  ps_ranks <- rank_names(ps)
  if (is_picrust2(ps)) {
    picrust_rank <- c("Picrust_trait", "Picrust_description")
    diff_rank <- setdiff(ps_ranks, picrust_rank)
    if (length(diff_rank)) {
      stop("ranks of picrust2 functional profile must be one of ",
           paste(picrust_rank, collapse = ", "))
    }
  } else {
    if (!all(ps_ranks %in% available_ranks)) {
      stop(
        "ranks of taxonimic profile must be one of ",
        paste(available_ranks, collapse = ", ")
      )
    }
  }
  
  invisible(ps)
}



check_taxa_rank <- function(ps, taxa_rank) {
  ranks <- rank_names(ps)
  all_taxa_rank <- c("all", "none", ranks)
  if (!taxa_rank %in% all_taxa_rank) {
    stop(
      "`taxa_rank` must be one of ",
      paste(all_taxa_rank, collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(ps)
}




is_picrust2 <- function(ps) {
  ps_ranks <- rank_names(ps)
  if ("Picrust_trait" %in% ps_ranks) TRUE else FALSE
}



# available taxonomic ranks, Summarize represents summarized tax
available_ranks <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)
available_ranks <- factor(
  available_ranks,
  levels = available_ranks
)

#usethis::use_data(available_ranks, internal = TRUE, overwrite = TRUE)

check_tax_summarize <- function(ps) {
  taxa <- row.names(otu_table(ps))
  # whether taxa is separated by `|`,
  # may be required to add extra separate strings in the future
  has_separate <- any(grepl("[|]", taxa))
  
  has_separate
}



keep_taxa_in_rows <- function(ps) {
  if (!taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  
  ps
}


phyloseq_qc <- function(ps) {
  prune_taxa(taxa_sums(ps) > 0, ps)
}


fix_duplicate_tax <- function(ps) {
  # convert na to Unknown first
  ps <- fix_na_tax(ps)
  
  tax <- tax_table(ps)
  if (ncol(tax) == 1) {
    return(ps)
  }
  
  for (i in 2:ncol(tax)) {
    tax_uniq <- unique(tax[, i])
    for (j in seq_along(tax_uniq)) {
      if (is.na(tax_uniq[j])) next
      ind <- which(tax[, i] == as.character(tax_uniq[j]))
      if (length(unique(tax[ind, i - 1])) > 1) {
        tax[ind, i] <- paste(tax[ind, i - 1], tax[ind, i], sep = "_")
      }
    }
  }
  
  tax_table(ps) <- tax
  
  ps
}


fix_na_tax <- function(ps) {
  tax <- as.data.frame(tax_table(ps))
  
  tax_fixed <- purrr::imap_dfc(
    tax,
    ~ ifelse(is.na(.x), NA, .x) # changed get_prefix(.y) to NA
  ) %>%
    as.matrix()
  row.names(tax_fixed) <- taxa_names(ps)
  tax_table(ps) <- tax_fixed
  
  ps
}


get_prefix <- function(ranks) {
  prefix <- substr(ranks, 1, 1) %>%
    tolower() %>%
    paste("__", sep = "")
  
  prefix
}


preprocess_ps <- function(ps) {
  zero_sample <- check_samples(ps)
  if (!is.null(zero_sample)) {
    warning(
      "The library size of sample(s): ",
      paste(zero_sample, collapse = ", "),
      " is/are zero, and will be removed in the subsequent analysis."
    )
    
    keep <- setdiff(sample_names(ps), zero_sample)
    ps <- prune_samples(keep, ps)
  }
  
  # keep taxa in rows
  ps <- keep_taxa_in_rows(ps)
  # filter the taxa whose abundance is zero
  ps <- phyloseq_qc(ps)
  # fix duplicated tax
  ps <- fix_duplicate_tax(ps)
  
  ps
}

# check samples in ps, make sure at least one non zero features in a sample
check_samples <- function(ps) {
  if (!taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  lib_size <- colSums(otu_table(ps))
  zero_ind <- which(lib_size == 0)
  
  if (length(zero_ind) == 0) {
    return(NULL)
  }
  
  return(sample_names(ps)[zero_ind])
}



lefse_format_grp <- function(sample_meta, group, subgroup = NULL) {
  groups <- sample_meta[[group]]
  group_nms <- unique(groups)
  
  if (is.null(subgroup)) {
    subgroups <- paste0(groups, "_subgrp")
  } else {
    subgroups <- paste(groups, sample_meta[[subgroup]], sep = "_")
  }
  
  group_hie <- split(subgroups, groups) %>%
    purrr::map(unique)
  
  return(list(group = groups, subgroup = subgroups, group_hie = group_hie))
}


add_missing_levels <- function(feature) {
  if (!taxa_are_rows(feature)) {
    feature <- t(feature)
  }
  feature_nms <- row.names(feature)
  feature <- feature@.Data %>% data.frame()
  
  # the missing feature names
  feature_nms2 <-
    strsplit(feature_nms, "|", fixed = TRUE) %>%
    purrr::map(
      ~ Reduce(
        function(x, y) paste(x, y, sep = "|"),
        .x,
        accumulate = TRUE
      )
    )
  
  unq_nms <- unlist(feature_nms2) %>% unique()
  missing_nms <- setdiff(unq_nms, feature_nms)
  
  if (length(missing_nms) == 0) {
    return(feature)
  }
  missing_nms_split <- strsplit(missing_nms, split = "|", fixed = TRUE)
  missing_mns_level <- lengths(missing_nms_split)
  missing_level_range <- range(missing_mns_level)
  
  # only sum the next level of tax, so we need first add the missing tax at
  # the most depth level
  for (i in missing_level_range[2]:missing_level_range[1]) {
    missing_nms_i <- missing_nms[missing_mns_level == i]
    taxs <- row.names(feature)
    
    indx <- purrr::map(
      missing_nms_i,
      ~ purrr::map_lgl(taxs, function(x) grepl(.x, x, fixed = TRUE))
    )
    
    # only sum the next level of tax
    feature_nms_level <- strsplit(taxs, split = "|", fixed = TRUE) %>%
      lengths()
    indx <- purrr::map(indx, ~ .x & feature_nms_level == (i + 1))
    abd <- purrr::map_df(
      feature,
      ~ purrr::map_dbl(indx, function(x) sum(.x[x]))
    )
    
    feature <- rbind(feature, abd)
    row.names(feature) <- c(taxs, missing_nms_i)
  }
  
  otu_table(feature, taxa_are_rows = TRUE)
}
  
summarize_taxa <- function(ps,
                           level = rank_names(ps)[1],
                           absolute = TRUE,
                           sep = "|") {
  
  # return ps if it has been summarized
  summarized <- check_tax_summarize(ps)
  if (summarized) {
    otu_summarized <- otu_table(ps) %>%
      add_missing_levels()
    tax_summarized <- row.names(otu_summarized) %>%
      matrix() %>%
      tax_table()
    row.names(tax_summarized) <- row.names(otu_summarized)
    return(phyloseq(otu_summarized, tax_summarized, sample_data(ps)))
  }
  
  if (!has_prefix(ps)) {
    ps <- add_prefix(ps)
    ps@tax_table[,2] <- str_replace(ps@tax_table[,2], "p__NA", replacement = "")
    ps@tax_table[,3] <- str_replace(ps@tax_table[,3], "c__NA", replacement = "")
    ps@tax_table[,4] <- str_replace(ps@tax_table[,4], "o__NA", replacement = "")
    ps@tax_table[,5] <- str_replace(ps@tax_table[,5], "f__NA", replacement = "")
    ps@tax_table[,6] <- str_replace(ps@tax_table[,6], "g__NA", replacement = "")
  }
  
  ps_ranks <- rank_names(ps)
  if (!level %in% ps_ranks) {
    stop("`level` must in the ranks of `ps` (rank_names(ps))")
  }
  
  ind <- match(level, ps_ranks)
  levels <- ps_ranks[ind:length(ps_ranks)]
  res <- purrr::map(
    levels,
    ~ .summarize_taxa_level(
      ps,
      rank = .x,
      absolute = absolute,
      sep = sep
    )
  )
  tax_nms <- purrr::map(res, row.names) %>% unlist()
  res <- bind_rows(res)
  row.names(res) <- tax_nms
  
  otu_summarized <- otu_table(res, taxa_are_rows = TRUE)
  tax_summarized <- row.names(otu_summarized) %>%
    matrix() %>%
    tax_table()
  row.names(tax_summarized) <- row.names(otu_summarized)
  row.names(otu_summarized) <- row.names(tax_summarized)
  
  # To ensure the rank of the summarized object is valid (one of  "domain"  
  # "kingdom" "phylum"  "class"   "order"   "family"  "genus"   "species"),  
  # set it (column names of tax_summarized) as the top level rank in the ps
  # object.
  #  
  # rank_prefix <- extract_prefix(ps_ranks)
  # colnames(tax_summarized) <- paste0(rank_prefix, collapse = sep)
  colnames(tax_summarized) <- ps_ranks[1]
  
  return(phyloseq(otu_summarized, tax_summarized, sample_data(ps)))
}

has_prefix <- function(ps) {
  sample_tax <- tax_table(ps)[1, 1]
  if (substr(sample_tax, 2, 3) == "__") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# add ranks prefix, e.g k__, p__, only worked for unsummarized data
add_prefix <- function(ps) {
  tax <- as(tax_table(ps), "matrix") %>%
    as.data.frame()
  lvl <- colnames(tax)
  prefix <- get_prefix(lvl)
  
  tax_new <- mapply(function(x, y) paste0(x, y), 
                    prefix, tax, SIMPLIFY = FALSE)
  tax_new <- do.call(cbind, tax_new)
  row.names(tax_new) <- row.names(tax)
  colnames(tax_new) <- lvl
  tax_table(ps) <- tax_new
  
  ps
}



pre_ps_taxa_rank <- function(ps, taxa_rank) {
  if (is_picrust2(ps)) {
    warning("para `taxa_rank` is not worked for picrust2 function profile",
            " and it will be ignored")
    return(ps)
  }
  
  ps <- check_taxa_rank(ps, taxa_rank)
  if (taxa_rank == "all") {
    ps_orig_summarized <- summarize_taxa(ps)
  } else if (taxa_rank == "none") {
    ps_orig_summarized <- extract_rank(ps, taxa_rank)
  } else {
    ps_orig_summarized <- aggregate_taxa(ps, taxa_rank) %>%
      extract_rank(taxa_rank)
  }
  
  return(ps_orig_summarized)
}



test_rep_wilcoxon <- function(subcls,
                              cls_hie,
                              feats_abd,
                              feats_name,
                              strict = 0,
                              wilcoxon_cutoff = 0.05,
                              multicls_strat = FALSE,
                              sample_min = 10,
                              only_same_subcls = FALSE,
                              curv = FALSE) {
  if (!strict %in% c(0, 1, 2)) {
    stop("`strict` must be 0, 1 or 2")
  }
  
  cls_nms <- names(cls_hie)
  pairs <- utils::combn(cls_nms, 2, simplify = FALSE)
  tot_ok <- 0
  
  all_diff <- list()
  
  for (pair in pairs) {
    dir_cmp <- "not_set"
    
    subcls1 <- cls_hie[[pair[1]]]
    subcls1_n <- length(subcls1)
    subcls2 <- cls_hie[[pair[2]]]
    subcls2_n <- length(subcls2)
    
    # multiple tests
    if (strict != 0) {
      wilcoxon_cutoff <- ifelse(
        strict == 2,
        wilcoxon_cutoff * subcls1_n * subcls2_n,
        1 - (1 - wilcoxon_cutoff)^(subcls1_n * subcls2_n)
      )
    }
    
    ok <- 0
    curv_sign <- 0
    first <- TRUE
    
    for (i in seq_along(subcls1)) {
      br <- FALSE
      for (j in seq_along(subcls2)) {
        if (only_same_subcls &&
            gsub(pair[1], "", subcls1[i]) !=
            gsub(pair[2], "", subcls2[j])) {
          ok <- ok + 1
          next
        }
        
        cls1_abd <- feats_abd[subcls == subcls1[i]]
        cls2_abd <- feats_abd[subcls == subcls2[j]]
        med_comp <- FALSE
        
        if (length(cls1_abd) < sample_min ||
            length(cls2_abd) < sample_min) {
          med_comp <- TRUE
        }
        
        sx <- stats::median(cls1_abd)
        sy <- stats::median(cls2_abd)
        
        if (cls1_abd[1] == cls2_abd[1] &&
            length(unique(cls1_abd)) == 1 &&
            length(unique(cls2_abd)) == 1
        ) {
          tres <- FALSE
          first <- FALSE
        } else if (!med_comp) {
          x <- c(cls1_abd, cls2_abd)
          y <- factor(
            c(
              rep(1, length(cls1_abd)),
              rep(2, length(cls2_abd))
            )
          )
          pv <- coin::wilcox_test(
            x ~ y,
            data = data.frame(x, y)
          )
          pv <- coin::pvalue(pv)
          tres <- pv < wilcoxon_cutoff * 2
        }
        
        if (first) {
          first <- FALSE
          
          if (!curv && (med_comp || tres)) {
            dir_cmp <- sx < sy
          } else if (curv) {
            dir_cmp <- NULL
            if (med_comp || tres) {
              curv_sign <- curv_sign + 1
              dir_cmp <- sx < sy
            }
          } else {
            br <- TRUE
          }
        } else if (!curv && med_comp) {
          if ((sx < sy) != dir_cmp || sx == sy) {
            br <- TRUE
          }
        } else if (curv) {
          if (tres && is.null(dir_cmp)) {
            curv_sign <- curv_sign + 1
            dir_cmp <- sx < sy
          }
          
          if (tres && dir_cmp != (sx < sy)) {
            br <- TRUE
            curv_sign <- curv_sign - 1
          }
        } else if (!tres || (sx < sy) != dir_cmp || sx == sy) {
          br <- TRUE
        }
        
        if (br) {
          break
        }
        ok <- ok + 1
      }
      if (br) {
        break
      }
    }
    
    diff <- ifelse(
      curv,
      curv_sign > 0,
      ok == subcls1_n * subcls2_n
    )
    
    if (diff) tot_ok <- tot_ok + 1
    if (!diff && multicls_strat) {
      return(FALSE)
    }
    if (diff && !multicls_strat) all_diff <- c(all_diff, pair)
  }
  
  if (!multicls_strat) {
    tot_k <- length(cls_hie)
    
    for (k in names(cls_hie)) {
      nk <- 0
      for (a in all_diff) {
        if (k %in% a) nk <- nk + 1
      }
      
      if (nk == tot_k - 1) {
        return(TRUE)
      }
    }
    
    return(FALSE)
  }
  
  return(TRUE)
}




get_feature_enrich_group <- function(class, feature) {
  feature$class <- class
  feature_mean <- group_by(feature, class) %>%
    group_modify(~ purrr::map_df(.x, mean)) %>%
    ungroup()
  
  feature_enrich_index <- select(feature_mean, -class) %>%
    purrr::map_dbl(which.max)
  feature_enrich_group <- feature_mean$class[feature_enrich_index]
  names(feature_enrich_group) <- names(feature)[names(feature) != "class"]
  feature_max_mean <- purrr::map2_dbl(
    select(feature_mean, -class),
    feature_enrich_index,
    ~ .x[.y]
  )
  feature_max_mean[feature_max_mean < 1] <- 1
  
  return(list(
    group = feature_enrich_group,
    log_max_mean = log10(feature_max_mean)
  ))
}


bootstap_lda <- function(feature_abundance,
                         boot_n,
                         class,
                         sample_fract,
                         seed = 2020) {
  # Bioconductor not allows set.seed
  ldas <- purrr::map(
    1:boot_n,\(i)
    bootstap_lda_one(
      feature_abundance,
      class,
      sample_fract
    )
  ) %>%
    purrr::transpose() %>%
    purrr::map(~ do.call(bind_rows, .x)) %>%
    bind_rows()
  
  mean_lds <- colMeans(ldas)
  mean_lds <- sign(mean_lds) * log10(1 + abs(mean_lds))
  
  mean_lds
}



bootstap_lda_one <- function(feature_abundance,
                             class,
                             sample_fract) {
  sample_groups <- unique(class)
  class_count <- table(class)
  feature_abundance$class <- class
  feature_abundance <- preprocess_feature_all(feature_abundance, class)
  
  sample_n <- nrow(feature_abundance)
  random_n <- floor(sample_n * sample_fract)
  class_n <- length(sample_groups)
  sample_min <- floor(
    min(class_count) * sample_fract * sample_fract * 0.5) %>%
    max(1)
  
  # class vs class
  pairs <- utils::combn(sample_groups, 2, simplify = FALSE) %>%
    purrr::map(sort, decreasing = TRUE)
  
  for (i in seq_len(1000)) {
    # random select samples using bootstrap method
    sample_indx <- sample(sample_n, random_n, replace = TRUE)
    
    is_checked <- check_bootstrap_sample(
      feature_abundance,
      sample_indx,
      sample_min,
      class
    )
    if (is_checked) {
      break
    }
  }
  
  if (!is_checked) {
    stop(
      "Too small samples in each class",
      " or the variance of feature abundances within a",
      " class too small (zero or near zero)",
      call. = FALSE
    )
  }
  
  lda <- purrr::map(
    pairs,
    ~ cal_pair_lda(feature_abundance, sample_indx, .x)
  )
  names(lda) <- purrr::map(pairs, paste, collapse = " -VS- ")
  
  lda
}




preprocess_feature_all <- function(x, class) {
  res <- group_by(x, class) %>%
    group_modify(~ purrr::map_df(.x, preprocess_feature)) %>%
    ungroup()
  
  res
}

preprocess_feature <- function(x) {
  if (length(unique(x)) <= max(length(x) * 0.5, 4)) {
    x <- purrr::map_dbl(x, ~ abs(.x + rnorm(1, 0, max(.x * 0.05, 0.01))))
  }
  
  x
}


check_bootstrap_sample <- function(feature_abundance,
                                   sample_indx,
                                   sample_min,
                                   class) {
  if ("class" %in% names(feature_abundance)) {
    feature_abundance$class <- NULL
  }
  feature_abundance <- feature_abundance[sample_indx, ]
  class_n <- length(unique(class))
  class <- class[sample_indx]
  
  if (length(unique(class)) < class_n) {
    return(FALSE)
  }
  
  for (cls in unique(class)) {
    if (sum(class == cls) < sample_min) {
      return(FALSE)
    }
    
    # sig feature smaller than min sample count
    cls_abundance <- feature_abundance[class == cls, ]
    
    for (i in seq_along(ncol(cls_abundance))) {
      unique_abd <- length(unique(cls_abundance[[i]]))
      
      if ((unique_abd <= sample_min && sample_min > 1) ||
          (unique_abd <= 1 && sample_min == 1)) {
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}



cal_pair_lda <- function(feature_abundance,
                         sample_indx,
                         pair) {
  sample_feature_abundance <- feature_abundance[sample_indx, ]
  
  # reference lefse.py in lefse
  lda_res <- suppressWarnings(
    MASS::lda(
      class ~ .,
      data = sample_feature_abundance,
      tol = 1.0e-10
    )
  )
  w <- lda_res$scaling[, 1]
  w_unit <- w / sqrt(sum(w^2))
  feature_remove_class <- sample_feature_abundance[-1]
  
  # not support subclass and subject argument in lefse
  
  ld <- as.matrix(feature_remove_class) %*% w_unit
  group1_indx <- sample_feature_abundance$class == pair[1]
  effect_size <- abs(mean(ld[group1_indx]) - mean(ld[-group1_indx]))
  wfinal <- w_unit * effect_size
  lda_means <- lda_res$means
  lda_row_nms <- row.names(lda_means)
  feature_n <- ncol(lda_means)
  coeff <- ifelse(is.nan(wfinal), 0, abs(wfinal))
  
  res <- purrr::map(
    pair,
    function(x) {
      if (x %in% lda_row_nms) {
        # fixes #7, Since `pair` is a level, and `lda_means[pair[i], ]`
        # corced pair[i]` to numeric rather than use the corresponding
        # level of pair[i]
        ind <- match(x, lda_row_nms)
        lda_means[ind, ]
      } else {
        rep(0, feature_n)
      }
    }
  )
  
  names(res) <- pair
  
  feature <- names(feature_remove_class)
  lda_score <- purrr::map_dbl(
    seq_along(feature),
    function(i) {
      gm <- abs(res[[1]][i] - res[[2]][i])
      return(gm + coeff[i] * 0.5)
    }
  )
  names(lda_score) <- feature
  
  lda_score
}


return_marker <- function(sig_feature, all_feature) {
  if (nrow(sig_feature)) {
    row.names(sig_feature) <- paste0("marker", seq_len(nrow(sig_feature)))
    marker <- marker_table(sig_feature)
    
  } else {
    warning("No marker was identified", call. = FALSE)
    marker <- NULL
  }
  
  marker
}


get_norm_method <- function(norm) {
  new_norm <- ifelse(
    is.numeric(norm),
    paste("per-sample normalized (sum of all taxa) to", norm),
    norm
  )
  
  new_norm
}


.summarize_taxa_level <- function(ps,
                                  rank_name,
                                  absolute = TRUE,
                                  sep = "|") {
  if (!absolute) {
    ps <- transform_sample_counts(ps, function(x) x / sum(x))
  }
  
  otus <- otu_table(ps)
  otus_extend <- slot(otus, ".Data") %>%
    tibble::as_tibble()
  
  taxas <- tax_table(ps)@.Data %>%
    tibble::as_tibble()
  
  ranks <- setdiff(available_ranks, "Summarize")
  rank_level <- match(rank_name, ranks)
  select_ranks <- intersect(ranks[seq_len(rank_level)], rank_names(ps))
  
  consensus <- taxas[, select_ranks] %>%
    purrr::pmap_chr(paste, sep = sep)
  otus_extend$consensus <- consensus
  
  taxa_summarized <- group_split(otus_extend, consensus) %>%
    purrr::map(.sum_consensus)
  taxa_summarized <- do.call(rbind, taxa_summarized)
  # filter taxa of which abundance is zero
  ind <- rowSums(taxa_summarized) != 0
  taxa_summarized <- taxa_summarized[ind, ]
  
  taxa_summarized
}


.sum_consensus <- function(x) {
  consensus <- unique(x$consensus)
  if (length(consensus) != 1) {
    stop("consensus in the same group muste be the same")
  }
  
  x$consensus <- NULL
  res <- as.data.frame(t(colSums(x)))
  row.names(res) <- consensus
  
  return(res)
}


run_lefse_rw <- function (ps, group, subgroup = NULL, taxa_rank = "all", transform = c("identity", 
                                                                       "log10", "log10p"), norm = "CPM", norm_para = list(), kw_cutoff = 0.05, 
          lda_cutoff = 2, bootstrap_n = 30, bootstrap_fraction = 2/3, 
          wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("0", 
                                                                     "1", "2"), sample_min = 10, only_same_subgrp = FALSE, 
          curv = FALSE) 
{
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be phyloseq object", call. = FALSE)
  }
  ps <- check_rank_names(ps)
  ps <- check_taxa_rank(ps, taxa_rank)
  transform <- match.arg(transform, c("identity", "log10", 
                                      "log10p"))
  strict <- match.arg(strict, c("0", "1", "2"))
  strict <- as.numeric(strict)
  summarized <- check_tax_summarize(ps)
  if (summarized && norm != "CPM") {
    stop("`norm` must be a 'CPM' or 'none' while `ps` has been summarized", 
         call. = FALSE)
  }
  ps <- preprocess_ps(ps)
  ps <- transform_abundances(ps, transform = transform)
  norm_para <- c(norm_para, method = norm, object = list(ps))
  ps_normed <- do.call(normalize, norm_para)
  sample_meta <- sample_data(ps_normed)
  grp_info <- lefse_format_grp(sample_meta, group, subgroup = subgroup)
  grp <- grp_info$group
  subgrps <- grp_info$subgroup
  grp_hie <- grp_info$group_hie
  ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
  otus <- abundances(ps_summarized, norm = TRUE)
  otus_test <- as.data.frame(t(otus), stringsAsFactors = FALSE)
  feature <- tax_table(ps_summarized)@.Data[, 1]
  names(otus_test) <- feature
  tax <- matrix(feature) %>% tax_table()
  row.names(tax) <- row.names(otus)
  kw_p <- purrr::map_dbl(otus_test, ~kruskal.test(.x, grp)$p.value)
  na_ind <- is.na(kw_p)
  if (sum(na_ind) >= 1) {
    otus_test <- otus_test[!na_ind]
    kw_p <- kw_p[!na_ind]
  }
  sig_ind <- kw_p <= kw_cutoff
  sig_otus <- otus_test[, sig_ind]
  features_nms <- names(sig_otus)
  wilcoxon_p <- purrr::map2_lgl(sig_otus, features_nms, ~test_rep_wilcoxon(subgrps, 
                                                                           grp_hie, .x, .y, wilcoxon_cutoff = wilcoxon_cutoff, multicls_strat = multigrp_strat, 
                                                                           strict = strict, sample_min = sample_min, only_same_subcls = only_same_subgrp, 
                                                                           curv = curv))
  sig_otus <- sig_otus[, wilcoxon_p, drop = FALSE]
  if (ncol(sig_otus) == 0) {
    warning("No marker was identified", call. = FALSE)
    mm <- microbiomeMarker(marker_table = NULL, norm_method = get_norm_method(norm), 
                           diff_method = "lefse", otu_table = otu_table(otus, 
                                                                        taxa_are_rows = TRUE), sam_data = sample_data(ps_normed), 
                           tax_table = tax)
    return(mm)
  }
  otus_enriched_group <- get_feature_enrich_group(grp, sig_otus)
  ldas <- bootstap_lda(sig_otus, boot_n = bootstrap_n, class = grp, 
                       sample_fract = bootstrap_fraction)
  lefse_res <- data.frame(feature = names(sig_otus), enrich_group = otus_enriched_group$group, 
                          ef_lda = ldas, pvalue = kw_p[sig_ind][wilcoxon_p], stringsAsFactors = FALSE)
  lefse_sig <- filter(lefse_res, .data$ef_lda >= lda_cutoff) %>% 
    arrange(.data$enrich_group, desc(.data$ef_lda))
  lefse_out <- return_marker(lefse_sig, lefse_res)
  lefse_out$padj <- lefse_out$pvalue
  mm <- microbiomeMarker(marker_table = lefse_out, norm_method = get_norm_method(norm), 
                         diff_method = "lefse", otu_table = otu_table(otus, taxa_are_rows = TRUE), 
                         sam_data = sample_data(ps_normed), tax_table = tax)
  mm
}









