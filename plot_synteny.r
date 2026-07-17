library(ggplot2)
library(dplyr)
library(patchwork)

rh_ru_block_file  <- "rhodnius2rubida.block.txt"
rh_sa_block_file  <- "rhodnius2sanguisuga.block.txt"
ru_sa_block_file  <- "rubida2sanguisuga.block.txt"
rh_sa_anchor_file <- "rhodnius2sanguisuga.anchors.txt"
ru_sa_anchor_file <- "rubida2sanguisuga.anchors.txt"

sp_rh <- "R. prolixus"
sp_ru <- "T. rubida"
sp_sa <- "T. sanguisuga"

TOP_N_RH <- 64
TOP_N_SA <- 181
TOP_N_RU <- 394

BIG_PAL <- c("#B71C1C","#D32F2F","#F44336","#EF9A9A","#880E4F","#C2185B","#E91E63","#F48FB1",
  "#AD1457","#F06292","#FF4081","#4A148C","#6A1B9A","#7B1FA2","#9C27B0","#CE93D8","#7C4DFF",
  "#651FFF","#0D47A1","#1565C0","#1976D2","#1E88E5","#42A5F5","#82B1FF","#2962FF","#004D40",
  "#00695C","#00897B","#26A69A","#80CBC4","#006064","#00838F","#00ACC1","#1B5E20","#2E7D32",
  "#388E3C","#43A047","#66BB6A","#33691E","#558B2F","#689F38","#9CCC65","#CCFF90","#827717",
  "#9E9D24","#AFB42B","#C6A700","#E65100","#EF6C00","#F57C00","#FB8C00","#FFA726","#FF6D00",
  "#FF9100","#3E2723","#4E342E","#5D4037","#6D4C41","#8D6E63","#A1887F","#795548","#212121",
  "#424242","#616161","#757575","#9E9E9E","#37474F","#546E7A","#607D8B")

pal_n <- function(n) {
  if (n <= 0) return(character(0))
  if (n <= length(BIG_PAL)) {
    idx <- round(seq(1, length(BIG_PAL), length.out = n))
    return(BIG_PAL[idx])
  }
  grDevices::colorRampPalette(BIG_PAL)(n)
}

pal_rep <- pal_n

pal_top <- function(n) {
  warm <- c("#B71C1C","#D32F2F","#F44336","#EF9A9A","#880E4F","#C2185B","#E91E63","#F48FB1",
    "#E65100","#EF6C00","#F57C00","#FFA726","#FF4081","#FF6D00","#FF9100","#F9A825","#3E2723",
    "#5D4037","#8D6E63","#795548")
  if (n <= 0) return(character(0))
  if (n <= length(warm)) return(warm[seq_len(n)])
  grDevices::colorRampPalette(warm)(n)
}

pal_bot <- function(n) {
  # Cool subset for bottom bar: blues, teals, greens, purples
  cool <- c("#0D47A1","#1565C0","#1976D2","#1E88E5","#42A5F5","#2962FF","#82B1FF","#004D40",
    "#00695C","#00897B","#26A69A","#006064","#00838F","#00ACC1","#80CBC4","#1B5E20","#2E7D32",
    "#388E3C","#66BB6A","#33691E","#558B2F","#9CCC65","#4A148C","#6A1B9A","#9C27B0","#CE93D8",
    "#37474F","#546E7A","#607D8B","#9E9E9E")
  if (n <= 0) return(character(0))
  if (n <= length(cool)) return(cool[seq_len(n)])
  grDevices::colorRampPalette(cool)(n)
}

top_chroms <- function(df, grp_col, end_col, n) {
  df %>%
    group_by(.data[[grp_col]]) %>%
    summarise(size = max(.data[[end_col]]), .groups = "drop") %>%
    arrange(desc(size)) %>%
    slice_head(n = n) %>%
    pull(.data[[grp_col]])
}

build_layout <- function(df, grp_col, end_col, keep_chroms) {
  df %>%
    filter(.data[[grp_col]] %in% keep_chroms) %>%
    group_by(.data[[grp_col]]) %>%
    summarise(size = max(.data[[end_col]]), .groups = "drop") %>%
    rename(chrom = 1) %>%
    mutate(chrom = factor(chrom, levels = keep_chroms)) %>%
    arrange(chrom) %>%
    mutate(
      offset = lag(cumsum(size), default = 0),
      mid    = offset + size / 2,
      label  = as.character(chrom)          # full actual name
    )
}

order_by_partner <- function(all_blocks, synt_blocks,
                             focal_col, focal_start_col, focal_end_col,
                             partner_col, partner_start_col, partner_end_col,
                             partner_layout, keep_n) {

  focal_sizes <- all_blocks %>%
    group_by(chrom = .data[[focal_col]]) %>%
    summarise(size = max(.data[[focal_end_col]]), .groups = "drop") %>%
    arrange(desc(size)) %>%
    slice_head(n = keep_n)
  
  focal_keep <- focal_sizes$chrom
  
  stats <- synt_blocks %>%
    filter(.data[[focal_col]] %in% focal_keep) %>%
    mutate(
      block_len   = .data[[partner_end_col]] - .data[[partner_start_col]],
      partner_mid = (.data[[partner_start_col]] + .data[[partner_end_col]]) / 2,
      focal_mid   = (.data[[focal_start_col]]   + .data[[focal_end_col]])   / 2
    ) %>%
    group_by(focal = .data[[focal_col]], partner = .data[[partner_col]]) %>%
    summarise(
      total_hits   = sum(hits),
      weighted_pos = sum(partner_mid * block_len) / sum(block_len),
      orientation  = ifelse(
        n() >= 3,
        cor(partner_mid, focal_mid, method = "spearman"),
        ifelse(
          n() == 2,
          sign(focal_mid[which.max(partner_mid)] -
                 focal_mid[which.min(partner_mid)]),
          NA_real_
        )
      ),
      .groups = "drop"
    )
  
  dominant <- stats %>%
    group_by(focal) %>%
    slice_max(total_hits, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(flipped = !is.na(orientation) & orientation < 0)
  
  partner_pos <- partner_layout %>%
    mutate(chrom = as.character(chrom)) %>%
    select(chrom, partner_offset = offset)
  
  ordered_df <- dominant %>%
    left_join(partner_pos, by = c("partner" = "chrom")) %>%
    mutate(global_pos = coalesce(partner_offset, 0) + weighted_pos) %>%
    arrange(global_pos) %>%
    select(focal, flipped)

  missing_chroms <- focal_sizes %>%
    filter(!chrom %in% ordered_df$focal) %>%
    arrange(desc(size)) %>%
    pull(chrom)
  
  missing <- data.frame(
    focal   = missing_chroms,
    flipped = rep(FALSE, length(missing_chroms))
  )
  
  dplyr::bind_rows(ordered_df, missing)
}

rh_ru <- read.table(rh_ru_block_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, comment.char = "")
rh_sa <- read.table(rh_sa_block_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, comment.char = "")
ru_sa <- read.table(ru_sa_block_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, comment.char = "")

for (df_name in c("rh_ru", "rh_sa", "ru_sa")) {
  df <- get(df_name)
  names(df)[names(df) == "X.hits"]  <- "hits"
  names(df)[names(df) == "X.gene1"] <- "gene1"
  names(df)[names(df) == "X.gene2"] <- "gene2"
  assign(df_name, df)
}

anc_rh_sa <- read.table(rh_sa_anchor_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
anc_ru_sa <- read.table(ru_sa_anchor_file,  header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

parse_anchor_chroms <- function(anc) {
  anc %>%
    mutate(
      grp1 = sub("^([^.]+)\\..*",            "\\1", block),
      grp2 = sub("^[^.]+\\.([^.]+)\\.\\d+$", "\\1", block)
    )
}
anc_rh_sa <- parse_anchor_chroms(anc_rh_sa)
anc_ru_sa <- parse_anchor_chroms(anc_ru_sa)

.sa_all  <- sort(unique(c(rh_sa$grp2, ru_sa$grp2)))
.ru_all  <- sort(unique(c(rh_ru$grp2, ru_sa$grp1)))
.rh_all  <- sort(unique(c(rh_sa$grp1, rh_ru$grp1)))

GLOBAL_SA_COLS <- setNames(pal_n(length(.sa_all)),   .sa_all)
GLOBAL_RU_COLS <- setNames(pal_n(length(.ru_all)),   .ru_all)
GLOBAL_RH_COLS <- setNames(pal_top(length(.rh_all)), .rh_all)


make_ribbon_plot <- function(blocks, sp1_label, sp2_label,
                             top_n1, top_n2,
                             mapping_file = NULL) {
  
  top1 <- top_chroms(blocks, "grp1", "end1", top_n1)
  lay1 <- build_layout(blocks, "grp1", "end1", top1)
  lay1$flipped <- FALSE

  top2_info <- order_by_partner(
    all_blocks        = blocks,
    synt_blocks       = blocks %>% filter(grp1 %in% top1),
    focal_col         = "grp2", focal_start_col   = "start2", focal_end_col   = "end2",
    partner_col       = "grp1", partner_start_col = "start1", partner_end_col = "end1",
    partner_layout    = lay1,   keep_n            = top_n2
  )
  top2 <- top2_info$focal
  lay2 <- build_layout(blocks, "grp2", "end2", top2) %>%
    left_join(top2_info %>% rename(chrom = focal), by = "chrom")

  lay1 <- lay1 %>% mutate(num = seq_len(n()), label = as.character(num))
  lay2 <- lay2 %>% mutate(num = seq_len(n()), label = as.character(num))
  if (!is.null(mapping_file)) {
    map_out <- rbind(
      data.frame(number = lay1$num, name = as.character(lay1$chrom), species = sp1_label),
      data.frame(number = lay2$num, name = as.character(lay2$chrom), species = sp2_label)
    )
    tryCatch(
      { write.table(map_out, mapping_file, sep = "\t", row.names = FALSE, quote = FALSE)
        message("Mapping saved: ", mapping_file) },
      error = function(e) warning("Could not write mapping file '", mapping_file,
                                  "': ", conditionMessage(e),
                                  "\nTip: set mapping_file to a full path you have write access to.")
    )
  }

  lay1_total <- sum(lay1$size)
  lay2_total <- sum(lay2$size)
  sc2 <- function(x) x / lay2_total * lay1_total

  lay2_sc <- lay2 %>%
    mutate(offset_sc = sc2(offset),
           size_sc   = sc2(size),
           mid_sc    = offset_sc + size_sc / 2)

  col_by_top <- nrow(lay1) <= nrow(lay2)
  
  b <- blocks %>%
    filter(grp1 %in% top1, grp2 %in% top2) %>%
    left_join(lay1 %>% select(chrom, off1 = offset, size1 = size, flip1 = flipped),
              by = c("grp1" = "chrom")) %>%
    left_join(lay2 %>% select(chrom, off2 = offset, size2 = size, flip2 = flipped),
              by = c("grp2" = "chrom")) %>%
    mutate(
      s1  = ifelse(flip1, size1 - end1,   start1),
      e1  = ifelse(flip1, size1 - start1, end1),
      s2  = ifelse(flip2, size2 - end2,   start2),
      e2  = ifelse(flip2, size2 - start2, end2),
      x1s = s1 + off1,       x1e = e1 + off1,
      x2s = sc2(s2 + off2),  x2e = sc2(e2 + off2),
      col_id   = if (col_by_top) grp1 else grp2,
      block_id = paste(grp1, grp2, block, sep = "_")
    )
  
  n_pts <- 30
  t_seq <- seq(0, 1, length.out = n_pts)
  
  bezier_ribbon <- function(xs, xe, ys, ye, block_id, col_id) {
    left_x  <- (1-t_seq)^3*xs + 3*(1-t_seq)^2*t_seq*xs +
      3*(1-t_seq)*t_seq^2*ys + t_seq^3*ys
    left_y  <- (1-t_seq)^3*1  + 3*(1-t_seq)^2*t_seq*0.65 +
      3*(1-t_seq)*t_seq^2*0.35 + t_seq^3*0
    right_x <- (1-t_seq)^3*xe + 3*(1-t_seq)^2*t_seq*xe +
      3*(1-t_seq)*t_seq^2*ye + t_seq^3*ye
    right_y <- left_y
    data.frame(x = c(left_x, rev(right_x)), y = c(left_y, rev(right_y)),
               block_id = block_id, col_id = col_id)
  }
  
  ribbon_df <- b %>%
    rowwise() %>%
    do(bezier_ribbon(.$x1s, .$x1e, .$x2s, .$x2e, .$block_id, .$col_id)) %>%
    ungroup()
  
  bar1 <- lay1 %>%
    mutate(xmin = offset, xmax = offset + size, ymin = 1.01, ymax = 1.06)
  bar2 <- lay2_sc %>%
    mutate(xmin = offset_sc, xmax = offset_sc + size_sc, ymin = -0.06, ymax = -0.01)

  bar_fill_col <- function(bar_df, sp_label) {
    bar_df %>% mutate(fill_col = dplyr::case_when(
      sp_label == sp_sa ~ as.character(chrom),
      sp_label == sp_rh ~ "#2C4A6E",          
      TRUE              ~ "#7B3F00"           
    ))
  }
  bar1 <- bar_fill_col(bar1, sp1_label)
  bar2 <- bar_fill_col(bar2, sp2_label)
  
  all_col_vals <- c(GLOBAL_RH_COLS, GLOBAL_SA_COLS, GLOBAL_RU_COLS,
                    c("#2C4A6E" = "#2C4A6E", "#7B3F00" = "#7B3F00"))
  
  lsz <- max(0.7, min(2.0, 55 / max(nrow(lay1), nrow(lay2))))
  
  ggplot() +
    geom_polygon(data = ribbon_df,
                 aes(x = x, y = y, group = block_id, fill = col_id),
                 alpha = 0.5, colour = NA) +
    geom_rect(data = bar1,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
              colour = "white", linewidth = 0.25) +
    geom_rect(data = bar2,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
              colour = "white", linewidth = 0.25) +
    geom_text(data = lay1,
              aes(x = mid, y = 1.08, label = label),
              size = lsz, angle = 0, hjust = 0.5, vjust = 0, colour = "grey10") +
    geom_text(data = lay2_sc,
              aes(x = mid_sc, y = -0.08, label = label),
              size = lsz, angle = 0, hjust = 0.5, vjust = 1, colour = "grey10") +
    annotate("text", x = -lay1_total * 0.008, y = 1.035,
             label = sp1_label, hjust = 1, size = 3.0, fontface = "italic") +
    annotate("text", x = -lay1_total * 0.008, y = -0.035,
             label = sp2_label, hjust = 1, size = 3.0, fontface = "italic") +
    scale_fill_manual(values = all_col_vals) +
    scale_x_continuous(limits = c(0, lay1_total), expand = c(0.10, 0)) +
    scale_y_continuous(limits = c(-0.65, 1.65)) +
    labs(title    = "Ribbon / Link Plot",
         subtitle = paste(sp1_label, "(top) vs", sp2_label, "(bottom)")) +
    theme_void(base_size = 9) +
    theme(
      legend.position = "none",
      plot.title    = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5, face = "italic"),
      plot.margin   = margin(t = 10, r = 5, b = 10, l = 80, unit = "pt"),
      plot.clip     = "off"
    )
}

p2a <- make_ribbon_plot(
  blocks       = rh_sa,
  sp1_label    = sp_rh,
  sp2_label    = sp_sa,
  top_n1       = TOP_N_RH,
  top_n2       = TOP_N_SA,
  mapping_file = "mapping_p2a_rh_sa.tsv"
)

p2b <- make_ribbon_plot(
  blocks       = rh_ru,
  sp1_label    = sp_rh,
  sp2_label    = sp_ru,
  top_n1       = TOP_N_RH,
  top_n2       = TOP_N_RU,
  mapping_file = "mapping_p2b_rh_ru.tsv"
)

ru_sa_swapped <- ru_sa %>%
  rename(grp1 = grp2, grp2 = grp1,
         start1 = start2, end1 = end2,
         start2 = start1, end2 = end1)
p2c <- make_ribbon_plot(
  blocks       = ru_sa_swapped,
  sp1_label    = sp_sa,
  sp2_label    = sp_ru,
  top_n1       = TOP_N_SA,
  top_n2       = TOP_N_RU,
  mapping_file = "mapping_p2c_sa_ru.tsv"
)

p2 <- p2a

make_multispecies_plot <- function(blocks_top_mid,
                                   blocks_mid_bot,
                                   sp_top, sp_mid, sp_bot,
                                   top_n_top, top_n_mid, top_n_bot,
                                   mapping_file_p3 = NULL) {

  top_mid <- top_chroms(blocks_top_mid, "grp2", "end2", top_n_mid)
  lay_mid <- build_layout(blocks_top_mid, "grp2", "end2", top_mid)
  lay_mid$flipped <- FALSE
  
  top_top_info <- order_by_partner(
    all_blocks        = blocks_top_mid,
    synt_blocks       = blocks_top_mid %>% filter(grp2 %in% top_mid),
    focal_col         = "grp1", focal_start_col   = "start1", focal_end_col   = "end1",
    partner_col       = "grp2", partner_start_col = "start2", partner_end_col = "end2",
    partner_layout    = lay_mid, keep_n           = top_n_top
  )
  top_top <- top_top_info$focal
  lay_top <- build_layout(blocks_top_mid, "grp1", "end1", top_top) %>%
    left_join(top_top_info %>% rename(chrom = focal), by = "chrom")
  
  top_bot_info <- order_by_partner(
    all_blocks        = blocks_mid_bot,
    synt_blocks       = blocks_mid_bot %>% filter(grp2 %in% top_mid),
    focal_col         = "grp1", focal_start_col   = "start1", focal_end_col   = "end1",
    partner_col       = "grp2", partner_start_col = "start2", partner_end_col = "end2",
    partner_layout    = lay_mid, keep_n           = top_n_bot
  )
  top_bot <- top_bot_info$focal
  lay_bot <- build_layout(blocks_mid_bot, "grp1", "end1", top_bot) %>%
    left_join(top_bot_info %>% rename(chrom = focal), by = "chrom")

  lay_mid <- lay_mid %>% mutate(num = seq_len(n()), label = as.character(num))
  lay_top <- lay_top %>% mutate(num = seq_len(n()), label = as.character(num))
  lay_bot <- lay_bot %>% mutate(num = seq_len(n()), label = as.character(num))
  if (!is.null(mapping_file_p3)) {
    map_out <- rbind(
      data.frame(number = lay_top$num, name = as.character(lay_top$chrom), species = sp_top),
      data.frame(number = lay_mid$num, name = as.character(lay_mid$chrom), species = sp_mid),
      data.frame(number = lay_bot$num, name = as.character(lay_bot$chrom), species = sp_bot)
    )
    tryCatch(
      { write.table(map_out, mapping_file_p3, sep = "\t", row.names = FALSE, quote = FALSE)
        message("Mapping saved: ", mapping_file_p3) },
      error = function(e) warning("Could not write mapping file '", mapping_file_p3,
                                  "': ", conditionMessage(e),
                                  "\nTip: set mapping_file_p3 to a full path you have write access to.")
    )
  }
  
  mid_total <- sum(lay_mid$size)
  top_total <- sum(lay_top$size)
  bot_total <- sum(lay_bot$size)
  
  scale_to_mid <- function(x, src) x / src * mid_total
  
  n_pts <- 25
  t_seq <- seq(0, 1, length.out = n_pts)
  
  make_ribbons <- function(b, lay1, lay2, off1_fn, off2_fn,
                           y_top, y_bot, ctrl_top, ctrl_bot,
                           col_key = "grp2") {
    b2 <- b %>%
      left_join(lay1 %>% select(chrom, off1 = offset, size1 = size, flip1 = flipped),
                by = c("grp1" = "chrom")) %>%
      left_join(lay2 %>% select(chrom, off2 = offset, size2 = size, flip2 = flipped),
                by = c("grp2" = "chrom")) %>%
      filter(!is.na(off1), !is.na(off2)) %>%
      mutate(
        s1 = ifelse(flip1, size1 - end1,   start1),
        e1 = ifelse(flip1, size1 - start1, end1),
        s2 = ifelse(flip2, size2 - end2,   start2),
        e2 = ifelse(flip2, size2 - start2, end2),
        x1s      = off1_fn(s1 + off1),
        x1e      = off1_fn(e1 + off1),
        x2s      = off2_fn(s2 + off2),
        x2e      = off2_fn(e2 + off2),
        block_id = paste(grp1, grp2, block, sep = "_"),
        col_id   = .data[[col_key]]
      )
    
    b2 %>% rowwise() %>%
      do({
        lx <- (1-t_seq)^3*.$x1s + 3*(1-t_seq)^2*t_seq*.$x1s +
          3*(1-t_seq)*t_seq^2*.$x2s + t_seq^3*.$x2s
        rx <- (1-t_seq)^3*.$x1e + 3*(1-t_seq)^2*t_seq*.$x1e +
          3*(1-t_seq)*t_seq^2*.$x2e + t_seq^3*.$x2e
        ly <- (1-t_seq)^3*y_top + 3*(1-t_seq)^2*t_seq*ctrl_top +
          3*(1-t_seq)*t_seq^2*ctrl_bot + t_seq^3*y_bot
        data.frame(x = c(lx, rev(rx)), y = c(ly, rev(ly)),
                   block_id = .$block_id, col_id = .$col_id)
      }) %>% ungroup()
  }

  rh_sa_filt <- blocks_top_mid %>% filter(grp1 %in% top_top, grp2 %in% top_mid)
  ribbons_tm <- make_ribbons(
    rh_sa_filt, lay_top, lay_mid,
    off1_fn   = function(x) scale_to_mid(x, top_total),
    off2_fn   = function(x) x,
    y_top = 2, y_bot = 1, ctrl_top = 1.6, ctrl_bot = 1.4
  )

  ru_sa_filt <- blocks_mid_bot %>% filter(grp1 %in% top_bot, grp2 %in% top_mid)
  ru_sa_swap <- ru_sa_filt %>%
    rename(grp1 = grp2, grp2 = grp1,
           start1 = start2, end1 = end2,
           start2 = start1, end2 = end1)
  ribbons_mb <- make_ribbons(
    ru_sa_swap, lay_mid, lay_bot,
    off1_fn   = function(x) x,
    off2_fn   = function(x) scale_to_mid(x, bot_total),
    y_top = 1, y_bot = 0, ctrl_top = 0.6, ctrl_bot = 0.4,
    col_key   = "grp1"
  )

  sa_col_map <- GLOBAL_SA_COLS

  bar_top <- lay_top %>%
    mutate(xmin = scale_to_mid(offset,        top_total),
           xmax = scale_to_mid(offset + size, top_total),
           ymin = 2.01, ymax = 2.06)
  bar_mid <- lay_mid %>%
    mutate(xmin = offset, xmax = offset + size,
           ymin = 0.97,   ymax = 1.03)
  bar_bot <- lay_bot %>%
    mutate(xmin = scale_to_mid(offset,        bot_total),
           xmax = scale_to_mid(offset + size, bot_total),
           ymin = -0.06, ymax = -0.01)
  
  ggplot() +
    geom_polygon(data = ribbons_tm,
                 aes(x = x, y = y, group = block_id, fill = col_id),
                 alpha = 0.4, colour = NA) +
    geom_polygon(data = ribbons_mb,
                 aes(x = x, y = y, group = block_id, fill = col_id),
                 alpha = 0.4, colour = NA) +
    geom_rect(data = bar_top,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#2C4A6E", colour = "white", linewidth = 0.2) +
    geom_rect(data = bar_mid,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = chrom),
              colour = "white", linewidth = 0.2, alpha = 0.9) +
    geom_rect(data = bar_bot,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#7B3F00", colour = "white", linewidth = 0.2) +
    # Numeric labels — angle=0, all shown
    geom_text(data = bar_top,
              aes(x = (xmin + xmax) / 2, y = 2.07, label = label),
              size = 1.6, angle = 0, hjust = 0.5, vjust = 0, colour = "grey25") +
    geom_text(data = bar_mid,
              aes(x = (xmin + xmax) / 2, y = 1.04, label = label),
              size = 1.6, angle = 0, hjust = 0.5, vjust = 0, colour = "grey10") +
    geom_text(data = bar_bot,
              aes(x = (xmin + xmax) / 2, y = -0.07, label = label),
              size = 1.6, angle = 0, hjust = 0.5, vjust = 1, colour = "grey25") +
    # Species name labels — italic on the left, at the vertical midpoint of each bar
    annotate("text", x = -mid_total * 0.008, y = 2.035,
             label = sp_top, hjust = 1, size = 3.2, fontface = "italic",
             colour = "#2C4A6E") +
    annotate("text", x = -mid_total * 0.008, y = 1.0,
             label = sp_mid, hjust = 1, size = 3.2, fontface = "italic",
             colour = "#8B1A1A") +
    annotate("text", x = -mid_total * 0.008, y = -0.035,
             label = sp_bot, hjust = 1, size = 3.2, fontface = "italic",
             colour = "#7B3F00") +
    scale_fill_manual(values = sa_col_map) +
    scale_x_continuous(expand = c(0.12, 0)) +
    # Extended y range to accommodate vertical labels
    scale_y_continuous(limits = c(-0.65, 2.65)) +
    labs(title    = "Multi-Species Linear Synteny",
         subtitle = paste(sp_top, "—", sp_mid, "—", sp_bot)) +
    theme_void(base_size = 9) +
    theme(
      legend.position = "none",
      plot.title    = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5, face = "italic"),
      plot.margin   = margin(t = 5, r = 5, b = 5, l = 90, unit = "pt"),
      plot.clip     = "off"
    )
}

p3 <- make_multispecies_plot(
  blocks_top_mid  = rh_sa,
  blocks_mid_bot  = ru_sa,
  sp_top          = sp_rh,
  sp_mid          = sp_sa,
  sp_bot          = sp_ru,
  top_n_top       = TOP_N_RH,
  top_n_mid       = TOP_N_SA,
  top_n_bot       = TOP_N_RU,
  mapping_file_p3 = "mapping_p3_multispecies.tsv"
)

combined <- (p1 | p2a) / p3 +
  plot_annotation(
    title    = "SyMAP Synteny Results",
    subtitle = paste(sp_rh, "·", sp_ru, "·", sp_sa),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40", face = "italic")
    )
  ) +
  plot_layout(heights = c(1, 0.9))

ggsave("symap_synteny.pdf", combined,
       width = 18, height = 13, device = cairo_pdf)
ggsave("symap_synteny.png", combined,
       width = 18, height = 13, dpi = 300, bg = "white")

ribbons_all <- p2a / p2b / p2c +
  plot_annotation(
    title    = "Pairwise Synteny Ribbon Plots",
    subtitle = paste(sp_rh, "·", sp_ru, "·", sp_sa),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40", face = "italic")
    )
  )

ggsave("symap_ribbons_all.pdf", ribbons_all,
       width = 14, height = 14, device = cairo_pdf)
ggsave("symap_ribbons_all.png", ribbons_all,
       width = 14, height = 14, dpi = 300, bg = "white")

ggsave("plot2a_ribbon_rh_sa.pdf",      p2a, width = 12, height = 5,  device = cairo_pdf)
ggsave("plot2b_ribbon_rh_ru.pdf",      p2b, width = 12, height = 5,  device = cairo_pdf)
ggsave("plot2c_ribbon_ru_sa.pdf",      p2c, width = 12, height = 5,  device = cairo_pdf)
ggsave("plot3_multispecies.pdf",       p3,  width = 14, height = 5,  device = cairo_pdf)
