library(ggtree)
library(treeio)
library(ggplot2)
library(tidyverse)

newick <- "(Acyrthosiphon_pisum,(Halyomorpha_halys,((Cimex_lectularius,nesidiocoris_tenuis),(Rhynocoris_fuscipes,((Triatoma_rubida,Triatoma_sanguisuga),Rhodnius_prolixus)))));"

raw <- read.delim("list.csv",
                  header      = TRUE,
                  sep         = ",",
                  strip.white = TRUE)

data <- raw %>%
  mutate(
    Branch   = strsplit(as.character(Branch), ",\\s*"),
    Category = str_trim(COG_Category)
  ) %>%
  unnest(Branch) %>%
  mutate(
    Branch   = str_trim(Branch),
    Category = str_squish(Category),
    Category = tolower (Category)
  )

node_map <- tribble(
  ~table_name, ~ggtree_node,
  "Node8",     "14",          # MRCA of rfusc + triatominae
  "Node4",     "12",          # MRCA of clect + reduviidae
  "Node10",    "15",          # MRCA of rprol + triatoma
  "Node11",    "16",          # MRCA of trubi + tsang
)

data <- data %>%
  left_join(node_map, by = c("Branch" = "table_name")) %>%
  mutate(Branch = if_else(!is.na(ggtree_node), ggtree_node, Branch)) %>%
  select(-ggtree_node)

print(unique(data$Category))

category_colors <- c(
  "metabolism and detoxification"   = "#993C1D",
  "proteases and digestion related" = "#D85A30",
  "sensory and signalling"          = "#1D9E75",
  "structural and regulatory"       = "#534AB7",
  "vesicle trafficking"             = "#378ADD",
  "unknown"                         = "#BA7517"
)

category_labels <- c(
  "metabolism and detoxification"   = "Metabolism & Detoxification Related",
  "proteases and digestion related" = "Proteases & Digestion Related",
  "sensory and signalling"          = "Sensory & Signalling Related",
  "structural and regulatory"       = "Structural & Regulatory Related",
  "vesicle trafficking"             = "Vesicle Trafficking Related",
  "unknown"                         = "Not Classified"
)

cat_order <- names(category_colors)

counts <- data %>%
  count(Branch, Category) %>%
  filter(Category %in% names(category_colors))

print(counts)

tree <- read.tree(text = newick)

p_base <- ggtree(
  tree,
  layout        = "rectangular",
  branch.length = "none",
  color         = "#5F5E5A",
  size          = 0.65
)

node_data <- p_base$data

p_check <- p_base +
  geom_text(
    aes(label = node),
    hjust = -0.3,
    size  = 2.8,
    color = "#D85A30"
  ) +
  geom_tiplab(
    size   = 3,
    offset = 0.1
  )

print(p_check)

cat_offset <- tibble(
  Category = cat_order,
  x_offset = seq(0.15, by = 0.40, length.out = length(cat_order))
)

tip_nodes <- node_data %>%
  filter(isTip == TRUE) %>%
  mutate(node_label = label) %>%
  select(node_label, x, y, node, isTip)

tip_plot <- tip_nodes %>%
  left_join(counts, by = c("node_label" = "Branch")) %>%
  filter(!is.na(Category)) %>%
  left_join(cat_offset, by = "Category") %>%
  mutate(
    x_dot    = x + x_offset,
    Category = factor(Category, levels = cat_order)
  )

internal_nodes <- node_data %>%
  filter(isTip == FALSE) %>%
  mutate(node_label = as.character(node)) %>%
  select(node_label, x, y, node, isTip)

internal_plot <- internal_nodes %>%
  left_join(counts, by = c("node_label" = "Branch")) %>%
  filter(!is.na(Category)) %>%
  left_join(cat_offset, by = "Category") %>%
  mutate(
    x_dot    = x + x_offset,
    Category = factor(Category, levels = cat_order)
  )

plot_data <- bind_rows(tip_plot, internal_plot)

show_node_labels <- TRUE   # set FALSE when done verifying

mrca_descriptions <- tribble(
  ~table_name, ~description,
  "Node8",     "MRCA rfusc + triatominae",
  "Node4",     "MRCA clect + reduviidae",
  "Node10",    "MRCA rprol + triatoma",
  "Node11",    "MRCA trubi + tsang"
)

internal_label_data <- node_data %>%
  filter(isTip == FALSE) %>%
  mutate(node_label = as.character(node)) %>%
  left_join(
    node_map %>% mutate(ggtree_node = as.character(ggtree_node)),
    by = c("node_label" = "ggtree_node")
  ) %>%
  left_join(mrca_descriptions, by = "table_name") %>%
  mutate(
    display_label = case_when(
      !is.na(table_name) ~ paste0(table_name, "\n", description),
      TRUE               ~ node_label
    )
  )

x_max <- max(node_data$x, na.rm = TRUE)

p_final <- p_base +
  
  # Tip labels
  geom_tiplab(
    aes(label = label),
    size     = 3.2,
    fontface = "italic",
    color    = "#2C2C2A",
    offset   = 0.12
  ) +
  
  {
    if (show_node_labels)
      geom_label(
        data          = internal_label_data %>% filter(!is.na(table_name)),
        aes(
          x     = x,
          y     = y,
          label = display_label
        ),
        size          = 1.9,
        color         = "#534AB7",
        fill          = "white",
        label.size    = 0.2,
        label.padding = unit(0.10, "lines"),
        vjust         = 1.8,
        inherit.aes   = FALSE
      )
  } +
  
  geom_point(
    data = plot_data,
    aes(
      x    = x_dot,
      y    = y,
      fill = Category,
      size = n
    ),
    shape       = 21,
    color       = "white",
    stroke      = 0.4,
    alpha       = 0.92,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = category_colors,
    labels = category_labels,
    name   = "Functional category",
    guide  = guide_legend(
      title.position = "top",
      override.aes   = list(size = 4),
      keywidth       = unit(0.5, "cm"),
      keyheight      = unit(0.5, "cm")
    )
  ) +
  
  scale_size_area(
    max_size = 10,
    breaks   = c(1, 3, 5, 10),
    labels   = c("1", "3", "5", "10"),
    name     = "No. of OGs",
    guide    = guide_legend(
      title.position = "top",
      override.aes   = list(
        fill  = "#888780",
        color = "white",
        shape = 21
      )
    )
  ) +
  
  xlim(NA, x_max + 5.5) +
  
  theme_tree2() +
  theme(
    legend.position   = "right",
    legend.box        = "vertical",
    legend.spacing.y  = unit(0.4, "cm"),
    legend.background = element_rect(
      fill      = "white",
      color     = "#D3D1C7",
      linewidth = 0.4
    ),
    legend.title = element_text(
      size  = 8.5,
      face  = "bold",
      color = "#2C2C2A"
    ),
    legend.text = element_text(
      size  = 8,
      color = "#444441"
    ),
    legend.key.size  = unit(0.45, "cm"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(15, 10, 15, 15)
  )

ggsave(
  "positive_selection_tree.pdf",
  plot   = p_final,
  width  = 12,
  height = 8,
  dpi    = 600,
  device = cairo_pdf
)

ggsave(
  "positive_selection_tree.png",
  plot   = p_final,
  width  = 11,
  height = 7,
  dpi    = 300
)
