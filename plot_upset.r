rm(list = ls())
library(UpSetR)
library(ggplot2)

data <- read.table("Orthogroups.GeneCount.tsv", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)

df <- data[, c("acryp", "clect", "hhaly", "ntenu", "rfusc", "rprol", "trubi", "tsang")]
df[df > 0] <- 1L
df[df == 0] <- 0L
df <- data.frame(df, check.names = FALSE)

api <- "Acyrthosiphon_pisum"
cim <- "Cimex_lectularius"
hha <- "Halyomorpha_halys"
nte <- "Nesidiocoris_tenuis"
rhy <- "Rhynocoris_fuscipes"
rho <- "Rhodnius_prolixus"
rub <- "Triatoma_rubida"
san <- "Triatoma_sanguisuga"

colnames(df) <- c(api, cim, hha, nte, rhy, rho, rub, san)


targets <- list(
  list(san),
  list(rub),
  list(rho),
  list(rhy),
  list(cim),

  list(api, hha, nte, cim, rhy, rho, rub, san),
  list(rhy, rho, rub, san), 
  list(cim, rho, rub, san), 
  list(rho, rub, san)       
)

cairo_pdf("Final_upset.pdf", width = 18, height = 10, fallback_resolution = 600)

upset(
  df,
  intersections = targets,
  sets          = c(san, rub, rho, rhy, cim, nte, hha, api),  # bottom → top
  keep.order    = TRUE,
  order.by      = c("freq", "degree"),
  decreasing    = c(FALSE, TRUE),
  mb.ratio      = c(0.65, 0.35),
  sets.bar.color = "#4A90E2",
  main.bar.color = "gray35",
  point.size    = 4.5,
  line.size     = 1.2,
  shade.color   = "gray85",
  shade.alpha   = 0.3,
  query.legend  = "top",
  text.scale    = c(1.5, 1.2, 1.2, 1, 1.5, 1.2),
      
  queries = list(
    list(
      query      = intersects,
      params     = list(api, hha, nte, cim, rhy, rho, rub, san),
      color      = "purple4",
      active     = TRUE,
      query.name = "Shared by All"
    ),
        
    list(
      query      = intersects,
      params     = list(rhy, rho, rub, san),
      color      = "darkblue",
      active     = TRUE,
      query.name = "Rhynocoris + Triatomines"
    ),
        
    list(
      query      = intersects,
      params     = list(cim, rho, rub, san),
      color      = "firebrick",
      active     = TRUE,
      query.name = "Cimex + Triatomines"
    ),
        
    list(
      query      = intersects,
      params     = list(rho, rub, san),
      color      = "darkgreen",
      active     = TRUE,
      query.name = "Triatomines"
    ),
        

    list(
      query      = intersects,
      params     = list(cim),
      color      = "firebrick",
      active     = TRUE,
      query.name = "Unique: Cimex"
    ),
        
    list(
      query      = intersects,
      params     = list(rhy),
      color      = "darkblue",
      active     = TRUE,
      query.name = "Unique: Rhynocoris"
    ),
        
    list(
      query      = intersects,
      params     = list(rho),
      color      = "darkgreen",
      active     = TRUE,
      query.name = "Unique: Rhodnius"
    ),
        
    list(
      query      = intersects,
      params     = list(rub),
      color      = "darkorange",
      active     = TRUE,
      query.name = "Unique: T. rubida"
    ),
        
    list(
      query      = intersects,
      params     = list(san),
      color      = "goldenrod3",
      active     = TRUE,
      query.name = "Unique: T. sanguisuga"
    )
  )
)

dev.off()


rownames(df) <- data[, 1]

get_orthogroups <- function(df, present_spp, all_spp) {
  absent_spp <- setdiff(all_spp, present_spp)
  present_filter <- rowSums(df[, present_spp, drop = FALSE] == 1) == length(present_spp)

  if (length(absent_spp) > 0) {
    absent_filter <- rowSums(df[, absent_spp, drop = FALSE] == 0) == length(absent_spp)
    return(rownames(df)[present_filter & absent_filter])
  } else {
    return(rownames(df)[present_filter])
  }
}

all_spp <- c(api, cim, hha, nte, rhy, rho, rub, san)

og_shared_all      <- get_orthogroups(df, all_spp,                all_spp)
og_rhy_triatomines <- get_orthogroups(df, c(rhy, rho, rub, san),  all_spp)
og_cim_triatomines <- get_orthogroups(df, c(cim, rho, rub, san),  all_spp)
og_triatomines     <- get_orthogroups(df, c(rho, rub, san),       all_spp)
og_unique_cim      <- get_orthogroups(df, cim,                    all_spp)
og_unique_rhy      <- get_orthogroups(df, rhy,                    all_spp)
og_unique_rho      <- get_orthogroups(df, rho,                    all_spp)
og_unique_rub      <- get_orthogroups(df, rub,                    all_spp)
og_unique_san      <- get_orthogroups(df, san,                    all_spp)

cat("Shared by all:           ", length(og_shared_all),      "\n")
cat("Rhynocoris + Triatomines:", length(og_rhy_triatomines), "\n")
cat("Cimex + Triatomines:     ", length(og_cim_triatomines), "\n")
cat("Triatomines only:        ", length(og_triatomines),     "\n")
cat("Unique Cimex:            ", length(og_unique_cim),      "\n")
cat("Unique Rhynocoris:       ", length(og_unique_rhy),      "\n")
cat("Unique Rhodnius:         ", length(og_unique_rho),      "\n")
cat("Unique T. rubida:        ", length(og_unique_rub),      "\n")
cat("Unique T. sanguisuga:    ", length(og_unique_san),      "\n")

og_table <- rbind(
  data.frame(Category = "Shared_by_all",           Orthogroup = og_shared_all),
  data.frame(Category = "Rhynocoris_Triatomines",  Orthogroup = og_rhy_triatomines),
  data.frame(Category = "Cimex_Triatomines",       Orthogroup = og_cim_triatomines),
  data.frame(Category = "Triatomines_only",        Orthogroup = og_triatomines),
  data.frame(Category = "Unique_Cimex",            Orthogroup = og_unique_cim),
  data.frame(Category = "Unique_Rhynocoris",       Orthogroup = og_unique_rhy),
  data.frame(Category = "Unique_Rhodnius",         Orthogroup = og_unique_rho),
  data.frame(Category = "Unique_T_rubida",         Orthogroup = og_unique_rub),
  data.frame(Category = "Unique_T_sanguisuga",     Orthogroup = og_unique_san)
)

write.table(
  og_table,
  "orthogroup_categories.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
