require(dplyr)
require(tidyr)
require(vegan)
require(ggbiplot)
require(dplyr)
require(tidyr)
require(ggplot2)
require(tidyverse)
require(readr)
require(scales)
require(car)
require(ARTool)
require(clinfun) 
require(emmeans)
require(multcomp)       
require(multcompView)
require(purrr)
require(tibble)
require(ggrepel)



setwd("insert path here")

# Hypothesis 1
####### PCA Legacy Vs Control ########

# 0) Optional: keep only Legacy phase if present
df_legacy <- if ("Phase" %in% names(df)) {
  df %>% filter(Phase == "Legacy")
} else {
  df
}

# 1) Abundance matrix: counts per Mesocosm × Treatment × Organism
abundance_matrix <- df_legacy %>%
  count(Mesocosm, Treatment, Organism, name = "Abundance") %>%
  pivot_wider(names_from = Organism,
              values_from = Abundance,
              values_fill = 0)

# Ensure Treatment factor order is stable/explicit
abundance_matrix <- abundance_matrix %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "U", "R", "UR", "4UR")))

# 2) PCA on Hellinger-transformed species matrix
species <- abundance_matrix %>% select(-Mesocosm, -Treatment)
species_hell <- decostand(as.matrix(species), method = "hellinger")

pca <- prcomp(species_hell, center = TRUE, scale. = FALSE)

# Percent variance explained
pc_var <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
pc12_lab <- paste0("PC1 (", round(pc_var[1], 1), "%) · PC2 (", round(pc_var[2], 1), "%)")

# 3) Plot PC1–PC2 by Treatment
ggbiplot(pca,
              groups    = abundance_matrix$Treatment,
              ellipse   = TRUE,
              circle    = TRUE,
              var.axes  = FALSE,
              obs.scale = 1,
              var.scale = 1) +
  scale_color_manual(
    name   = "Treatment",
    values = c("Control" = "black",
               "U"       = "purple",
               "R"       = "yellow",
               "UR"      = "lightgreen",
               "4UR"     = "pink"),
    drop = FALSE
  ) +
  scale_fill_manual(
    name   = "Treatment",
    values = c("Control" = "black",
               "U"       = "purple",
               "R"       = "yellow",
               "UR"      = "lightgreen",
               "4UR"     = "pink"),
    drop = FALSE
  ) +
  labs(title = "PCA of Legacy Community Composition",
       subtitle = pc12_lab,
       x = "PC1", y = "PC2") +
  theme_minimal()

# Euclidean distance on Hellinger-transformed abundances
dist <- vegdist(species_hell, method = "euclidean")

# PERMANOVA: test community differences by Treatment
set.seed(1)
adonis2(dist ~ Treatment, data = abundance_matrix, permutations = 9999)
####### PCA 4c Vs 4UR #########

setwd("INSERT_PATHWAY_HERE")

df <- read.csv("4cv4ur.csv") %>% filter(Dosing == "Legacy")

ab <- df %>%
  count(Mesocosm, Treatment, Organism, name = "Abundance") %>%
  pivot_wider(names_from = Organism, values_from = Abundance, values_fill = 0) %>%
  mutate(Treatment = factor(Treatment, levels = c("4c","4UR")))

sp <- ab %>% select(-Mesocosm, -Treatment) %>% as.matrix()
sp_hel <- decostand(sp, method = "hellinger")
pca <- prcomp(sp_hel)

p <- ggbiplot(pca,
              groups = ab$Treatment,
              ellipse = TRUE, circle = TRUE,
              var.axes = FALSE) +
  scale_color_manual(values = c("4c"="#FFD580","4UR"="pink")) +
  scale_fill_manual(values = c("4c"="#FFD580","4UR"="pink")) +
  theme_minimal()

ggsave("PCA_Legacy_4v4.png", p, width = 8, height = 6, dpi = 600, bg = "white")

adonis2(vegdist(sp_hel, "euclidean") ~ Treatment, data = ab, permutations = 9999)
####### PCA All 2025 Samples #########
setwd("INSERT_PATHWAY_HERE")

df <- read.csv("final_data_legacy.csv")
df <- df[-c(26852:26863), ]
df <- df[order(df$Mesocosm), ]

ab <- df %>% 
  count(Mesocosm, Treatment, Organism, name = "Abundance") %>%
  pivot_wider(names_from = Organism, values_from = Abundance, values_fill = 0)

sp <- ab %>% select(-Mesocosm, -Treatment) %>% as.matrix()
sp_hel <- decostand(sp, method = "hellinger")
pca <- prcomp(sp_hel)

p <- ggbiplot(pca, groups = ab$Treatment, ellipse = FALSE, circle = FALSE, var.axes = FALSE) +
  theme_minimal()
ggsave("PCA_AllData_2025.png", p, width = 8, height = 6, dpi = 600, bg = "white")

set.seed(1)
adonis2(vegdist(sp_hel, "euclidean") ~ Treatment, data = ab, permutations = 9999)

####### RDA Warming ##########
setwd("INSERT_PATHWAY_HERE")

df <- read.csv("temp_spring_both.csv")
df <- df[order(df$Mesocosm), ]

ab <- df %>%
  count(Mesocosm, Treatment, Organism, Season, Year, name = "Abundance") %>%
  pivot_wider(names_from = Organism, values_from = Abundance, values_fill = 0)

sp <- ab %>% select(-Mesocosm, -Treatment, -Season, -Year) %>% as.matrix()
sp_hel <- decostand(sp, method = "hellinger")

ab$Mesocosm <- factor(ab$Mesocosm)
ab$Year     <- factor(ab$Year, levels = c("2019","2025"))

mod <- rda(sp_hel ~ Treatment * Year + Condition(Mesocosm), data = ab)
anova(mod, by = "terms", permutations = 1000)

site <- scores(mod, display = "sites", scaling = 2)
spcs <- scores(mod, display = "species", scaling = 2)
top6 <- order(rowSums(spcs^2), decreasing = TRUE)[1:6]
sp_top <- spcs[top6, ]

eigp <- 100 * mod$CCA$eig / sum(mod$CCA$eig)
xlab <- paste0("RDA1 (", round(eigp[1], 1), "%)")
ylab <- paste0("RDA2 (", round(eigp[2], 1), "%)")

pal <- c("black", "#00008B", "#4169E1", "#6495ED",
         "#FF8C00", "#FF4500", "#DC143C", "#8B0000")
tnum <- as.integer(as.character(ab$Treatment))
cols <- pal[tnum + 1L]
pchv <- ifelse(ab$Year == "2019", 16, 17)

png("RDA_warming.png", width = 8, height = 6, units = "in", res = 600, bg = "white")
plot(site[,1], site[,2], type = "n", xlab = xlab, ylab = ylab, xlim = c(-1,1), ylim = c(-1,1))
abline(h = 0, v = 0, lty = 3)
points(site[,1], site[,2], pch = pchv, col = cols, cex = 1.2)
arrows(0, 0, sp_top[,1], sp_top[,2], length = 0.07)
text(sp_top[,1]*1.05, sp_top[,2]*1.05, labels = rownames(sp_top), cex = 0.8)
legend("topright", legend = c("2019", "2025"), pch = c(16,17), bty = "n", title = "Year")
legend("bottomright", legend = 0:7, col = pal, pch = 15, bty = "n", title = "Treatment")
dev.off()

####### RDA BACI Pesticides #######
setwd("INSERT_PATHWAY_HERE")

ab <- read.csv("abundance_matrix_chem.csv")

sp <- ab %>% select(-Mesocosm, -Dosing, -Treatment) %>% as.matrix()
sp_hel <- decostand(sp, method = "hellinger")

mod <- rda(sp_hel ~ Treatment * Dosing + Condition(Mesocosm), data = ab)
anova(mod, by = "terms", permutations = 1000)

site <- scores(mod, display = "sites", scaling = 2)
spcs <- scores(mod, display = "species", scaling = 2)
top6 <- order(rowSums(spcs^2), decreasing = TRUE)[1:6]
sp_top <- spcs[top6, ]

eigp <- 100 * mod$CCA$eig / sum(mod$CCA$eig)
xlab <- paste0("RDA1 (", round(eigp[1], 1), "%)")
ylab <- paste0("RDA2 (", round(eigp[2], 1), "%)")

ab$Dosing <- factor(ab$Dosing, levels = c("Pre","Post","Legacy"))
cols <- c(Pre="forestgreen", Post="red", Legacy="gold")[ab$Dosing]

tf <- factor(ab$Treatment)
pch_map <- setNames(c(16,15,17,18,8)[seq_along(levels(tf))], levels(tf))
pchv <- unname(pch_map[tf])

png("RDA_chem.png", width = 8, height = 6, units = "in", res = 600, bg = "white")
plot(site[,1], site[,2], type = "n", xlab = xlab, ylab = ylab, xlim = c(-1,1), ylim = c(-1,1))
abline(h=0, v=0, lty=3)
points(site[,1], site[,2], pch = pchv, col = cols, cex = 1.2)
arrows(0, 0, sp_top[,1]*0.7, sp_top[,2]*0.7, length = 0.07)
text(sp_top[,1]*0.735, sp_top[,2]*0.735, labels = rownames(sp_top), cex = 0.8)
legend("topright", title = "Dosing", legend = levels(ab$Dosing),
       col = c("forestgreen","red","gold"), pch = 16, bty = "n")
legend("bottomright", title = "Treatment", legend = names(pch_map),
       pch = unname(pch_map), bty = "n")
dev.off()


#Hypothesis 2
####### Biomass Vs Time #######
setwd("INSERT_PATHWAY_HERE")

# --- Load data ---
df <- read.csv("temp_spring_both.csv")
df <- df[order(df$Mesocosm), ]
weights <- read.csv("weightmg.csv")

# --- Prep (counts, join weights, per-L scaling) ---
VOL_L <- 38.4
counts <- df %>%
  count(Mesocosm, Organism, Year, Treatment, name = "Abundance") %>%
  mutate(Organism = ifelse(Organism == "Daphnia.", "Daphnia", Organism))

dat <- counts %>%
  left_join(weights, by = "Organism") %>%
  mutate(
    Abundance_L = Abundance / VOL_L,
    Biomass_mg_L = (Abundance * Weight.mg.) / VOL_L
  )

# --- Community totals per mesocosm (biomass & mean individual size) ---
historical_total <- dat %>%
  group_by(Year, Mesocosm, Treatment) %>%
  summarise(
    Total_biomass_mg_L = sum(Biomass_mg_L, na.rm = TRUE),
    Total_abundance_L  = sum(Abundance_L,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_indiv_biomass_mg = ifelse(Total_abundance_L > 0,
                                   Total_biomass_mg_L / Total_abundance_L, NA_real_),
    Year = factor(Year),
    Treatment = factor(Treatment)
  )

# --- ART: total biomass (Year × Treatment) ---
art_total <- art(Total_biomass_mg_L ~ Year * Treatment,
                 data = historical_total %>% filter(is.finite(Total_biomass_mg_L)))
anova(art_total)

# --- ART: mean individual body mass (Year × Treatment) + post-hoc when needed ---
art_size <- art(Mean_indiv_biomass_mg ~ Year * Treatment,
                data = historical_total %>% filter(is.finite(Mean_indiv_biomass_mg)))
anova(art_size)

mod_year <- artlm(art_size, "Year")
emmeans(mod_year, ~ Year)
pairs(emmeans(mod_year, ~ Year), adjust = "tukey")

mod_trt <- artlm(art_size, "Treatment")
emm_trt <- emmeans(mod_trt, ~ Treatment)
pairs(emm_trt, adjust = "tukey")
cld(emm_trt, adjust = "tukey")

# --- Per-taxon mean individual mass per mesocosm (ART by taxon) ---
mass_taxon <- dat %>%
  group_by(Year, Mesocosm, Treatment, Organism) %>%
  summarise(
    total_biomass_mg = sum(Abundance * Weight.mg., na.rm = TRUE),
    total_abund      = sum(Abundance,               na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_mass_mg = ifelse(total_abund > 0, total_biomass_mg / total_abund, NA_real_),
    Year = factor(Year), Treatment = factor(Treatment), Organism = factor(Organism)
  ) %>%
  filter(is.finite(mean_mass_mg))

run_art_taxon <- function(d) {
  m <- art(mean_mass_mg ~ Year * Treatment, data = d)
  a <- as.data.frame(anova(m)); a$Effect <- rownames(a); rownames(a) <- NULL
  dplyr::transmute(a, Organism = as.character(unique(d$Organism)),
                   Effect, Df, Df.res, F = `F value`, p = `Pr(>F)`)
}
taxa_results <- mass_taxon %>%
  group_split(Organism) %>%
  purrr::map_df(run_art_taxon)

# --- Warmed-only deltas: Δabundance vs Δmean size (figure) ---
collapsed <- dat %>%
  filter(Treatment != 0) %>%
  group_by(Year, Organism) %>%
  summarise(
    total_abundance_L  = sum(Abundance_L,  na.rm = TRUE),
    total_biomass_mg_L = sum(Biomass_mg_L, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_indiv_biomass_mg = ifelse(total_abundance_L > 0,
                                   total_biomass_mg_L / total_abundance_L, NA_real_),
    mean_indiv_biomass_ug = mean_indiv_biomass_mg * 1000
  )

delta_tbl <- collapsed %>%
  tidyr::pivot_wider(
    names_from = Year,
    values_from = c(total_abundance_L, mean_indiv_biomass_ug),
    names_sep = "_"
  ) %>%
  mutate(
    delta_abundance = total_abundance_L_2025 - total_abundance_L_2019,
    delta_size_ug   = mean_indiv_biomass_ug_2025 - mean_indiv_biomass_ug_2019
  ) %>%
  filter(is.finite(delta_abundance), is.finite(delta_size_ug))

p <- ggplot(delta_tbl, aes(x = delta_abundance, y = delta_size_ug, label = Organism)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(x = expression(Delta~abundance~"(ind L"^{-1}*", 2025 − 2019)"),
       y = expression(Delta~mean~individual~biomass~"(µg, 2025 − 2019)")) +
  theme_classic(base_size = 11)
ggplot2::ggsave("delta_abundance_vs_size_warmed.png", p, width = 7, height = 5, dpi = 600, bg = "white")

# --- Controls-only check (Year-to-Year change in Treatment==0) ---
ctrl <- historical_total %>%
  dplyr::filter(Treatment == 0) %>%
  dplyr::select(Mesocosm, Year, Mean_indiv_biomass_mg) %>%
  tidyr::pivot_wider(names_from = Year, values_from = Mean_indiv_biomass_mg, names_prefix = "Y_")

eps <- 1e-9
ctrl <- ctrl %>%
  mutate(delta = Y_2025 - Y_2019,
         log_ratio = log((Y_2025 + eps) / (Y_2019 + eps)))

wilcox.test(ctrl$delta,     mu = 0, alternative = "less")
wilcox.test(ctrl$log_ratio, mu = 0, alternative = "less")


#Hypothesis 3
####### Biomass Pel vs Benthic#######
setwd("INSERT_PATHWAY_HERE")

# --- Pelagic biomass (from counts × weights) ---
pelagic  <- read.csv("final_data_legacy.csv")
weights  <- read.csv("weightmg.csv")
pelagic  <- pelagic[-c(26852:26863), ]

pel_abund <- aggregate(list(Abundance = pelagic$Organism),
                       by = list(Mesocosm = pelagic$Mesocosm,
                                 Treatment = pelagic$Treatment,
                                 Organism = pelagic$Organism),
                       FUN = length)
pel_abund <- subset(pel_abund, !Treatment %in% c("U","R","UR"))

pel_biomass <- merge(pel_abund, weights, by = "Organism", all.x = TRUE)
pel_biomass$Group_biomass_mg <- pel_biomass$Abundance * pel_biomass$Weight.mg.

pel_totals <- aggregate(list(Total_biomass_mg = pel_biomass$Group_biomass_mg),
                        by = list(Mesocosm = pel_biomass$Mesocosm,
                                  Treatment = pel_biomass$Treatment),
                        FUN = sum, na.rm = TRUE)
pel_totals$Biomass_mg_L <- pel_totals$Total_biomass_mg / 38.4
pel_totals$Biomass_g_L  <- pel_totals$Biomass_mg_L / 1000
pel_totals$Biomass_g_mesocosm <- pel_totals$Biomass_g_L * 2000

write.csv(pel_totals, "pelagic_raw_biomass.csv", row.names = FALSE)

# --- Benthic biomass (scale cores to whole mesocosm top 1 cm) ---
benthic <- read.csv("benthic.csv")
core_volume_cm3    <- 453
n_cores            <- 4
total_sampled_cm3  <- core_volume_cm3 * n_cores
mesocosm_area_cm2  <- pi * (80^2)
sediment_depth_cm  <- 1
mesocosm_volume_cm3 <- mesocosm_area_cm2 * sediment_depth_cm

ben_totals <- aggregate(list(Total_biomass_g = benthic$Drymass_g),
                        by = list(Mesocosm = benthic$Mesocosm,
                                  Treatment = benthic$Treatment),
                        FUN = sum, na.rm = TRUE)
ben_totals$Biomass_g_per_cm3  <- ben_totals$Total_biomass_g / total_sampled_cm3
ben_totals$Biomass_g_mesocosm <- ben_totals$Biomass_g_per_cm3 * mesocosm_volume_cm3

write.csv(ben_totals, "benthic_raw_biomass.csv", row.names = FALSE)

# --- Merge & model ---
pel <- read.csv("pelagic_raw_biomass.csv")
ben <- read.csv("benthic_raw_biomass.csv")

pel$Mesocosm <- as.character(pel$Mesocosm)
ben$Mesocosm <- as.character(ben$Mesocosm)

dat <- merge(pel[, c("Mesocosm","Treatment","Biomass_g_mesocosm")],
             ben[, c("Mesocosm","Biomass_g_mesocosm")],
             by = "Mesocosm", all = FALSE)
names(dat) <- c("Mesocosm","Treatment","biomass_g_meso_pel","biomass_g_meso_benthic")

m_simple <- lm(biomass_g_meso_benthic ~ biomass_g_meso_pel, data = dat)
summary(m_simple)

m_int <- lm(biomass_g_meso_benthic ~ biomass_g_meso_pel * Treatment, data = dat)
summary(m_int)

# --- Plot ---
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 needed")
library(ggplot2)
p <- ggplot(dat, aes(biomass_g_meso_pel, biomass_g_meso_benthic, colour = Treatment)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c(
    "Control"="black","1c"="#4575b4","2c"="#74add1","3c"="#abd9e9",
    "4c"="#fdae61","5c"="#f46d43","6c"="#d73027","7c"="#a50026"
  )) +
  labs(x = expression(Pelagic~biomass~"(g L"^{-1}*")"),
       y = expression(Benthic~biomass~"(g cm"^{-3}*")"),
       colour = "Treatment") +
  theme_minimal(base_size = 13)

ggplot2::ggsave("fig_benthic_vs_pelagic_by_treatment_biomass.png",
                p, width = 180, height = 120, units = "mm", dpi = 600, bg = "white")

####### Met Rate ######
# Raw biomass coupling (whole-mesocosm scaling already computed elsewhere) 
pel <- read.csv("pelagic_raw_biomass.csv")      # expects: Mesocosm, Treatment, Biomass_g_mesocosm
ben <- read.csv("benthic_raw_biomass.csv")      # expects: Mesocosm, Biomass_g_mesocosm

pel$Mesocosm <- as.character(pel$Mesocosm)
ben$Mesocosm <- as.character(ben$Mesocosm)

dat_biomass <- merge(
  pel[, c("Mesocosm","Treatment","Biomass_g_mesocosm")],
  ben[, c("Mesocosm","Biomass_g_mesocosm")],
  by = "Mesocosm"
)
names(dat_biomass) <- c("Mesocosm","Treatment","Pelagic_g","Benthic_g")

m_biomass <- lm(Benthic_g ~ Pelagic_g * Treatment, data = dat_biomass)
summary(m_biomass)

library(ggplot2)
p_biomass <- ggplot(dat_biomass, aes(Pelagic_g, Benthic_g, colour = Treatment)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(x = expression(Pelagic~biomass~"(g per mesocosm)"),
       y = expression(Benthic~biomass~"(g per mesocosm)"),
       colour = "Treatment") +
  theme_minimal(base_size = 13)
ggplot2::ggsave("scatter_benthic_vs_pelagic_biomass.png", p_biomass, width = 7, height = 5, dpi = 600, bg = "white")

#Metabolic capacity (MCI) and coupling tests 
# Inputs with individual dry masses:
pel_i <- read.csv("pelagic_with_weights.csv")   # expects: Mesocosm, Treatment, Weight_g (per individual)
ben_i <- read.csv("benthic_metrate.csv")        # expects: Mesocosm, Treatment, Drymass_g (per individual)

E <- 0.63; k <- 8.617e-5; ambient_C <- 20
calc_MCI <- function(df, mass_col) {
  df$Temp_C <- ambient_C + df$Treatment
  aggregate((df[[mass_col]]^0.75), by = list(Mesocosm = df$Mesocosm, Treatment = df$Treatment, Temp_C = df$Temp_C), FUN = sum) |>
    transform(MCI = x * exp(-E / (k * (Temp_C + 273.15)))) |>
    subset(select = c(Mesocosm, Treatment, Temp_C, MCI))
}

mci_ben <- calc_MCI(ben_i, "Drymass_g"); names(mci_ben)[4] <- "MCI_benthic"
mci_pel <- calc_MCI(pel_i, "Weight_g");  names(mci_pel)[4] <- "MCI_pelagic"

mci <- merge(mci_ben, mci_pel, by = c("Mesocosm","Treatment","Temp_C"))

# Temperature effects (per methods i–ii)
m_pel_T <- lm(MCI_pelagic  ~ Temp_C, data = mci)
m_ben_T <- lm(MCI_benthic  ~ Temp_C, data = mci)
summary(m_pel_T); summary(m_ben_T)

# Coupling: benthic ~ pelagic (+ temperature; polynomial and interaction; AIC-based)
m_add <- lm(MCI_benthic ~ poly(MCI_pelagic, 2, raw = TRUE) + Temp_C, data = mci)
m_int <- lm(MCI_benthic ~ poly(MCI_pelagic, 2, raw = TRUE) * Temp_C, data = mci)
AIC(m_add, m_int); anova(m_add, m_int); summary(m_add); summary(m_int)

# Optional non-linear check (GAM)
if (requireNamespace("mgcv", quietly = TRUE)) {
  library(mgcv)
  m_gam <- gam(MCI_benthic ~ s(MCI_pelagic, k = 4) + Temp_C, data = mci, method = "REML")
  summary(m_gam); AIC(m_gam)
}
