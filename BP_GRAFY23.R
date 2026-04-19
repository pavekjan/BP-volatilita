# =====================================================================
# BAKALÁŘSKÁ PRÁCE
# Analýza intradenních vzorců volatility a likvidity akciových trhů
# NASDAQ & NYSE · 2024 | RK Bartlett estimátor
#
# SKRIPT 02: Finální analytické grafy
# ---------------------------------------------------------------------
# Skript načítá výstupní tabulky ze Skriptu 01 a produkuje 8 grafů:
#
#   Graf 01 — Intradenní U-tvar volatility (všechny sektory, zoom od 09:35)
#   Graf 02 — U-tvar: Technologie vs. Defenzivní spotřeba (dvoupanel, 95% CI)
#   Graf 03 — Hloubka U-tvaru vs. roční volatilita (scatter + OLS regrese)
#   Graf 04 — Test MDH hypotézy: r(vol,liq) vs. r(vol,int)
#   Graf 05 — Srovnání volatility, intenzity a likvidity (Cleveland dot plot)
#   Graf 06 — Koeficient determinace R² + signifikance koeficientů (OLS+HAC)
#   Graf 07 — Heatmapa relativní volatility: sektor × čas
#   Graf 08 — Roční anualizovaná RK Bartlett volatilita po sektorech
#
# Předpoklad: Skript 01 byl úspěšně spuštěn a tabulky existují v
#   BP_srovnani_vystupy/tabulky/
# =====================================================================

library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)
library(gridExtra)
library(grid)

# ── Cesta k výstupům ──────────────────────────────────────────────────
OUT_DIR <- file.path(getwd(), "BP_srovnani_vystupy")
dir.create(file.path(OUT_DIR, "grafy_final"), recursive = TRUE,
           showWarnings = FALSE)

# ── Načtení dat ───────────────────────────────────────────────────────
souhrn  <- fread(file.path(OUT_DIR, "tabulky", "01_souhrn_akcie.csv"))
profily <- fread(file.path(OUT_DIR, "tabulky", "02_sezonni_profily.csv"))
testy   <- fread(file.path(OUT_DIR, "tabulky", "03_statisticke_testy.csv"))
ols     <- fread(file.path(OUT_DIR, "tabulky", "04_regrese_ols_hac.csv"))

# ── Čistá jména symbolů (bez přípony .O / .N) ────────────────────────
for (dt in list(souhrn, profily, testy, ols)) {
  if ("symbol" %in% names(dt))
    dt[, sym_clean := gsub("\\.[ON]$", "", symbol)]
}

# ── Převod datových typů ──────────────────────────────────────────────
num_souhrn <- c("rano_0930","zaverecna_1555","u_hloubka","rkb_ann_pct",
                "r_vol_int","r_vol_liq","avg_ticks_den")
souhrn[, (num_souhrn) := lapply(.SD, as.numeric), .SDcols = num_souhrn]
profily[, rkb_mean          := as.numeric(rkb_mean)]
profily[, int_mean          := as.numeric(int_mean)]
profily[, liq_mean          := as.numeric(liq_mean)]
profily[, minutes_from_open := as.integer(minutes_from_open)]
testy[, r_vol_int := as.numeric(r_vol_int)]
testy[, r_vol_liq := as.numeric(r_vol_liq)]
ols[, c("beta_int","beta_int_se_hac","beta_liq","beta_liq_se_hac",
        "r_squared","beta_int_signif","beta_liq_signif") :=
      list(as.numeric(beta_int), as.numeric(beta_int_se_hac),
           as.numeric(beta_liq),  as.numeric(beta_liq_se_hac),
           as.numeric(r_squared),  as.character(beta_int_signif),
           as.character(beta_liq_signif))]

# ── Barvy sektorů ─────────────────────────────────────────────────────
# Konzistentní barevné schéma pro všechny grafy.
SCOLS <- c(
  "Technologie"         = "#1B6CA8",
  "Finance"             = "#C0392B",
  "Energie"             = "#E07B39",
  "Zdravotnictvi"       = "#27AE60",
  "Defenzivni spotreba" = "#8E44AD",
  "Cyklicka spotreba"   = "#D4AC0D",
  "Prumysl"             = "#2C3E50",
  "Telekomunikace"      = "#16A085"
)

# ── Jednotné grafické téma ────────────────────────────────────────────
TH <- theme_bw(base_size = 12) %+replace% theme(
  plot.title       = element_text(size = 13, face = "bold",
                                  margin = margin(b = 4)),
  plot.subtitle    = element_text(size = 10, color = "grey25",
                                  lineheight = 1.3, margin = margin(b = 6)),
  plot.caption     = element_text(size = 8, color = "grey50", hjust = 0,
                                  margin = margin(t = 6)),
  axis.title       = element_text(size = 11),
  axis.text        = element_text(size = 10),
  legend.position  = "bottom",
  legend.title     = element_text(size = 10, face = "bold"),
  legend.text      = element_text(size = 9),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey93", linewidth = 0.35),
  strip.background = element_rect(fill = "grey96", color = "grey80"),
  strip.text       = element_text(size = 10, face = "bold"),
  plot.margin      = margin(10, 14, 8, 10)
)

# ── Pomocná funkce pro uložení grafů ─────────────────────────────────
uloz <- function(g, n, w = 15, h = 8) {
  ggsave(file.path(OUT_DIR, "grafy_final", n),
         plot = g, width = w, height = h, dpi = 300, bg = "white")
  message("  \u2713 ", n)
}

# Rovnoměrné rozdělení = 1/78 ≈ 1,28 % (referenční hodnota pro grafy)
EQ <- 100 / 78


# ====================================================================
# GRAF 01 — Intradenní U-tvar volatility (všechny sektory)
# ====================================================================
# Hlavní výsledkový graf práce. Zobrazuje průměrný sezónní podíl
# RK Bartlett volatility v každém 5minutovém intervalu, zprůměrovaný
# přes všechny akcie v sektoru a přes 252 obchodních dní.
#
# Graf je záměrně zobrazen od 09:35 — ranní spike v 09:30 dosahuje
# 14–23 % denní volatility a při jeho zahrnutí by byl zbytek dne
# nečitelný. Hodnota 09:30 je uvedena v popisku grafu.
#
# Kruskal-Wallisův test: p < 0,001 pro všech 23 akcií → U-tvar je
# statisticky průkazný a není způsoben několika málo extrémy.
# ====================================================================

# Průměrný profil za každý sektor
profil_s <- profily[, .(rkb_mean = mean(rkb_mean, na.rm = TRUE)),
                    by = .(sektor, minutes_from_open)]
setorder(profil_s, sektor, minutes_from_open)

profil_s_zoom <- profil_s[minutes_from_open >= 5]

g1 <- ggplot(profil_s_zoom,
             aes(x = minutes_from_open, y = rkb_mean,
                 color = sektor, group = sektor)) +
  annotate("rect", xmin = 5,   xmax = 60,  ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.55) +
  annotate("rect", xmin = 120, xmax = 240, ymin = -Inf, ymax = Inf,
           fill = "#D1ECF1", alpha = 0.45) +
  annotate("rect", xmin = 330, xmax = 390, ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.55) +
  annotate("text", x = 32,  y = 4.45, label = "Otevírací\nfáze",
           size = 3.0, color = "#7D6608", fontface = "bold", hjust = 0.5) +
  annotate("text", x = 180, y = 4.45, label = "Polední\nklid",
           size = 3.0, color = "#0B5563", fontface = "bold", hjust = 0.5) +
  annotate("text", x = 360, y = 4.45, label = "Závěrečný\nnárůst",
           size = 3.0, color = "#7D6608", fontface = "bold", hjust = 0.5) +
  # Šipky znázorňující tvar U-křivky
  annotate("segment", x = 30, xend = 180, y = 3.8, yend = 0.6,
           arrow = arrow(length = unit(0.2,"cm"), type = "closed"),
           color = "grey50", linewidth = 0.5) +
  annotate("segment", x = 330, xend = 200, y = 3.8, yend = 0.6,
           arrow = arrow(length = unit(0.2,"cm"), type = "closed"),
           color = "grey50", linewidth = 0.5) +
  annotate("text", x = 195, y = 0.35, label = "Polední\nminimum",
           size = 2.7, color = "grey40", hjust = 0.5, fontface = "italic") +
  geom_hline(yintercept = EQ, linetype = "dashed",
             color = "grey45", linewidth = 0.7) +
  annotate("text", x = 280, y = EQ + 0.12,
           label = sprintf("Rovnoměrné = %.2f %%", EQ),
           size = 2.8, color = "grey40") +
  geom_line(linewidth = 1.3, alpha = 0.92) +
  scale_color_manual(values = SCOLS, name = "Sektor") +
  scale_x_continuous(
    breaks = seq(5, 390, 60),
    labels = c("09:35","10:35","11:35","12:35","13:35","14:35","15:35"),
    limits = c(0, 395)) +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), " %"),
    limits = c(0, 5),
    expand = expansion(mult = c(0.01, 0.06))) +
  labs(
    title    = "Intradenní U-tvar volatility (RK Bartlett) — NASDAQ & NYSE · 2024",
    subtitle = paste0(
      "Průměrný sezónní podíl RKB(5min)/RKB(den). Graf zobrazuje od 09:35 ",
      "(ranní spike 09:30 dosahuje 14–23 %, odříznut pro přehlednost).\n",
      "Kruskal-Wallis test sezónnosti: p < 0,001 pro všech 23 akcií. ",
      "Andersen & Bollerslev (1997)."),
    x = "Čas obchodování (Eastern Time)",
    y = "Sezónní podíl na denní volatilitě [%]",
    caption = paste0(
      "Pozn.: ranní spike při 09:30 dosahuje průměrně 14–23 % denní volatility. ",
      "RK Bartlett: Barndorff-Nielsen et al. (2008). N = 23 akcií, 252 dní.")
  ) +
  TH + guides(color = guide_legend(nrow = 2))

uloz(g1, "Graf01_Utvar_zoom.png", 16, 8)


# ====================================================================
# GRAF 02 — U-tvar: Technologie vs. Defenzivní spotřeba (dvoupanel)
# ====================================================================
# Srovnání dvou sektorů s nejodlišnějším chováním:
#   Panel A: celý den od 09:30 — viditelný rozdíl ranního spike
#   Panel B: zoom od 09:35 — viditelný rozdíl poledního minima
#
# Stínovaná pásma = 95% intervaly spolehlivosti (±1,96 × SE).
# Kde se pásma nepřekrývají, je rozdíl statisticky průkazný.
#
# Klíčové zjištění: defenzivní spotřeba má paradoxně vyšší ranní
# spike (22,6 %) než technologie (17,6 %), ale technologie drží
# vyšší volatilitu v poledních hodinách → mělčí U-tvar.
# ====================================================================

pr2_full <- profily[sektor %in% c("Technologie", "Defenzivni spotreba"),
                    .(mean_r = mean(rkb_mean, na.rm = TRUE),
                      se_r   = sd(rkb_mean,  na.rm = TRUE) / sqrt(.N)),
                    by = .(sektor, minutes_from_open, hour_minute)]
pr2_full[, ci_lo := mean_r - 1.96 * se_r]
pr2_full[, ci_hi := mean_r + 1.96 * se_r]
setorder(pr2_full, sektor, minutes_from_open)

pr2_zoom <- pr2_full[minutes_from_open >= 5]

barvy  <- c("Technologie"         = "#1B6CA8",
            "Defenzivni spotreba" = "#8E44AD")
labely <- c("Technologie"         = "Technologie (NVDA, AMD, MSFT, AAPL, …)",
            "Defenzivni spotreba" = "Defenzivní spotřeba (KO, PG)")

# Panel A: celý den (viditelný ranní spike)
gA <- ggplot(pr2_full,
             aes(x = minutes_from_open, y = mean_r,
                 color = sektor, fill = sektor, group = sektor)) +
  annotate("rect", xmin = 0,   xmax = 60,  ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.5) +
  annotate("rect", xmin = 120, xmax = 240, ymin = -Inf, ymax = Inf,
           fill = "#D1ECF1", alpha = 0.4) +
  annotate("rect", xmin = 330, xmax = 390, ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.5) +
  annotate("text", x = 30,  y = 21, label = "Otevírací fáze",
           size = 3.0, color = "#7D6608", fontface = "bold", hjust = 0.5) +
  annotate("text", x = 180, y = 21, label = "Polední klid",
           size = 3.0, color = "#0B5563", fontface = "bold", hjust = 0.5) +
  annotate("text", x = 360, y = 21, label = "Závěrečný nárůst",
           size = 3.0, color = "#7D6608", fontface = "bold", hjust = 0.5) +
  annotate("text", x = 3, y = 23.5,
           label = sprintf("Def.: %.1f %%",
                           pr2_full[sektor == "Defenzivni spotreba" &
                                      minutes_from_open == 0, mean_r]),
           size = 3.0, color = "#8E44AD", hjust = 0, fontface = "bold") +
  annotate("text", x = 3, y = 18,
           label = sprintf("Tech.: %.1f %%",
                           pr2_full[sektor == "Technologie" &
                                      minutes_from_open == 0, mean_r]),
           size = 3.0, color = "#1B6CA8", hjust = 0, fontface = "bold") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, color = NA) +
  geom_hline(yintercept = EQ, linetype = "dashed",
             color = "grey45", linewidth = 0.6) +
  geom_line(linewidth = 1.6) +
  scale_color_manual(values = barvy, labels = labely, name = NULL) +
  scale_fill_manual(values  = barvy, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 60, 120, 180, 240, 300, 360, 390),
    labels = c("09:30","10:30","11:30","12:30","13:30","14:30","15:30","16:00"),
    limits = c(0, 395)) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), " %"),
                     limits = c(0, 25)) +
  labs(title    = "A) Celý obchodní den — ranní spike ukazuje klíčový rozdíl",
       subtitle = "Defenzivní spotřeba dosahuje při 09:30 vyššího spike (~23 %) než technologie (~17 %) → hlubší U-tvar.",
       x = NULL, y = "Sezónní podíl na denní volatilitě [%]") +
  TH %+replace% theme(legend.position = "none",
                      plot.title = element_text(size = 11, face = "bold"))

# Panel B: zoom od 09:35 (viditelný rozdíl poledního minima)
tech_pol <- pr2_zoom[sektor == "Technologie" &
                       minutes_from_open %in% 120:235, min(mean_r)]
def_pol  <- pr2_zoom[sektor == "Defenzivni spotreba" &
                       minutes_from_open %in% 120:235, min(mean_r)]

gB <- ggplot(pr2_zoom,
             aes(x = minutes_from_open, y = mean_r,
                 color = sektor, fill = sektor, group = sektor)) +
  annotate("rect", xmin = 5,   xmax = 60,  ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.5) +
  annotate("rect", xmin = 120, xmax = 240, ymin = -Inf, ymax = Inf,
           fill = "#D1ECF1", alpha = 0.4) +
  annotate("rect", xmin = 330, xmax = 390, ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.5) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, color = NA) +
  geom_hline(yintercept = EQ, linetype = "dashed",
             color = "grey45", linewidth = 0.6) +
  annotate("text", x = 270, y = EQ + 0.12,
           label = sprintf("Rovnoměrné = %.2f %%", EQ),
           size = 2.8, color = "grey40") +
  geom_line(linewidth = 1.6) +
  # Šipka znázorňující rozdíl v poledním minimu
  annotate("segment", x = 185, xend = 185,
           y = tech_pol + 0.04, yend = def_pol - 0.04,
           arrow = arrow(ends = "both", length = unit(0.15,"cm"),
                         type = "closed"),
           color = "grey30", linewidth = 0.7) +
  annotate("text", x = 200, y = (tech_pol + def_pol) / 2,
           label = "Tech. výše\nv poledne\n→ mělčí U",
           size = 2.6, color = "grey25", hjust = 0, fontface = "italic") +
  scale_color_manual(values = barvy, labels = labely, name = "Sektor") +
  scale_fill_manual(values  = barvy, guide = "none") +
  scale_x_continuous(
    breaks = c(5, 65, 125, 185, 245, 305, 365, 390),
    labels = c("09:35","10:35","11:35","12:35","13:35","14:35","15:35","16:00"),
    limits = c(0, 395)) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), " %"),
                     limits = c(0, 5),
                     expand = expansion(mult = c(0.01, 0.06))) +
  labs(title    = "B) Zoom od 09:35 — technologie drží v poledne vyšší volatilitu",
       subtitle = "Technologie neklesají tak hluboko do poledního minima → volatilita je rovnoměrněji rozložena přes celý den.",
       x = "Čas obchodování (Eastern Time)",
       y = "Sezónní podíl na denní volatilitě [%]",
       caption = paste0("RK Bartlett: Barndorff-Nielsen et al. (2008). ",
                        "95% CI = ±1,96 × SE průměru přes akcie v sektoru. ",
                        "N = 8 technologických a 2 defenzivní akcie, 252 obchodních dní.")) +
  TH %+replace% theme(legend.position = "bottom",
                      plot.title = element_text(size = 11, face = "bold"))

# Složení do jednoho grafu
nadpis <- textGrob(
  "U-tvar volatility: Technologie vs. Defenzivní spotřeba · 2024",
  gp = gpar(fontsize = 13, fontface = "bold"))

g2_final <- arrangeGrob(gA, gB, nrow = 2, top = nadpis, heights = c(1, 1))

ggsave(file.path(OUT_DIR, "grafy_final", "Graf02_Dvoupanel.png"),
       plot = g2_final, width = 15, height = 12, dpi = 300, bg = "white")
message("  \u2713 Graf02_Dvoupanel.png")


# ====================================================================
# GRAF 03 — Hloubka U-tvaru vs. roční volatilita (scatter)
# ====================================================================
# Testuje hypotézu: liší se hloubka U-tvaru systematicky podle
# celkové úrovně volatility akcie?
#
# Hloubka U-tvaru = ranní spike (09:30) / polední minimum.
# Vyšší hodnota = výraznější U-tvar = silnější intradenní sezónnost.
#
# Přerušovaná čára = OLS regresní přímka s 95% intervalem spolehlivosti.
# Šedý pás = 95% CI odhadu.
#
# Výsledek: negativní vztah — technologické akcie s vysokou roční
# volatilitou mají mělčí U-tvar, defenzivní akcie s nízkou roční
# volatilitou mají U-tvar nejhlubší. Interpretace: technologické akcie
# reagují na informace průběžně celý den, defenzivní se rozhoupe
# hlavně při otevření a pak zklidní.
# ====================================================================

g3 <- ggplot(souhrn,
             aes(x = rkb_ann_pct, y = u_hloubka,
                 color = sektor, label = sym_clean)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              fill = "grey88", linewidth = 0.9,
              linetype = "dashed", fullrange = FALSE) +
  # Přerušované čáry značící průměry (kvadranty)
  geom_vline(xintercept = mean(souhrn$rkb_ann_pct),
             linetype = "dotted", color = "grey60", linewidth = 0.6) +
  geom_hline(yintercept = mean(souhrn$u_hloubka, na.rm = TRUE),
             linetype = "dotted", color = "grey60", linewidth = 0.6) +
  annotate("text", x = 14, y = 50,
           label = "Nízká volatilita\nHluboký U-tvar",
           size = 2.9, color = "grey45", fontface = "italic", hjust = 0.5) +
  annotate("text", x = 52, y = 50,
           label = "Vysoká volatilita\nHluboký U-tvar",
           size = 2.9, color = "grey45", fontface = "italic", hjust = 0.5) +
  annotate("text", x = 14, y = 30.5,
           label = "Nízká volatilita\nMělký U-tvar",
           size = 2.9, color = "grey45", fontface = "italic", hjust = 0.5) +
  annotate("text", x = 52, y = 30.5,
           label = "Vysoká volatilita\nMělký U-tvar",
           size = 2.9, color = "grey45", fontface = "italic", hjust = 0.5) +
  geom_point(size = 5.0, alpha = 0.90) +
  geom_text_repel(size = 3.2, fontface = "bold",
                  box.padding = 0.45, show.legend = FALSE,
                  segment.color = "grey65", segment.size = 0.35,
                  max.overlaps = 30) +
  scale_color_manual(values = SCOLS, name = "Sektor") +
  scale_x_continuous(labels = function(x) paste0(x, " %"),
                     breaks = seq(10, 60, 10)) +
  scale_y_continuous(labels = function(x) paste0(round(x), "\u00d7"),
                     breaks = seq(30, 55, 5)) +
  labs(
    title    = "Hloubka U-tvaru vs. roční volatilita — 23 akcií · 2024",
    subtitle = paste0(
      "Negativní vztah: akcie s vyšší roční volatilitou (technologie) mají MĚLČÍ U-tvar.\n",
      "Defenzivní akcie mají volatilitu více soustředěnou v otevírací fázi \u2192 hlubší U-tvar. ",
      "Šedý pás = 95% CI OLS regrese."),
    x = "Roční anualizovaná RK Bartlett volatilita [% p.a.]",
    y = "Hloubka U-tvaru (ranní spike \u00f7 polední minimum)",
    caption = paste0("RK Bartlett: Barndorff-Nielsen et al. (2008). ",
                     "Sezónní podíly: Andersen & Bollerslev (1997). N = 23 akcií.")
  ) +
  TH

uloz(g3, "Graf03_Hloubka_vs_Volatilita_scatter.png", 13, 9)


# ====================================================================
# GRAF 04 — Test MDH hypotézy: r(vol,liq) vs. r(vol,int)
# ====================================================================
# Mixture of Distributions Hypothesis (Clark, 1973) předpokládá,
# že volatilitu a objem pohání společný latentní informační tok.
# Pokud MDH platí, korelace r(volatilita, objem) > r(volatilita, počet
# transakcí), protože objem lépe zachycuje informační tok.
#
# Osa x: korelace volatility s intenzitou (počtem transakcí)
# Osa y: korelace volatility s likviditou (objemem)
# Přerušovaná diagonála: situace kdy jsou obě korelace stejné.
#
# Výsledek: všech 23 akcií leží nad diagonálou → MDH podporováno.
# ====================================================================

ref <- unique(souhrn[, .(sym_clean, sektor)])
t2  <- merge(testy, ref, by = "sym_clean", all.x = TRUE)
if ("sektor.x" %in% names(t2)) {
  t2[, sektor := fifelse(!is.na(sektor.y), sektor.y, sektor.x)]
  t2[, c("sektor.x","sektor.y") := NULL]
}

g4 <- ggplot(t2, aes(x = r_vol_int, y = r_vol_liq,
                     color = sektor, label = sym_clean)) +
  annotate("polygon",
           x = c(-0.05, 0.95, 0.95, -0.05),
           y = c(-0.05, 0.95, -0.05, -0.05),
           fill = "grey96", color = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.9) +
  annotate("text", x = 0.70, y = 0.62, label = "r(liq) = r(int)",
           color = "grey50", size = 3.1, angle = 40, fontface = "italic") +
  annotate("text", x = 0.08, y = 0.925,
           label = "\u2190 Likvidita silnější (MDH podporováno)",
           color = "#16A085", size = 3.1, fontface = "bold.italic") +
  geom_point(size = 4.8, alpha = 0.88) +
  geom_text_repel(size = 3.0, fontface = "bold", box.padding = 0.4,
                  show.legend = FALSE, max.overlaps = 30,
                  segment.color = "grey70", segment.size = 0.4) +
  scale_color_manual(values = SCOLS, name = "Sektor") +
  scale_x_continuous(limits = c(-0.05, 0.95),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(limits = c(0.45, 0.95),
                     labels = function(x) sprintf("%.1f", x)) +
  labs(
    title    = "Test MDH: Pearsonovy korelace r(vol, liq) vs. r(vol, int) — 2024",
    subtitle = paste0(
      "Každý bod = jedna akcie. Všech 23 akcií leží nad přerušovanou diagonálou:\n",
      "r(volatilita, likvidita) > r(volatilita, intenzita) \u2192 ",
      "konzistentní s Mixture of Distributions Hypothesis (Clark, 1973)."),
    x = "Pearsonova korelace r(volatilita, intenzita obchodování)",
    y = "Pearsonova korelace r(volatilita, likvidita)",
    caption = paste0("Korelace spočítány na průměrných sezónních podílech (78 intervalů). ",
                     "MDH: Clark (1973). RK Bartlett: Barndorff-Nielsen et al. (2008).")
  ) +
  TH

uloz(g4, "Graf04_MDH_korelace.png", 13, 9)


# ====================================================================
# GRAF 05 — Srovnání volatility, intenzity a likvidity (Cleveland)
# ====================================================================
# Cleveland dot plot srovnává tři klíčové metriky najednou pro každou
# akcii. Všechny hodnoty jsou normalizovány na škálu 0–100, aby byly
# srovnatelné přes různé jednotky.
#
# Normalizace (min-max): x_norm = (x - min) / (max - min) × 100
#
# Metriky:
#   ▲ Intenzita = průměrný počet ticků za den
#   ■ Volatilita = RK Bartlett anualizovaná [% p.a.]
#   ● Likvidita = průměrný denní objem (počet zobchodovaných akcií)
#
# Poznámka k NVDA: všechny tři metriky dosahují hodnoty 100 —
# NVDA je absolutní maximum vzorku ve všech dimenzích. Kolečko
# (likvidita) překrývá čtverec i trojúhelník, protože hodnoty jsou
# shodné (100). Částečně jde o efekt akciového splitu 10:1 v červnu
# 2024, který dramaticky zvýšil počet zobchodovaných kusů.
# ====================================================================

# Načti avg_volume_den z cache souborů
cache_dir <- file.path(OUT_DIR, "cache")
rds_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)

vol_list <- vector("list", length(rds_files))
for (i in seq_along(rds_files)) {
  tryCatch({
    res     <- readRDS(rds_files[i])
    if (is.null(res) || nrow(res) == 0) next
    vol_col <- intersect(c("volume_sum","VOLUME_SUM","vol"), names(res))
    if (length(vol_col) == 0) next
    sym     <- res$symbol[1]
    avg_vol <- res[, .(vol_den = sum(get(vol_col[1]), na.rm = TRUE)),
                   by = date_ny][, mean(vol_den, na.rm = TRUE)]
    vol_list[[i]] <- data.table(symbol = sym,
                                avg_volume_den = round(avg_vol))
  }, error = function(e) NULL)
}
vol_dt <- rbindlist(Filter(Negate(is.null), vol_list))

if ("avg_volume_den" %in% names(souhrn)) souhrn[, avg_volume_den := NULL]
souhrn <- merge(souhrn, vol_dt, by = "symbol", all.x = TRUE)

# Min-max normalizace
souhrn[, vol_norm := round((rkb_ann_pct    - min(rkb_ann_pct,    na.rm=TRUE)) /
                             (max(rkb_ann_pct,    na.rm=TRUE) -
                                min(rkb_ann_pct,    na.rm=TRUE)) * 100, 1)]
souhrn[, int_norm := round((avg_ticks_den  - min(avg_ticks_den,  na.rm=TRUE)) /
                             (max(avg_ticks_den,  na.rm=TRUE) -
                                min(avg_ticks_den,  na.rm=TRUE)) * 100, 1)]
souhrn[, liq_norm := round((avg_volume_den - min(avg_volume_den, na.rm=TRUE)) /
                             (max(avg_volume_den, na.rm=TRUE) -
                                min(avg_volume_den, na.rm=TRUE)) * 100, 1)]

setorder(souhrn, sektor, vol_norm)
poradi <- souhrn$sym_clean

dtl <- rbindlist(list(
  souhrn[, .(sym_clean, sektor, x = int_norm,
             metrika = "\u25b2  Intenzita (avg. tick/den)")],
  souhrn[, .(sym_clean, sektor, x = vol_norm,
             metrika = "\u25a0  Volatilita (RK Bartlett, % p.a.)")],
  souhrn[, .(sym_clean, sektor, x = liq_norm,
             metrika = "\u25cf  Likvidita (avg. denní objem, ks)")]
))
dtl[, sym_clean := factor(sym_clean, levels = poradi)]
dtl[, metrika   := factor(metrika, levels = c(
  "\u25b2  Intenzita (avg. tick/den)",
  "\u25a0  Volatilita (RK Bartlett, % p.a.)",
  "\u25cf  Likvidita (avg. denní objem, ks)"))]

rozsah5 <- dtl[, .(xmin = min(x, na.rm=TRUE), xmax = max(x, na.rm=TRUE)),
               by = .(sym_clean, sektor)]

g5 <- ggplot(dtl, aes(x = x, y = sym_clean)) +
  geom_segment(data = rozsah5,
               aes(x = xmin, xend = xmax, y = sym_clean, yend = sym_clean,
                   color = sektor),
               linewidth = 0.55, alpha = 0.28, inherit.aes = FALSE) +
  geom_point(aes(shape = metrika, color = sektor), size = 3.6, stroke = 0.9,
             na.rm = TRUE) +
  geom_vline(xintercept = 50, linetype = "dotted",
             color = "grey60", linewidth = 0.5) +
  scale_shape_manual(
    values = c("\u25b2  Intenzita (avg. tick/den)"        = 17,
               "\u25a0  Volatilita (RK Bartlett, % p.a.)"  = 15,
               "\u25cf  Likvidita (avg. denní objem, ks)"   = 16),
    name = "Metrika") +
  scale_color_manual(values = SCOLS, name = "Sektor") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0\n(min)", "25", "50", "75", "100\n(max)"),
                     limits = c(-3, 113)) +
  facet_grid(sektor ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "Srovnání volatility, intenzity a likvidity — 23 akcií · 2024",
    subtitle = paste0(
      "Hodnoty normalizovány na škálu 0–100 (0 = minimum, 100 = maximum celého vzorku).\n",
      "\u25b2 Intenzita = průměrný počet ticků/den   ",
      "\u25a0 Volatilita = RK Bartlett anualizovaná [% p.a.]   ",
      "\u25cf Likvidita = průměrný denní objem (počet zobchodovaných akcií)"),
    x = "Normalizovaná hodnota (0 = nejnižší, 100 = nejvyšší v celém vzorku)",
    y = NULL,
    caption = paste0("Zdroj: vlastní výpočty z high-frequency dat. ",
                     "RK Bartlett: Barndorff-Nielsen et al. (2008). N = 23 akcií.")
  ) +
  TH %+replace% theme(
    strip.text.y       = element_text(angle = 0, size = 9.5, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90")
  ) +
  guides(
    shape = guide_legend(order = 1,
                         override.aes = list(size = 4.5, color = "grey30")),
    color = guide_legend(order = 2, nrow = 2)
  )

uloz(g5, "Graf05_3metriky.png", 15, 13)


# ====================================================================
# GRAF 06 — Koeficient determinace R² (OLS model + HAC)
# ====================================================================
# Zobrazuje R² pro každou akcii zvlášť.
# Model: RKB_podíl = α + β_int × Intenzita + β_liq × Likvidita
# Standardní chyby: HAC (Newey-West), lag = 2.
#
# Barva sloupce ukazuje signifikanci koeficientů:
#   Zelená  = oba β signifikantní (p < 0,05)
#   Modrá   = jen β_liq signifikantní
#   Šedá    = oba nesignifikantní
#
# Interpretace: vysoké R² u technologií (AMD 0,929) ukazuje, že
# intenzita a likvidita spolehlivě kopírují tvar sezónního profilu.
# Nízké R² u defenzivních akcií (JNJ 0,332, KO 0,368) signalizuje
# slabší propojení likvidity s volatilitou — tyto akcie mají sezónní
# profil méně závislý na aktivitě obchodování.
# ====================================================================

ols2 <- merge(ols, ref, by = "sym_clean", all.x = TRUE)
if ("sektor.x" %in% names(ols2)) {
  ols2[, sektor := fifelse(!is.na(sektor.y), sektor.y, sektor.x)]
  ols2[, c("sektor.x","sektor.y") := NULL]
}
setorder(ols2, r_squared)
ols2[, sym_f3 := factor(sym_clean, levels = unique(sym_clean))]

ols2[, kat := fcase(
  beta_liq_signif %in% c("*","**","***") & beta_int_signif %in% c("*","**","***"),
  "Oba \u03b2 sign.",
  beta_liq_signif %in% c("*","**","***") & !(beta_int_signif %in% c("*","**","***")),
  "Jen \u03b2_liq sign.",
  !(beta_liq_signif %in% c("*","**","***")) & beta_int_signif %in% c("*","**","***"),
  "Jen \u03b2_int sign.",
  default = "Oba nesign."
)]
ols2[, kat := factor(kat, levels = c("Oba \u03b2 sign.", "Jen \u03b2_liq sign.",
                                     "Jen \u03b2_int sign.", "Oba nesign."))]
prumer_r2 <- mean(ols2$r_squared)

g6 <- ggplot(ols2, aes(x = r_squared, y = sym_f3, fill = kat)) +
  geom_col(width = 0.74, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.3f", r_squared)),
            hjust = -0.12, size = 3.1, color = "grey20") +
  geom_vline(xintercept = prumer_r2, linetype = "dashed",
             color = "grey30", linewidth = 0.9) +
  annotate("text", x = prumer_r2 + 0.013, y = 3,
           hjust = 0, size = 3.0, color = "grey25",
           label = sprintf("Průměr\nR\u00b2 = %.3f", prumer_r2)) +
  scale_fill_manual(
    values = c("Oba \u03b2 sign."     = "#27AE60",
               "Jen \u03b2_liq sign." = "#2980B9",
               "Jen \u03b2_int sign." = "#E67E22",
               "Oba nesign."         = "#BDC3C7"),
    name = "Signifikance koeficientů (HAC, p < 0,05)") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = function(x) sprintf("%.2f", x),
                     limits = c(0, 1.05)) +
  labs(
    title    = "Koeficient determinace R\u00b2 — OLS model sezónního profilu · 2024",
    subtitle = paste0(
      "Model: RKB_podíl = \u03b1 + \u03b2_int×Intenzita + \u03b2_liq×Likvidita  (HAC std. chyby).\n",
      "Průměrné R\u00b2 = ", sprintf("%.3f", prumer_r2),
      ". Nejsilnější modely: AMD (0,929), GS (0,861), GE (0,838). ",
      "Nejslabší: JNJ (0,332), KO (0,368), T (0,357)."),
    x = "Koeficient determinace R\u00b2",
    y = NULL,
    caption = "HAC: Newey & West (1987). N = 78 sezónních intervalů na akcii."
  ) +
  TH %+replace% theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom"
  )

uloz(g6, "Graf06_R2_signifikance.png", 13, 10)


# ====================================================================
# GRAF 07 — Heatmapa relativní volatility: sektor × čas
# ====================================================================
# Syntetický pohled na celý intradenní vzorec najednou.
# Index = průměrný RK Bartlett podíl (sektor × 30min blok) /
#         průměr celého vzorku.
#
# Hodnota 1,00 = průměr. Nad 1,00 = nadprůměrná volatilita.
# Pod 1,00 = podprůměrná volatilita.
#
# Barevná škála: červená = nadprůměr, bílá = průměr, modrá = podprůměr.
#
# Klíčové pozorování: sloupec 09:30 je sytě červený pro všechny
# sektory (index 4,67–5,72), což potvrzuje universálnost ranního spike.
# Telekomunikace a defenzivní spotřeba mají ranní index vyšší než
# technologie — konzistentní se zjištěním o hloubce U-tvaru (Graf 03).
# ====================================================================

profily[, blok30 := {
  h   <- 9L + (30L + minutes_from_open) %/% 60L
  m   <- (30L + minutes_from_open) %% 60L
  m30 <- (m %/% 30L) * 30L
  sprintf("%02d:%02d", h, m30)
}]

blk <- c("09:30","10:00","10:30","11:00","11:30","12:00","12:30",
         "13:00","13:30","14:00","14:30","15:00","15:30")

heat <- profily[blok30 %in% blk,
                .(rkb = mean(rkb_mean, na.rm = TRUE)),
                by = .(sektor, blok30)]
gm   <- mean(heat$rkb, na.rm = TRUE)
heat[, idx    := round(rkb / gm, 2)]
heat[, blok30 := factor(blok30, levels = blk)]

# Seřazení sektorů od nejvyšší ranní volatility
poradi_sektoru <- heat[blok30 == "09:30"][order(-idx), sektor]
heat[, sektor  := factor(sektor, levels = rev(poradi_sektoru))]

# Barva textu: bílá na tmavém pozadí, tmavá na světlém
heat[, txt_color := ifelse(idx > 2.0 | idx < 0.60, "white", "grey15")]

g7 <- ggplot(heat, aes(x = blok30, y = sektor, fill = idx)) +
  geom_tile(color = "white", linewidth = 1.0) +
  geom_text(aes(label = sprintf("%.2f", idx), color = txt_color),
            size = 3.2, fontface = "bold") +
  # Svislé čáry oddělující fáze dne
  geom_vline(xintercept = 1.5, color = "grey20",
             linewidth = 1.4, linetype = "solid") +
  geom_vline(xintercept = 2.5, color = "grey45",
             linewidth = 0.7, linetype = "dashed") +
  geom_vline(xintercept = 12.5, color = "grey45",
             linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradientn(
    colours = c("#2166AC","#92C5DE","#F7F7F7",
                "#FDDBC7","#F4A582","#D6604D","#B2182B"),
    values  = rescale(c(0.30, 0.70, 1.00, 1.30, 1.80, 2.50, 5.80)),
    limits  = c(0.30, 5.80),
    name    = "Index relativní volatility  (1,00 = průměr celého vzorku)",
    guide   = guide_colorbar(barwidth = 22, barheight = 0.85,
                             title.position = "top", title.hjust = 0.5,
                             ticks.colour = "grey40",
                             frame.colour = "grey40")) +
  scale_color_identity() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title    = "Heatmapa relativní volatility — sektory \u00d7 čas · 2024",
    subtitle = paste0(
      "Index = RK Bartlett podíl / průměr vzorku  |  ",
      "Bílá = průměr (1,00)  |  Červená = nadprůměrná volatilita  |  Modrá = podprůměrná\n",
      "Silná svislá čára = konec ranního spike (po 09:30)  |  ",
      "Přerušované čáry = přechod do poledního klidu (10:30) a závěrečné fáze (15:30)"),
    x = "30minutový blok (začátek, Eastern Time)",
    y = NULL,
    caption = paste0("Průměr přes všechny akcie daného sektoru. ",
                     "RK Bartlett: Barndorff-Nielsen et al. (2008). ",
                     "N = 23 akcií, 252 obchodních dní.")
  ) +
  TH %+replace% theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    plot.margin     = margin(8, 14, 8, 10),
    axis.title.x    = element_text(margin = margin(t = 15))
  )

uloz(g7, "Graf07_Heatmap.png", 18, 7)


# ====================================================================
# GRAF 08 — Roční anualizovaná volatilita po sektorech
# ====================================================================
# Přehledový graf charakterizující výběrový vzorek.
# σ_ann = √(252 × průměrná denní RKB)
#
# Graf ukazuje variabilitu volatility napříč sektory a akciemi —
# klíčový kontext pro interpretaci ostatních výsledků.
#
# Přerušovaná čára = průměr celého vzorku (23 akcií).
# ====================================================================

setorder(souhrn, sektor, rkb_ann_pct)
souhrn[, sym_f8 := factor(sym_clean, levels = unique(sym_clean))]
prumer_vol <- mean(souhrn$rkb_ann_pct)

g8 <- ggplot(souhrn, aes(x = rkb_ann_pct, y = sym_f8, color = sektor)) +
  geom_vline(xintercept = prumer_vol, linetype = "dashed",
             color = "grey35", linewidth = 0.8) +
  geom_segment(aes(x = 0, xend = rkb_ann_pct, y = sym_f8, yend = sym_f8),
               color = "grey82", linewidth = 0.85) +
  geom_point(size = 5.0, alpha = 0.92) +
  geom_text(aes(label = paste0(rkb_ann_pct, " %")),
            hjust = -0.35, size = 3.2, fontface = "bold") +
  scale_color_manual(values = SCOLS, name = "Sektor") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22)),
                     labels = function(x) paste0(x, " %")) +
  facet_grid(sektor ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  labs(
    title    = "Roční anualizovaná RK Bartlett volatilita — 23 akcií · 2024",
    subtitle = paste0(
      "\u03c3_ann = \u221a(252 \u00d7 RKB_denn\u00ed). Průměr vzorku = ",
      sprintf("%.1f %%", prumer_vol),
      ". Svislá přerušovaná čára = průměr celého vzorku.\n",
      "Technologie: NVDA (58,5 %), AMD (42,9 %). ",
      "Defenzivní spotřeba nejnižší: KO (13,5 %), PG (14,1 %). Rozdíl 4,3\u00d7."),
    x = "Anualizovaná RK Bartlett volatilita [% p.a.]",
    y = NULL,
    caption = paste0("\u03c3_ann = \u221a(252 \u00d7 pr\u016fm\u011brn\u00e1 denn\u00ed RKB). ",
                     "RK Bartlett: Barndorff-Nielsen et al. (2008). N = 252 obchodních dní.")
  ) +
  TH %+replace% theme(
    legend.position = "none",
    strip.text.y    = element_text(angle = 0, size = 9.5)
  )

uloz(g8, "Graf08_Rocni_volatilita.png", 13, 11)


# ====================================================================
message("\n", strrep("=", 62))
message("  HOTOVO \u2014 v\u0161ech 8 graf\u016f ulo\u017eeno.")
message("  Slo\u017eka: ", file.path(OUT_DIR, "grafy_final"))
message(strrep("=", 62))
cat(paste0(
  "\n  Graf01  U-tvar: v\u0161echny sektory (zoom od 09:35)\n",
  "  Graf02  U-tvar: Technologie vs. Defenzivn\u00ed (dvoupanel, 95% CI)\n",
  "  Graf03  Hloubka U-tvaru vs. ro\u010dn\u00ed volatilita (scatter)\n",
  "  Graf04  MDH test: r(vol,liq) vs. r(vol,int)\n",
  "  Graf05  Cleveland dot plot: volatilita, intenzita, likvidita\n",
  "  Graf06  R\u00b2 a signifikance OLS modelu (HAC)\n",
  "  Graf07  Heatmapa: sektor \u00d7 \u010das\n",
  "  Graf08  Ro\u010dn\u00ed volatilita po sektorech\n"
))