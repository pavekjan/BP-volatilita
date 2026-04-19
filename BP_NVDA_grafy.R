# =====================================================================
# BAKALÁŘSKÁ PRÁCE
# Analýza intradenních vzorců volatility a likvidity
#
# Metodologie:
#   Primární estimátor: Realized Kernel s Bartlettovým jádrem (RK).
#   Srovnávací estimátor: Pre-averaging (PA).
#   Oba filtrují mikrostrukturní šum, který systematicky nafukuje
#   naivní realizovanou varianci (RV).
#
# Výstupy (složka grafy/):
#   Graf01 — Volatility Signature Plot     Graf07 — Kvartální stabilita
#   Graf02 — U-shape volatility (zoom)     Graf08 — Box plot distribucí
#   Graf03 — Volatilita, intenzita, liq.   Graf09 — Hustotní grafy
#   Graf04 — RK Bartlett vs. naivní RV     Graf10 — Denní časová řada
#   Graf05 — Normalizované U-křivky        Graf11 — Měsíční panely
#   Graf06 — Heatmapa den × čas            Graf12 — Scatter plot
# =====================================================================


# =====================================================================
# ČÁST A — KNIHOVNY A KONFIGURACE
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(highfrequency)   # rKernelCov() = RK estimátor, rMRCov() = PA
  library(xts)             # Formát časových řad vyžadovaný highfrequency
  library(zoo)             # na.locf(): forward-fill mezer v gridu
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(viridis)
})

# ── Parametry analýzy ────────────────────────────────────────────────
SYMBOL        <- "NVDA.O"
ROK           <- 2024
TZ_MKT        <- "America/New_York"
RTH_MIN_START <- 570L          # 09:30 ET = 570 min od půlnoci
RTH_MIN_END   <- 960L          # 16:00 ET
INTERVAL_MIN  <- 5L            # Délka okna v minutách
N_BUCKETS     <- 78L           # Počet 5-min oken za obchodní den
MIN_TRADES    <- 20L           # Min. počet ticků pro výpočet RK/PA
MIN_TICKS_DAY <- 500L          # Min. ticků za den (filtr výpadků)

BASE_DATA_DIR <- "/storage/brno2/home/janpavek03/ondemand/data/sys/dashboard/batch_connect/sys/bc_meta_rstudio/output/49fc93f0-b6db-4442-aca3-f85cbe220cdf/"

OUT_DIR <- file.path(getwd(), "BP_NVDA_vystupy")
dir.create(file.path(OUT_DIR, "grafy"), recursive = TRUE, showWarnings = FALSE)

# ── Barevná paleta (konzistentní napříč všemi grafy) ─────────────────
# Každá veličina má pevně přiřazenou barvu — čtenář nemusí číst
# legendu znova u každého grafu.
COL_RKB <- "#1B6CA8"   # RK Bartlett   — modrá
COL_RV  <- "#E07B39"   # Naivní RV     — oranžová
COL_PA  <- "#9B59B6"   # Pre-averaging — fialová
COL_INT <- "#C0392B"   # Intenzita     — červená
COL_LIQ <- "#27AE60"   # Likvidita     — zelená

# ── Grafický styl ────────────────────────────────────────────────────
# theme_bw() = bílé pozadí, šedá mřížka → vhodné pro tisk v práci.
THEME_BP <- theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(size = 12, face = "bold",
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 9, color = "grey30",
                                    margin = margin(b = 8), lineheight = 1.4),
    plot.caption     = element_text(size = 8, color = "grey45",
                                    hjust = 0, margin = margin(t = 8),
                                    lineheight = 1.4),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8.5),
    legend.key.width = unit(1.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    strip.background = element_rect(fill = "grey96"),
    strip.text       = element_text(size = 9, face = "bold"),
    plot.margin      = margin(8, 12, 6, 8)
  )

uloz_graf <- function(graf, nazev, sirka = 13, vyska = 7) {
  ggsave(file.path(OUT_DIR, "grafy", nazev),
         plot = graf, width = sirka, height = vyska, dpi = 300, bg = "white")
  message("  ✓ ", nazev)
}

# Referenční hodnota rovnoměrného rozdělení:
# pokud by volatilita nezávisela na čase, každý ze 78 intervalů by
# obsahoval přesně 1/78 ≈ 1.28 % denní volatility.
EQ_SHARE <- 100 / N_BUCKETS


# =====================================================================
# ČÁST B — NAČTENÍ A ČIŠTĚNÍ TICK DAT
# =====================================================================
# Tick data = nejpodrobnější forma tržních dat: každý záznam je jeden
# realizovaný obchod s cenou, objemem a přesným časem.
# NVDA 2024: průměrně ~45 000 ticků za obchodní den.
#
# Čištění probíhá v krocích:
#   1) Načtení .rda souborů (1 soubor = 1 obchodní den)
#   2) Filtrace symbolu NVDA.O
#   3) Převod UTC → Eastern Time (ET)
#   4) Ořez na obchodní hodiny 09:30–16:00 ET
#   5) Odstranění neplatných hodnot (záporné ceny, nulové objemy)
#   6) Deduplikace: dva obchody ve stejné sekundě → ponecháme pozdější

nacti_rda <- function(cesta) {
  tryCatch({
    env  <- new.env()
    nazvy <- load(cesta, envir = env)
    if (length(nazvy) == 0) return(NULL)
    as.data.table(env[[nazvy[1]]])
  }, error = function(e) NULL)
}

cistit_ticky <- function(dt, symbol_cil) {
  if (is.null(dt) || !is.data.table(dt) || nrow(dt) == 0) return(data.table())

  setnames(dt, toupper(names(dt)))
  potrebne <- c("TIMESTAMP", "VALUE", "VOLUME", "SYMBOL")
  if (!all(potrebne %in% names(dt))) return(data.table())

  dt <- dt[SYMBOL == symbol_cil]
  if (nrow(dt) == 0) return(data.table())

  # UTC → Eastern Time
  dt[, ts_utc := ymd_hms(TIMESTAMP, tz = "UTC", quiet = TRUE)]
  dt <- dt[!is.na(ts_utc)]
  dt[, ts_ny   := with_tz(ts_utc, TZ_MKT)]
  dt[, date_ny := as.Date(ts_ny, tz = TZ_MKT)]

  # Ořez na pravidelné obchodní hodiny
  dt[, min_od_pulnoci := hour(ts_ny) * 60L + minute(ts_ny)]
  dt <- dt[min_od_pulnoci >= RTH_MIN_START & min_od_pulnoci < RTH_MIN_END]
  if (nrow(dt) == 0) return(data.table())

  dt[, VALUE  := as.numeric(VALUE)]
  dt[, VOLUME := as.numeric(VOLUME)]
  dt <- dt[is.finite(VALUE) & VALUE > 0 & is.finite(VOLUME) & VOLUME >= 0]
  if (nrow(dt) == 0) return(data.table())

  # Deduplikace + log-cena pro výpočet výnosů
  setorder(dt, ts_ny)
  dt <- dt[, .SD[.N], by = ts_ny]
  dt[, logp := log(VALUE)]
  dt[, .(ts_ny, date_ny, VALUE, VOLUME, logp)]
}

message("\n[B] Načítám tick data...")
data_slozka <- file.path(BASE_DATA_DIR, paste0("NASDAQ_NVDA_", ROK))
if (!dir.exists(data_slozka))
  stop("Složka nenalezena: ", data_slozka)

soubory <- list.files(data_slozka, pattern = "\\.rda$", full.names = TRUE)
message("  Souborů: ", length(soubory))

vsechny_ticky <- vector("list", length(soubory))
for (i in seq_along(soubory)) {
  dt_raw <- nacti_rda(soubory[i])
  dt     <- cistit_ticky(dt_raw, SYMBOL)
  if (is.data.table(dt) && nrow(dt) > 0) vsechny_ticky[[i]] <- dt
}

ticky <- rbindlist(Filter(Negate(is.null), vsechny_ticky),
                   use.names = TRUE, fill = TRUE)
setorder(ticky, ts_ny)
rm(vsechny_ticky)

# Vyřazení dní s málo daty (svátky, výpadky burzy)
platne_dny <- ticky[, .(n = .N), by = date_ny][n >= MIN_TICKS_DAY, date_ny]
ticky      <- ticky[date_ny %in% platne_dny]
message("  Platných dní: ", length(platne_dny),
        " | ticků: ", format(nrow(ticky), big.mark = " "))


# =====================================================================
# ČÁST C — VÝPOČET ESTIMÁTORŮ VOLATILITY
# =====================================================================
# Pro každé 5-minutové okno počítáme:
#
#   rv_tick     — Naivní RV = Σ r²ᵢ (log-výnosy tick po ticku).
#                 Trpí mikrostrukturním šumem (bid-ask bounce,
#                 zaokrouhlování cen). Nadhodnocuje skutečnou volatilitu,
#                 zvláště při vysoké obchodní aktivitě ráno.
#
#   rk_bartlett — Realized Kernel s Bartlettovým jádrem:
#                 RK = Σₕ k(h/H) · Σₜ rₜ rₜ₋ₕ
#                 kde k() = Bartlettovo jádro, H = šířka pásma.
#                 Filtruje šum vahou autokovarancí a dává konzistentní
#                 odhad IV bez ohledu na vzorkovací interval.
#                 → highfrequency::rKernelCov()
#
#   rk_parzen   — Totéž s Parzen jádrem (pro sensitivní analýzu).
#
#   pa          — Pre-averaging: výnosy se průměrují přes sousední
#                 pozorování (blok kₙ ≈ θ·√n) před umocněním.
#                 Tím se ruší vliv i.i.d. šumu.
#                 → highfrequency::rMRCov()
#
# Okna s méně než MIN_TRADES ticků dostávají NA — estimátor by měl
# příliš velký rozptyl (nedostatek informace).

jako_skalar <- function(obj) {
  # highfrequency vrací matici nebo seznam → chceme skalár
  if (is.null(obj))   return(NA_real_)
  if (is.matrix(obj)) return(as.numeric(obj[1L, 1L]))
  if (is.numeric(obj)) return(as.numeric(obj[1L]))
  if (is.list(obj))   return(jako_skalar(obj[[1L]]))
  NA_real_
}

vypocti_RK <- function(okno) {
  if (nrow(okno) < MIN_TRADES) return(list(rk_bartlett = NA_real_,
                                           rk_parzen   = NA_real_))
  px <- xts(okno$VALUE, order.by = okno$ts_ny)
  colnames(px) <- "price"
  # H = min(10, n/3): standardní heuristika pro bias/variance kompromis
  H <- min(10L, floor(nrow(okno) / 3L))

  rk_b <- tryCatch(
    jako_skalar(rKernelCov(px, makeReturns = TRUE,
                            kernelType = "bartlett",
                            kernelParam = H, kernelDOFadj = TRUE)),
    error = function(e) NA_real_)
  rk_p <- tryCatch(
    jako_skalar(rKernelCov(px, makeReturns = TRUE,
                            kernelType = "parzen",
                            kernelParam = H, kernelDOFadj = TRUE)),
    error = function(e) NA_real_)

  list(rk_bartlett = ifelse(is.finite(rk_b) & rk_b > 0, rk_b, NA_real_),
       rk_parzen   = ifelse(is.finite(rk_p) & rk_p > 0, rk_p, NA_real_))
}

vypocti_PA <- function(okno) {
  if (nrow(okno) < MIN_TRADES * 2L) return(NA_real_)
  # PA potřebuje rovnoměrný raster → grid 5 s, prázdná místa forward-fill
  t0   <- okno$ts_ny[1L]
  t1   <- okno$ts_ny[nrow(okno)]
  grid <- data.table(ts_ny = seq(t0, t1, by = 5L))
  merged <- merge(grid, okno[, .(ts_ny, VALUE)], by = "ts_ny", all.x = TRUE)
  merged[, VALUE := zoo::na.locf(VALUE, na.rm = FALSE)]
  merged <- merged[!is.na(VALUE)]
  if (nrow(merged) < 10L) return(NA_real_)

  px <- xts(merged$VALUE, order.by = merged$ts_ny)
  colnames(px) <- "PRICE"
  # theta = 0.8: doporučená hodnota dle Jacod et al. (2009)
  pa <- tryCatch(
    jako_skalar(rMRCov(list(asset = px), makePsd = TRUE, theta = 0.8)),
    error = function(e) NA_real_)
  ifelse(is.finite(pa) & pa > 0, pa, NA_real_)
}

vypocti_okno <- function(okno) {
  n    <- nrow(okno)
  rk   <- vypocti_RK(okno)
  pa   <- vypocti_PA(okno)
  list(n_trades    = n,
       volume_sum  = sum(okno$VOLUME, na.rm = TRUE),
       price_last  = okno$VALUE[n],
       rv_tick     = sum(diff(okno$logp)^2, na.rm = TRUE),
       rk_bartlett = rk$rk_bartlett,
       rk_parzen   = rk$rk_parzen,
       pa          = pa)
}

message("\n[C] Počítám estimátory volatility...")
message("    (Výpočet RK a PA je náročný — trvá několik minut.)")

ticky[, bucket := floor_date(ts_ny, unit = paste0(INTERVAL_MIN, " minutes"))]
ticky <- ticky[
  bucket >= as.POSIXct(paste0(date_ny, " 09:30:00"), tz = TZ_MKT) &
    bucket <  as.POSIXct(paste0(date_ny, " 16:00:00"), tz = TZ_MKT)]
setorder(ticky, date_ny, bucket, ts_ny)

results <- ticky[, vypocti_okno(.SD), by = .(date_ny, bucket)]

results[, minutes_from_open := as.integer(difftime(
  bucket,
  as.POSIXct(paste0(date_ny, " 09:30:00"), tz = TZ_MKT),
  units = "mins"))]
results[, `:=`(hour_minute = format(bucket, "%H:%M"), symbol = SYMBOL)]

message("  Oken celkem:     ", format(nrow(results), big.mark = " "))
message("  Oken s RK:       ",
        format(sum(!is.na(results$rk_bartlett)), big.mark = " "),
        " (", round(mean(!is.na(results$rk_bartlett)) * 100, 1), " %)")


# =====================================================================
# ČÁST D — SEZÓNNÍ PODÍLY A AGREGACE
# =====================================================================
# Sezónní podíl eliminuje denní úroveň volatility a izoluje čistě
# intradenní časový vzorec (metodologie dle Andersen & Bollerslev 1997):
#
#   s_{i,t} = RK_{i,t} / RK_{denní,t}
#
# Průměr s̄_i = (1/T) Σₜ s_{i,t} je intradenní sezónní faktor
# pro interval i. Součet s̄_i přes všechny intervaly = 1 (validace).

# RK ≤ 0 je matematicky neplatné (RK je semidefinitně pozitivní)
results[, rk_imp := fifelse(is.na(rk_bartlett) | rk_bartlett <= 0,
                            NA_real_, rk_bartlett)]

# Denní součet RK s korekcí pro chybějící okna:
# Přeškálování: RK_denní = sum(validní) × (78 / počet_validních)
# Předpoklad: chybějící okna mají stejnou průměrnou volatilitu jako
# validní. Tento předpoklad je rozumný — chybějící okna jsou náhodná
# (nízká aktivita v poledne), nikoli systematicky odlehlá.
daily_corr <- results[, {
  n_valid <- sum(!is.na(rk_imp))
  rk_sum  <- sum(rk_imp, na.rm = TRUE)
  .(rkb_day = if (n_valid > 0L) rk_sum * (N_BUCKETS / n_valid) else NA_real_)
}, by = date_ny]

results <- merge(results, daily_corr, by = "date_ny", all.x = TRUE)

daily_trading <- results[, .(total_trades = sum(n_trades, na.rm = TRUE),
                             total_volume = sum(volume_sum, na.rm = TRUE),
                             rv_day       = sum(rv_tick, na.rm = TRUE)),
                         by = date_ny]
results <- merge(results, daily_trading, by = "date_ny", all.x = TRUE)

results[, `:=`(
  rkb_share    = rk_imp / rkb_day,
  rv_share     = rv_tick / rv_day,
  trades_share = n_trades   / total_trades,
  volume_share = volume_sum / total_volume
)]

# Validace: průměr podílů musí ≈ 1/78 (odchylka > 1 % = chyba ve výpočtu)
avg_s <- mean(results$rkb_share, na.rm = TRUE)
message("\n[D] Validace: průměrný podíl = ", round(avg_s, 5),
        " | očekáváno 1/78 = ", round(1/N_BUCKETS, 5),
        " | odchylka = ",
        round(abs(avg_s - 1/N_BUCKETS) / (1/N_BUCKETS) * 100, 3), " %")

# Pomocné proměnné pro grafy
results[, wday_cz := factor(
  wday(as.Date(date_ny), week_start = 1L), 1:5,
  c("Pondělí", "Úterý", "Středa", "Čtvrtek", "Pátek"))]
results[, month_num   := month(as.Date(date_ny))]
results[, quarter_lbl := paste0("Q", quarter(as.Date(date_ny)))]
N_DAYS <- uniqueN(results$date_ny)
message("  Platných obchodních dní: ", N_DAYS)

# Agregace průměrného sezónního profilu (78 intervalů)
pct <- function(x, p) quantile(x, p, na.rm = TRUE) * 100

ushape <- results[, .(
  rkb_mean = mean(rkb_share, na.rm = TRUE) * 100,
  rkb_p10  = pct(rkb_share, 0.10),
  rkb_p25  = pct(rkb_share, 0.25),
  rkb_p50  = pct(rkb_share, 0.50),
  rkb_p75  = pct(rkb_share, 0.75),
  rkb_p90  = pct(rkb_share, 0.90),
  rv_mean  = mean(rv_share,     na.rm = TRUE) * 100,
  int_mean = mean(trades_share, na.rm = TRUE) * 100,
  int_p25  = pct(trades_share,  0.25),
  int_p75  = pct(trades_share,  0.75),
  liq_mean = mean(volume_share, na.rm = TRUE) * 100,
  liq_p25  = pct(volume_share,  0.25),
  liq_p75  = pct(volume_share,  0.75)
), by = .(hour_minute, minutes_from_open)]

setorder(ushape, minutes_from_open)
ushape[, time_f := factor(hour_minute, levels = unique(hour_minute))]

# Tiky osy x každých 30 min (každý 6. ze 78 intervalů)
x_ticky <- ushape$time_f[seq(1L, nrow(ushape), 6L)]


# =====================================================================
# ČÁST E — GRAFY
# =====================================================================

message("\n[E] Generuji grafy...")


# ── Graf 01: Volatility Signature Plot ───────────────────────────────
# Zdůvodňuje volbu RK Bartlett. Při vysoké frekvenci mikrostrukturní
# šum (bid-ask bounce, zaokrouhlování cen) nadhodnocuje naivní RV.
# Čím kratší interval, tím větší nadhodnocení → rostoucí křivka vlevo.

price_serie <- results[order(date_ny, minutes_from_open),
                       .(date_ny, minutes_from_open, price_last)]

rv_ze_sub <- function(dt, k) {
  sub <- dt[minutes_from_open %% (k * 5L) == 0L]
  setorder(sub, date_ny, minutes_from_open)
  sub[, r := c(NA_real_, diff(log(price_last + 1e-12))), by = date_ny]
  daily_rv <- sub[!is.na(r), .(rv = sum(r^2, na.rm = TRUE)), by = date_ny]
  sqrt(252 * mean(daily_rv$rv, na.rm = TRUE)) * 100
}

rv_tick_ann <- sqrt(252 * mean(
  results[, .(rv = sum(rv_tick, na.rm = TRUE)), by = date_ny]$rv,
  na.rm = TRUE)) * 100
rkb_ann <- sqrt(252 * mean(daily_corr$rkb_day, na.rm = TRUE)) * 100

kroky_k   <- c(1L, 2L, 3L, 6L, 12L)
kroky_min <- kroky_k * 5L
subsamp   <- sapply(kroky_k, function(k) rv_ze_sub(price_serie, k))

avg_tick_s <- 23400 / mean(
  results[, .(n = sum(n_trades, na.rm = TRUE)), by = date_ny]$n, na.rm = TRUE)

IV_day <- (rkb_ann / 100)^2 / 252
xi2    <- max(0, ((rv_tick_ann / 100)^2 / 252 - IV_day) /
              (2 * 23400 / avg_tick_s))

delta_th <- exp(seq(log(0.2), log(3600), length.out = 300))
ann_th   <- sqrt(252 * (IV_day + 2 * (23400 / delta_th) * xi2)) * 100
th_dt    <- data.table(delta_s = delta_th, ann_v = ann_th)

emp_dt <- data.table(
  delta_s = c(avg_tick_s, kroky_k * 5L * 60L),
  ann_v   = c(rv_tick_ann, subsamp),
  typ     = c("tick", rep("sub", 5L)),
  popis   = c(sprintf("Tick (~%.1f s)", avg_tick_s), paste0(kroky_min, " min")))
bias_pct <- round((rv_tick_ann / rkb_ann - 1) * 100)

g1 <- ggplot() +
  annotate("rect", xmin = 0.15, xmax = 4200,
           ymin = rkb_ann * 0.95, ymax = rkb_ann * 1.05,
           fill = COL_RKB, alpha = 0.10) +
  geom_hline(yintercept = rkb_ann, linetype = "dashed",
             color = COL_RKB, linewidth = 1.1) +
  geom_line(data = th_dt, aes(x = delta_s, y = ann_v),
            color = "#5B7BAC", linewidth = 2.2, alpha = 0.75) +
  geom_point(data = emp_dt[typ == "sub"], aes(x = delta_s, y = ann_v),
             color = "#1B3A5C", size = 5, shape = 19) +
  geom_text(data = emp_dt[typ == "sub"],
            aes(x = delta_s, y = ann_v, label = popis),
            nudge_y = diff(range(emp_dt$ann_v)) * 0.04,
            size = 3.1, color = "#1B3A5C", fontface = "bold") +
  geom_point(data = emp_dt[typ == "tick"], aes(x = delta_s, y = ann_v),
             color = COL_RV, size = 7, shape = 17) +
  geom_text(data = emp_dt[typ == "tick"],
            aes(x = delta_s, y = ann_v, label = popis),
            nudge_y = diff(range(emp_dt$ann_v)) * 0.04,
            size = 3.1, color = COL_RV, fontface = "bold") +
  annotate("segment",
           x = avg_tick_s * 0.88, xend = avg_tick_s * 0.88,
           y = rv_tick_ann * 0.97, yend = rkb_ann * 1.07,
           arrow = arrow(length = unit(0.25, "cm"), ends = "both",
                         type = "closed"),
           color = "#C0392B", linewidth = 1.0) +
  annotate("label", x = avg_tick_s * 3.5, y = (rv_tick_ann + rkb_ann) / 2,
           label = paste0("+", bias_pct, " % bias\n(mikrostrukturní šum)"),
           size = 3.2, color = "#C0392B", fill = "#FFF2F0",
           label.size = 0.3, hjust = 0, fontface = "bold") +
  annotate("label", x = 3200, y = rkb_ann,
           label = paste0("RK Bartlett: ", round(rkb_ann, 1),
                          " % p.a.  (robustní estimátor)"),
           size = 3.2, color = COL_RKB, fill = "white",
           label.size = 0.3, hjust = 1, fontface = "bold") +
  geom_vline(xintercept = 300, linetype = "dotted",
             color = COL_LIQ, linewidth = 1.0, alpha = 0.75) +
  annotate("text", x = 320, y = rkb_ann * 0.88,
           label = "5 min\n(zvolený\ninterval)",
           size = 2.8, color = COL_LIQ, hjust = 0, lineheight = 1.2) +
  scale_x_log10(
    breaks = c(0.2,0.5,1,2,5,10,30,60,120,300,600,1800,3600),
    labels = c("0.2s","0.5s","1s","2s","5s","10s","30s",
               "1min","2min","5min","10min","30min","60min"),
    expand = expansion(mult = c(0.04, 0.06))) +
  scale_y_continuous(labels = function(x) paste0(round(x), " %"),
                     expand = expansion(mult = c(0.05, 0.10))) +
  labs(
    title    = "Volatility Signature Plot — NVDA · NASDAQ · 2024",
    subtitle = paste0(
      "Naivní tick-by-tick RV nadhodnocuje volatilitu při vysoké frekvenci ",
      "(mikrostrukturní šum: zaokrouhlení cen, bid-ask bounce).\n",
      "RK Bartlett je robustní — hodnota nezávisí na vzorkovacím intervalu. ",
      "Volba 5-min oken zdůvodněna konvergencí k robustnímu odhadu."),
    x = "Vzorkovací interval Δ   ←  vysoká frekvence  ·  nízká frekvence  →",
    y = "Anualizovaná realizovaná volatilita [% p.a.]",
    caption = paste0(
      "Šumový model: E[RV(Δ)] = IV + 2·(T/Δ)·ξ²  (Aït-Sahalia et al. 2011). ",
      "RK Bartlett: Barndorff-Nielsen et al. (2008). N = ", N_DAYS, " dní.")
  ) + THEME_BP + theme(axis.text.x = element_text(angle = 40, hjust = 1))

uloz_graf(g1, "Graf01_Signature_Plot.png", 13, 7.5)


# ── Graf 02: U-shape intradenní sezónnosti volatility ────────────────
# Hlavní výsledkový graf. Tři fáze obchodního dne:
#   Ráno: nejvyšší volatilita (price discovery — trh zpracovává
#         informace nashromážděné přes noc).
#   Poledne: minimum (~0.63 % denní volatility).
#   Závěr: mírný nárůst (~3.93 % při 15:55) → MOC příkazy.
#
# Osa začíná od 09:35 (interval 09:30 = 18.85 % je o řád výše).
# Bez přiblížení by závěrečný nárůst splynul s nulou.

val_0930 <- ushape[minutes_from_open == 0L, rkb_mean]
y_zoom   <- max(ushape[minutes_from_open > 20 & minutes_from_open < 360,
                       rkb_mean] * 2.3, EQ_SHARE * 5)
ush_z    <- ushape[minutes_from_open >= 5L]
ush_z[, tf := factor(hour_minute, levels = unique(hour_minute))]
x_z    <- ush_z$tf[seq(1L, nrow(ush_z), 6L)]
x_1000 <- which(ush_z$hour_minute == "10:00")
x_1430 <- which(ush_z$hour_minute == "14:30")

min_idx   <- which.min(ush_z$rkb_mean)
min_value <- ush_z$rkb_mean[min_idx]
close_max <- ush_z[hour_minute >= "15:00"]
close_idx <- which(ush_z$hour_minute ==
                     close_max[which.max(rkb_mean), hour_minute])

g2 <- ggplot(ush_z, aes(x = tf, y = rkb_mean, group = 1)) +
  geom_ribbon(aes(ymin = rkb_p10, ymax = rkb_p90),
              fill = COL_RKB, alpha = 0.09, color = NA) +
  geom_ribbon(aes(ymin = rkb_p25, ymax = rkb_p75),
              fill = COL_RKB, alpha = 0.22, color = NA) +
  geom_hline(yintercept = EQ_SHARE, linetype = "longdash",
             color = "grey55", linewidth = 0.55) +
  annotate("text", x = round(nrow(ush_z) * 0.50),
           y = EQ_SHARE + y_zoom * 0.022,
           label = sprintf("Rovnoměrné rozdělení = %.2f %%", EQ_SHARE),
           size = 2.9, color = "grey45", hjust = 0.5) +
  geom_vline(xintercept = x_1000, linetype = "dotdash",
             color = "grey58", linewidth = 0.45) +
  geom_vline(xintercept = x_1430, linetype = "dotdash",
             color = "grey58", linewidth = 0.45) +
  annotate("text", x = x_1000 + 0.4, y = y_zoom * 0.80,
           label = "10:00\n(makro)", size = 2.7, color = "grey40",
           hjust = 0, lineheight = 1.2) +
  annotate("text", x = x_1430 + 0.4, y = y_zoom * 0.80,
           label = "14:30\n(Fed)", size = 2.7, color = "grey40",
           hjust = 0, lineheight = 1.2) +
  geom_line(color = COL_RKB, linewidth = 1.3) +
  geom_line(aes(y = rkb_p50), color = COL_RKB, linewidth = 0.9,
            linetype = "dashed", alpha = 0.65) +
  annotate("text", x = 1.5, y = y_zoom * 0.96,
           label = sprintf("Interval 09:30: %.2f %%\n(mimo zobrazenou osu)",
                           val_0930),
           hjust = 0, vjust = 1, size = 2.9, color = "grey35",
           fontface = "italic") +
  annotate("text", x = min_idx + 1, y = min_value + y_zoom * 0.08,
           label = sprintf("Polední minimum\n~%.2f %%", min_value),
           size = 2.8, color = "grey35", hjust = 0, lineheight = 1.2) +
  annotate("text",
           x = close_idx - 2,
           y = close_max[which.max(rkb_mean), rkb_mean] + y_zoom * 0.06,
           label = sprintf("Závěrečný nárůst\n~%.2f %%\n(MOC příkazy)",
                           close_max[which.max(rkb_mean), rkb_mean]),
           size = 2.8, color = COL_RKB, hjust = 1,
           lineheight = 1.2, fontface = "bold") +
  scale_x_discrete(breaks = x_z) +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), " %"),
                     limits = c(0, y_zoom),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title    = "U-shape intradenní sezónnosti volatility — NVDA · 2024 (RK Bartlett)",
    subtitle = paste0(
      "s̄ᵢ = (1/T) Σₜ (RKᵢₜ / RK_denní). Přiblížená osa od 09:35. ",
      "Tmavší pás = P25–P75 (IQR). Světlý pás = P10–P90. N = ", N_DAYS, " dní."),
    x = "Čas (Eastern Standard Time)",
    y = "RK Bartlett podíl na denní volatilitě [%]",
    caption = paste0("Sezónní podíly: Andersen & Bollerslev (1997). ",
                     "RK Bartlett: Barndorff-Nielsen et al. (2008).")
  ) + THEME_BP

uloz_graf(g2, "Graf02_Ushape_Zoom.png")


# ── Graf 03: Porovnání intradenních vzorců ────────────────────────────
# Všechny tři veličiny sdílejí U-tvar, ale s důležitými rozdíly:
#   Volatilita: nejhlubší polední minimum.
#   Likvidita: dvojité maximum (velké blokové obchody ráno → Kyle 1985;
#              institucionální MOC příkazy večer).
#   Intenzita: největší nárůst až při závěru, ne ráno → průměrná
#              velikost obchodu ráno >> večer (Admati & Pfleiderer 1988).

il_long <- rbindlist(list(
  ushape[, .(time_f, mean = rkb_mean, p25 = rkb_p25, p75 = rkb_p75,
             serie = "Volatilita (RK Bartlett)")],
  ushape[, .(time_f, mean = int_mean, p25 = int_p25, p75 = int_p75,
             serie = "Intenzita (počet obchodů)")],
  ushape[, .(time_f, mean = liq_mean, p25 = liq_p25, p75 = liq_p75,
             serie = "Likvidita (objem)")]
))

g3 <- ggplot(il_long, aes(x = time_f, y = mean, group = serie,
                          color = serie)) +
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = serie),
              alpha = 0.10, color = NA, show.legend = FALSE) +
  geom_hline(yintercept = EQ_SHARE, linetype = "longdash",
             color = "grey65", linewidth = 0.5) +
  annotate("text", x = round(nrow(ushape) * 0.50), y = EQ_SHARE + 0.08,
           label = sprintf("Rovnoměrné rozdělení = %.2f %%", EQ_SHARE),
           size = 2.7, color = "grey45", hjust = 0.5) +
  geom_vline(xintercept = which(ushape$hour_minute == "10:00"),
             linetype = "dotdash", color = "grey65", linewidth = 0.4) +
  geom_vline(xintercept = which(ushape$hour_minute == "14:30"),
             linetype = "dotdash", color = "grey65", linewidth = 0.4) +
  annotate("text", x = which(ushape$hour_minute == "10:00") + 0.4,
           y = 3.30, label = "10:00", size = 2.6, color = "grey45", hjust = 0) +
  annotate("text", x = which(ushape$hour_minute == "14:30") + 0.4,
           y = 3.30, label = "14:30", size = 2.6, color = "grey45", hjust = 0) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    values = c("Volatilita (RK Bartlett)"  = COL_RKB,
               "Intenzita (počet obchodů)" = COL_INT,
               "Likvidita (objem)"         = COL_LIQ), name = NULL) +
  scale_fill_manual(
    values = c("Volatilita (RK Bartlett)"  = COL_RKB,
               "Intenzita (počet obchodů)" = COL_INT,
               "Likvidita (objem)"         = COL_LIQ)) +
  scale_x_discrete(breaks = x_ticky) +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), " %"),
                     expand = expansion(mult = c(0.01, 0.08))) +
  labs(
    title    = "Porovnání intradenních vzorců — Volatilita, Intenzita, Likvidita — NVDA · 2024",
    subtitle = paste0("Průměrný podíl každé veličiny na denním součtu (09:35–15:55). ",
                      "Pásy = P25–P75 (IQR). N = ", N_DAYS, " dní."),
    x = "Čas (Eastern Standard Time)", y = "Podíl na denním součtu [%]",
    caption = paste0(
      "Volatilita: nejhlubší polední minimum. Likvidita: vrcholí při 09:30 i 15:55 (Kyle 1985). ",
      "Intenzita: nejsilnější nárůst při závěru — průměrná velikost obchodu ráno > večer ",
      "(Admati & Pfleiderer 1988).\nSezónní podíly: Andersen & Bollerslev (1997).")
  ) + THEME_BP

uloz_graf(g3, "Graf03_Porovnani_Vzorcu.png")


# ── Graf 04: Srovnání estimátorů RK Bartlett vs. naivní RV ───────────
# Šrafovaná plocha = intradenní průběh mikrostrukturního biasu.
# Bias je největší ráno (nejvíce ticků za nejkratší čas → šum dominuje).
# Ke konci dopoledne se křivky sbližují (klesající tick-frekvence).

g4 <- ggplot(ushape, aes(x = time_f, group = 1)) +
  geom_ribbon(aes(ymin = rkb_mean, ymax = rv_mean),
              fill = COL_RV, alpha = 0.18, color = NA) +
  geom_hline(yintercept = EQ_SHARE, linetype = "longdash",
             color = "grey65", linewidth = 0.5) +
  geom_vline(xintercept = which(ushape$hour_minute == "10:00"),
             linetype = "dotdash", color = "grey65", linewidth = 0.4) +
  annotate("text", x = which(ushape$hour_minute == "10:00") + 0.5,
           y = max(ushape$rv_mean, na.rm = TRUE) * 0.88,
           label = "10:00", size = 2.7, color = "grey45", hjust = 0) +
  annotate("label", x = 3, y = max(ushape$rv_mean, na.rm = TRUE) * 0.90,
           label = sprintf("Bias při 09:30:\nRV nadhodnocuje o ~%d %%",
                           bias_pct),
           size = 3.0, color = COL_RV, fill = "#FFF5F0",
           label.size = 0.3, hjust = 0, fontface = "bold") +
  annotate("text", x = which(ushape$hour_minute == "12:30"),
           y = EQ_SHARE + 0.12, label = "Konvergence\nv poledne",
           size = 2.6, color = "grey40", hjust = 0.5, lineheight = 1.2) +
  geom_line(aes(y = rv_mean,  color = "Naivní RV (tick-by-tick)"),
            linewidth = 1.5) +
  geom_line(aes(y = rkb_mean, color = "RK Bartlett (robustní)"),
            linewidth = 1.5) +
  scale_color_manual(
    values = c("Naivní RV (tick-by-tick)" = COL_RV,
               "RK Bartlett (robustní)"   = COL_RKB), name = "Estimátor") +
  scale_x_discrete(breaks = x_ticky) +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), " %"),
                     expand = expansion(mult = c(0.01, 0.08))) +
  labs(
    title    = "Srovnání estimátorů volatility — RK Bartlett vs. Naivní RV — NVDA · 2024",
    subtitle = paste0(
      "Průměrný podíl na denní volatilitě pro každý 5-min interval (09:35–15:55). ",
      "Šrafovaná plocha = mikrostrukturní bias = RV − RK. N = ", N_DAYS, " dní."),
    x = "Čas (Eastern Standard Time)", y = "Podíl na denní volatilitě [%]",
    caption = paste0(
      "Bias je největší při 09:30 (+", bias_pct, " %) — nejvyšší tick-frekvence (bid-ask bounce). ",
      "RK Bartlett a PA filtrují šum a konvergují ke stejné hodnotě.\n",
      "RK Bartlett: Barndorff-Nielsen et al. (2008). Sezónní podíly: Andersen & Bollerslev (1997).")
  ) + THEME_BP +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5)))

uloz_graf(g4, "Graf04_Srovnani_Estimatoru.png")


# ── Graf 05: Normalizované srovnání tvarů U-křivek ───────────────────
# Každá série normalizována na vlastní maximum = 100 %, aby bylo možné
# srovnat tvary bez zkreslení absolutními hodnotami.
# Klíčový výsledek: volatilita klesá nejhlouběji (~7 % maxima v poledne).
# Intenzita se při 15:55 vrací na 100 %, volatilita jen na ~18 % —
# závěrečné obchodování je aktivní, ale informačně nenáročné.

nm <- function(x) x / max(x, na.rm = TRUE) * 100

norm_dt <- rbindlist(list(
  ushape[, .(time_f, val = nm(rkb_mean),
             serie = "Volatilita — RK Bartlett  (max: 09:30)")],
  ushape[, .(time_f, val = nm(int_mean),
             serie = "Intenzita — počet obchodů  (max: 15:55)")],
  ushape[, .(time_f, val = nm(liq_mean),
             serie = "Likvidita — objem  (max: 09:30)")]
))

g5 <- ggplot(norm_dt, aes(x = time_f, y = val, group = serie,
                          color = serie, linetype = serie)) +
  geom_hline(yintercept = 100 / 78 / max(ushape$rkb_mean) * 100,
             linetype = "longdash", color = "grey65", linewidth = 0.5) +
  geom_vline(xintercept = which(ushape$hour_minute %in% c("10:00","14:30")),
             linetype = "dotdash", color = "grey72", linewidth = 0.35) +
  annotate("text", x = which(ushape$hour_minute == "12:30"), y = 12,
           label = "Polední\nminimum", size = 2.7, color = "grey38",
           hjust = 0.5, lineheight = 1.2) +
  annotate("text", x = 1.5, y = 98, label = "max\n09:30",
           size = 2.6, color = COL_LIQ, hjust = 0, lineheight = 1.2) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(
    values = c("Volatilita — RK Bartlett  (max: 09:30)"  = COL_RKB,
               "Intenzita — počet obchodů  (max: 15:55)" = COL_INT,
               "Likvidita — objem  (max: 09:30)"          = COL_LIQ),
    name = NULL) +
  scale_linetype_manual(values = c("solid","dashed","dotdash"), name = NULL) +
  scale_x_discrete(breaks = x_ticky) +
  scale_y_continuous(labels = function(x) paste0(round(x), " %"),
                     expand = expansion(mult = c(0.01, 0.05))) +
  labs(
    title    = "Normalizované srovnání tvarů U-křivek — NVDA · 2024",
    subtitle = paste0(
      "Každá série normalizována samostatně: maximum = 100 %.\n",
      "Volatilita a likvidita vrcholí při 09:30; intenzita při 15:55 ",
      "— silný závěrečný spike MOC příkazů."),
    x = "Čas (Eastern Standard Time)",
    y = "Normalizovaná hodnota [% maxima dané série]",
    caption = paste0(
      "Normalizace: y → y / max(y) × 100. Volatilita klesá nejhlouběji ",
      "(~7 % maxima) → nejsilnější sezónní vzorec ze všech tří veličin.\n",
      "Andersen & Bollerslev (1997). RK Bartlett: Barndorff-Nielsen et al. (2008). ",
      "N = ", N_DAYS, " dní.")
  ) + THEME_BP

uloz_graf(g5, "Graf05_Normalizovane_Srovnani.png")


# ── Graf 06: Heatmapa relativní volatility: den × čas ────────────────
# Index = průměrný RK Bartlett podíl (den × 30-min blok) / celkový průměr.
# Hodnota 1.0 = roční průměr. Viditelné:
#   - Pondělní a čtvrteční otevření: výrazně nadprůměrné (víkendové info,
#     čtvrtletní výsledky).
#   - Středa: konzistentně nejklidnější den (midweek liquidity concentration).

results[, abs_minute := 570L + as.integer(minutes_from_open)]
results[, blok30 := paste0(
  sprintf("%02d", abs_minute %/% 60L), ":",
  sprintf("%02d", (abs_minute %% 60L) %/% 30L * 30L))]

blk_lvl <- c("09:30","10:00","10:30","11:00","11:30","12:00","12:30",
             "13:00","13:30","14:00","14:30","15:00","15:30")

heat_dt <- results[!is.na(rkb_share),
                   .(rkb_rel = mean(rkb_share, na.rm = TRUE)),
                   by = .(wday_cz, blok30)]
gm <- mean(heat_dt$rkb_rel, na.rm = TRUE)
heat_dt[, idx := rkb_rel / gm]
heat_dt[, blok30 := factor(blok30, levels = blk_lvl)]
heat_dt <- heat_dt[!is.na(blok30)]

g6 <- ggplot(heat_dt, aes(x = blok30, y = wday_cz, fill = idx)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", idx),
                color = ifelse(idx > 1.7 | idx < 0.55, "white", "#1a1a1a")),
            size = 2.8, fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#313695","#4575b4","#74add1","#abd9e9","#f7f7f7",
                "#fee090","#fdae61","#f46d43","#d73027","#a50026"),
    values  = scales::rescale(c(0.2,0.45,0.65,0.85,1.0,1.15,1.35,1.65,2.2,3.5)),
    name    = "Index relativní volatility\n(1.0 = průměr)",
    guide   = guide_colorbar(barwidth = 14, barheight = 0.7,
                             title.position = "top", title.hjust = 0.5)) +
  scale_color_identity() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev(levels(heat_dt$wday_cz))) +
  labs(
    title    = "Heatmapa relativní volatility — NVDA · 2024 (RK Bartlett, 30-min bloky)",
    subtitle = paste0(
      "Index = průměrný RK Bartlett podíl (den × blok) / celkový průměr. ",
      "1.0 = průměr.  > 1.0 = nadprůměr;  < 1.0 = podprůměr.  N = ", N_DAYS, " dní."),
    x = "30minutový blok (začátek, ET)", y = NULL,
    caption = paste0(
      "Páteční otevírací spike: kombinace NFP/CPI makrodat (8:30 ET) a weekly options expiry.   ",
      "Středa: nejnižší průměrná volatilita.\n",
      "RK Bartlett: Barndorff-Nielsen et al. (2008).   ",
      "Sezónní podíly: Andersen & Bollerslev (1997).")
  ) + THEME_BP +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

uloz_graf(g6, "Graf06_Heatmapa_Den_Cas.png", 14, 5.5)


# ── Graf 07: Kvartální stabilita U-shape vzorce ───────────────────────
# Pokud je U-tvar přítomen ve všech čtyřech kvartálech, jde o trvalou
# strukturální vlastnost trhu, ne o náhodný artefakt pár dní.
# Otevírací hodnoty: Q4 > Q2 > Q1 > Q3 (year-end rebalancing vs. letní klid).

ush_q <- results[, .(rkb_mean = mean(rkb_share, na.rm = TRUE) * 100),
                 by = .(hour_minute, minutes_from_open, quarter_lbl)]
setorder(ush_q, quarter_lbl, minutes_from_open)
ush_q[, tf := factor(hour_minute, levels = unique(ushape$hour_minute))]
q_vals <- ush_q[minutes_from_open == 0L, .(quarter_lbl, rkb_mean)]
Q_COLS <- c("Q1"="#1B6CA8","Q2"="#27AE60","Q3"="#E07B39","Q4"="#9B59B6")

pA <- ggplot(ush_q, aes(x = tf, y = rkb_mean,
                        group = quarter_lbl, color = quarter_lbl)) +
  geom_hline(yintercept = EQ_SHARE, linetype = "longdash",
             color = "grey65", linewidth = 0.5) +
  geom_line(linewidth = 0.95, alpha = 0.90) +
  annotate("text", x = 1.5, y = 19,
           label = paste0("09:30 oříznut\n(",
             paste(q_vals[order(quarter_lbl),
                          sprintf("%s=%.1f %%", quarter_lbl, rkb_mean)],
                   collapse = ", "), ")"),
           size = 2.5, color = "grey35", hjust = 0, fontface = "italic") +
  scale_color_manual(values = Q_COLS, name = "Kvartál") +
  scale_x_discrete(breaks = x_ticky) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), " %"),
                     expand = expansion(mult = c(0.01, 0.08))) +
  labs(title = "A  Plná škála (09:30–15:55)", x = NULL,
       y = "RK Bartlett podíl na denní volatilitě [%]") +
  THEME_BP + theme(plot.title = element_text(size = 10), legend.position = "none")

ush_qz <- ush_q[minutes_from_open >= 5L]
ush_qz[, tf := factor(hour_minute, levels = unique(hour_minute))]
x_qz <- ush_qz$tf[seq(1L, nrow(ush_qz[quarter_lbl == "Q1"]), 6L)]

pB <- ggplot(ush_qz, aes(x = tf, y = rkb_mean,
                         group = quarter_lbl, color = quarter_lbl)) +
  geom_hline(yintercept = EQ_SHARE, linetype = "longdash",
             color = "grey65", linewidth = 0.5) +
  annotate("text",
           x = round(nrow(ush_qz[quarter_lbl == "Q1"]) * 0.50),
           y = EQ_SHARE + 0.12,
           label = sprintf("Rovnoměrné = %.2f %%", EQ_SHARE),
           size = 2.6, color = "grey50") +
  annotate("text", x = nrow(ush_qz[quarter_lbl == "Q1"]) - 1, y = 1.9,
           label = "Q4 year-end\nrebalancing", size = 2.4,
           color = Q_COLS["Q4"], hjust = 1, lineheight = 1.2) +
  geom_line(linewidth = 1.1, alpha = 0.90) +
  scale_color_manual(values = Q_COLS, name = "Kvartál") +
  scale_x_discrete(breaks = x_qz) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), " %"),
                     limits = c(0, 9),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(title = "B  Zoom od 09:35 — viditelný U-tvar a kvartální rozdíly",
       x = "Čas (Eastern Standard Time)",
       y = "RK Bartlett podíl na denní volatilitě [%]") +
  THEME_BP + theme(plot.title = element_text(size = 10), legend.position = "bottom")

g7 <- pA + pB +
  plot_layout(widths = c(1, 1.8), guides = "collect") &
  theme(legend.position = "bottom") &
  plot_annotation(
    title    = "Stabilita U-shape vzorce podle kvartálu — NVDA · 2024 (RK Bartlett)",
    subtitle = paste0(
      "Průměrný podíl RK volatility pro každý 5-min slot, zvlášť pro Q1–Q4. ",
      "Panel A = plná škála. Panel B = zoom (09:35–15:55).\n",
      "Konzistence U-tvaru napříč kvartály potvrzuje systematický (sezónní) charakter vzorce."),
    caption  = paste0(
      "NVDA earnings: Q1=21. 2., Q2=22. 5., Q3=28. 8., Q4=20. 11. 2024. ",
      "Q4 silnější závěrečný nárůst = year-end portfolio rebalancing.\n",
      "RK Bartlett: Barndorff-Nielsen et al. (2008). N = 63 dní/kvartál."),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey30", lineheight = 1.4),
      plot.caption  = element_text(size = 8, color = "grey45",
                                   hjust = 0, lineheight = 1.4)
    )
  )

uloz_graf(g7, "Graf07_Kvartal_final.png", 18, 7)


# ── Graf 08: Box plot distribucí pro vybrané 5-min intervaly ─────────
# U-tvar přítomný v mediánech (ne jen průměrech) = vzorec není efektem
# pár extrémních dní

casy_A <- c("09:30","09:35","09:40","09:45","09:50","09:55","10:00")
casy_B <- c("10:30","11:00","12:00","13:00","14:00",
            "15:30","15:35","15:40","15:45","15:50","15:55")

box_raw <- results[hour_minute %in% c(casy_A,casy_B) & !is.na(rkb_share),
                   .(hour_minute, val = rkb_share * 100)]
box_raw[, faze := fcase(
  hour_minute %in% c("09:30","09:35","09:40","09:45",
                     "09:50","09:55","10:00"), "Otevírací (09:30–10:00)",
  hour_minute %in% c("10:30","11:00"),          "Dopoledne (10:30–11:00)",
  hour_minute %in% c("12:00","13:00"),          "Polední (12:00–13:00)",
  hour_minute %in% c("14:00"),                  "Odpoledne (14:00)",
  default                                     = "Závěr (15:00–15:55)")]

FAZE_COLS <- c("Otevírací (09:30–10:00)"="#1B6CA8",
               "Dopoledne (10:30–11:00)"="#74B3D8",
               "Polední (12:00–13:00)"  ="#27AE60",
               "Odpoledne (14:00)"      ="#E07B39",
               "Závěr (15:00–15:55)"    ="#C0392B")
med_dt <- box_raw[, .(med = median(val, na.rm = TRUE)), by = hour_minute]

dat_A <- box_raw[hour_minute %in% casy_A]
dat_A[, hm := factor(hour_minute, levels = casy_A)]
med_A <- med_dt[hour_minute %in% casy_A]
med_A[, hm := factor(hour_minute, levels = casy_A)]
y_A   <- quantile(dat_A$val, 0.98, na.rm = TRUE) * 1.12

bA <- ggplot(dat_A, aes(x = hm, y = val, fill = faze)) +
  geom_hline(yintercept = EQ_SHARE, linetype = "dashed",
             color = "grey55", linewidth = 0.6) +
  geom_boxplot(outlier.size = 1.2, outlier.alpha = 0.35,
               outlier.color = "grey45", width = 0.65,
               color = "grey20", linewidth = 0.4) +
  geom_text(data = med_A,
            aes(x = hm, y = med, label = paste0(round(med, 2), " %")),
            inherit.aes = FALSE, vjust = -0.55,
            size = 2.8, fontface = "bold", color = "grey15") +
  scale_fill_manual(values = FAZE_COLS) +
  coord_cartesian(ylim = c(0, y_A)) +
  scale_y_continuous(labels = function(x) paste0(round(x), " %")) +
  labs(title    = "A  Otevírací hodina (09:30–10:00)",
       subtitle = "Vlastní škála — detail price discovery",
       x = NULL, y = "RK Bartlett podíl na denní volatilitě [%]") +
  THEME_BP + theme(legend.position = "none")

dat_B <- box_raw[hour_minute %in% casy_B]
dat_B[, hm := factor(hour_minute, levels = casy_B)]
med_B <- med_dt[hour_minute %in% casy_B]
med_B[, hm := factor(hour_minute, levels = casy_B)]
y_B   <- quantile(dat_B$val, 0.995, na.rm = TRUE) * 1.15

bB <- ggplot(dat_B, aes(x = hm, y = val, fill = faze)) +
  annotate("rect",
           xmin = which(casy_B == "15:30") - 0.5,
           xmax = length(casy_B) + 0.5,
           ymin = 0, ymax = y_B, fill = "#C0392B", alpha = 0.05) +
  geom_vline(xintercept = which(casy_B == "15:30") - 0.5,
             linetype = "dotdash", color = "#C0392B",
             linewidth = 0.6, alpha = 0.65) +
  annotate("text", x = which(casy_B == "15:30") + 0.1, y = y_B * 0.91,
           label = "Závěrečná\nhodina", size = 2.6, color = "#C0392B",
           hjust = 0, lineheight = 1.2) +
  geom_hline(yintercept = EQ_SHARE, linetype = "dashed",
             color = "grey55", linewidth = 0.6) +
  geom_boxplot(outlier.size = 1.0, outlier.alpha = 0.30,
               outlier.color = "grey45", width = 0.68,
               color = "grey20", linewidth = 0.4) +
  geom_text(data = med_B,
            aes(x = hm, y = med, label = paste0(round(med, 2), " %")),
            inherit.aes = FALSE, vjust = -0.5,
            size = 2.4, fontface = "bold", color = "grey15") +
  scale_fill_manual(values = FAZE_COLS) +
  coord_cartesian(ylim = c(0, y_B)) +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), " %")) +
  labs(title    = "B  Zbytek dne (10:30–15:55)",
       subtitle = "Zoom — viditelný U-tvar: polední minimum + závěrečný nárůst",
       x = "5-minutový interval (ET)",
       y = "RK Bartlett podíl na denní volatilitě [%]") +
  THEME_BP + theme(legend.position = "none")

g8 <- bA + bB +
  plot_layout(widths = c(1.1, 2.8)) &
  plot_annotation(
    title    = "Distribuce RK volatility pro vybrané 5-min intervaly — NVDA · 2024",
    subtitle = paste0(
      "Každý box = mezidenní distribuce podílu RK Bartlett volatility přes ",
      N_DAYS, " dní. Číslo = medián. Box = IQR. Vousy = 1.5×IQR.\n",
      "Tečky = extrémní dny (earnings, makrošoky). ",
      "Přerušovaná čára = rovnoměrné rozdělení (",
      round(EQ_SHARE, 2), " %)."),
    caption  = paste0(
      "Panel A: 09:30 má medián ~19 % a obrovský IQR — extrémní variabilita ",
      "způsobená earnings dny.\n",
      "RK Bartlett: Barndorff-Nielsen et al. (2008). N = ", N_DAYS, " dní."),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey30", lineheight = 1.4),
      plot.caption  = element_text(size = 8, color = "grey45",
                                   hjust = 0, lineheight = 1.4)
    )
  )

uloz_graf(g8, "Graf08_BoxPlot_Distribuce.png", 18, 7.5)


# ── Graf 09: Hustotní grafy distribuce podle fáze dne ────────────────
# KDE (Kernel Density Estimation) místo histogramu: plynulejší tvar
# lépe odhaluje pravý tail (otevírací hodina = earnings efekt).
# Tři panely: otevření (pravostranně zešikmené) / poledne (symetrické,
# úzké) / závěr (posunuto doprava, bez výrazného tailu).

dens_dt <- results[!is.na(rkb_share),
                   .(minutes_from_open, val = rkb_share * 100)]
dens_dt[, faze := fcase(
  minutes_from_open >= 0   & minutes_from_open < 60,
  "A  Otevírací hodina (09:30–10:29)",
  minutes_from_open >= 120 & minutes_from_open < 240,
  "B  Polední hodiny (11:30–13:29)",
  minutes_from_open >= 330 & minutes_from_open <= 385,
  "C  Závěrečná hodina (15:00–15:55)")]
dens_dt <- dens_dt[!is.na(faze)]

DENS_COLS <- c("A  Otevírací hodina (09:30–10:29)" = COL_RKB,
               "B  Polední hodiny (11:30–13:29)"   = COL_LIQ,
               "C  Závěrečná hodina (15:00–15:55)" = COL_INT)
med_dens <- dens_dt[, .(med = median(val, na.rm = TRUE)), by = faze]

g9 <- ggplot(dens_dt, aes(x = val, fill = faze, color = faze)) +
  geom_density(alpha = 0.30, linewidth = 0.9, bw = "SJ") +
  geom_vline(data = med_dens, aes(xintercept = med, color = faze),
             linetype = "dashed", linewidth = 0.9, show.legend = FALSE) +
  geom_vline(xintercept = EQ_SHARE, linetype = "dotdash",
             color = "grey35", linewidth = 0.7) +
  annotate("text", x = EQ_SHARE + 0.05, y = Inf,
           label = "Rovnoměrné\nrozdělení", vjust = 1.4, size = 2.6,
           color = "grey40", hjust = 0) +
  facet_wrap(~ faze, scales = "free", nrow = 1) +
  scale_fill_manual(values  = DENS_COLS, guide = "none") +
  scale_color_manual(values = DENS_COLS, guide = "none") +
  scale_x_continuous(labels = function(x) paste0(round(x, 1), " %")) +
  labs(
    title    = "Hustotní grafy: distribuce RK volatility podle fáze dne — NVDA · 2024",
    subtitle = paste0(
      "Každý panel = distribuce podílu RK Bartlett volatility přes ", N_DAYS,
      " dní zvlášť pro tři fáze. Přerušovaná čára = medián. ",
      "Tečkovaná = rovnoměrné rozdělení (1.28 %)."),
    x = "RK Bartlett podíl na denní volatilitě [%]", y = "Hustota",
    caption = paste0(
      "Panel A: pravostranně zkosená — velká variabilita (normální vs. earnings dny). ",
      "Panel B: úzká, symetrická — polední klid.\n",
      "Panel C: posunuta doprava — závěrečný rebalancing (MOC příkazy). ",
      "N = ", N_DAYS, " dní. RK Bartlett: Barndorff-Nielsen et al. (2008).")
  ) + THEME_BP +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.text  = element_text(size = 8.5))

uloz_graf(g9, "Graf09_Hustota_FazeDne.png", 14, 6.5)


# ── Graf 10: Denní časová řada realizované volatility ────────────────
# σ_ann = √(252 × RK_denní): anualizovaná denní volatilita.
# Klouzavý průměr 20 dní = volatility clustering (shlukování volatility).
# Osa y oříznutá na 95. percentil — earnings spiky nedeformují škálu,
# ale jsou zachyceny jako zkrácené sloupce.

daily_ts <- results[, .(
  date   = unique(as.Date(date_ny)),
  rv_ann  = sqrt(252 * sum(rv_tick, na.rm = TRUE)) * 100,
  rkb_ann = sqrt(252 * sum(rk_imp,  na.rm = TRUE)) * 100
), by = date_ny]
daily_ts[, date := as.Date(date_ny)]
setorder(daily_ts, date)
daily_ts[, rkb_roll20 := frollmean(rkb_ann, 20L, align = "right", na.rm = TRUE)]
avg_rkb <- mean(daily_ts$rkb_ann, na.rm = TRUE)
y_cap   <- quantile(daily_ts$rv_ann, 0.95, na.rm = TRUE) * 1.05

earnings_dt <- data.table(
  datum = as.Date(c("2024-02-21","2024-05-22","2024-08-28","2024-11-20")),
  popis = c("Earnings\n21. Feb","Earnings\n22. May",
            "Earnings\n28. Aug","Earnings\n20. Nov"))

g10 <- ggplot(daily_ts, aes(x = date)) +
  geom_col(aes(y = pmin(rv_ann, y_cap)), fill = "grey78", alpha = 0.65, width = 1) +
  geom_hline(yintercept = avg_rkb, linetype = "dashed",
             color = "grey35", linewidth = 0.65) +
  annotate("text", x = min(daily_ts$date) + 7,
           y = avg_rkb + y_cap * 0.028,
           label = paste0("Roční průměr RK: ", round(avg_rkb, 1), " % p.a."),
           size = 2.9, color = "grey30", hjust = 0) +
  geom_vline(data = earnings_dt, aes(xintercept = datum),
             linetype = "dashed", color = "#C0392B",
             linewidth = 1.0, alpha = 0.75) +
  geom_text(data = earnings_dt,
            aes(x = datum, y = y_cap * 0.92, label = popis),
            size = 2.7, color = "#C0392B", hjust = -0.08, lineheight = 1.2) +
  geom_line(aes(y = pmin(rkb_ann, y_cap)),
            color = COL_RKB, linewidth = 0.5, alpha = 0.80) +
  geom_line(aes(y = rkb_roll20),
            color = COL_RV, linewidth = 1.8, na.rm = TRUE) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = function(x) paste0(round(x), " %"),
                     limits = c(0, y_cap),
                     expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title    = "Denní anualizovaná realizovaná volatilita — NVDA · NASDAQ · 2024",
    subtitle = paste0(
      "Sloupce: tick-by-tick RV (šedá). Modrá: RK Bartlett (robustní, nikdy záporný). ",
      "Oranžová: klouzavý průměr 20 dní.  σ_ann = √(252 × RV_denní).  ",
      "Roční průměr RK: ", round(avg_rkb, 1), " % p.a."),
    x = NULL, y = "Anualizovaná realizovaná volatilita [% p.a.]",
    caption = paste0(
      "Hodnoty nad 99. percentilem oříznuty — spiky při earnings zachyceny jako zkrácené sloupce. ",
      "RK Bartlett je vždy ≥ 0 (matematická vlastnost estimátoru).\n",
      "Červené čáry: NVDA earnings (21. Feb, 22. May, 28. Aug, 20. Nov 2024). ",
      "N = ", N_DAYS, " dní.  RK Bartlett: Barndorff-Nielsen et al. (2008).")
  ) + THEME_BP +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

uloz_graf(g10, "Graf10_DenniRV_TimeSeries.png", 15, 6)


# ── Graf 11: Měsíční panely intradenního profilu ──────────────────────
# Dva účely:
#   (1) Ověřit konzistenci U-tvaru po celý rok (12 panelů).
#   (2) Earnings efekt: únor, květen, srpen, listopad mají výrazně
#       vyšší ranní hodnoty oproti ostatním měsícům.

ush_m <- results[, .(rkb_mean = mean(rkb_share, na.rm = TRUE) * 100),
                 by = .(hour_minute, minutes_from_open, month_num)]
setorder(ush_m, month_num, minutes_from_open)
ush_m[, tf    := factor(hour_minute, levels = unique(ushape$hour_minute))]
ush_m[, mesic := factor(month.abb[month_num], levels = month.abb)]

g11 <- ggplot(ush_m, aes(x = tf, y = rkb_mean, group = 1)) +
  geom_area(fill = COL_RKB, alpha = 0.15, color = NA) +
  geom_line(color = COL_RKB, linewidth = 0.65) +
  geom_hline(yintercept = EQ_SHARE, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  facet_wrap(~ mesic, nrow = 3, ncol = 4) +
  scale_x_discrete(breaks = x_ticky) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), " %"),
                     expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title    = "Intradenní vzorec volatility per měsíc — NVDA · 2024 (RK Bartlett)",
    subtitle = paste0(
      "Průměrný podíl RK volatility pro každý 5-min slot zvlášť pro každý měsíc.\n",
      "Přerušovaná čára = rovnoměrné rozdělení (",
      round(EQ_SHARE, 2), " %). Vyplněná plocha zvýrazňuje oblast pod křivkou."),
    x = "Čas (Eastern Standard Time)", y = "RK Bartlett [%]",
    caption = paste0(
      "NVDA earnings dny (Feb 21, May 22, Aug 28, Nov 20) způsobují výrazně vyšší ranní hodnoty. ",
      "N ≈ 21 dní/měsíc.")
  ) + THEME_BP +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7.5),
        strip.text  = element_text(size = 10, face = "bold"))

uloz_graf(g11, "Graf11_Mesicny_Panely.png", 16, 10)


# ── Graf 12: Scatter plot — volatilita vs. intenzita a likvidita ──────
# Každý bod = jeden průměrný 5-min interval (přes N_DAYS dní).
# Barva kóduje čas od otevření → vidíme, které intervaly táhnou korelaci.
# Pearsonovo r ≈ 0.58–0.59, ale je taženo především ranními body
# (09:30–10:30), kde jsou všechny tři veličiny současně vysoké.
# Odpolední shluk: nízká volatilita, nízká likvidita, ale stoupající
# intenzita → počet obchodů roste, ale cenová pohyblivost ne.

scatter_dt <- ushape[, .(time_f, minutes_from_open, rkb_mean,
                         int_mean, liq_mean, hour_minute)]

r_int <- cor(scatter_dt$rkb_mean, scatter_dt$int_mean, use = "complete.obs")
r_liq <- cor(scatter_dt$rkb_mean, scatter_dt$liq_mean, use = "complete.obs")

pSa <- ggplot(scatter_dt, aes(x = int_mean, y = rkb_mean,
                              color = minutes_from_open)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linetype = "dashed", linewidth = 0.8, fill = "grey85") +
  geom_point(size = 3.5, alpha = 0.90) +
  geom_text(
    data = scatter_dt[hour_minute %in% c("09:30","09:40","09:50",
                                          "11:00","14:00","15:55")],
    aes(label = hour_minute), nudge_y = 0.12, size = 2.5,
    color = "grey20", fontface = "bold") +
  annotate("text",
           x = min(scatter_dt$int_mean, na.rm = TRUE),
           y = max(scatter_dt$rkb_mean, na.rm = TRUE) * 0.93,
           label = sprintf("r = %.3f", r_int),
           size = 4.5, color = COL_LIQ, fontface = "bold", hjust = 0) +
  scale_color_viridis_c(
    name = "Čas od otevření\n(minuty)", option = "plasma", direction = -1,
    labels = function(x) paste0(
      sprintf("%02d", (570L + as.integer(x)) %/% 60L), ":",
      sprintf("%02d", (570L + as.integer(x)) %% 60L))) +
  scale_x_continuous(labels = function(x) paste0(round(x, 2), " %")) +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), " %")) +
  labs(title = "a)  Volatilita vs. Intenzita",
       x = "Intenzita (podíl počtu obchodů) [%]",
       y = "RK Bartlett volatilita [%]") +
  THEME_BP + theme(legend.position = "right",
                   legend.key.height = unit(1.2, "cm"))

pSb <- ggplot(scatter_dt, aes(x = liq_mean, y = rkb_mean,
                              color = minutes_from_open)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linetype = "dashed", linewidth = 0.8, fill = "grey85") +
  geom_point(size = 3.5, alpha = 0.90) +
  geom_text(
    data = scatter_dt[hour_minute %in% c("09:30","09:40","09:50",
                                          "11:00","14:00","15:55")],
    aes(label = hour_minute), nudge_y = 0.12, size = 2.5,
    color = "grey20", fontface = "bold") +
  annotate("text",
           x = min(scatter_dt$liq_mean, na.rm = TRUE),
           y = max(scatter_dt$rkb_mean, na.rm = TRUE) * 0.93,
           label = sprintf("r = %.3f", r_liq),
           size = 4.5, color = COL_LIQ, fontface = "bold", hjust = 0) +
  scale_color_viridis_c(
    name = "Čas od otevření\n(minuty)", option = "plasma", direction = -1,
    labels = function(x) paste0(
      sprintf("%02d", (570L + as.integer(x)) %/% 60L), ":",
      sprintf("%02d", (570L + as.integer(x)) %% 60L))) +
  scale_x_continuous(labels = function(x) paste0(round(x, 2), " %")) +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), " %")) +
  labs(title = "b)  Volatilita vs. Likvidita",
       x = "Likvidita (podíl objemu) [%]",
       y = "RK Bartlett volatilita [%]") +
  THEME_BP + theme(legend.position = "right",
                   legend.key.height = unit(1.2, "cm"))

g12 <- pSa + pSb +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  plot_annotation(
    title    = "Vztah volatility, intenzity a likvidity — NVDA · 2024",
    subtitle = paste0(
      "Každý bod = průměrný 5-min interval (N = ", N_DAYS, " dní). ",
      "Barva = čas od otevření. OLS přímka + 95% CI."),
    caption  = paste0(
      "RK Bartlett: Barndorff-Nielsen et al. (2008). ",
      "MDH (Mixture of Distributions Hypothesis): Clark (1973).\n",
      "Korelace jsou taženy převážně ranními intervaly (09:30–10:30), ",
      "kde jsou všechny tři veličiny současně vysoké."),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      plot.caption  = element_text(size = 8, color = "grey45",
                                   hjust = 0, lineheight = 1.4)
    )
  )

uloz_graf(g12, "Graf12_Scatter_final.png", 16, 7.5)


# =====================================================================
# ZÁVĚR
# =====================================================================

message("\n", strrep("=", 65))
message("✓ HOTOVO — grafy uloženy do: ", file.path(OUT_DIR, "grafy"))
message(strrep("=", 65))
cat(paste0(
  "\n  Graf01_Signature_Plot.png       Zdůvodnění volby estimátoru\n",
  "  Graf02_Ushape_Zoom.png          Hlavní výsledek — U-shape volatility\n",
  "  Graf03_Porovnani_Vzorcu.png     Volatilita, intenzita, likvidita\n",
  "  Graf04_Srovnani_Estimatoru.png  RK Bartlett vs. naivní RV\n",
  "  Graf05_Normalizovane_Srovnani.png Normalizované U-křivky\n",
  "  Graf06_Heatmapa_Den_Cas.png     Heatmapa den × čas\n",
  "  Graf07_Kvartal_final.png        Kvartální stabilita Q1–Q4\n",
  "  Graf08_BoxPlot_Distribuce.png   Box plot distribucí\n",
  "  Graf09_Hustota_FazeDne.png      Hustotní grafy fází dne\n",
  "  Graf10_DenniRV_TimeSeries.png   Denní časová řada RV\n",
  "  Graf11_Mesicny_Panely.png       Měsíční panely\n",
  "  Graf12_Scatter_final.png        Scatter plot vztahů\n"
))
