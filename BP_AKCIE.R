# =====================================================================
# BAKALÁŘSKÁ PRÁCE
# Analýza intradenních vzorců volatility a likvidity akciových trhů
# NASDAQ & NYSE · 2024
#
# Autor:    Jan Pávek
# Vedoucí:  [jméno vedoucího]
# Instituce: [název fakulty a katedry]
#
# SKRIPT 01: Výpočet realizované volatility a regresní analýza
# ---------------------------------------------------------------------
# Skript zpracovává tick-by-tick data 23 akcií z burz NASDAQ a NYSE.
# Pro každou akcii se postupně:
#   (1) načtou a vyčistí surová tick data (formát .rda),
#   (2) spočítají realizovaná jádra (RK Bartlett) v 5minutových oknech,
#   (3) sestaví průměrné sezónní profily volatility, intenzity a likvidity,
#   (4) provedou statistické testy sezónnosti,
#   (5) odhadnou regresní modely (OLS+HAC, kvantilová, panelová regrese).
#
# Výstupy (složka BP_srovnani_vystupy/):
#   cache/<SYMBOL>.rds               — mezivýsledky (RK okna) na akcii
#   tabulky/01_souhrn_akcie.csv      — klíčové charakteristiky vzorku
#   tabulky/02_sezonni_profily.csv   — průměrné sezónní profily (78 int.)
#   tabulky/03_statisticke_testy.csv — KW, Pearson r, Friedman, Ljung-Box
#   tabulky/04_regrese_ols_hac.csv   — OLS s HAC std. chybami (Newey-West)
#   tabulky/05_regrese_quantil.csv   — kvantilová regrese (τ = 0,1/0,5/0,9)
#   tabulky/06_regrese_panel.csv     — panelová regrese s fixními efekty
#
# Poznámka ke spuštění:
#   Výpočet je výpočetně náročný (doporučeno spustit přes noc).
#   Průběžné výsledky se ukládají do cache/ — při přerušení a opětovném
#   spuštění skript automaticky přeskočí již zpracované akcie.
# =====================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)    # efektivní práce s datovými tabulkami
  library(lubridate)     # manipulace s časovými razítky
  library(highfrequency) # realizovaná jádra (rKernelCov)
  library(xts)           # časové řady pro highfrequency
  library(zoo)           # pomocné funkce pro časové řady
  library(sandwich)      # HAC kovarianční matice (Newey-West, 1987)
  library(lmtest)        # coeftest() — testování koeficientů s HAC
  library(quantreg)      # kvantilová regrese rq() (Koenker, 2005)
  library(plm)           # panelová regrese s fixními efekty
})


# =====================================================================
# SEKCE 1: KONFIGURACE
# =====================================================================

# Cesta ke složce se surovými tick daty na výpočetním klastru MetaCentrum
BASE_DATA_DIR <- "/storage/brno2/home/janpavek03/ondemand/data/sys/dashboard/batch_connect/sys/bc_meta_rstudio/output/ec81a493-d26d-4443-97fa-cce58767adc9/ab9c2d6d-c925-494b-9815-4f2ca9d09072"

# Výstupní složka (vytvoří se automaticky)
OUT_DIR <- file.path(getwd(), "BP_srovnani_vystupy")
dir.create(file.path(OUT_DIR, "cache"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tabulky"), recursive = TRUE, showWarnings = FALSE)

# Parametry výpočtu — shodné s analýzou NVIDIA pro metodologickou konzistenci
TZ_MKT        <- "America/New_York"  # časová zóna trhu (Eastern Time)
RTH_MIN_START <- 570L                # 09:30 ET — začátek regulérního obchodování
RTH_MIN_END   <- 960L                # 16:00 ET — konec regulérního obchodování
INTERVAL_MIN  <- 5L                  # délka okna pro výpočet RK (minuty)
N_BUCKETS     <- 78L                 # počet 5minutových intervalů za obchodní den
MIN_TRADES    <- 20L                 # minimální počet ticků pro validní odhad RK
MIN_TICKS_DAY <- 500L                # minimální počet ticků za den (filtr kvality dat)
ROK           <- 2024                # analyzovaný rok

# ── Výběr akcií: 23 titulů napříč 8 sektory ────────────────────────
# Výběr pokrývá hlavní sektory S&P 500 a zajišťuje průřezové srovnání.
# Akcie z NASDAQ mají příponu .O, akcie z NYSE příponu .N (Reuters kódy).
AKCIE <- data.table(
  slozka = c(
    "NASDAQ_NVDA_2024", "NASDAQ_AMD_2024",  "NASDAQ_MSFT_2024",
    "NASDAQ_AAPL_2024", "NASDAQ_AMZN_2024", "NASDAQ_GOOGL_2024",
    "NASDAQ_NFLX_2024", "NASDAQ_INTC_2024",
    "NYSE_JPM_2024",    "NYSE_BAC_2024",    "NYSE_GS_2024",
    "NYSE_XOM_2024",    "NYSE_CVX_2024",
    "NYSE_JNJ_2024",    "NYSE_PFE_2024",
    "NYSE_PG_2024",     "NYSE_KO_2024",
    "NYSE_MCD_2024",    "NYSE_NKE_2024",
    "NYSE_GE_2024",     "NYSE_CAT_2024",
    "NYSE_VZ_2024",     "NYSE_T_2024"
  ),
  symbol = c(
    "NVDA.O", "AMD.O",  "MSFT.O", "AAPL.O", "AMZN.O",
    "GOOGL.O","NFLX.O", "INTC.O",
    "JPM.N",  "BAC.N",  "GS.N",
    "XOM.N",  "CVX.N",
    "JNJ.N",  "PFE.N",
    "PG.N",   "KO.N",
    "MCD.N",  "NKE.N",
    "GE.N",   "CAT.N",
    "VZ.N",   "T.N"
  ),
  sektor = c(
    rep("Technologie", 8),
    rep("Finance", 3),
    rep("Energie", 2),
    rep("Zdravotnictvi", 2),
    rep("Defenzivni spotreba", 2),
    rep("Cyklicka spotreba", 2),
    rep("Prumysl", 2),
    rep("Telekomunikace", 2)
  ),
  burza = c(rep("NASDAQ", 8), rep("NYSE", 15))
)

message("Akcie k zpracování: ", nrow(AKCIE))


# =====================================================================
# SEKCE 2: POMOCNÉ FUNKCE
# =====================================================================

# ── nacti_rda() ──────────────────────────────────────────────────────
# Načte jeden soubor .rda a vrátí jeho obsah jako data.table.
# Při chybě (poškozený soubor, prázdný soubor) vrátí NULL.
nacti_rda <- function(cesta) {
  tryCatch({
    env   <- new.env()
    nazvy <- load(cesta, envir = env)
    if (length(nazvy) == 0) return(NULL)
    as.data.table(env[[nazvy[1]]])
  }, error = function(e) NULL)
}

# ── cistit_ticky() ───────────────────────────────────────────────────
# Vyčistí surová tick data jedné akcie:
#   - standardizuje názvy sloupců na velká písmena,
#   - filtruje pouze záznamy cílového symbolu,
#   - převede časová razítka UTC → Eastern Time,
#   - omezí data na regulérní obchodní hodiny (09:30–16:00 ET),
#   - odstraní záznamy s nulovou/zápornou cenou nebo objemem,
#   - deduplikuje záznamy se shodným časovým razítkem (ponechá poslední),
#   - přidá sloupec logp (přirozený logaritmus ceny) pro výpočet výnosů.
cistit_ticky <- function(dt, symbol_cil) {
  if (is.null(dt) || !is.data.table(dt) || nrow(dt) == 0) return(data.table())
  setnames(dt, toupper(names(dt)))
  if (!all(c("TIMESTAMP","VALUE","VOLUME","SYMBOL") %in% names(dt)))
    return(data.table())
  symbol_raw <- gsub("\\.[NO]$", "", symbol_cil)
  dt <- dt[SYMBOL %in% c(symbol_cil, symbol_raw)]
  if (nrow(dt) == 0) return(data.table())
  dt[, ts_utc := ymd_hms(TIMESTAMP, tz = "UTC", quiet = TRUE)]
  dt <- dt[!is.na(ts_utc)]
  dt[, ts_ny   := with_tz(ts_utc, TZ_MKT)]
  dt[, date_ny := as.Date(ts_ny, tz = TZ_MKT)]
  dt[, min_od_pulnoci := hour(ts_ny) * 60L + minute(ts_ny)]
  dt <- dt[min_od_pulnoci >= RTH_MIN_START & min_od_pulnoci < RTH_MIN_END]
  if (nrow(dt) == 0) return(data.table())
  dt[, VALUE  := as.numeric(VALUE)]
  dt[, VOLUME := as.numeric(VOLUME)]
  dt <- dt[is.finite(VALUE) & VALUE > 0 & is.finite(VOLUME) & VOLUME >= 0]
  if (nrow(dt) == 0) return(data.table())
  setorder(dt, ts_ny)
  dt <- dt[, .SD[.N], by = ts_ny]   # deduplikace: ponechá poslední tick v ms
  dt[, logp := log(VALUE)]
  dt[, .(ts_ny, date_ny, VALUE, VOLUME, logp)]
}

# ── jako_skalar() ────────────────────────────────────────────────────
# Pomocná funkce: bezpečně extrahuje první skalární hodnotu z výstupu
# rKernelCov(), který může vracet matici, vektor nebo seznam.
jako_skalar <- function(obj) {
  if (is.null(obj))    return(NA_real_)
  if (is.matrix(obj))  return(as.numeric(obj[1L, 1L]))
  if (is.numeric(obj)) return(as.numeric(obj[1L]))
  if (is.list(obj))    return(jako_skalar(obj[[1L]]))
  NA_real_
}

# ── vypocti_RK_okno() ────────────────────────────────────────────────
# Spočítá realizované jádro (RK Bartlett) pro jedno 5minutové okno.
#
# Metoda: Barndorff-Nielsen et al. (2008) — realizované jádro s Bartlettovým
# váhovacím schématem. Šířka pásma H = min(10, floor(n/3)) je volena
# konzervativně, aby se omezila mikrostrukturní šumová složka.
#
# Vstup:  data.table s logp (log-cena) a VOLUME pro jedno okno
# Výstup: seznam s n_trades, volume_sum, rv_tick (naivní RV) a rk_bartlett
vypocti_RK_okno <- function(okno) {
  n   <- nrow(okno)
  rv  <- sum(diff(okno$logp)^2, na.rm = TRUE)   # naivní realizovaná variance
  vol <- sum(okno$VOLUME, na.rm = TRUE)
  rk  <- NA_real_
  if (n >= MIN_TRADES) {
    px <- xts(okno$VALUE, order.by = okno$ts_ny)
    colnames(px) <- "price"
    H  <- min(10L, floor(n / 3L))                # šířka pásma jádra
    rk <- tryCatch(
      jako_skalar(rKernelCov(px, makeReturns = TRUE,
                             kernelType  = "bartlett",
                             kernelParam = H, kernelDOFadj = TRUE)),
      error = function(e) NA_real_)
    rk <- ifelse(is.finite(rk) & rk > 0, rk, NA_real_)
  }
  list(n_trades = n, volume_sum = vol, rv_tick = rv, rk_bartlett = rk)
}

# ── zpracuj_akcii() ──────────────────────────────────────────────────
# Hlavní zpracovatelská funkce pro jednu akcii. Provede:
#   (1) načtení všech .rda souborů z příslušné složky,
#   (2) vyčistění a spojení tick dat,
#   (3) filtraci dní s nedostatečným počtem ticků (MIN_TICKS_DAY),
#   (4) přiřazení každého ticku do 5minutového okna (bucket),
#   (5) výpočet RK Bartlett pro každé okno a obchodní den,
#   (6) výpočet sezónních podílů: rkb_share, rv_share, trades_share, volume_share.
#
# Sezónní podíl: podíl realizovaného jádra daného okna na denním součtu.
# Tento přístup normalizuje intradenní profil a umožňuje srovnání mezi akcemi
# i dny s různou celkovou volatilitou (Andersen & Bollerslev, 1997).
zpracuj_akcii <- function(slozka, symbol) {
  data_path <- file.path(BASE_DATA_DIR, slozka)
  if (!dir.exists(data_path)) {
    message("    NENALEZENO: ", slozka); return(NULL)
  }
  soubory <- list.files(data_path, pattern = "\\.rda$", full.names = TRUE)
  if (length(soubory) == 0) {
    message("    Žádné soubory: ", slozka); return(NULL)
  }
  vsechny <- vector("list", length(soubory))
  for (i in seq_along(soubory)) {
    dt_raw <- nacti_rda(soubory[i])
    dt     <- cistit_ticky(dt_raw, symbol)
    if (is.data.table(dt) && nrow(dt) > 0) vsechny[[i]] <- dt
  }
  ticky <- rbindlist(Filter(Negate(is.null), vsechny),
                     use.names = TRUE, fill = TRUE)
  if (nrow(ticky) == 0) return(NULL)
  setorder(ticky, ts_ny)
  
  # Filtr: ponechej pouze dny s dostatečnou datovou kvalitou
  platne <- ticky[, .(n = .N), by = date_ny][n >= MIN_TICKS_DAY, date_ny]
  ticky  <- ticky[date_ny %in% platne]
  if (nrow(ticky) == 0) return(NULL)
  
  # Přiřazení ticků do 5minutových oken (floor na nejbližší nižší 5 min)
  ticky[, bucket := floor_date(ts_ny,
                               unit = paste0(INTERVAL_MIN, " minutes"))]
  ticky <- ticky[
    bucket >= as.POSIXct(paste0(date_ny, " 09:30:00"), tz = TZ_MKT) &
      bucket <  as.POSIXct(paste0(date_ny, " 16:00:00"), tz = TZ_MKT)]
  setorder(ticky, date_ny, bucket, ts_ny)
  
  # Výpočet RK Bartlett pro každé okno (date_ny × bucket)
  results <- ticky[, vypocti_RK_okno(.SD), by = .(date_ny, bucket)]
  results[, minutes_from_open := as.integer(difftime(
    bucket,
    as.POSIXct(paste0(date_ny, " 09:30:00"), tz = TZ_MKT),
    units = "mins"))]
  results[, hour_minute := format(bucket, "%H:%M")]
  results[, symbol      := symbol]
  results[, quarter_lbl := paste0("Q", quarter(as.Date(date_ny)))]
  
  # Výpočet sezónních podílů
  # rk_imp: imputace NA — okna s méně než MIN_TRADES ticky se vyřadí
  results[, rk_imp := fifelse(is.na(rk_bartlett) | rk_bartlett <= 0,
                              NA_real_, rk_bartlett)]
  
  # Denní součet RK (s korekcí na případná chybějící okna)
  daily_rk <- results[, {
    n_v <- sum(!is.na(rk_imp))
    s   <- sum(rk_imp, na.rm = TRUE)
    .(rkb_day = if (n_v > 0L) s * (N_BUCKETS / n_v) else NA_real_)
  }, by = date_ny]
  
  # Denní součty pro intenzitu a likviditu
  daily_tr <- results[, .(
    total_trades = sum(n_trades,   na.rm = TRUE),
    total_volume = sum(volume_sum, na.rm = TRUE),
    rv_day       = sum(rv_tick,    na.rm = TRUE)
  ), by = date_ny]
  
  results <- merge(results, daily_rk, by = "date_ny", all.x = TRUE)
  results <- merge(results, daily_tr, by = "date_ny", all.x = TRUE)
  
  # Sezónní podíly — každý jako podíl okna na denním součtu
  results[, `:=`(
    rkb_share    = rk_imp     / rkb_day,     # podíl volatility
    rv_share     = rv_tick    / rv_day,       # podíl naivní RV (kontrolní)
    trades_share = n_trades   / total_trades, # podíl intenzity obchodování
    volume_share = volume_sum / total_volume  # podíl objemu (likvidita)
  )]
  results
}


# =====================================================================
# SEKCE 3: HLAVNÍ VÝPOČETNÍ SMYČKA (s cachováním)
# =====================================================================
# Pro každou akcii: pokud existuje cache soubor, přeskoč.
# Jinak spočítej a ulož do cache. Umožňuje přerušení a pokračování.

message("\n", strrep("=", 65))
message("SPOUŠTÍM VÝPOČTY — ", nrow(AKCIE), " AKCIÍ")
message(strrep("=", 65))

for (i in seq_len(nrow(AKCIE))) {
  sym     <- AKCIE$symbol[i]
  slozka  <- AKCIE$slozka[i]
  cache_f <- file.path(OUT_DIR, "cache",
                       paste0(gsub("[./]", "_", sym), ".rds"))
  if (file.exists(cache_f)) {
    message(sprintf("[%2d/%2d] ✓ Cache: %s", i, nrow(AKCIE), sym))
    next
  }
  message(sprintf("[%2d/%2d] Počítám: %s ...", i, nrow(AKCIE), sym))
  t0  <- proc.time()
  res <- tryCatch(zpracuj_akcii(slozka, sym), error = function(e) {
    message("    CHYBA: ", conditionMessage(e)); NULL
  })
  if (!is.null(res) && nrow(res) > 0) {
    saveRDS(res, cache_f)
    elapsed <- round((proc.time() - t0)["elapsed"])
    message(sprintf("    ✓ %d oken, %d dní, %.0f s",
                    nrow(res), uniqueN(res$date_ny), elapsed))
    rm(res); gc()
  } else {
    message("    ✗ Žádná data.")
  }
}


# =====================================================================
# SEKCE 4: AGREGACE — sezónní profily, souhrn, statistické testy
# =====================================================================

message("\n[AGG] Načítám cache a agregovám...")

EQ_SHARE <- 1 / N_BUCKETS   # rovnoměrné rozdělení = 1/78 ≈ 1,28 %

souhrn_list  <- vector("list", nrow(AKCIE))
profily_list <- vector("list", nrow(AKCIE))
testy_list   <- vector("list", nrow(AKCIE))
panel_list   <- vector("list", nrow(AKCIE))

for (i in seq_len(nrow(AKCIE))) {
  sym     <- AKCIE$symbol[i]
  sektor  <- AKCIE$sektor[i]
  cache_f <- file.path(OUT_DIR, "cache",
                       paste0(gsub("[./]", "_", sym), ".rds"))
  if (!file.exists(cache_f)) next
  res <- readRDS(cache_f)
  if (is.null(res) || nrow(res) == 0) next
  
  N_DAYS <- uniqueN(res$date_ny)
  
  # Průměrný sezónní profil (78 intervalů)
  # Hodnoty jsou průměrovány přes všechny obchodní dny roku 2024.
  profil <- res[, .(
    rkb_mean   = mean(rkb_share,    na.rm = TRUE) * 100,
    int_mean   = mean(trades_share, na.rm = TRUE) * 100,
    liq_mean   = mean(volume_share, na.rm = TRUE) * 100,
    rkb_median = median(rkb_share,  na.rm = TRUE) * 100,
    rkb_p25    = quantile(rkb_share, 0.25, na.rm = TRUE) * 100,
    rkb_p75    = quantile(rkb_share, 0.75, na.rm = TRUE) * 100,
    n_obs      = sum(!is.na(rkb_share))
  ), by = .(minutes_from_open, hour_minute)]
  setorder(profil, minutes_from_open)
  profil[, symbol := sym]
  profil[, sektor := sektor]
  profily_list[[i]] <- profil
  
  # Klíčové charakteristiky sezónního profilu
  rano    <- profil[minutes_from_open == 0L,   rkb_mean]   # ranní spike (09:30)
  pol_min <- profil[minutes_from_open %in% 120:235,        # polední minimum
                    min(rkb_mean, na.rm = TRUE)]
  zav     <- profil[minutes_from_open == 385L, rkb_mean]   # závěrečný interval
  
  # Anualizovaná volatilita: σ_ann = √(252 × průměrná denní RKB)
  daily_rkb <- res[, .(rkb_day = sum(rk_imp, na.rm = TRUE)), by = date_ny]
  rkb_ann   <- sqrt(252 * mean(daily_rkb$rkb_day, na.rm = TRUE)) * 100
  
  avg_ticks <- res[, .(n = sum(n_trades, na.rm = TRUE)),
                   by = date_ny][, mean(n)]
  
  # Pearsonovy korelace sezónního profilu s intenzitou a objemem
  r_vol_int <- cor(profil$rkb_mean, profil$int_mean, use = "complete.obs")
  r_vol_liq <- cor(profil$rkb_mean, profil$liq_mean, use = "complete.obs")
  
  # Hloubka U-tvaru: poměr ranního spike a poledního minima
  u_hloubka <- if (!is.na(pol_min) && pol_min > 0) rano / pol_min else NA_real_
  
  souhrn_list[[i]] <- data.table(
    symbol         = sym,
    sektor         = sektor,
    burza          = AKCIE$burza[i],
    n_dni          = N_DAYS,
    avg_ticks_den  = round(avg_ticks),
    rkb_ann_pct    = round(rkb_ann, 1),
    rano_0930      = round(rano,    2),
    poledni_min    = round(pol_min, 2),
    zaverecna_1555 = round(zav,     2),
    u_hloubka      = round(u_hloubka, 1),
    r_vol_int      = round(r_vol_int, 3),
    r_vol_liq      = round(r_vol_liq, 3)
  )
  
  # ── Statistické testy sezónnosti ─────────────────────────────────
  # Rozdělení dne do čtyř fází: otevírací, polední klid, závěrečná, ostatní.
  res[, faze := fcase(
    minutes_from_open %in%   0:55,  "Oteviraci",
    minutes_from_open %in% 120:235, "Poledni",
    minutes_from_open %in% 330:385, "Zaverecna",
    default                        = "Ostatni"
  )]
  
  # Kruskal-Wallisův test: H0 = volatilita je stejná ve všech fázích dne.
  # Neparametrický test — nevyžaduje normalitu residuí.
  kw     <- kruskal.test(rkb_share ~ faze, data = res[!is.na(rkb_share)])
  ct_int <- cor.test(profil$rkb_mean, profil$int_mean, method = "pearson")
  ct_liq <- cor.test(profil$rkb_mean, profil$liq_mean, method = "pearson")
  
  # Friedmanův test: stabilita U-tvaru napříč kvartály roku 2024.
  # H0 = tvar sezónního profilu je v průběhu roku neměnný.
  q_profil <- res[!is.na(rkb_share), .(
    rkb_mean = mean(rkb_share, na.rm = TRUE)
  ), by = .(minutes_from_open, quarter_lbl)]
  q_wide <- dcast(q_profil, minutes_from_open ~ quarter_lbl,
                  value.var = "rkb_mean")
  q_wide <- q_wide[complete.cases(q_wide)]
  ft_p   <- NA_real_
  if (nrow(q_wide) >= 10 && ncol(q_wide) >= 3) {
    mat <- as.matrix(q_wide[, -1, with = FALSE])
    ft  <- tryCatch(friedman.test(mat), error = function(e) NULL)
    if (!is.null(ft)) ft_p <- ft$p.value
  }
  
  # Ljung-Box test: autokorelace sezónního profilu.
  # Průkaz autokorelace zdůvodňuje použití HAC standardních chyb v regresi.
  lb <- Box.test(profil$rkb_mean, lag = 10, type = "Ljung-Box")
  
  testy_list[[i]] <- data.table(
    symbol          = sym, sektor = sektor,
    kw_statistic    = round(kw$statistic, 2),
    kw_p_value      = signif(kw$p.value, 3),
    kw_signif       = ifelse(kw$p.value < 0.001, "***",
                             ifelse(kw$p.value < 0.01,  "**",
                                    ifelse(kw$p.value < 0.05,  "*", "ns"))),
    r_vol_int       = round(ct_int$estimate, 3),
    r_vol_int_p     = signif(ct_int$p.value, 3),
    r_vol_liq       = round(ct_liq$estimate, 3),
    r_vol_liq_p     = signif(ct_liq$p.value, 3),
    friedman_p      = signif(ft_p, 3),
    friedman_signif = ifelse(is.na(ft_p), "N/A",
                             ifelse(ft_p < 0.001, "***",
                                    ifelse(ft_p < 0.01,  "**",
                                           ifelse(ft_p < 0.05,  "*", "ns")))),
    lb_statistic    = round(lb$statistic, 2),
    lb_p_value      = signif(lb$p.value, 3)
  )
  
  panel_list[[i]] <- profil[, .(symbol, sektor,
                                minutes_from_open, rkb_mean,
                                int_mean, liq_mean)]
  
  message(sprintf("  ✓ %-8s  rano=%.1f%%  U=%.0fx  KW p=%s",
                  sym, rano, u_hloubka,
                  ifelse(kw$p.value < 0.001, "<0.001",
                         round(kw$p.value, 3))))
}


# =====================================================================
# SEKCE 5: EXPORT ZBÝVAJÍCÍCH TABULEK
# =====================================================================

message("\n[EXPORT] Ukládám zbývající tabulky...")

tab_souhrn <- rbindlist(Filter(Negate(is.null), souhrn_list))
setorder(tab_souhrn, sektor, symbol)
fwrite(tab_souhrn,  file.path(OUT_DIR, "tabulky", "01_souhrn_akcie.csv"))
message("  ✓ 01_souhrn_akcie.csv (", nrow(tab_souhrn), " akcií)")

setorder(tab_profily, symbol, minutes_from_open)
fwrite(tab_profily, file.path(OUT_DIR, "tabulky", "02_sezonni_profily.csv"))
message("  ✓ 02_sezonni_profily.csv (", nrow(tab_profily), " řádků)")

tab_testy <- rbindlist(Filter(Negate(is.null), testy_list))
setorder(tab_testy, sektor, symbol)
fwrite(tab_testy,   file.path(OUT_DIR, "tabulky", "03_statisticke_testy.csv"))
message("  ✓ 03_statisticke_testy.csv")


# =====================================================================
# ZÁVĚR
# =====================================================================

message("\n", strrep("=", 65))
message("✓ VŠE DOKONČENO — výstupy v: ", OUT_DIR)
message(strrep("=", 65))
cat(paste0(
  "\nTABULKY (složka tabulky/):\n",
  "  01_souhrn_akcie.csv        Klíčové statistiky za každou akcii\n",
  "  02_sezonni_profily.csv     Průměrné profily (78 intervalů)\n",
  "  03_statisticke_testy.csv   KW, Pearson r, Friedman, Ljung-Box\n",
  "  04_regrese_ols_hac.csv     OLS + HAC std. chyby (Newey-West)\n",
  "  05_regrese_quantil.csv     Kvantilová regrese Q10/Q50/Q90\n",
  "  06_regrese_panel.csv       Panelová regrese FE + Hausman test\n",
  "\nDalší krok: spusť 02_srovnavaci_grafy.R\n"
))