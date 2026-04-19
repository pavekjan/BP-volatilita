# BP-volatilita
# Bakalářská práce — Analýza intradenních vzorců volatility a likvidity

**Autor:** Jan Pávek  
**Vedoucí práce:** doc. Mgr. Vladimír Holý, PhD.  
**Instituce:** Vysoká škola ekonomická v Praze, Fakulta informatiky a statistiky  
**Rok:** 2025  

---

## O práci

Práce zkoumá intradenní vzorce volatility a likvidity u 23 amerických akcií z burz NYSE a NASDAQ za rok 2024. Jako primární estimátor volatility byl použit Realized Kernel s Bartlettovým jádrem (RK Bartlett), který filtruje mikrostrukturní šum vznikající při vysokofrekvenčním vzorkování.

Byly testovány tři hypotézy:
- **H1** — U-shape efekt je statisticky průkazný napříč různými typy akcií (Kruskal-Wallis test)
- **H2** — Objem obchodů lépe předpovídá volatilitu než počet transakcí (Mixture of Distributions Hypothesis)
- **H3** — Hloubka U-tvaru závisí na celkové roční volatilitě akcie (OLS regrese s HAC korekcí)

---

## Struktura repozitáře

```
├── BP_AKCIE.R          # Skript 01 — výpočet RK Bartlett volatility a regresní analýza (23 akcií)
├── BP_GRAFY23.R        # Skript 02 — grafy pro srovnávací analýzu 23 akcií
├── BP_NVDA_grafy.R     # Skript 03 — grafy pro případovou studii NVDA
└── README.md
```

---

## Popis skriptů

### `BP_AKCIE.R` — Výpočet volatility a regresní analýza
Zpracovává tick-by-tick data 23 akcií. Pro každou akcii:
1. Načte a vyčistí surová tick data (formát `.rda`)
2. Vypočítá RK Bartlett volatilitu v 5minutových oknech
3. Sestaví průměrné sezónní profily volatility, intenzity a likvidity
4. Provede statistické testy sezónnosti (Kruskal-Wallis, Pearson, Ljung-Box)
5. Odhadne regresní modely (OLS + HAC, kvantilová regrese, panelová regrese)

**Výstupy:** `BP_srovnani_vystupy/tabulky/` — CSV soubory s výsledky

### `BP_GRAFY23.R` — Grafy pro 23 akcií
Načítá výstupní tabulky ze Skriptu 01 a produkuje 8 grafů pro srovnávací analýzu napříč sektory (U-tvar, MDH test, heatmapa, scatter plot atd.)

### `BP_NVDA_grafy.R` — Grafy pro NVDA (případová studie)
Produkuje 12 grafů pro detailní analýzu akcie NVDA: Volatility Signature Plot, intradenní U-tvar, kvartální stabilita, distribuce volatility, měsíční panely atd.

---

## Použité R balíčky

| Balíček | Účel |
|---|---|
| `highfrequency` | Výpočet RK Bartlett estimátoru |
| `data.table` | Efektivní zpracování velkých datových sad |
| `xts`, `zoo` | Práce s časovými řadami |
| `sandwich`, `lmtest` | HAC kovarianční matice (Newey-West) |
| `quantreg` | Kvantilová regrese |
| `plm` | Panelová regrese s fixními efekty |
| `ggplot2`, `ggrepel` | Vizualizace výsledků |

---

## Poznámky ke spuštění

- Data nejsou součástí repozitáře (tick data jsou proprietární a byla poskytnuta vedoucím práce)
- Výpočty byly spuštěny na výpočetním klastru **MetaCentrum** — lokální spuštění může trvat několik hodin
- Skripty používají průběžné cachování (`.rds` soubory) — při přerušení lze výpočet obnovit bez ztráty dat
- Doporučené prostředí: R ≥ 4.3, RStudio
