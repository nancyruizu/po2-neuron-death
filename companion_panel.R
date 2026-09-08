# ---------------------------------------------------------------
# Companion panel for Figure d: within-animal centred pO2,
# surviving vs disappearing neurons, by genotype.
#
# Summary statistics are ANIMAL-WEIGHTED (mean of per-animal means),
# so that the figure matches the within-animal permutation test.
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(systemfonts)   # font lookup
  library(ragg)          # PNG device that honours system fonts
})

# ---- font -------------------------------------------------------
# Arial is used when installed (Mac/Windows); otherwise fall back to a
# metric-compatible substitute so the layout is identical either way.
pick_font <- function(prefs = c("Arial", "Liberation Sans", "Helvetica",
                                "Nimbus Sans", "Arimo", "sans")) {
  fams <- unique(systemfonts::system_fonts()$family)
  hit  <- prefs[prefs %in% fams]
  if (length(hit)) hit[1] else "sans"
}
FONT <- pick_font()
message("Font in use: ", FONT)

COL_SURV <- "#BEBEBE"   # grey  = surviving   (matches panel d)
COL_DIE  <- "#E8442D"   # red   = disappearing
COL_INK  <- "black"

set.seed(20260904)

read_one <- function(path, genotype) {
  read_excel(path, sheet = 1) |>
    mutate(Genotype = genotype, CaseID = as.character(CaseID),
           isDying = as.logical(isDying))
}

dat <- bind_rows(
  read_one("data/Tg_results.xlsx", "rTg4510"),
  read_one("data/WT_results.xlsx", "WT")
) |>
  group_by(CaseID) |>
  mutate(pO2c = pO2 - mean(pO2)) |>
  ungroup() |>
  mutate(
    Genotype = factor(Genotype, levels = c("rTg4510", "WT")),
    Status   = factor(isDying, levels = c(FALSE, TRUE),
                      labels = c("Surviving", "Disappearing"))
  )

# ---- animal-level means (the experimental unit is the animal) ----
animal_means <- dat |>
  group_by(Genotype, CaseID, Status) |>
  summarise(m = mean(pO2c), n_neurons = n(), .groups = "drop")

# animal-weighted group summary
grp <- animal_means |>
  group_by(Genotype, Status) |>
  summarise(
    n_animals = n(),
    mean      = mean(m),
    sd        = sd(m),
    sem       = sd(m) / sqrt(n()),
    .groups   = "drop"
  )
print(as.data.frame(grp), digits = 3)

# neuron-pooled equivalent, for contrast (NOT plotted)
cat("\nNeuron-pooled means (differ from animal-weighted -- see note):\n")
dat |> group_by(Genotype, Status) |>
  summarise(n = n(), pooled_mean = mean(pO2c), .groups = "drop") |>
  as.data.frame() |> print(digits = 3)

# ---- annotation ------------------------------------------------
YMIN <- -30; YMAX <- 46
BR_Y <- 34; TICK <- 2

ns_lab <- tibble(
  Genotype = factor(c("rTg4510", "WT"), levels = c("rTg4510", "WT")),
  label    = "n.s."
)

n_lab <- dat |>
  count(Genotype, Status) |>
  mutate(label = paste0("n = ", format(n, big.mark = ",", trim = TRUE)))

bracket <- list(
  annotate("segment", x = 1, xend = 2, y = BR_Y, yend = BR_Y,
           linewidth = 0.4, colour = COL_INK),
  annotate("segment", x = 1, xend = 1, y = BR_Y, yend = BR_Y - TICK,
           linewidth = 0.4, colour = COL_INK),
  annotate("segment", x = 2, xend = 2, y = BR_Y, yend = BR_Y - TICK,
           linewidth = 0.4, colour = COL_INK)
)

p <- ggplot(dat, aes(Status, pO2c)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey60") +
  # individual neurons
  geom_jitter(aes(colour = Status), width = 0.20, height = 0,
              size = 0.55, alpha = 0.35, stroke = 0) +
  # per-animal means, paired within animal
  geom_line(data = animal_means, aes(x = Status, y = m, group = CaseID),
            colour = "grey30", linewidth = 0.3, alpha = 0.9) +
  geom_point(data = animal_means, aes(x = Status, y = m),
             colour = COL_INK, fill = "white", shape = 21,
             size = 1.8, stroke = 0.45) +
  # animal-weighted group mean (matches the permutation test)
  geom_errorbar(data = grp,
                aes(x = Status, ymin = mean - sem, ymax = mean + sem),
                inherit.aes = FALSE, width = 0.12,
                linewidth = 0.45, colour = COL_INK) +
  geom_segment(data = grp,
               aes(x = as.numeric(Status) - 0.20, xend = as.numeric(Status) + 0.20,
                   y = mean, yend = mean),
               inherit.aes = FALSE, linewidth = 0.55, colour = COL_INK) +
  bracket +
  geom_text(data = ns_lab, aes(x = 1.5, y = BR_Y + 3.5, label = label),
            inherit.aes = FALSE, size = 3.3, fontface = "italic",
            family = FONT) +
  geom_text(data = n_lab, aes(x = Status, y = YMIN + 2.5, label = label),
            inherit.aes = FALSE, size = 2.5, colour = "grey35",
            family = FONT) +
  scale_colour_manual(values = c(Surviving = COL_SURV, Disappearing = COL_DIE)) +
  scale_x_discrete(labels = c("Surviving" = "Surviving",
                              "Disappearing" = "Disappearing")) +
  scale_y_continuous(limits = c(YMIN, YMAX), breaks = seq(-30, 40, 10),
                     expand = c(0, 0)) +
  facet_wrap(~ Genotype) +
  # Axis label: plain string, NOT expression(). plotmath draws Delta from a
  # symbol font (embedding a second typeface) and cannot render Unicode at all.
  # A plain string is drawn entirely in FONT. If your Arial build lacks the
  # subscript-two glyph, swap the label for "\u0394pO2 from animal mean (mmHg)".
  labs(x = NULL, y = "\u0394pO\u2082 from animal mean (mmHg)") +
  theme_classic(base_size = 10, base_family = FONT) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.x      = element_text(size = 8, colour = COL_INK,
                                    angle = 30, hjust = 1),
    axis.text.y      = element_text(size = 8.5, colour = COL_INK),
    axis.title.y     = element_text(size = 9),
    axis.line        = element_line(colour = COL_INK, linewidth = 0.4),
    axis.ticks       = element_line(colour = COL_INK, linewidth = 0.4),
    panel.spacing    = unit(10, "pt"),
    legend.position  = "none",
    plot.margin      = margin(6, 8, 2, 6)
  )

# ragg + cairo_pdf both resolve system fonts and embed them in the output
ggsave("figures/panel_e.png", p, width = 3.9, height = 3.3, dpi = 400,
       bg = "white", device = ragg::agg_png)
ggsave("figures/panel_e.pdf", p, width = 3.9, height = 3.3,
       device = cairo_pdf)   # embeds the font; required by most journals

cat("\nWrote figures/panel_e.png and figures/panel_e.pdf\n")
