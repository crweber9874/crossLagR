library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(magick)
library(ggplot2)
library(cowplot)
library(showtext)

font_add_google("Inter", "inter")
showtext_auto()

# 1. Render the CLPM diagram as a transparent PNG.
clpm <- grViz("
digraph CLPM {
  graph [layout = neato, bgcolor = transparent, splines = true, margin = 0.15]
  node [shape = plaintext, fontsize = 18, fontcolor = '#fafafa', fontname = 'Helvetica-Bold']
  edge [color = '#d4d4d8', penwidth = 0.9, arrowsize = 0.4]

  Xt1 [pos = '-1.0,0.7!', label = <X<sub>1</sub>>]
  Xt2 [pos = '1.0,0.7!',  label = <X<sub>2</sub>>]
  Yt1 [pos = '-1.0,-0.7!', label = <Y<sub>1</sub>>]
  Yt2 [pos = '1.0,-0.7!',  label = <Y<sub>2</sub>>]

  Xt1 -> Xt2
  Yt1 -> Yt2
  Xt1 -> Yt2
  Yt1 -> Xt2
}
")

svg_path <- tempfile(fileext = ".svg")
png_path <- tempfile(fileext = ".png")
writeLines(export_svg(clpm), svg_path)
rsvg_png(svg_path, png_path, width = 800)
diagram_img <- image_read(png_path) |> image_trim()

# 2. Build the hex sticker entirely in ggplot2.
# Pointy-top hex: width = sqrt(3), height = 2. Center at (0,0).
hex_pts <- data.frame(
  x = c(0, sqrt(3)/2, sqrt(3)/2, 0, -sqrt(3)/2, -sqrt(3)/2),
  y = c(1, 0.5,        -0.5,      -1, -0.5,        0.5)
)

# Light grey band across the upper portion, clipped to hex shape.
band_y <- 0.42
band_pts <- data.frame(
  x = c(0, sqrt(3)/2, sqrt(3)/2, -sqrt(3)/2, -sqrt(3)/2),
  y = c(1, 0.5,        band_y,    band_y,     0.5)
)

p <- ggplot() +
  geom_polygon(data = hex_pts,  aes(x, y), fill = "#0a0a0a") +
  geom_polygon(data = band_pts, aes(x, y), fill = "#d4d4d8") +
  geom_polygon(data = hex_pts,  aes(x, y),
               fill = NA, color = "#a1a1aa", linewidth = 1.6) +
  annotate("text", x = 0, y = 0.78, label = "crossLagR",
           family = "inter", fontface = "bold",
           size = 38, color = "#0a0a0a") +
  annotate("text", x = 0, y = 0.52,
           label = "Panel Based Structural Equation Models in R",
           family = "inter", size = 13, color = "#27272a") +
  annotate("text", x = 0, y = -0.92,
           label = "crweber9874.github.io/crossLagR",
           family = "inter", size = 10, color = "#71717a") +
  coord_fixed(xlim = c(-sqrt(3)/2 - 0.03, sqrt(3)/2 + 0.03),
              ylim = c(-1.03, 1.03), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))

# Composite the CLPM diagram in the lower (dark) section.
final <- ggdraw(p) +
  draw_image(diagram_img,
             x = 0.5, y = 0.32,
             width = 0.55,
             hjust = 0.5, vjust = 0.5)

out_dir  <- "man/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_path <- file.path(out_dir, "logo.png")

ggsave(out_path, final,
       width = 4, height = 4.62, units = "in",
       dpi = 600, bg = "transparent")

cat("Saved:", out_path, "\n")
