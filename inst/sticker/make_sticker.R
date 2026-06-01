library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(hexSticker)
library(showtext)
library(magick)

font_add_google("Inter", "inter")
showtext_auto()

clpm <- grViz("
digraph CLPM {
  graph [layout = neato, bgcolor = transparent, splines = true, margin = 0.15]
  node [shape = plaintext, fontsize = 28, fontcolor = '#fafafa', fontname = 'Helvetica-Bold']
  edge [color = '#d4d4d8', penwidth = 1.6, arrowsize = 0.7]

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

# Trim then pad into a square canvas so hexSticker's subplot scaling
# is predictable. Width = max(w,h)*1.15.
img    <- image_read(png_path) |> image_trim()
info   <- image_info(img)
canvas <- round(max(info$width, info$height) * 1.45)
img    <- image_extent(img, geometry = sprintf("%dx%d", canvas, canvas),
                       color = "transparent", gravity = "center")
trimmed_png <- tempfile(fileext = ".png")
image_write(img, trimmed_png)

out_dir <- "man/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

sticker(
  subplot   = trimmed_png,
  package   = "crossLagR",
  s_x       = 1,
  s_y       = 0.82,
  s_width   = 0.90,
  s_height  = 0.90,
  p_x       = 1,
  p_y       = 1.50,
  p_size    = 26,
  p_color   = "#fafafa",
  p_family  = "inter",
  h_fill    = "#0a0a0a",
  h_color   = "#a1a1aa",
  h_size    = 1.4,
  url       = "crweber9874.github.io/crossLagR",
  u_color   = "#71717a",
  u_size    = 4.5,
  u_family  = "inter",
  white_around_sticker = FALSE,
  filename  = file.path(out_dir, "logo.png"),
  dpi       = 600
)

cat("Saved:", file.path(out_dir, "logo.png"), "\n")
