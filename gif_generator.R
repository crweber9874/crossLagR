library(dagitty)
library(ggdag)
library(dplyr)
library(lavaan)
library(crossLagR)
library(ggplot2)
library(tibble)
library(ggridges)
library(gganimate)
library(knitr)

# Define the nodes for the SEM with three waves (two rows and three columns)
nodes <- data.frame(
  node = 1:6,
  label = c("X1", "Y1", "X2", "Y2", "X3", "Y3"),
  x = c(1, 1, 2, 2, 3, 3),
  y = c(2, 1, 2, 1, 2, 1)
)

# Define the edges for the SEM
edges <- data.frame(
  from = c(1, 1, 2, 2, 3, 3, 4, 4),
  to =   c(3, 4, 3, 4, 5, 6, 5, 6)
)

# Function to offset the arrowheads
offset_arrow <- function(x1, y1, x2, y2, offset = 0.10) {
  angle <- atan2(y2 - y1, x2 - x1)
  x1_new <- x1 + offset * cos(angle)
  y1_new <- y1 + offset * sin(angle)
  x2_new <- x2 - offset * cos(angle)
  y2_new <- y2 - offset * sin(angle)
  return(data.frame(x1 = x1_new, y1 = y1_new, x2 = x2_new, y2 = y2_new))
}

# Apply the offset to the edges
edges_offset <- do.call(rbind, apply(edges, 1, function(row) {
  offset_arrow(nodes$x[row["from"]], nodes$y[row["from"]],
               nodes$x[row["to"]], nodes$y[row["to"]])
}))

# Define the paths for the balls through the SEM
path1 <- c(1, 3, 5) # X1 to X2 to Y3
path2 <- c(1, 3, 6) # X1 to X2 to X3
path3 <- c(1, 4, 5) # X1 to X2 to Y3
path4 <- c(1, 4, 6) # X1 to X2 to X3

path5 <- c(2, 3, 5) # X1 to X2 to Y3
path6 <- c(2, 3, 6) # X1 to X2 to X3
path7 <- c(2, 4, 5) # X1 to X2 to Y3
path8 <- c(2, 4, 6) # X1 to X2 to X3


path_df1 <- data.frame(node = path1, time = 1:length(path1), ball = "Ball 1", color = "blue")
path_df2 <- data.frame(node = path2, time = 1:length(path2), ball = "Ball 2", color = "blue")
path_df3 <- data.frame(node = path3, time = 1:length(path3), ball = "Ball 3", color = "red")
path_df4 <- data.frame(node = path4, time = 1:length(path4), ball = "Ball 4", color = "red")

path_df5 <- data.frame(node = path5, time = 1:length(path5), ball = "Ball 5", color = "pink")
path_df6 <- data.frame(node = path6, time = 1:length(path6), ball = "Ball 6", color = "pink")
path_df7 <- data.frame(node = path7, time = 1:length(path7), ball = "Ball 7", color = "lightblue")
path_df8 <- data.frame(node = path8, time = 1:length(path8), ball = "Ball 8", color = "lightblue")


path_df <- rbind(path_df1, path_df2, path_df3, path_df4, path_df5, path_df6, path_df7, path_df8)
path_df <- merge(path_df, nodes, by = "node")

# Create a ggplot2 plot of the SEM
p <- ggplot() +
  geom_segment(data = edges_offset, aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.4, "cm")), size = 1, lineend = "round") +
  geom_point(data = nodes, aes(x = x, y = y), size = 12, shape = 21, fill = "white", color = "darkblue", stroke = 1.5) +
  geom_text(data = nodes, aes(x = x, y = y, label = label), vjust = 0.5, hjust = 0.5) +
  geom_point(data = path_df, aes(x = x, y = y, frame = time, color = color), size = 15, alpha = 0.35) +
  scale_color_manual(values = c("blue", "red", "pink", "lightblue")) +
  theme_void() + # Remove all axes and background
  theme(panel.background = element_rect(fill = "white", color = NA)) +
  # no legend
  theme(legend.position = "none")

# Animate the balls moving through the SEM
anim <- p + transition_time(time) + ease_aes('linear')
anim
# save this as a gif
anim_save("~/Dropbox/github_repos/crossLag_p/crossLagR/vignettes/presentations/clpm.gif", width = 800, height = 600, res = 96)


library(gganimate)
# Define the nodes for the SEM with latent factors
nodes <- data.frame(
  node = 1:12,
  label = c("X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4", "X5", "Y5", "X", "Y"),
  x = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 3, 3),  # Center Factor X and Factor Y above X3 and Y3
  y = c(2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 3, 0)   # Adjust y positions for factors
)

# Define the edges for the SEM, including latent factors and cross-lagged paths
edges <- data.frame(
  from = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12),
  to =   c(3, 4, 3, 4, 5, 6, 5, 6, 7, 8, 7, 8, 9, 10, 9, 10, 1, 3, 5, 7, 9, 2, 4, 6, 8, 10) # Added connections for X and Y
)

# Function to offset the arrowheads
offset_arrow <- function(x1, y1, x2, y2, offset = 0.10) {
  angle <- atan2(y2 - y1, x2 - x1)
  x1_new <- x1 + offset * cos(angle)
  y1_new <- y1 + offset * sin(angle)
  x2_new <- x2 - offset * cos(angle)
  y2_new <- y2 - offset * sin(angle)
  return(data.frame(x1 = x1_new, y1 = y1_new, x2 = x2_new, y2 = y2_new))
}

# Apply the offset to the edges
edges_offset <- do.call(rbind, apply(edges, 1, function(row) {
  offset_arrow(nodes$x[row["from"]], nodes$y[row["from"]],
               nodes$x[row["to"]], nodes$y[row["to"]])
}))

# Define the paths for the balls through the SEM
# Create paths originating from all nodes, including cross-lagged paths
path1 <- c(11, 1, 3, 5, 7, 9) # Factor X to X1 to X2 to X3 to X4 to X5
path2 <- c(12, 2, 4, 6, 8, 10) # Factor Y to Y1 to Y2 to Y3 to Y4 to Y5

path3 <- c(11, 1, 4, 5, 8, 9) # Factor X to X1 to X2 to X3 to X4 to X5
path4 <- c(12, 2, 3, 6, 7, 10) # Factor Y to Y1 to Y2 to Y3 to Y4 to Y5

path5 <- c(11, 3, 5, 7, 9) # Factor X to X1 to X2 to X3 to X4 to X5
path6 <- c(12, 4, 6, 8, 10) # Factor Y to Y1 to Y2 to Y3 to Y4 to Y5

path7 <- c(11, 5, 7, 9) # Factor X to X1 to X2 to X3 to X4 to X5
path8 <- c(12, 6, 8, 10) # Factor Y to Y1 to Y2 to Y3 to Y4 to Y5

path9 <- c(11, 7, 9) # Factor X to X1 to X2 to X3 to X4 to X5
path10 <- c(12,8, 10) # Factor Y to Y1 to Y2 to Y3 to Y4 to Y5


# Assign timing to ensure animations originate from all nodes
path_df1 <- data.frame(node = path1, time = 1:length(path1), ball = "Ball 1", color = "red")
path_df2 <- data.frame(node = path2, time = 1:length(path2), ball = "Ball 2", color = "blue")

path_df3 <- data.frame(node = path3, time = 1:length(path3), ball = "Ball 3", color = "red")
path_df4 <- data.frame(node = path4, time = 1:length(path4), ball = "Ball 4", color = "blue")

path_df5 <- data.frame(node = path5, time = 1:length(path5), ball = "Ball 5", color = "red")
path_df6 <- data.frame(node = path6, time = 1:length(path6), ball = "Ball 6", color = "blue")

path_df7 <- data.frame(node = path7, time = 1:length(path7), ball = "Ball 7", color = "red")
path_df8 <- data.frame(node = path8, time = 1:length(path8), ball = "Ball 8", color = "blue")

path_df9 <- data.frame(node = path9, time = 1:length(path9), ball = "Ball 9", color = "red")
path_df10 <- data.frame(node = path10, time = 1:length(path10), ball = "Ball 10", color = "blue")

path_df <- rbind(path_df1, path_df2, path_df3, path_df4, path_df5, path_df6, path_df7, path_df8,
                 path_df9, path_df10)
path_df <- merge(path_df, nodes, by = "node")

# Create a ggplot2 plot of the SEM
p <- ggplot() +
  geom_segment(data = edges_offset, aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.4, "cm")), size = 1, lineend = "round") +
  geom_point(data = nodes, aes(x = x, y = y), size = 12, shape = 21, stroke = 1.5, fill = "white") +
  geom_text(data = nodes, aes(x = x, y = y, label = label), vjust = 0.5, hjust = 0.5) +
  geom_point(data = path_df, aes(x = x, y = y, frame = time, color = color), size = 12, alpha = 0.55) +
  scale_color_identity() + # Use the colors defined in the path_df
  theme_void() + # Remove all axes and background
  theme(panel.background = element_rect(fill = "white", color = NA)) +
  theme(legend.position = "none") # No legend

# Animate the balls moving through the SEM
anim <- p + transition_time(time) + ease_aes('linear')
anim
anim_save("~/Dropbox/github_repos/crossLag_p/crossLagR/vignettes/presentations/riclpm.gif")



library(ggplot2)
library(gganimate)

# Define the nodes for the SEM with latent factors
nodes <- data.frame(
  node = 1:12,
  label = c("X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4", "X5", "Y5", "X", "Y"),
  x = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 3, 3),
  y = c(2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 3, 0)
)

# Define the edges for the SEM, including latent factors and cross-lagged paths
edges <- data.frame(
  from = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12),
  to =    c(3, 4, 3, 4, 5, 6, 5, 6, 7, 8, 7, 8, 9, 10, 9, 10, 1, 3, 5, 7, 9, 2, 4, 6, 8, 10)
)

# Function to offset the arrowheads
offset_arrow <- function(x1, y1, x2, y2, offset = 0.10) {
  angle <- atan2(y2 - y1, x2 - x1)
  x1_new <- x1 + offset * cos(angle)
  y1_new <- y1 + offset * sin(angle)
  x2_new <- x2 - offset * cos(angle)
  y2_new <- y2 - offset * sin(angle)
  return(data.frame(x1 = x1_new, y1 = y1_new, x2 = x2_new, y2 = y2_new))
}

# Apply the offset to the edges
edges_offset <- do.call(rbind, apply(edges, 1, function(row) {
  offset_arrow(nodes$x[row["from"]], nodes$y[row["from"]],
               nodes$x[row["to"]], nodes$y[row["to"]])
}))

# Define the paths for the balls through the SEM
path1 <- c(1, 3, 5, 7, 9) # X1 -> X2 -> X3 -> X4 -> X5
path2 <- c(2, 4, 6, 8, 10) # Y1 -> Y2 -> Y3 -> Y4 -> Y5
path3 <- c(1, 4, 5, 8, 9) # X1 -> X2 -> X3 -> X4 -> X5
path4 <- c(2, 3, 6, 7, 10) # Y1 -> Y2 -> Y3 -> Y4 -> Y5

# Assign timing to ensure animations originate from all nodes
path_df1 <- data.frame(node = path1, time = 1:length(path1), ball = "Ball 1", color = "red")
path_df2 <- data.frame(node = path2, time = 1:length(path2), ball = "Ball 2", color = "blue")
path_df3 <- data.frame(node = path3, time = 1:length(path3), ball = "Ball 1", color = "red")
path_df4 <- data.frame(node = path4, time = 1:length(path4), ball = "Ball 2", color = "blue")

# Create grey balls for X1-X5 and Y1-Y5 (static)
grey_balls_x <- data.frame(
  node = 1:5,
  time = 10,
  ball = paste0("Grey X", 1:5),
  color = "darkgreen"
)
grey_balls_y <- data.frame(
  node = 6:10,
  time = 10,
  ball = paste0("Grey Y", 1:5),
  color = "darkgreen"
)

path_df <- rbind(path_df1, path_df2, path_df3, path_df4,  grey_balls_x, grey_balls_y)
path_df <- merge(path_df, nodes, by = "node")

# Create a ggplot2 plot of the SEM
p <- ggplot() +
  geom_segment(data = edges_offset, aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.4, "cm")), size = 1, lineend = "round") +
  geom_point(data = nodes, aes(x = x, y = y), size = 12, shape = 21, stroke = 1.5, fill = "white") +
  geom_text(data = nodes, aes(x = x, y = y, label = label), vjust = 0.5, hjust = 0.5) +
  geom_point(data = path_df, aes(x = x, y = y, frame = time, color = color), size = 12, alpha = 0.55) +
  # Add static grey points
  geom_point(data = nodes[1:10, ], aes(x = x, y = y), size = 12, color = "darkgreen", alpha = 0.25) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA)) +
  theme(legend.position = "none")

# Animate the balls moving through the SEM
anim <- p + transition_time(time) + ease_aes('linear')
anim
# Render the animation with a larger number of frames (e.g., 200)
animate(anim, nframes = 200 ,loop = TRUE)

anim_save("~/Dropbox/github_repos/crossLag_p/crossLagR/vignettes/presentations/riclpm.gif")
