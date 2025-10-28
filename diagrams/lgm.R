library(DiagrammeR)
# Linear Growth
lgm_diagram <- grViz("
digraph LGM {

  # Graph attributes
  graph [layout = neato, rankdir = TB, bgcolor = white]

  # Node attributes
  node [shape = box, style = filled, fillcolor = white, fontsize = 11]

  # Growth factors for X
  I_X [pos = '-3,2.75!', shape = circle, fillcolor = white, label = <I<sub>x</sub>>, width = 0.5, height = 0.5]
  S_X [pos = '-2.75,2!', shape = circle, fillcolor = white, label = <S<sub>x</sub>>, width = 0.5, height = 0.5]

  # Observed X variables
  X1 [pos = '-2,3!', label = 'x1', width = 0.6, height = 0.3]
  X2 [pos = '-1,3!', label = 'x2', width = 0.6, height = 0.3]
  X3 [pos = '0,3!', label = 'x3', width = 0.6, height = 0.3]
  X4 [pos = '1,3!', label = 'x4', width = 0.6, height = 0.3]
  X5 [pos = '2,3!', label = 'x5', width = 0.6, height = 0.3]

  # Latent X variables (optional - can be removed if direct paths from I_X/S_X to observed)
  wX1 [pos = '-2,1!', shape = circle, fillcolor = white, label = 'q1']
  wX2 [pos = '-1,1!', shape = circle, fillcolor = white, label = 'q2']
  wX3 [pos = '0,1!', shape = circle, fillcolor = white, label = 'q3']
  wX4 [pos = '1,1!', shape = circle, fillcolor = white, label = 'q4']
  wX5 [pos = '2,1!', shape = circle, fillcolor = white, label = 'q5']

  # Latent Y variables
  wY1 [pos = '-2,-1!', shape = circle, fillcolor = white, label = 'p1']
  wY2 [pos = '-1,-1!', shape = circle, fillcolor = white, label = 'p2']
  wY3 [pos = '0,-1!', shape = circle, fillcolor = white, label = 'p3']
  wY4 [pos = '1,-1!', shape = circle, fillcolor = white, label = 'p4']
  wY5 [pos = '2,-1!', shape = circle, fillcolor = white, label = 'p5']

  # Observed Y variables
  Y1 [pos = '-2,-3!', label = 'y1', width = 0.6, height = 0.3]
  Y2 [pos = '-1,-3!', label = 'y2', width = 0.6, height = 0.3]
  Y3 [pos = '0,-3!', label = 'y3', width = 0.6, height = 0.3]
  Y4 [pos = '1,-3!', label = 'y4', width = 0.6, height = 0.3]
  Y5 [pos = '2,-3!', label = 'y5', width = 0.6, height = 0.3]

  # Growth factors for Y
  I_Y [pos = '-3,-2!', shape = circle, fillcolor = white, label = <I<sub>y</sub>>, width = 0.5, height = 0.5]
  S_Y [pos = '-2.75,-2.5!', shape = circle, fillcolor = white, label = <S<sub>y</sub>>, width = 0.5, height = 0.5]

  # INTERCEPT FACTOR PATHS (all fixed to 1)
  I_X -> wX1 [label = '1', arrowsize = 0.5]
  I_X -> wX2 [label = '1', arrowsize = 0.5]
  I_X -> wX3 [label = '1', arrowsize = 0.5]
  I_X -> wX4 [label = '1', arrowsize = 0.5]
  I_X -> wX5 [label = '1', arrowsize = 0.5]

  I_Y -> wY1 [label = '1', arrowsize = 0.5]
  I_Y -> wY2 [label = '1', arrowsize = 0.5]
  I_Y -> wY3 [label = '1', arrowsize = 0.5]
  I_Y -> wY4 [label = '1', arrowsize = 0.5]
  I_Y -> wY5 [label = '1', arrowsize = 0.5]

  # SLOPE FACTOR PATHS (linear growth: 0, 1, 2, 3, 4)
  S_X -> wX1 [label = '0', arrowsize = 0.5]
  S_X -> wX2 [label = '1', arrowsize = 0.5]
  S_X -> wX3 [label = '2', arrowsize = 0.5]
  S_X -> wX4 [label = '3', arrowsize = 0.5]
  S_X -> wX5 [label = '4', arrowsize = 0.5]

  S_Y -> wY1 [label = '0', arrowsize = 0.5]
  S_Y -> wY2 [label = '1', arrowsize = 0.5]
  S_Y -> wY3 [label = '2', arrowsize = 0.5]
  S_Y -> wY4 [label = '3', arrowsize = 0.5]
  S_Y -> wY5 [label = '4', arrowsize = 0.5]

  # MEASUREMENT MODEL (latent to observed)
  wX1 -> X1 [label = '', arrowsize = 0.5]
  wX2 -> X2 [label = '', arrowsize = 0.5]
  wX3 -> X3 [label = '', arrowsize = 0.5]
  wX4 -> X4 [label = '', arrowsize = 0.5]
  wX5 -> X5 [label = '', arrowsize = 0.5]

  wY1 -> Y1 [label = '', arrowsize = 0.5]
  wY2 -> Y2 [label = '', arrowsize = 0.5]
  wY3 -> Y3 [label = '', arrowsize = 0.5]
  wY4 -> Y4 [label = '', arrowsize = 0.5]
  wY5 -> Y5 [label = '', arrowsize = 0.5]

  # GROWTH FACTOR CORRELATIONS
  I_X -> I_Y [dir = both, color = black, style = dashed, arrowsize = 0.5]
  I_X -> S_Y [dir = both, color = black, style = dashed, arrowsize = 0.5]
  S_X -> I_Y [dir = both, color = black, style = dashed, arrowsize = 0.5]
  S_X -> S_Y [dir = both, color = black, style = dashed, arrowsize = 0.5]

  # RESIDUALS ON OBSERVED VARIABLES (measurement error)
  X1 -> X1 [dir = both, color = black, style = dashed, label = '', tailport = n, headport = n, arrowsize = 0.5]
  X2 -> X2 [dir = both, color = black, style = dashed, label = '', tailport = n, headport = n, arrowsize = 0.5]
  X3 -> X3 [dir = both, color = black, style = dashed, label = '', tailport = n, headport = n, arrowsize = 0.5]
  X4 -> X4 [dir = both, color = black, style = dashed, label = '', tailport = n, headport = n, arrowsize = 0.5]
  X5 -> X5 [dir = both, color = black, style = dashed, label = '', tailport = n, headport = n, arrowsize = 0.5]

  Y1 -> Y1 [dir = both, color = black, style = dashed, label = '', tailport = s, headport = s, arrowsize = 0.5]
  Y2 -> Y2 [dir = both, color = black, style = dashed, label =  '', tailport = s, headport = s, arrowsize = 0.5]
  Y3 -> Y3 [dir = both, color = black, style = dashed, label = '', tailport = s, headport = s, arrowsize = 0.5]
  Y4 -> Y4 [dir = both, color = black, style = dashed, label = '', tailport = s, headport = s, arrowsize = 0.5]
  Y5 -> Y5 [dir = both, color = black, style = dashed, label = '', tailport = s, headport = s, arrowsize = 0.5]

  # RESIDUALS ON GROWTH FACTORS (factor variances)
  I_X -> I_X [dir = both, color = black, style = dashed, label = '', labelangle = 135, labeldistance = 1.5, tailport = w, headport = s, arrowsize = 0.5]
  S_X -> S_X [dir = both, color = black, style = dashed, label = '', labelangle = 135, labeldistance = 1.5, tailport = w, headport = s, arrowsize = 0.5]
  I_Y -> I_Y [dir = both, color = black, style = dashed, label = '', labelangle = 135, labeldistance = 1.5, tailport = w, headport = n, arrowsize = 0.5]
  S_Y -> S_Y [dir = both, color = black, style = dashed, label = '', labelangle = 135, labeldistance = 1.5, tailport = w, headport = n, arrowsize = 0.5]

}
")
# Display the diagram
print("Dual Latent Growth Model - 5 Wave")
lgm_diagram
