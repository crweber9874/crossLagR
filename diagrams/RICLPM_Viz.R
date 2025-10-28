library(DiagrammeR)
# RI-CLPM with proper residuals as curved self-loops and correlated errors - 5 WAVE
riclpm_diagram <- grViz("
digraph RICLPM {

  # Graph attributes
  graph [layout = neato, rankdir = TB, bgcolor = white]

  # Node attributes
  node [shape = box, style = filled, fillcolor = white, fontsize = 11]

  RI_X [pos = '-3,2!', shape = circle, fillcolor = white, label = <ξ<sub>x</sub>>,  width = 0.5, height = 0.5]
  X1 [pos = '-2,2!', label = 'x1', width = 0.6, height = 0.3]
  X2 [pos = '-1,2!', label = 'x2', width = 0.6, height = 0.3]
  X3 [pos = '0,2!', label = 'x3', width = 0.6, height = 0.3]
  X4 [pos = '1,2!', label = 'x4', width = 0.6, height = 0.3]
  X5 [pos = '2,2!', label = 'x5', width = 0.6, height = 0.3]

  # Row 2: Within-person X latents (middle-upper)
  wX1 [pos = '-2,1!', shape = circle, fillcolor = white, label = 'q1']
  wX2 [pos = '-1,1!', shape = circle, fillcolor = white, label = 'q2']
  wX3 [pos = '0,1!', shape = circle, fillcolor = white, label = 'q3']
  wX4 [pos = '1,1!', shape = circle, fillcolor = white, label = 'q4']
  wX5 [pos = '2,1!', shape = circle, fillcolor = white, label = 'q5']

  # Row 3: Within-person Y latents (middle-lower)
  wY1 [pos = '-2,-1!', shape = circle, fillcolor = white, label = 'p1']
  wY2 [pos = '-1,-1!', shape = circle, fillcolor = white, label = 'p2']
  wY3 [pos = '0,-1!', shape = circle, fillcolor = white, label = 'p3']
  wY4 [pos = '1,-1!', shape = circle, fillcolor = white, label = 'p4']
  wY5 [pos = '2,-1!', shape = circle, fillcolor = white, label = 'p5']

  # Row 4: observed Y variables and RI_Y (bottom)
  Y1 [pos = '-2,-2!', label = 'y1', width = 0.6, height = 0.3]
  Y2 [pos = '-1,-2!', label = 'y2', width = 0.6, height = 0.3]
  Y3 [pos = '0,-2!', label = 'y3', width = 0.6, height = 0.3]
  Y4 [pos = '1,-2!', label = 'y4', width = 0.6, height = 0.3]
  Y5 [pos = '2,-2!', label = 'y5', width = 0.6, height = 0.3]
  RI_Y [pos = '-3,-2!', shape = circle, fillcolor = white, label = <ξ<sub>y</sub>>, width = 0.5, height = 0.5]


  # Measurement model paths (RI to observed)
  RI_X -> wX1 [label = '', arrowsize = 0.5]
  RI_X -> wX2 [label = '', arrowsize = 0.5]
  RI_X -> wX3 [label = '',  arrowsize = 0.5]
  RI_X -> wX4 [label = '', arrowsize = 0.5]
  RI_X -> wX5 [label = '', arrowsize = 0.5]
  RI_Y -> wY1 [label = '', arrowsize = 0.5]
  RI_Y -> wY2 [label = '', arrowsize = 0.5]
  RI_Y -> wY3 [label = '', arrowsize = 0.5]
  RI_Y -> wY4 [label = '', arrowsize = 0.5]
  RI_Y -> wY5 [label = '', arrowsize = 0.5]

  # Measurement model paths (w to observed)
  wX1 -> X1 [label = '1', arrowsize = 0.5]
  wX2 -> X2 [label = '1', arrowsize = 0.5]
  wX3 -> X3 [label = '1', arrowsize = 0.5]
  wX4 -> X4 [label = '1', arrowsize = 0.5]
  wX5 -> X5 [label = '1', arrowsize = 0.5]
  wY1 -> Y1 [label = '1', arrowsize = 0.5]
  wY2 -> Y2 [label = '1', arrowsize = 0.5]
  wY3 -> Y3 [label = '1', arrowsize = 0.5]
  wY4 -> Y4 [label = '1', arrowsize = 0.5]
  wY5 -> Y5 [label = '1', arrowsize = 0.5]

  # Autoregressive paths
  wX1 -> wX2 [color = black, penwidth = 2, arrowsize = 0.5]
  wX2 -> wX3 [color = black, penwidth = 2, arrowsize = 0.5]
  wX3 -> wX4 [color = black, penwidth = 2, arrowsize = 0.5]
  wX4 -> wX5 [color = black, penwidth = 2, arrowsize = 0.5]
  wY1 -> wY2 [color = black, penwidth = 2, arrowsize = 0.5]
  wY2 -> wY3 [color = black, penwidth = 2, arrowsize = 0.5]
  wY3 -> wY4 [color = black, penwidth = 2, arrowsize = 0.5]
  wY4 -> wY5 [color = black, penwidth = 2, arrowsize = 0.5]

  # Cross-lagged paths with offset gamma labels and subscripts
  wX1 -> wY2 [color = black, penwidth = 2, arrowsize = 0.5]
  wX2 -> wY3 [color = black, penwidth = 2, arrowsize = 0.5]
  wX3 -> wY4 [color = black, penwidth = 2, arrowsize = 0.5]
  wX4 -> wY5 [color = black, penwidth = 2, arrowsize = 0.5]
  wY1 -> wX2 [color = black, penwidth = 2, arrowsize = 0.5]
  wY2 -> wX3 [color = black, penwidth = 2, arrowsize = 0.5]
  wY3 -> wX4 [color = black, penwidth = 2, arrowsize = 0.5]
  wY4 -> wX5 [color = black, penwidth = 2, arrowsize = 0.5]

  # Within-time correlations - ONLY these are slightly rounded
  wX1 -> wY1 [dir = both, color = black, style = dashed, splines = curved, constraint = false, arrowsize = 0.5 ]
  wX2 -> wY2 [dir = both, color = black, style = dashed, splines = curved, constraint = false, arrowsize = 0.5 ]
  wX3 -> wY3 [dir = both, color = black, style = dashed, splines = curved, constraint = false, arrowsize = 0.5 ]
  wX4 -> wY4 [dir = both, color = black, style = dashed, splines = curved, constraint = false, arrowsize = 0.5 ]
  wX5 -> wY5 [dir = both, color = black, style = dashed, splines = curved, constraint = false, arrowsize = 0.5 ]

  # Random intercept correlation - ONLY this is slightly rounded
  RI_X -> RI_Y [dir = both, color = black, style = dashed, penwidth = 1, splines = curved, constraint = false,  arrowsize = 0.5]


  # Self-loops connected to bottom of each node using ports
  wX2 -> wX2 [dir = both, color = black, style = dashed, label = '', tailport = sw, headport = s , arrowsize = 0.5]
  wX3 -> wX3 [dir = both, color = black, style = dashed, label = '', tailport = sw, headport = s, arrowsize = 0.5]
  wX4 -> wX4 [dir = both, color = black, style = dashed, label = '', tailport = sw, headport = s, arrowsize = 0.5]
  wX5 -> wX5 [dir = both, color = black, style = dashed, label = '', tailport = sw, headport = s, arrowsize = 0.5]
  wY2 -> wY2 [dir = both, color = black, style = dashed, label = '', tailport = nw, headport = n, arrowsize = 0.5]
  wY3 -> wY3 [dir = both, color = black, style = dashed, label = '', tailport = nw, headport = n, arrowsize = 0.5]
  wY4 -> wY4 [dir = both, color = black, style = dashed, label = '', tailport = nw, headport = n, arrowsize = 0.5]
  wY5 -> wY5 [dir = both, color = black, style = dashed, label = '', tailport = nw, headport = n, arrowsize = 0.5]

  # Self-loops connected to bottom of each node using ports
  X1 -> X1 [dir = both, color = black, style = dashed, label = '0', tailport = n, headport = n , arrowsize = 0.5]
  X2 -> X2 [dir = both, color = black, style = dashed, label = '0', tailport = n, headport = n , arrowsize = 0.5]
  X3 -> X3 [dir = both, color = black, style = dashed, label = '0', tailport = n, headport = n , arrowsize = 0.5]
  X4 -> X4 [dir = both, color = black, style = dashed, label = '0', tailport = n, headport = n , arrowsize = 0.5]
  X5 -> X5 [dir = both, color = black, style = dashed, label = '0', tailport = n, headport = n , arrowsize = 0.5]

  Y1 -> Y1 [dir = both, color = black, style = dashed, label = '0', tailport = s, headport = s , arrowsize = 0.5]
  Y2 -> Y2 [dir = both, color = black, style = dashed, label = '0', tailport = s, headport = s , arrowsize = 0.5]
  Y3 -> Y3 [dir = both, color = black, style = dashed, label = '0', tailport = s, headport = s , arrowsize = 0.5]
  Y4 -> Y4 [dir = both, color = black, style = dashed, label = '0', tailport = s, headport = s , arrowsize = 0.5]
  Y5 -> Y5 [dir = both, color = black, style = dashed, label = '0', tailport = s, headport = s , arrowsize = 0.5]


  # Residuals for random intercepts (epsilon)
  RI_X -> RI_X [dir = both, color = black, style = dashed, label = 'ε', labelangle = 135, labeldistance = 1.5, tailport = w, headport = s, arrowsize = 0.5, width = 0.6, height = 0.3]
  RI_Y -> RI_Y [dir = both, color = black, style = dashed, label = 'ε', labelangle = 135, labeldistance = 1.5, tailport = w, headport = n, arrowsize = 0.5, width = 0.6, height = 0.3]

}
")
# Display the diagram
print("5-Wave RI-CLPM with proper curved residual self-loops")
riclpm_diagram
