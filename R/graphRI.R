#' @title riclpm_graph
#' @description This function draws a tidySEM graph of a RI-CLPM model
#'
#' Create a simulated dataset of a cross-lagged model with a specified number of waves and structural parameters.
#'
#' @param fit A fitted lavaan model object.
#' @param waves Number of waves in the fitted model.
#'
#' @return RI-CLPM SEM plot
#'
#' @export
#'
#'
#'
riclpm_graph = function(model = fit, waves = c(3,4,5,6)){
  if (waves <3 | waves >6){
    print('Please choose a wave length of between 3 and 6. Otherwise, the graph gets muddled')
  }
  if(waves == 6){
    layout = get_layout(
      ""  ,     ""  ,   "",  "bx", "",   "" ,      "", "", "",
      ""  ,    ""  ,  "x1", "x2", "x3", "x4" ,   "x5" , "x6" ,      "",
      ""  ,    ""  ,   "p1",  "p2", "p3", "p4" , "p5" ,  "p6" ,      "",
      ""  ,    ""  ,   "",  "", "", "" ,     "", "", "",
      ""  ,    ""  ,   "q1",  "q2", "q3", "q4" , "q5" , "q6",     "",
      ""  ,    ""  ,  "y1", "y2", "y3", "y4" ,   "y5" ,   "y6",   "",
      ""  ,     ""  ,   "",  "by", "",   "" ,      "",   "", "" ,
      rows =7)
  }
  else if(waves == 5){
    layout = get_layout(
      ""  ,     ""  ,   "",  "bx", "",   "" ,      "", "",
      ""  ,    ""  ,  "x1", "x2", "x3", "x4" ,   "x5" ,      "",
      ""  ,    ""  ,   "p1",  "p2", "p3", "p4" , "p5" ,      "",
      ""  ,    ""  ,   "",  "", "", "" ,     "", "",
      ""  ,    ""  ,   "q1",  "q2", "q3", "q4" , "q5" ,      "",
      ""  ,    ""  ,  "y1", "y2", "y3", "y4" ,   "y5" ,      "",
      ""  ,     ""  ,   "",  "by", "",   "" ,      "",   "",
      rows =7)
  }
  else if(waves == 4){
    layout = get_layout(
      ""  ,     ""  ,   "",  "bx", "",   "" ,      "",
      ""  ,    ""  ,  "x1", "x2", "x3", "x4" ,        "",
      ""  ,    ""  ,   "p1",  "p2", "p3", "p4" ,      "",
      ""  ,    ""  ,   "",  "", "", "" ,     "",
      ""  ,     ""  ,   "q1",  "q2", "q3", "q4" ,     "",
      ""  ,     ""  ,   "y1",  "y2", "y3", "y4" ,     "",
      ""  ,     ""  ,   "",  "by", "",   "" ,      "",
      rows =7)
  }
  else if(waves == 3){
    layout = get_layout(
      ""  ,     ""  ,   "",  "bx", "",
      ""  ,    ""  ,  "x1", "x2", "x3",        "",
      ""  ,    ""  ,   "p1",  "p2", "p3",      "",
      ""  ,    ""  ,   "",  "", "", "" ,     "",
      ""  ,     ""  ,   "q1",  "q2", "q3",      "",
      ""  ,     ""  ,   "y1",  "y2", "y3",      "",
      ""  ,     ""  ,   "",  "by", "",        "",
      rows =7)
  }
  graphData = tidySEM::prepare_graph(model, layout = layout)

  tidySEM::edges(graphData) %>%
    mutate(
      label = ifelse(op == "~~", "", "")) -> tidySEM::edges(graphData)

  graphData %>%
    tidySEM::edit_graph({
      label_color <- "black"  # Set label color to blue
    })  %>% plot()

}

# # Not run
#   simulate_riclpm(waves = 4, sample.nobs = 15000)$data
#      %>%   dplyr::select(-kappa, -omega, -U))
#
## This is the way to structure the plot
# lavaan::lavaan(model_syntax_cl  pm(waves = 4, model_type = "ri-clpm"),
#                      simulate_riclpm(waves = 4, sample.nobs = 15000)$data
#                      %>%   dplyr::select(-U)) %>% ri_graph(waves = 4)

