
bound = function(x, lower, upper) {
  pmin(pmax(lower, x), upper)
}

get.stochasticM = function(gstarM_astar, Y_preds1, Y_preds0) {
  Y_preds1*gstarM_astar + Y_preds0*(1 - gstarM_astar)
}

