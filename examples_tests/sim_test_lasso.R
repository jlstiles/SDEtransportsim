res1 = res
res2 = res
res = append(res, res1)
res = append(res, res2)
length(res)

test_zero = lapply(res, FUN = function(x) x$res[1] == 0)
sum(unlist(test_zero))

res = res[!unlist(test_zero)]
length(res)
data = res[[40]]$data
set.seed(res[[40]]$p)

data = data.frame(W1 = fucked$W1, W2 = fucked$W2, A = fucked$A, Z = fucked$Z, 
                  M = fucked$M, Y = fucked$Y, weights = fucked$weights)

test_zero = SDE_tmle_lasso(data, formsNT, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                           transport = FALSE, pooled = FALSE, gstar_S = 0, truth = truth_NT)

test_zero = SDE_tmle_lasso(data, formsNP, RCT = 0.5, Wnames = Wnames, 
                           Wnamesalways = Wnamesalways, transport = TRUE, 
                           pooled = FALSE, gstar_S = 0, truth = truth)

debug(get_gstarM_lasso)
debug(SDE_tmle_lasso)

!unlist(test_zero)
results = lapply(res[!unlist(test_zero)], FUN = function(x) unlist(x[1]))
results = lapply(res, FUN = function(x) x$res)
results = do.call(rbind, results)

mean(results[,2] <= results[,4] & results[,4] <= results[,3])
mean(results[,6] <= results[,8] & results[,8] <= results[,7])

test_zero = lapply(res, FUN = function(x) x[[1]] == 0)

results1 = lapply(res[!unlist(test_zero)], FUN = function(x) unlist(x[1:8]))
results1 = do.call(rbind, results1)

results1[1:5,]
mean(results1[,2] <= results1[,4] & results1[,4] <= results1[,3])
mean(results1[,6] <= results1[,8] & results1[,8] <= results1[,7])

test_zero = lapply(res, FUN = function(x) x[[1]] == 0)

results2 = lapply(res[!unlist(test_zero)], FUN = function(x) unlist(x[1:8]))
results2 = do.call(rbind, results2)

results2[1:5,]
mean(results2[,2] <= results2[,4] & results2[,4] <= results2[,3])
mean(results2[,6] <= results2[,8] & results2[,8] <= results2[,7])

results = rbind(results, results1, results2)
nrow(results)
mean(results[,2] <= results[,4] & results[,4] <= results[,3])
mean(results[,6] <= results[,8] & results[,8] <= results[,7])
