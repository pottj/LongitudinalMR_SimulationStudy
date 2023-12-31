mrest_me = function(mrob, max_iter = 100, no_ini = 1){
  # all credits to Steve and his collegue! This function was copy-pasted from github (check with Steve for reference!)
  bxhat = as.matrix(mrob@betaX)
  byhat = mrob@betaY
  sebx = as.matrix(mrob@betaXse)
  seby = mrob@betaYse
  p = length(byhat)
  K = dim(bxhat)[2]
  S = diag(seby^-2)
  SigX = lapply(1:p, function(j){diag(sebx[j, ]^2, length(sebx[j, ]))})
  
  l = matrix(nrow = max_iter, ncol = no_ini)
  thest = matrix(nrow = K, ncol = no_ini)
  for (k in 1:no_ini){
    bxtilde = t(matrix(sapply(1:p, function(j){mvrnorm(1, bxhat[j, ], SigX[[j]])}), ncol = p))
    for (i in 1:100){
      thest[, k] = solve(t(bxtilde) %*% S %*% bxtilde, t(bxtilde) %*% S %*% byhat)
      l[i, k] = -0.5 * sum(sapply(1:p, function(j){
        (byhat[j] - t(bxhat[j, ]) %*% thest[, k])^2 / (seby[j]^2 + t(thest[, k]) %*% SigX[[j]] %*% thest[, k])
      }))
      for (j in 1:p){
        bxtilde[j, ] = t(solve(thest[, k] %*% t(thest[, k]) / seby[j]^2 + solve(SigX[[j]]), byhat[j] * thest[, k] / seby[j]^2 + solve(SigX[[j]], bxhat[j, ])))
      }
      if (i > 1){
        if (abs(l[i] - l[(i-1)]) < 1e-4) {break}
      }
    }
  }
  k0 = which.max(apply(as.matrix(l[is.na(l[, 1]) == F, ]), 2, max))
  th = thest[, k0]
  
  v = sapply(1:p, function(j){seby[j]^2 + t(th) %*% SigX[[j]] %*% th})
  e = sapply(1:p, function(j){byhat[j] - bxhat[j, ] %*% th})
  t = sapply(1:p, function(j){e[j] / sqrt(v[j])})
  
  dt = matrix(sapply(1:p, function(j){(-v[j] * bxhat[j, ] - e[j] * SigX[[j]] %*% th) / v[j]^(3/2)}), ncol = p)
  B = dt %*% t(dt)
  
  dt2 = vector(length = p, mode = "list")
  for (j in 1:p){
    dt2[[j]] = matrix(nrow = K, ncol = K)
    S = SigX[[j]] %*% th
    for (k in 1:K){
      for (l in 1:K){
        dt2[[j]][k, l] = v[j]^(-3/2) * (-2 * S[l] * bxhat[j, k] + S[k] * bxhat[j, l] - e[j] * SigX[[j]][k, l]) +
          3 * v[j]^(-5/2) *(v[j] * bxhat[j, k] + e[j] * S[k]) * S[l]
      }
    }
  }
  a = Reduce('+', lapply(1:p, function(j){c(t[j]) * dt2[[j]]}))
  A = (dt %*% t(dt) + a)
  
  Var = solve(A, B) %*% t(solve(A))
  return(list("thest" = th, "l" = l, "Var" = Var))
}
