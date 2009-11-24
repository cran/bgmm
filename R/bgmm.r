predict.mModel <- function(object, X, knowns=NULL, B=NULL, P=NULL, ...) {
  # densities for all components
  if (is.null(dim(X)) || is.data.frame(X)) X = as.matrix(X)
  if (!is.null(knowns)) {
      if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
      X = rbind(knowns, X)         # change
  }
  if ((is.null(B) & is.null(P) & !is.null(knowns)) | ({!is.null(B) | !is.null(P)} & is.null(knowns))) {
     stop("If knowns are specified there should be also B or P specified as well!")
  }
  fik <- matrix(0, nrow(X), object$k)
  rownames(fik) = rownames(X)
  for (i in 1:object$k) {
      if (object$d > 1) {
         ss = svd(object$cvar[i,,])
         matc = t(ss$u) %*% diag(ss$d^(-1/2)) %*% ss$v
         tx = apply(X, 1, get("-"), object$mu[i,,drop=F])
         fik[,i] <- exp(-colSums((matc %*% tx)^2)/2)/sqrt(prod(2*pi*ss$d))
       } else {
         tx = apply(X, 1, get("-"), object$mu[i,,drop=F])
         fik[,i] <- exp(-(tx^2)/(2*object$cvar[i,,]))/sqrt(prod(2*pi*object$cvar[i,,]))
       }
  }
  b.pi <- repeat.rows(object$pi, nrow(X))
  if (!is.null(B))          # change
        b.pi[1:nrow(knowns),] = B
#        b.pi[nrow(X) - (nrow(knowns):1) +1,] = B
  if (!is.null(P))          # change 
        b.pi[1:nrow(knowns),] = P * b.pi[1:nrow(knowns),]
#        b.pi[nrow(X) - (nrow(knowns):1) +1,] = P * b.pi[nrow(X) - (nrow(knowns):0),]

  tik =  t(apply(fik * b.pi, 1, normalize))
  list(tij = tik, class=get.labels.from.beliefs(tik))
}


supervised <- function(knowns, class=NULL, k=length(unique(class)), B=NULL, P=NULL,clusterAssigment=TRUE, model.structure=getModelStructure()) {
  if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
  if (is.null(class)) {
    if (!is.null(B)) {
      class = apply(B, 1, which.max)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else if (!is.null(P)){
      class = apply(P, 1, which.max)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else 
      stop("Argument class need to be specified")
  }
  result = init.model.params.knowns(knowns, class, k, ncol(knowns)) 

  # new means 
  if  (model.structure$mean!="D") {
    result$mu = repeat.rows(colMeans(knowns), k)
  }
  # are variances equal?
  if (model.structure$between=="E") {
      # averaging among clusters
     ncvar = matrix(0, result$d, result$d)
     for (i in 1:k) 
        ncvar = ncvar + result$cvar[i, , ] * result$pi[i]
     for (i in 1:k) 
        result$cvar[i, , ] = ncvar
  }
  if (model.structure$within=="E" && ncol(knowns)>1) {
      # averaging among variables
     for (i in 1:k) {
        ndiag = sum(diag(result$cvar[i, , ]))
        sdiag = ndiag/(ncol(knowns))
        result$cvar[i, , ] = min(sdiag, (sum(result$cvar[i, , ])-ndiag)/(ncol(knowns)*(ncol(knowns)-1)))
        diag(result$cvar[i, , ]) = sdiag
     }
  }
  # are covariances equal to 0?
  if (model.structure$cov=="0") {
   # covariances are equal to 0
   for (i in 1:k) 
      result$cvar[i, , ] = diag(diag(result$cvar[i, , ]), nrow=ncol(knowns))
  }  
  result$m = nrow(knowns)
  result$n = nrow(knowns)
  result$k = k
  result$d = ncol(knowns)
  result$knowns = knowns
  result$class = class
  class(result) = c("supervisedModel", "mModel")
  result
}

semisupervised <- function(X, knowns, class=NULL, k=ifelse(!is.null(class),length(unique(class)),ifelse(!is.null(B),ncol(B),ncol(P))),B=NULL,P=NULL, ..., clusterAssigment=TRUE, all.possible.permutations=FALSE) {
  if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
  if (is.null(dim(X)) || is.data.frame(X)) X = as.matrix(X)
  if (is.null(class)) {
    if (!is.null(B)) {
      class = apply(B, 1, which.max)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else if (!is.null(P)){
      class = apply(B, 1, which.max)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else 
      stop("Argument class need to be specified")
  }
  if (ncol(X) != ncol(knowns))  
      stop("number of columns in X and knowns must agree")
  result = soft(X, knowns, get.simple.beliefs(class, k=k, b.min=0), k=k, ..., all.possible.permutations=all.possible.permutations) 
  result$X = X
  result$knowns = knowns
  result$class = class
  class(result) = c("semisupervisedModel", "mModel")
  result
}

belief <- function(X, knowns, B=NULL, k=ifelse(!is.null(B),ncol(B),ifelse(!is.null(P),ncol(P),length(unique(class)))), P=NULL, class=map(B), init.params=init.model.params(X, knowns, B=B, P=P, class=class, k=k), model.structure=getModelStructure(), stop.likelihood.change=10^-5, stop.max.nsteps=100, trace=FALSE, b.min=0.025, clusterAssigment=TRUE, all.possible.permutations=FALSE) {
  if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
  if (is.null(dim(X)) || is.data.frame(X)) X = as.matrix(X)
  if (is.null(B)) {
    if (!is.null(P)) {
      B=P
      if (is.null(k) | k < 2)
          k=ncol(B)
    } else if (!is.null(class)){
      B = get.simple.beliefs(class, b.min=b.min)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else 
      stop("Argument B need to be specified")
  }
  if (k > ncol(B))
    B = cbind(B, matrix(0,nrow(B),k - ncol(B)))
  if (ncol(X) != ncol(knowns))  
      stop("number of columns in X and knowns must agree")
  init.params$B = B
  init.params$m = nrow(knowns)
  init.params$n = nrow(knowns) + nrow(X)
  init.params$k = k
  init.params$d = ncol(X)
  result = bgmm.internal(internal.funct=belief.internal, X=rbind(knowns, X), init.params=init.params, model.structure=model.structure, stop.likelihood.change=stop.likelihood.change, stop.max.nsteps=stop.max.nsteps, trace=trace, all.possible.permutations=all.possible.permutations)
  result$X = X
  result$knowns = knowns
  result$B = B
#  result$likelihood = loglikelihood.mModel(result, X)
  result$model.structure = model.structure
  class(result) = c("beliefModel", "mModel")
  result
}

#
# do we need to consider all possible permutations?
#
bgmm.internal <- function(internal.funct=belief.internal, init.params, ..., all.possible.permutations=all.possible.permutations) {
  if (all.possible.permutations) {
    require(combinat)
    resmax = NULL
    likemax = -Inf
    lpermut = permn(seq_along(init.params$mu))
    for (perm in lpermut) {
       tmodel.params = init.params
       tmodel.params$mu = tmodel.params$mu[perm,,drop=F]
       tmodel.params$cvar = tmodel.params$cvar[perm,,,drop=F]
       tmodel.params$pi = tmodel.params$pi[perm,drop=F]
       tmpr = internal.funct(..., model.params=tmodel.params)
       if (tmpr$likelihood > likemax) {
            likemax = tmpr$likelihood
            resmax = tmpr
       }
    }
    return(resmax)
  } 
  return(internal.funct(..., model.params=init.params))
}
 

# bgmm.internal.call

belief.internal <- function(X, model.params, model.structure, stop.likelihood.change=10^-5, stop.max.nsteps=100, trace=F) {
  prev.likelihood = -Inf
  n.steps = 0
  stopP = FALSE
  repeat {
    n.steps = n.steps +1
    tmp = bgmm.e.step(X, model.params) 
    model.params = bgmm.m.step(X, model.params, model.structure, tmp$tik)
    if (stopP)
          break
    if (trace) {
      cat("step:          ", n.steps, "\n likelihood:   ", tmp$log.likelihood, "\n change:       ", tmp$log.likelihood - prev.likelihood, "\n\n")
    }
    if ((abs(tmp$log.likelihood - prev.likelihood)/ifelse(is.infinite(prev.likelihood), 1,  (1+abs(prev.likelihood))) < stop.likelihood.change) || 
        (n.steps >= stop.max.nsteps)) {
          model.params$likelihood = tmp$log.likelihood
          stopP = TRUE
      }
    prev.likelihood = tmp$log.likelihood
  }
  model.params$likelihood = prev.likelihood
  model.params$n.steps = n.steps
  model.params
}

#
# usefull functions
#

normalize <- function(x) x/sum(x)

repeat.rows <- function(x, k) matrix(x, k, length(x), byrow=T)

get.labels.from.beliefs <- function(B) {apply(B, 1, function(x) which.max(x)[1])}
map = get.labels.from.beliefs


determinant.numeric <- function (x, logarithm = TRUE, ...) {list(modulus = ifelse(logarithm,log(x),x), sign=sign(x))}

loglikelihood.mModel <- function(model, X) {
  # densities for all components
  lfik <- matrix(0, nrow(X), model$k)
  for (i in 1:model$k) {
      if (model$d > 1) {
         ss = svd(model$cvar[i,,])
         matc = t(ss$u) %*% diag(ss$d^(-1/2)) %*% ss$v
         tx = apply(X, 1, get("-"), model$mu[i,,drop=F])
         lfik[,i] <-  -colSums((matc %*% tx)^2)/2 - sum(log(2*pi*ss$d))/2
       } else {
        lfik[,i] <-   dnorm(X, model$mu[i,,drop=F], sqrt(model$cvar[i,,]), log=T)
       }
  }
  fb.ik = exp(lfik + repeat.rows(log(model$pi), nrow(X)))
  sum(log(rowSums(fb.ik)))
}



getGIC <- function(model, p=2) {
  if (p=="AIC") p = 2
  if (p=="BIC") p = log(max(model$n - model$m,1))
 -2*loglikelihood.mModel(model, model$X) + getDF(model)*p
}

# degress of fredom
getDF <- function(model) {
  k = model$k
  d = model$d
  ifelse(model$model.structure$mean=="D", k*d, d) + # mean value
                    ifelse(model$model.structure$between=="D", k, 1) * # variance between clusters
                    (ifelse(model$model.structure$within=="D", d, 1)+ # diagonal
                     ifelse(model$model.structure$cov=="0", 0, 1)*ifelse(model$model.structure$between=="D", d*(d-1)/2, 1) ) # covariances
}


chooseModels <- function(models, kList = NULL, struct = NULL) {
   if (!is.null(struct)) {  
    indStr = unlist(lapply(struct, function(x) {grep(x, models$names)}))
   } else {
    indStr = seq_along(models$names)
   }
   if (!is.null(kList)) {  
    indList = unlist(lapply(paste("k=",kList," ",sep=""), function(x) {grep(x, models$names)}))
   } else {
    indList = seq_along(models$names)
    kList = models$kList
   }
   indLS = intersect(indList, indStr)
   models2 = models
   models2$kList = kList
   models2$loglikelihoods = models$loglikelihoods[indLS]
   models2$names = models$names[indLS]
   models2$params = models$params[indLS]
   models2$models = models$models[indLS]
   models2
}


chooseOptimal <- function(models, penalty=2) {
   values = sapply(models$models, getGIC, p=penalty)
   models$models[[which.max(values)[1]]] 
}


