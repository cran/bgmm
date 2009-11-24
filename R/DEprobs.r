
DEprobs <- function(model, verbose=FALSE){
  k=model$k
  if (!k%in%c(2,3))
    stop("DEprobs: Wrong number of fitted model components. Only support for 2 or 3-component models is provided")

  ## Get the right component to correspond to the differentially expressed genes -- which looks like the differential one according to the posterior probabilities
  ## For models of size 2 -- should be the one with a broader range (encompassing the other, or the one with higher maximum deviation from 0 (we assume the data are centered around 0 and are both positive and negative)

  d=model$X
  if (!is.null(model$knowns))
    d=rbind(d, model$knowns)
  preds=predict(model, d)

  if (k==2)
    resp=check2res(d, preds, verbose)
  else
    resp=check3res(d, model, verbose)

  return(resp)
}


check2res <- function(d, preds,verbose){
  class=preds$class
  d1 <- d[class==1]
  d2 <- d[class==2]
  r1 <- range(d1)
  r2 <- range(d2)
  m1=max(abs(r1))
  m2=max(abs(r2))
  diff=2
  if (r1[1]<r2[1]){ ##1 sticks out on the left
    if(r1[2]>r2[2]){ ##1 sticks out on the right
      diff=1 
    }else{
      if (m1>m2){ ##1 sticks out more from 0 than the 2
        diff=1
        warning(paste("DEprobs: accoding to the data, the differentially expressed component cannot be identified as the one which is broader than the unchanged on both sides. Choosing the component with a larger dev from 0"))
      }
    }
  }## diff either encompasses or has higher dev from 0
  
  p <- preds$tij[,diff]
  downs <- (d<0)
  p[downs] <- p[downs]*(-1)
  if (verbose)
    print(paste("The differential component number:",diff))
  return(list(diff.p=p, diff.c=c(diff)))
}


check3res <- function(d, model, verbose){
  ##order by means
  ord=order(model$mu)
  preds=predict(model, d)
  tij=preds$tij[,ord]

  if (verbose){
    print(paste("the downregulated component number:", ord[1]))
    print(paste("the upregulated component number:", ord[3]))
  }

  p <- -1*tij[,1] ##probability of down-regulation
  ups <- (d>0) 
  p[ups] <- tij[ups,3]

  return(list(diff.p=p, diff.c=c(ord[1],ord[3])))
}
