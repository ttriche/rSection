# strawman dSection port to R and/or C++ via Rcpp
# (novelty comes from using a hierarchical prior)
#
# Y is a PxN matrix of beta values from N samples
# p0 is a TxN matrix of prior class probabilities
# W0 is a precision parameter for p0 ~ Dirichlet
# proposal defines the size of the proposal step
# X is a column matrix of treatment indicators
# burnin is the amount of burn-in to perform (>0)
# samples is the number of samples to compute > 0
# posterior is boolean (whether to sample from it)
#
MCMC.comp <- function(Y, p0, W0, proposal, X, burnin, samples, posterior=T){#{{{

  require('impute')
  require('methylumiCSV')
  if( max(Y) <= 1 && min(Y) >= 0 ) Y <- beta.probit(Y) # a hack
  if( anyMissing(Y) ) Y <- impute.knn(Y)$data # another hack (lm hates NAs)

  # this can and probably should be a design matrix, derp
  lvl <- levels(factor(as.factor(X))) # treatment levels not a matrix... yet
  unkl <- c(lvl, 'unknown') # we'll come back to this 
  nT <- dim(p0)[1]
  Tvec <- 1:nT
  UNKvec <- 1:(nT+1)
  i <- dim(Y)[1]
  j <- dim(Y)[2]

  # OLS solutions for each treatment level -- this is butt ugly
  xLS <- lapply(lvl, function(e) lm.fit(p0[which(X==e),],Y[which(X==e),])$coeff)
  names(xLS) <- lvl

  # residual variance, to estimate the precision parameter of the Dirichlet
  V <- rowMeans(data.matrix(as.data.frame(lapply(lvl, function(e) {
    rowVars(t(Y[which(X==e),])-t(xLS[[e]])%*%t(p0[which(X==e),]), na.rm=T)
  }))), na.rm=T)

  # estimate precision
  lambdaLS = V**(-1)
  require('MASS')
  params <- fitdistr(lambda_LS, 'gamma')[[1]] # could just get the MLE...
  alpha0 <- params[1]
  beta0 <- 1/params[2]
  maxiter <- (burnin+samples+1)
  MCData <- list(list(x=xLS, lambda=lambdaLS, p=beta.scale(p0, n=0.999)))
  nuLS <- alpha0/beta0

  # Hui, this is all you ;-)
  stop("This port of DSection to R is so very not finished. Please do not use!")
    
  # sample from chain
  for (idx in 2:maxiter) {

    cat(idx-1, '/', maxiter, "\n")
    if(idx > burnin + 1) { 
      iter <- idx-nBurnIn
      MCData[[iter]] <- MCData[[(iter-1)]]
    } else iter = 1 # still burning in
    MC <- MCData[[iter]]
    
    if(posterior) {
      for(m in 1:j) { # for each probe
        MC$p[,m] <- samplepvec( Y[,m],  # wtf is this??!?
                                MC$lambda, 
                                MC$x[,,Treatment(j)], 
                                MC$p(,m),
                                p0[,j],
                                W0,
                                W_proposal,
                                T )
      }
    }

    for( k in 1:i ) {
      DJ = matrix(0, 1, j)
      for( e in 1:E ) {
        DJ[ which(X==e) ] <- t(MC$x[,k,e]) %*% MC$p[,which(X==e)]
      }
            
      MC$lambda[k] = gamrnd(alpha0+J/2,1/(beta0+1/2*(Y(k,)-DJ)*t(Y(k,)-DJ)))
        
      for(e in lvl) {
        for(t in 1:nT) {
          #Helper data for sampling cell type expression
          Z <- t(MC$p[Tvec == t,Treatment == e])*MC$x[Tvec = t,k,e]
          P <- MC$lambda[k]*(Y[k,which(X==e)]*t(MC$p[t,which(X==e)])-MC$p[t,which(X==e)]*Z)+xLS[t,k,e]*nuLS
          Q <- MC$lambda[k]*sum(MC$p[t,which(X==e)]**2)+nuLS
                
          #Sample cell type expression
          MC$x[t,k,e] <- (P/Q)+(1/sqrt(Q))*runif(1,0,1)

        }
      }   
    }
    
    #Erase the tail of the print-out
    b = 1 + numel(num2str(idx-1)) + numel(num2str(nBurnIn+nSamples))
    str = repmat('\b',c(1,b))
    fprintf(str)
  }

  #Remove the first sample from the struct, it's not part of the chain
  MCData = MCData(2:end)

  #Stop timer
  fprintf(' done!\n')
  t = toc
  fprintf('Time elapsed: %fs\n',t)

  #Function for sampling mixing proportions using M-H
  p_new <- function(y, lambda, x, p_old, p0, W1, W2, T) {
    samplepvec(Y,lambda,x,p_old,p0,W1,W2,T)
  }

  #Propose a new p vector from the current Dirichlet distribution
  p_prop = zeros(T,1)
  for(t in 1:T)
      p_prop(t) = gamrnd(W2*p_old(t),1)
  p_prop = p_prop/sum(p_prop)

  #Push the proposal away from the boundaries, if it touches them
  p_prop(p_prop > 1-eps) = 1-eps
  p_prop(p_prop < eps) = eps

  #A fraction product of Gamma functions
  gamprod_frac = 1
  for( t in 1:T )
     gamprod_frac = gamprod_frac*gamma(W2*p_prop(t))/gamma(W2*p_old(t))

  #Calculate the jump factor
  alpha = exp(-1/2*(t(lambda)*((t(x)*p_prop)^2-2*Y*(t(x)*p_prop)))+t(W2*p_old-1)*log(p_prop)+t(W1*p0-1)*log(p_prop)+1/2*(t(lambda)*((t(x)*p_old)^2-2*Y*(t(x)*p_old)))-t(W2*p_prop-1)*log(p_old)-t(W1*p0-1)*log(p_old))*gamprod_frac

  #Make a decision whether to jump to a new state or to stay at the old one
  if(runif(1,0,1) < min(1,alpha))
    p_new = p_prop
  else
    p_new = p_old

}
