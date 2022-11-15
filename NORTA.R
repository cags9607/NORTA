
library('pracma')    ## For numerical integration and root finding
library('Matrix')    ## For matrix operations
library("MASS")     ## To simulate a MVN random vector

# Let's code a family of functions to obtain the mean and variance of the internal R distributions

meannorm = function(mean, sd) mean
varnorm = function(mean, sd) sd^2
meanpois = function(lambda) lambda
varpois = function(lambda) lambda
meanbinom = function(size, prob) size*prob
varbinom = function(size, prob) size*prob*(1 - prob)
meanbeta = function(shape1, shape2, ncp = NULL) {
  if (is.null(ncp)) {
    return (shape1/(shape1 + shape2))
  }
  else {
    return(mean(rbeta(1e6,shape1,shape2,ncp)))
    
  }
}
varbeta = function(shape1, shape2, ncp = NULL) {
  if (is.null(ncp)) {
    return ((shape1*shape2) / ((shape1 + shape2)^2)*(shape1 +shape2 + 1))
  }
  else {
    return(var(rbeta(1e6,shape1,shape2,ncp)))
    
  }
}
meanchisq = function(df, ncp) df + ncp
varchisq = function(df, ncp) 2*(df + 2*ncp)
meanexp = function(rate) 1/rate
varexp = function(rate) 1/rate^2
meanf = function(df1, df2, ncp) (ncp + df1)*df2/(df1*(df2 - 2))
varf = function(df1, df2, ncp) df2^2 * (ncp^2 + (2*df1 + 4)*ncp + df1*(df1 + 2))/(df1^2*(df2-4)(df2-2)^2)
meangamma = function(shape, scale) shape*scale
vargamma = function(shape, scale) shape*scale^2
meangeom = function(p) 1/p
vargeom = function(p) (1-p)/p^2
meanhyper = function(m, n, k) k*m/(m+n)
varhyper = function(m, n, k) (k*m/(m+n)) * (1 - m/(m+n))
meanlnorm = function(meanlog, sdlog) exp(meanlog + sdlog^2/2)
varlnorm = function(meanlog, sdlog) exp(2*meanlog + sdlog^2)*(exp(sdlog^2)-1)
meannbinom = function(size, prob) size*(1-p)/p
varnbinom = function(size, prob) size*(1-p)/p^2
meant = function(df, ncp) ncp*(1 - 3/(4*df -1))^-1
vart = function(df, ncp) df*(1 + ncp^2)/(df-2) - ncp^2*(1 - 3/(4*df -1))^-2
meanunif = function(min, max) (min+max)/2
varunif = function(min, max) (max - min)^2/12
meanweibull = function(shape, scale) scale*gamma(1 + 1/shape)
varweibull = function(shape, scale) scale^2*(gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)

NORTA = function(N, marginals, Sigma_X, parms, warning = T, message = T, xmin = -8, xmax = 8, ymin = -8, ymax = 8, 
                  tol = 1e-12){
  n = length(marginals)
  # A list containing the quantile functions of the marginals
  quantiles = sapply(1:n, function(k) function(x) (do.call(paste("q",marginals[k],sep=""),append(x,parms[[k]]))))

  means = sapply(1:n, function(k) do.call(paste("mean",marginals[k],sep=""),parms[[k]]))
  variances = sapply(1:n, function(k) do.call(paste("var",marginals[k],sep=""),parms[[k]]))

  
  
  Sigma_Z <- matrix(nrow=n,ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      if (i > j) { 
        
        ## These are part of the expected value integrand
        
        q_i <- function(z_i) quantiles[[i]](pnorm(z_i)) 
        q_j <- function(z_j) quantiles[[j]](pnorm(z_j))
        
        ## The pdf of the binomial standard with correlation pij
        
        f_ij <- function(z_i,z_j,pij) (1/(2*pi*sqrt(1-pij^2))) * exp(-(2*(1-pij^2))^-1 *(z_i^2-2*pij*z_i*z_j +z_j^2))
        
        ## The expected value
        
        E_ij <- function(z_i,z_j,pij) function(z_i, z_j) q_i(z_i)*q_j(z_j)*f_ij(z_i,z_j,pij)
        
        ## The correlation function
        
        c_ij <- function(pZ_ij,parms) (integral2(fun = E_ij(z_i,z_j,pij=pZ_ij[1]),
                                             xmin = xmin,xmax = xmax, ymin= ymin, ymax = ymax,
                                             vectorized = F, reltol = 1e-10)$Q - means[[i]]*means[[j]])/sqrt(variances[[i]]*variances[[j]])  - parms[1]
        
        # Calculate the bounds of the feasible interval
        
        pij_max = round(c_ij(0.99,0),5)
        
        pij_min = round(c_ij(-0.99,0),5)
        
        if (message) cat(paste("The feasible interval for marginal",j," and marginal", i, " is: \n [", pij_min,",", pij_max,"]\n"))
        
        # Find the root using Newton's method
        
        pZ_ij = newtonRaphson(function(pij) c_ij(pij,Sigma_X[i,j]),x0 = Sigma_X[i,j], tol = tol)
        Sigma_Z[i,j] = pZ_ij$root
  
      }
    }
  }
  
  diag(Sigma_Z) <- 1 # set the diagonal to ones
  Sigma_Z <- forceSymmetric(Sigma_Z,uplo="L") # as the correlation matrix is a sym. matrix.
  Sigma_Z = matrix(Sigma_Z,n,n)
  
    if (!isposdef(Sigma_Z)){ # to evaluate positive definiteness 
    Sigma_Z <- nearPD(Sigma_Z)
    if (isTRUE(warning)) print('Warning!!! Non-positive definite matrix found. A correction was made to make it Positive Definite')
  }
  
  M <- chol(Sigma_Z) # Cholesky Decomposition
  #W <- matrix(rnorm(n*N),nrow=N,ncol=n) # A matrix of standard normal values
  #Z <- W %*% t(M)          # To apply step 4 of CarioNelson pseudocode
  
  # Generate a multivariate random normal vector with correlation matrix SIgma_X
  
  Z = mvrnorm(N, rep(0,n),Sigma_Z)
  
  output <- matrix(nrow = nrow(Z),ncol = ncol(Z))
  
  # Apply the inverse transformation to the uniform distributed vector U = pnorm(Z)
  
  for (i in 1:ncol(Z)) {
    output[,i] <- sapply((pnorm(Z[,i])),quantiles[[i]])
  }
  
  return (list(data = output, Sigma_Z =Sigma_Z, M = M))
  
}
