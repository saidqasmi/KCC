# Define spline environment with n nodes
# Inputs : need a weighting function w to compute locally adaptative splines
base_spline = function(n,w=1,rho=1){

  # Base spline :
  x = 1:n	# knots
  rn = n+2;	# base spline dimension
  X =matrix(data = 0,nrow = rn,ncol = 1)

  X[1,1] = x[1]
  X[2,1] = mean(x[1:2])
  X[3:(rn-2),1] = x[2:(n-1)]
  X[rn-1,1] = mean(x[(n-1):n])
  X[rn,1] = x[n]

  Y0 =diag(x=1,nrow = rn,ncol = rn)
  Y1 = Y2 = Y3 = matrix(data=0,nrow=rn,ncol=rn)

  for (i in 1:rn){
    f <- splinefun(X,Y0[,i],method = "fmm")## Compute the spline function

    # Compute (f0), 1st derivative (f1), 2nd derivative (f2) et third (f3) at Xi
    Y1[,i] = f(X,deriv=1);
    Y2[,i] = f(X,deriv=2);
    Y3[,i] = f(X,deriv=3);

  }


  ## Checks on w
  if (length(w)==1) { w=rep(1,n) }
  nw = length(w)
  if (nw != n){ stop("Length of w() is not correct") }
  # Construct the final weights ww
  w_interv = apply(cbind(w[1:(n-1)],w[2:n]),1,mean)
  ww = w_interv[c(1,1:(n-1),n-1)]


  ## Matrix B :
  Zn = Y0[c(1, 3:(rn-2), rn),]
  ## Matrix G :
  Gn = G0=matrix(data = 0,nrow=dim(t(Zn)%*%Zn)[1],ncol = dim(t(Zn)%*%Zn)[2])
  interv = X[2:rn] - X[1:(rn-1)]
  interv =matrix(data = interv,nrow = length(interv),ncol = 1)
  Interv = interv%*% matrix(data = 1, nrow=1, ncol=length(Y2[1,]) )

  for(s1 in 1:rn){
    f2 = Y2[1:(dim(Y2)[1]-1),s1]%*% matrix(data = 1, nrow =1,ncol =length(Y2[1,])  )
    f3 = Y3[1:(dim(Y3)[1]-1),s1]%*% matrix(data = 1, nrow =1,ncol =length(Y2[1,])  )
    g2 = Y2[1:(dim(Y2)[1]-1),]
    g3 = Y3[1:(dim(Y3)[1]-1),]
    Int = f2 * g2 *Interv + (f2 * g3 + g2 * f3)*Interv^2/2 + f3 * g3 *Interv^3/3
    Gn[s1,] =apply(ww*Int, 2, sum)
  }

  Zn <<- Zn
  Gn <<- Gn

}



