#' Decomposition of simulated time series into natural and anthropogenic
#' responses
#'
#' \code{x_fit} fits a Generalized Additive Model (GAM, see equation 7 in Qasmi
#' and Ribes, 2021) to each climate model ensemble mean to estimate the natural
#' and anthropogenic responses. For each response, an ensemble of realisations
#' is returned accounting for uncertainty related to internal variability.
#'
#' @param Xd a 2-D array of dimension \code{[length(year),length(model)]}
#'     containing the time series associated with the model ensemble means. The
#'     first dimension has the time series of years, \code{year}, as names. The
#'     second dimension, \code{model}, must be a vector containing the names of
#'     all models.
#' @param Enat a 2-D array of dimension \code{[length(year),Nres]}
#'     containing the response to natural forcings.
#' @param Sigma the covariance matrix accounting for uncertainty in internal
#'     variability
#' @param x_df the number of degrees of freedom for the smoothing by splines
#' @param Nres the number of realisations in the gaussian sample of
#'     \code{Enat}
#' @param ant a logical value indicating whether the anthropogenic response
#'     should be returned
#'
#' @return a 4-D array of dimension
#'     \code{[length(year), Nres+1, length(forcing), length(model)]}, containing
#'     the response to several forcings. \code{forcing} is a character vector
#'     containing the types of external forcing: \code{nat} for natural,
#'     \code{all} for all forcings, \code{ant} for anthropogenic (if flagged at
#'     the function call). A best-estimate is provided with \code{Nres}
#'     realisations sampling uncertainty in each model.
#'
#' @importFrom stats lm
#' @importFrom stats rnorm
#' @importFrom gam gam
#' @importFrom gam s
#'
#' @examples
#'
#' @export
x_fit = function(Xd,Enat,Sigma,x_df,Nres,ant) {

  # Extract dimensions
  year_str = dimnames(Xd)$year
  year = as.numeric(year_str)
  ny = length(year_str)
  Models = dimnames(Xd)$model
  Nmod = length(Models)
  sample_str = c("be",paste0("nres",1:Nres))

  # Checks
  if (dim(Enat)[1]!=ny) {message("Error in x_fit_da.R: Xd and Enat have different ny"); return}
  if (dim(Enat)[2]!=Nres+1) {message("Error in x_fit_da.R: Enat have wrong Nres (ie dim[2])"); return}

  # Resize Sigma -- if needed
  # Sigma can be eihter one single covariance matrix (the same for all models), or a collection of matrices (one for each model).
  if (length(dim(Sigma))==2) {
    v_mod = array(1,dim=Nmod,dimnames=list(model=Models))
    Sigma = Sigma %o% v_mod
  }


  # Preliminary spline calculation
  #--------------------------------
  HatM = hatm(x_df,year)
  Cov_spline = 0*Sigma
  Cov_spline_12 = 0*Sigma
  for (mod in Models) {
    Cov_spline[,,mod] = t(HatM) %*% Sigma[,,mod] %*% HatM
    Cov_spline_12[,,mod] = matrix_sqrt(Cov_spline[,,mod])
  }
  # For extention: compute a Hat matrix (possibly with different methods), then apply to Xd...

  # Design X
  #----------
  year_str = dimnames(Xd)$year
  ny = length(year_str)
  Models = dimnames(Xd)$model
  Nmod = length(Models)
  sample_str = c("be",paste0("nres",1:Nres))

  X = array(NA,dim=c(ny,Nres+1,3,Nmod),#
            dimnames=list(year=year_str,#
                          sample=sample_str,#
                          forc=c("nat","ant","all"),#
                          model=Models))

  beta_nat = array(NA,dim=c(Nres+1,Nmod),#
                   dimnames=list(sample=sample_str,#
                                 model=Models))


  # Fit gam, then apply smoothing splines
  #---------------------------------------
  # x = 1 + Enat + f(t) + e
  # By default the intercept "1" is put together with Enat, so x_nat = 1 + Enat
  for (mod in Models){
    message(c("		",mod))
    x = Xd[,mod]
    # Best estimate
    gam_be = gam(x ~ s(year,df=x_df) + Enat[,1])		# Preliminary be
    #save(gam_be,file="toto.Rdata")
    beta_nat["be",mod] = gam_be$coefficients[3]					# beta_nat
    x_nat_brut = beta_nat["be",mod]*Enat[,1]						# Enat
    x_nat_removed_be = x - x_nat_brut								# x - Enat
    x_ant_uncentered_be = HatM %*% x_nat_removed_be				# mu + f(t)
    x_ant_be = x_ant_uncentered_be - x_ant_uncentered_be[1]	# f(t)
    X[,"be","ant",mod] = x_ant_be
    X[,"be","all",mod] = mean(x) + x_nat_brut - mean(x_nat_brut) + x_ant_be - mean(x_ant_be)	# = 1 + Enat + f(t)
    # Rough estimate of beta_nat variance (Gaussian distribution is assumed in resampling) -- OLS formula
    sd_beta_nat = sqrt(t(Enat[,1])%*%Sigma[,,mod]%*%Enat[,1] / (t(Enat[,1])%*%Enat[,1])^2)
    # Resampling
    for (i in 1:Nres){
      gam_re = gam(x ~ s(year,df=x_df) + Enat[,i+1])
      beta_nat[i+1,mod] = gam_re$coefficients[3]+sd_beta_nat*rnorm(1)	# beta_nat; Possibly sd_beta_nat should be calculated for each Enat_resampled...
      x_nat_brut_re = beta_nat[i+1,mod]*Enat[,i+1]					# Enat
      x_nat_removed_re = x - x_nat_brut_re										# x - Enat
      x_ant_uncentered_re = HatM %*% x_nat_removed_re + Cov_spline_12[,,mod] %*% rnorm(ny)	# mu + f(t)
      X[,i+1,"ant",mod] = x_ant_uncentered_re - x_ant_uncentered_re[1]				# f(t)
      X[,i+1,"all",mod] = mean(x) + x_nat_brut_re - mean(x_nat_brut_re) + X[,i+1,"ant",mod] - mean(X[,i+1,"ant",mod])	# = 1 + Enat + f(t)
    }
    X[,,"nat",mod] = X[,,"all",mod] - X[,,"ant",mod]	# X_nat = 1 + Enat
  }

  if (ant) {
    return(X)
  } else {
    return(X[,,c("all","nat"),])
  }

}


hatm = function(df,year) {
	base_spline(length(year),ones(year))
	HatM = proj_dl(df)
	return(HatM)
}

# Adjusting the degree of freedom by dichotomy
proj_dl = function(Deglib) {

  # Initialise
  rhomax=10^7
  rhomin=10^-2
  tol=10^-2

  Rho = c(rhomin,rhomax)
  while (Rho[2]/Rho[1]> (1+tol) ) {
    rho_new = mean(Rho)
    ## H, Gamma :
    Hn = t(Zn)%*%Zn + rho_new * Gn
    Gamma_n = Zn %*% solve(Hn) %*% t(Zn)
    dl_new = sum(diag(Gamma_n))

    if (dl_new>Deglib) {
  	 Rho[1] = rho_new
    } else {
  	 Rho[2] = rho_new
    }
  }

  #return(Rho[1])
  Hn = t(Zn)%*%Zn + rho_new * Gn
  Gamma_n = Zn %*% solve(Hn) %*% t(Zn)
  return(Gamma_n)

}

#' Remove trend from time series
#'
#' @param x a vector of time series
#'
#' @return a vector from which the linear trend in \code{x} is removed
#'
#' @importFrom stats lm
#'
#' @examples
#'
#' @export
pictl_detrend = function(x) {
  n_full = length(x)
  n_ok = sum(!is.na(x))
  year = 1:n_ok
  lmo = lm(x[1:n_ok]~year)
  y_res = c( lmo$residuals, rep(NA,n_full-n_ok) )
  return(y_res)
}

#' Replace each element of a numeric vector by 1
#'
#' @param x a numeric vector
#'
#' @return a vector of the same length as \code{x} but with each element equal
#'     to 1
#'
#' @examples
#'
#' @export
ones = function(x) {
  z = x-x+1
  z[is.na(z)] = 1
  return(z)
}

