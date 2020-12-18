#' Same as plot_cx, but inputs (X, X_Cons) include several scenarios, and are grouped into the same (big) array
#'
#' More detailed description
#'
#' @param Xd A piControl time series
#'
#' @return A detrended time series
#'
#' @examples
#'
#' @export
plot_scx = function(X_Cons,obs_x,ofile,ref_plot=NULL,sample=0,event=NA, name_multi="multi",ylim=NULL){
	pdf(ofile)
		plot_cx(X_Cons[,,"uncons"],#
				  X_Cons[,,"cons"],#
				  obs_x,ofile,ref_plot,sample, event, name_multi, ylim, add=T)
	dev.off()
}


#' Illustrate the observational constraint on covariate X
#'
#' More detailed description
#'
#' @param Xd A piControl time series
#'
#' @return A detrended time series
#'
#' @examples
#'
#' @export
plot_cx = function(X,X_Cons,obs_x,ofile,ref_plot=NULL,sample=0,event=NA,name_multi="multi",ylim=NULL,add=F){
	year = as.numeric(dimnames(X)$year)
	year_obs_x = as.numeric(names(obs_x))

	if ("model" %in% nnames(X)) {
		x			= X[,,name_multi]
		x_cons	= X_Cons[,,name_multi]
	} else {
		x			= X
		x_cons	= X_Cons
	}

	if (!is.null(ref_plot)) {
		# Check that ref_plot is included in year and year_obs_x
		if (prod(ref_plot %in% year_obs_x) & prod(ref_plot %in% year)) {
			x = x - ones(x[,1]) %o% apply(x[year %in% ref_plot,], 2, mean)
			x_cons = x_cons - ones(x_cons[,1]) %o% apply(x_cons[year %in% ref_plot,], 2, mean)
			obs_x = obs_x	- mean(obs_x[year_obs_x %in% ref_plot])
		}
	}
	
	if (is.na(event)) {
		event_year = NA
	} else {
		event_year = event$year
	}

	if (!add) {
		pdf(ofile)
	}
	
	par(font.lab=2,font.axis=2,cex.lab=1.2,mar=c(4,4,1,1),mgp=c(2.5,.7,0))
	x_q95 = apply(x[,-1],1,quantile,.95)
	x_q05 = apply(x[,-1],1,quantile,.05)
	xc_q95 = apply(x_cons[,-1],1,quantile,.95)
	xc_q05 = apply(x_cons[,-1],1,quantile,.05)
	if (is.null(ylim)) {
		ylim=range(obs_x,x_q05,x_q95,xc_q05,xc_q95)
	} else {
	}
	
	#plot(year_obs_x, obs_x, xlim=range(year), ylim=ylim, type="p", pch=16, cex=.8, xlab=TeX("\\textbf{Year}"), ylab=TeX("\\textbf{Temperature $$ $$ ($^o C$)}"), panel.first=abline(v=event_year,col="gray"))
	plot(year_obs_x, obs_x, xlim=range(year), ylim=ylim, type="p", pch=16, cex=.8, xlab="Year", ylab="Temperature (°C)", panel.first=abline(v=event_year,col="gray"))
	yaxp = par("yaxp")
	yticks = seq( yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3] )
	abline(h=yticks,lty=3)
	polygon(c(year,rev(year)), c(x_q95,rev(x_q05)),border=NA,col=rgb(1,0,0,alpha=.2))
	polygon(c(year,rev(year)), c(xc_q95,rev(xc_q05)),border=NA,col=rgb(1,0,0,alpha=.5))
	colb = col2rgb("brown")/255
	lines(year,x[,1]		,lwd=1.5,col=rgb(colb[1],colb[2],colb[3],alpha=.5))
	lines(year,x_cons[,1],lwd=2,col=rgb(colb[1],colb[2],colb[3],alpha=1))
	
	if (sample>0) {
		X_sample = matrix(0,ny,sample)
		X_Cons_sample = matrix(0,ny,sample)
		for (i in 1:sample){ 
			lines(year, X_sample[,i])
			X_Cons_sample[,i] = A_full %*% (mu_post + Sigma_post_sqrt %*% rnorm(mu_post))
			lines(year, X_Cons_sample[,i], col="red")
		}
	}
	
	if (!add) {
		dev.off()
	}
}

