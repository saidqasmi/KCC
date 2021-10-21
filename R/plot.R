#' Illustrate the observational constraint
#'
#' @param X_Cons a 3-D array of dimension
#'     \code{[length(year),Nres+1,c("uncons","cons")]} containing time series
#'     sampling the prior and posterior distributions. The first dimension has
#'     the time series of years, \code{year}, as names. The second dimension
#'     corresponds to the number of realisations sampling the distribution ; the
#'     first realisation is the best-estimate. For the third dimension,
#'     \code{"uncons"} (\code{"cons"}) refers to the unconstrained (constrained)
#'     distribution.
#' @param obs_x a vector of observations.
#' @param ofile a character giving the name of the file in which graphics are
#'     saved
#' @param ref_plot a vector containing the years corresponding to the reference
#'     period to plot anomalies. If \code{NULL}, raw values are plotted.
#' @param event a numeric value indicating a year corresponding to an event
#'     plotted as a vertical line
#' @param ylim a vector of two numeric values indicating the limits along the
#'     y-axis. If \code{NULL} a range is computed automatically.
#'
#' @return A detrended time series
#'
#' @examples
#'
#' @export
plot_scx = function(X_Cons,obs_x,ofile,ref_plot=NULL,event=NA, ylim=NULL){
	pdf(ofile)
		plot_cx(X_Cons[,,"uncons"],#
				  X_Cons[,,"cons"],#
				  obs_x,ofile,ref_plot, event, ylim)
	dev.off()
}


plot_cx = function(X,X_Cons,obs_x,ofile,ref_plot=NULL,event=NA,ylim=NULL){
	year = as.numeric(dimnames(X)$year)
	year_obs_x = as.numeric(names(obs_x))

	x = X
	x_cons = X_Cons

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

	pdf(ofile)

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
	plot(year_obs_x, obs_x, xlim=range(year), ylim=ylim, type="p", pch=16, cex=.8, xlab="Year", ylab="Temperature (K)", panel.first=abline(v=event_year,col="gray"))
	yaxp = par("yaxp")
	yticks = seq( yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3] )
	abline(h=yticks,lty=3)
	polygon(c(year,rev(year)), c(x_q95,rev(x_q05)),border=NA,col=rgb(1,0,0,alpha=.2))
	polygon(c(year,rev(year)), c(xc_q95,rev(xc_q05)),border=NA,col=rgb(1,0,0,alpha=.5))
	colb = col2rgb("brown")/255
	lines(year,x[,1]		,lwd=1.5,col=rgb(colb[1],colb[2],colb[3],alpha=.5))
	lines(year,x_cons[,1],lwd=2,col=rgb(colb[1],colb[2],colb[3],alpha=1))

	dev.off()

}
