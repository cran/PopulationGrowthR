#' Produces plot of the fitted spline function after adjusting for number of Specimens
#'
#' @param fit a "lagphase" fit object to plot
#' @param ylim vector of size 2 - limits of the $y$-axis
#' @param xlab Label for the $x$-axis
#' @param ylab Label for the $y$-axis
#' @param main Title of the plot
#' @param ... (optional) parameters to pass to plot
#'
#' @return Produces a plot of the fit with confidence bands
#' @examples
#' Species = unique(fdata$Species) #List of all species
#' fit1 = lagfit(fdata, yeardata, species=Species[1])
#' growthplot(fit1$fit)
#' @export
growthplot <- function(fit, ylim=NULL, xlab="Year", ylab="Adjusted Frequency", main = fit$name, ...)
{
  fits <- predict(fit, se.fit=TRUE)

  #Specimens <- model.matrix(fit)[,"Specimens"]
  #adjfits <- fits$fit - coef(fit)["Specimens"]*(Specimens - mean(Specimens))
  ladjfits <- fits$fit - fit$offset + mean(fit$offset,na.rm=TRUE)
  adjfits <- exp(ladjfits)
  up <- exp(ladjfits + 2*fits$se.fit)
  lo <- exp(ladjfits - 2*fits$se.fit)
  if(is.null(ylim))
    ylim <- range(lo,pmin(up,3*adjfits))

  j <- (fit$data$Specimens > 0)

  plot(fit$data$Year[j],adjfits, type="n", xlab=xlab,ylab=ylab,
       main=main,...)
  polygon(c((fit$data$Year[j]),rev(fit$data$Year[j])),c(lo,rev(up)),border=FALSE,col="gray")
  lines(fit$data$Year[j],adjfits)
  if(!is.null(fit$knots))
    abline(v=fit$knots[1],col="gray")
  rug(fit$Year[fit$data$Frequency > 0])
}
