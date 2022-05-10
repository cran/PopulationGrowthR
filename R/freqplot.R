#' Frequency plot for a lagphase fit
#'
#' @param fit1,fit2,fit3,fit4  "lagphase" fit objects to plot
#' @param xlab Label for the $x$-axis
#' @param ylab Label for the $y$-axis
#' @param main Title of the plot
#' @param cols Colors to be used to draw the lines
#' @param ... (optional) parameters to pass to plot
#'
#' @return Produces a plot of observed and predicted frequencies for the species against year
#' @examples
#' Species = unique(fdata$Species) #List of all species
#' fit1 = lagfit(fdata, yeardata, species=Species[1])
#' freqplot(fit1$fit)
#' @export

freqplot <- function(fit1, fit2=NULL, fit3=NULL, fit4=NULL,
                     xlab="Year", ylab="Frequency", main=fit1$name, cols=2:5, ...)
{
  if(is.element("data",names(fit1)))
    data <- fit1$data
  else
  {
    data <- fit1
    fit1 <- NULL
  }

  plot(Frequency ~ Year, data=data, xlab=xlab, ylab=ylab, main=main,
       ylim=range(0,data$Frequency,na.rm=TRUE),...)
  j <- (data$Specimens > 0)
  if(!is.null(fit1))
    lines(data$Year[j],fitted(fit1),col=cols[1])
  if(!is.null(fit2))
    lines(data$Year[j],fitted(fit2),col=cols[2])
  if(!is.null(fit3))
    lines(data$Year[j],fitted(fit3),col=cols[3])
  if(!is.null(fit4))
    lines(data$Year[j],fitted(fit4),col=cols[4])
}

