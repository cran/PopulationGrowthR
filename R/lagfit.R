#' Fits a piecewise glm model with lags
#'
#' This function fits a piecewise poisson model to the frequency data of different Species. It assumes that the data
#' contains columns Year, Frequency and Specimens.
#'
#' @param data a dataframe containing the columns Species (optional), Year, Frequency and Specimens.
#' @param yeardata a dataframe containing the columns Year and Specimens giving the total number of Specimens for each Year.
#' @param species list of species for which the model is to be fitted. Default is NULL, which fits the model for all species in the data.
#' @param knots a list of knots to be used for the piecewise model. Default is NULL, which chooses the optimal model with 0-4 knots.
#' @param zeros logical. Specifies whether missing year for the species will be filled with zeros. Default is TRUE.
#' @param plotlag logical. If TRUE a plot of the fitted model will be produced for each species.
#' @param plotfreq logical. If TRUE frquency plots will be created for each species.
#'
#' @return If the model is fit for a single species following are returned as a list
#' \itemize{
#' \item Species - Species name
#' \item Scene - Different scenario of the fit between the knots. A sequence of 0, + or - is returned. A 0 indicates constant, + indicates increasing and a - indicates decreasing.
#' \item Lag - Logical. Is there a lag present or not.
#' \item Laglength - Length of the first lag. Position of the First Knot - the first year for that species
#' \item FirstYear - The first year for that species for which data is available.
#' \item EndYear - The first knot position.
#' \item fit - the fitted model.
#' }
#'
#' @return If the number of species is more than one, then a list is returned with following items:
#' \itemize{
#' \item  fitdata - dataframe is returned with the items in the above list except for the fitted model.
#' \item  fitcoefs - list of coefficients for the piecewise fits for each Species
#'}
#'@examples
#'#Run lagfit for 1 species only
#'Species = unique(fdata$Species) #List of all species
#'
#'fit1 = lagfit(fdata, yeardata, species=Species[1])

#'#Run lagfit for multiple species
#'fit2 = lagfit(fdata, yeardata, species=Species[1:3])
#'fitdata = fit2$fitdata  #Dataframe containing fits
#'fitcoefs = fit2$fitcoefs #List containing slopes of the fitted splines
#'
#'\dontrun{
#'#Run lagfit for the whole dataset
#'fitall = lagfit(fdata, yeardata)
#'}
#' @importFrom graphics abline lines polygon rug
#' @importFrom stats fitted glm na.omit optim poisson predict quantile vcov
#' @export
lagfit = function(data, yeardata, species = NULL, knots=NULL, zeros=TRUE, plotlag=FALSE, plotfreq=FALSE){


  # plotlag - TRUE if you want plot.lagphase for each species
  # plotfreq - TRUE if you want frequency plot for each species

  outspecies<-NULL
  outscene <- NULL
  outslopes <- NULL
  outlag <- NULL
  outlaglength <- NULL
  outfirstyear <- NULL
  outendyear <- NULL

  if(is.null(species)){
    Species = unique(as.character(data$Species))
  }
  else Species = species

  # Run the main loop
  fitlist=list()

  for( i in 1: length(Species)){

    sdata = get.species(data,yeardata, species = Species[i], zeros = zeros)
    #print(sdata$Species)
    #sdata$Island=island
    if(length(sdata$Year) >= 2) { #Insufficient data check to fit models

      #outisland = c(outisland,island)

      outspecies = c(outspecies, Species[i])
      fit0 = lagphase(sdata, knots, zeros=zeros, order=1)
      # Check for the scenario

      if(!fit0$lagphase){
        endyear = NA
        if(fit0$scene=="constant"){
          scene = '0'
          slopes ='0'
        }
        if(fit0$scene=="linear"){
          if(fit0$coef[2] >0) scene = '+'
          if(fit0$coef[2] <0) scene = '-'
          slopes = as.character(round(fit0$coef[2],3))
        }
        coefl = fit0$coef
      }
      else{
        nknots = length(fit0$knots)
        vcoef = vcov(fit0)
        scene = NULL
        slopes = NULL
        coefl = fit0$coef[1]
        for(iknot in 1:(nknots+1)){
          beta = sum(fit0$coef[2:(iknot+1)])
          coefl = c(coefl,beta)

          varbeta = 0.0
          for(j in 2:(iknot+1)){
            varbeta = varbeta + sum(vcoef[j, 2:(iknot+1)])
          }
          sebeta = sqrt(varbeta)
          if(is.na(beta)) {
            tstat=0.0
          }
           else tstat = beta/sebeta

          if(tstat > 2.0) {
            scene = paste0(scene,'+')
            slopes = paste(slopes, round(beta,3), ';')
          }
          else {
            if(tstat < -2.0) {
              scene = paste0(scene,'-')
              slopes = paste(slopes, round(beta,3), ';')
            }
            else {
              scene = paste0(scene,'0')
              slopes = paste(slopes, '0', ';')
            }
          }
        }

      }
      firstscene = substr(scene,1,1)
      firstyear = fit0$Year[1]

      if(firstscene != '0'){
        fit0$lagphase = FALSE
        fit0$lengthlag = NA
      }
      if(fit0$lagphase){
        endyear = fit0$knots[1]
      }
      else endyear = NA
      outlag = c(outlag, fit0$lagphase)
      outlaglength = c(outlaglength,fit0$lengthlag)
      outfirstyear = c(outfirstyear, firstyear)
      outendyear = c(outendyear, endyear)

      outscene = c(outscene, scene)
      outslopes = c(outslopes, slopes)

      if(plotlag) growthplot(fit0)
      if(plotfreq) freqplot(fit0)
      sname= Species[i]
      coeflist = list(coefl)
      names(coeflist)=sname
      fitlist = append(fitlist, coeflist)
    }
  }

  if (length(Species) > 1){
    fitdata = data.frame(Species = outspecies, Scene = outscene, Lag =outlag, Laglength = outlaglength, FirstYear=outfirstyear, EndYear=outendyear, Slopes = outslopes)
    out = list(fitdata=fitdata, fitcoefs = fitlist)
  }
  else{
    if(length(sdata$Year) >= 2) {
      out = list(Species = outspecies, Scene = outscene, Lag =outlag, Laglength = outlaglength, FirstYear=outfirstyear, EndYear=outendyear)
      out$fit = fit0
    }
    else {
      print('There is insufficient data')
      return()
    }
  }
  out
}

# Extract Frequency data for given island and Species from data
# If zeros=TRUE, include zeros in returned data

get.species <- function(x, y, species, zeros=TRUE)
{

  if ("Species" %in% colnames(x)) out <- subset(x, x$Species==species)
  else out = x

  out <- out[,c("Year","Frequency", "Specimens")]

  # Sort the data by Year
  yorder = order(out$Year)
  out = out[yorder,]

  indx = which(diff(out$Year)==0)
  if(length(indx)>0){
    out$Frequency[indx]=out$Frequency[indx]+out$Frequency[indx+1]
    out = out[-(indx+1),]
  }

  if(zeros)
  {
    yrs <- min(out$Year):max(out$Year)
    zeros <- as.data.frame(matrix(0,nrow=length(yrs),ncol=3))
    colnames(zeros) <- colnames(out)
    zeros[,"Year"] <- yrs
    j <- is.element(zeros[,"Year"],out[,"Year"])
    zeros[j,"Frequency"] <- out[,"Frequency"]
    j1 <- is.element(y[,"Year"],zeros[,"Year"])
    j2 <- is.element(zeros[,"Year"], y[,"Year"])
    zeros[j2,"Specimens"] <- y[j1,"Specimens"]
    out <- zeros
  }
  # Either way, Frequency is missing if Specimens=0
  #out$Frequency[out$Specimens==0] <- NA
  out <- as.list(out)
  out$Species <- species
  #out$Island <- island
  return(out)
}

# Main function. Give it a set of data where the columns include
#  Year
#  Frequency
#  Specimens
# It will find appropriate knots if not specified
# It will choose an appropriate order if not specified
# Just don't give it knots but no order
# If gam=TRUE, it will return a gam model instead.

lagphase <- function(data, knots=NULL, order=1, gam=FALSE,zeros=TRUE)
{
  # Set zeros to missing
  if(!zeros)
    data$Frequency[data$Frequency==0] <- NA
  else
    data$Frequency[data$Specimens==0] <- NA

    # Fit gam
  if(gam)
  {
    gamfit <- gam(Frequency ~ s(Year), offset=log(data$Specimens), data=data, family=poisson)
    gamfit$Year <- data$Year
    gamfit$name <- data$Species
    return(gamfit)
  }

  # Otherwise fit a glm
  # Check if knots==0 meaning no knots to be included
  if(!is.null(knots))
  {
    if(length(knots)==1)
    {
      if(knots==0) # i.e., no knots to be included
      {
        # Fit model with no knots
        suppressWarnings(fit <- glm(Frequency ~ 1, offset=log(data$specimens), data=data, family=poisson, na.action=na.omit))
        fit$Year <- data$Year

        fit$name <- data$Species

        fit$data <- data
        fit$lengthlag <- NA
        fit$lagphase <- FALSE
        class(fit) <- c("lagphase","glm","lm")
        return(fit)
      }
    }
  }

  # Otherwise proceed
  # Choose order if not provided
  if(is.null(order))
  {
    if(!is.null(knots))
      stop("Not implemented. If you specify the knots, you need to specify the order.")
    fit1 <- lagphase(data, order=1)
    fit3 <- lagphase(data, order=3)
    bestfit <- fit1
    if(AICc(fit3) < AICc(bestfit))
      bestfit <- fit3
    return(bestfit)
  }
  # Otherwise proceed with specified order
  if(!is.null(knots))
  {
    return(lagphase.knots(knots, data, order))
  }
  # Otherwise order specified but knots unspecified

  # Choose some initial knots
  knots <- quantile(data$Year,prob=c(0.2,0.4,0.6,0.8))
  names(knots) <- NULL

  # Fit best 4, 3, 2, 1 and 0 knot models

  fit4 <- optim(knots, tryknots, data=data, order=order)
  fit3 <- optim(knots[2:4], tryknots, data=data, order=order)
  fit2 <- optim(knots[c(2,4)], tryknots, data=data, order=order)
  fit1 <- optim(knots[2], tryknots, data=data, order=order, method="Brent",
                lower=min(data$Year), upper=max(data$Year))
  suppressWarnings(fit0 <- glm(Frequency ~ 1, offset=log(data$Specimens), family=poisson, data=data, na.action=na.omit))
  fitl <- glm(Frequency ~ Year, offset=log(data$Specimens), family=poisson, data=data, na.action=na.omit)

  # Find best of these models:
  bestfit <- fit4
  if(fit3$value < bestfit$value)
    bestfit <- fit3
  if(fit2$value < bestfit$value)
    bestfit <- fit2
  if(fit1$value < bestfit$value)
    bestfit <- fit1
  if(AICc(fit0) < bestfit$value)
  {
    bestfit <- fit0
    bestfit$Year <- data$Year
    bestfit$name <- data$Species

    bestfit$data <- data
    bestfit$scene <- "constant"
  }
  else if(AICc(fitl) < bestfit$value)
  {
    bestfit <- fitl
    bestfit$Year <- data$Year
    bestfit$name <- data$Species
    bestfit$data <- data
    bestfit$scene <- "linear"

  }
  else { # Refit best model
    bestfit$par = round(bestfit$par)
    bestfit <- lagphase.knots(bestfit$par, data=data, order=order)
  }
  if(is.null(bestfit$knots))
  {
    bestfit$lagphase <- FALSE
    bestfit$lengthlag <- NA
  }


  return(bestfit)
}


# Fit model with lag phase followed by growth
# where knots and order are specified
lagphase.knots <- function(knots, data, order=1)
{
  x <- matrix(NA,ncol=length(knots)+1,nrow=length(data$Year))
  x[,1] <- data$Year
  for(i in 1:length(knots))
    x[,i+1] <- pmax((data$Year-knots[i])^order,0)

  suppressWarnings(fit <- glm(Frequency ~ x, offset=log(data$Specimens), family=poisson, data=data, na.action=na.omit))
  fit$knots <- knots
  names(fit$knots) <- paste("K",1:length(knots),sep="")
  fit$Year <- data$Year
  fit$order <- order
  fit$name <- data$Species
  fit$data <- data

  # Check if there is a lag phase and record it
  fit$lengthlag <- NA
  if(length(fit$knots) > 0)
  {
    if(fit$coef[2] != 0)
      fit$lengthlag <- fit$knots[1] - min(data$Year)
  }
  fit$lagphase <- !is.na(fit$lengthlag)

  class(fit) <- c("lagphase","glm","lm")
  return(fit)
}

# Check that the specified knots make sense.
# Then use lagphase.knots to fit the model
# Returns AIC of fitted model
tryknots <- function(knots, data, order)
{
  # Knots must be interior to the data
  knots = round(knots)
  if(min(knots) < min(data$Year))
    return(1e50)
  if(max(knots) > max(data$Year))
    return(1e50)
  # Knots must be five Years apart and ordered
  if(length(knots) > 1)
  {
    dk <- diff(knots)
    if(min(diff(knots)) < 5)
      return(1e50)
  }

  # Number of points in between each knots must at least 2
  if(length(knots) > 0)
  {

    npoint = tapply(data$Frequency, cut(data$Year, breaks=unique(c(min(data$Year), knots, max(data$Year)))), function(x) sum(x!=0))
    if(all(is.na(npoint))) return(1e50)
    if (min(npoint, na.rm=TRUE) < 2) return(1e50)
  }
  # OK. Now fit the model
  fit <- lagphase.knots(knots, data, order)
  # Return the AICc of the fitted model
  return(AICc(fit))
}

# Function to return corrected AIC from a fitted object
AICc <- function(object)
{
  aic <- object$aic
  k <- length(object$coefficients)
  n <- object$df.residual+k
  aicc <- aic + 2*k*(k+1)/(n-k-1)
  if(aicc == Inf) aicc=1e50
  return(aicc)
}


