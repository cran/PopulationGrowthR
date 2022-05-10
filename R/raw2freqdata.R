#' Extract Frequency and Specimen data from the raw data
#'
#' @param rawdata a dataframe containing species, year
#' @param species name of the column containing species names
#' @param year name of the column containing year
#'
#' @return Retirns a list of two dataframes
#' \itemize{
#' \item data - a dataframe conatining Species, Year, Freuqency and Specimens
#' \item yeardata - a dataframe containing Year and Specimens
#' }
#' @examples
#' cleandata = raw2freqdata(rawdata)
#' fdata = cleandata$data
#' yeardata = cleandata$yeardata
#' @export
raw2freqdata = function(rawdata, species='species', year='year'){
  df = dplyr::group_by(rawdata, species, year)
  df = dplyr::summarise(df,Frequency = dplyr::n(), .groups = 'keep')
  df2 = dplyr::group_by(df, year)
  data = dplyr::mutate(df2, Specimens = sum(df2$Frequency))
  yeardata = dplyr::summarise(df2, Specimens = sum(df2$Frequency))
  names(data) = c("Species", "Year", "Frequency", "Specimens")
  names(yeardata) = c("Year", "Specimens")
  list(data=as.data.frame(data), yeardata=as.data.frame(yeardata))
}


