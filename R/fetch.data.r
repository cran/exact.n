#' Download a target power library
#'
#' Function downloads one of 85 power libraries from chrislloyd.com.au.
#' The libraries are all four column matrices with roughly 20 million
#' rows and will be around 1.4Gb within R. Download should take roughly
#' a minute. The object will have a name of the form LIB.alpha.delta.Rdata.
#'
#' @aliases fetch.data
#' @param delta value of clinically relevant difference
#' @param alpha value of nominal size of test
#' @param prin If TRUE, error and/or warning messages are printed.
#' @return No return value. A library object will appear in the global environment.
#' @author Chris J. Lloyd
#' @importFrom utils download.file
#' @importFrom httr GET user_agent write_disk
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to size and power constraints. \doi{10.13140/RG.2.2.11828.94085}
#' @examples
#' #'
#' \dontrun{
#' fetch.data(alpha=0.05,delta=0.10)
#' # Object LIB.a050.d10 will appear with 83,289,960 elements
#' }
#'
#' @export fetch.data

fetch.data=function(alpha,delta, prin = FALSE){
  delta.vals=(5:20)/100
  if(!is.element(delta,delta.vals)){
    if (prin) cat("Allowable values of delta are 0.05, 0.06, ... , 0.19, 0.20","\n")
    stop()
  }
  alpha.vals=c(0.01,0.02,0.025,0.05,0.10)
  if(!is.element(alpha,alpha.vals)){
    if (prin) cat("Allowable values of alpha are 0.01, 0.02, 0.025, 0.05, 0.10","\n")
    stop()
  }
  #
  # Local function
  LIB.name=function(alpha,delta){
    digit1=floor(alpha*10)
    digit2=floor(alpha*100)-digit1*10
    digit3=floor(alpha*1000)-digit1*100-digit2*10
    char.alpha=paste(digit1,digit2,digit3,sep="")
    digit1=floor(delta*10)
    digit2=floor(delta*100)-digit1*10
    char.delta=paste(digit1,digit2,sep="")
    object.name=paste("LIB.a",char.alpha,".d",char.delta,sep="")
    return(object.name)
  }
  URL.name=paste("https://chrislloyd.com.au/wp-content/themes/chrislloydblog/software/",LIB.name(alpha,delta),".Rdata",sep="")
  temp_RData = tempfile(fileext = '.Rdata')
#
  if (prin) cat("This download may take a minute or two.....","\n")
  if (prin) cat("A large matrix starting with LIB should appear in your global environment.","\n")
  if (prin) cat("It will have more than 20 million rows so it is best","\n")
  if (prin) cat("not to bother trying to view it.","\n")
  if (prin) cat("\n")
#
  res_downl = httr::GET(URL.name,
                        httr::user_agent("Mozilla/5.0"),
                        httr::write_disk(temp_RData, overwrite = TRUE))
  load(file = temp_RData,envir = .GlobalEnv)
  rm(temp_RData,URL.name)
  NULL
}

