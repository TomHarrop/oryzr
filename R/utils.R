#' dateF
#' 
#' date formatted for console messages
dateF <- function() {
    return(format(Sys.time(), "[%a %b %d %H:%M:%S %Y]"))
} 
