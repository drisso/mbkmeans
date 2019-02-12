#' blocksize
#'
#' Return the maximum number of rows to use based on the amount of ram memory.
#'
#' @param data matrix-like object.
#' @param ram the max amount of ram (in bytes) to use.
#'
#' @importFrom benchmarkme get_ram
#'
#' @return  Numeric value of the maximum number of rows.
#'
#'@export
#'
blocksize<-function(data, ram = get_ram()){
  result<-min(floor(as.numeric(ram)/(2*8*nrow(data))), ncol(data))
  result
}
