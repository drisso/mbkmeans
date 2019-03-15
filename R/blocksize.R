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
#' @examples
#' data <- matrix(NA, nrow = 100, ncol=1000)
#' blocksize(data, ram=1e6)
#' @export
blocksize<-function(data, ram = get_ram()){
    if(is.na(ram)) {
        warning("Unable to determine available RAM. Setting it to 1GB.")
        ram <- 1e9
    }
    result<-min(floor(as.numeric(ram)/(2*8*nrow(data))), ncol(data))
    result
}
