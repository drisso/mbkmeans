#' @export
#' @rdname MbkmeansParam-class
#' @import methods
#' @importClassesFrom bluster BlusterParam
setClass("MbkmeansParam", 
         slots=c(centers="integer_OR_function", 
                 extra.args="list"),
         contains="BlusterParam")

