 order <- function ( ...) {
   UseMethod("order", ...)
 }

order.default <- function(...){
  base::order(...)
}
