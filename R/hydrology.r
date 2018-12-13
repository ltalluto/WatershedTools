#' @export
discharge <- function(ws, Q, x, y) {
	Q <- data.frame(x = x, y = y, Q = Q)
	sp::coordinates(Q) <- c(1,2)
	Q <- snapToStream(Q, )
}