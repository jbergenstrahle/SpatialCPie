#' Parse spot detector output
#'
#' Parses the output from the ST spot detector tool for use with SpatialCPie.
#'
#' @param file spot file
#' @return `data.frame` with columns "x" and "y" specifying the pixel
#' coordinates of each spot
#' @export
#' @importFrom utils read.table
#' @examples
#' ## Create spot file
#' data <- matrix(
#'     c(
#'         c(7, 18, 7.00, 18.07, 563.2, 947.0),
#'         c(8, 11, 8.00, 11.04, 612.5, 627.7)
#'     ),
#'     byrow = TRUE,
#'     nrow = 2
#' )
#' filename <- tempfile()
#' write.table(
#'     data,
#'     file = filename,
#'     sep = '\t',
#'     quote = FALSE,
#'     col.names = c("x", "y", "new_x", "new_y", "pixel_x", "pixel_y")
#' )
#'
#' ## Parse spot file
#' parseSpotFile(filename)
#'
#' ## Delete spot file
#' unlink(filename)
parseSpotFile <- function(file){
    spots <- read.table(file, header = TRUE)
    xcoord <- as.numeric(spots$pixel_x)
    ycoord <- as.numeric(spots$pixel_y)
    coords <- as.data.frame(cbind(x = xcoord, y = ycoord))
    rownames(coords) <- paste(spots$x, spots$y, sep = "x")
    coords
}
