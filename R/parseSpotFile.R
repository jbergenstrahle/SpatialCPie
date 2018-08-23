#' Parsing of spot detector output
#'
#' A function that parses the output from the ST spot detector tool to receive pixel coordinates for the spatial areas.
#'
#' @param file A spot file that contains the pixel coordinates for each spatial area after alignment to the tissue image via the
#' @return FIXME
#' ST spot detector tool (https://github.com/SpatialTranscriptomicsResearch/st_spot_detector).
#' @keywords spotdetector
#' @export
#' @examples
#' ## Create spot file
#' data <- matrix(
#'   c(
#'     c(7, 18, 7.00, 18.07, 563.2, 947.0),
#'     c(8, 11, 8.00, 11.04, 612.5, 627.7)
#'   ),
#'   byrow = T,
#'   nrow = 2
#' )
#' filename <- tempfile()
#' write.table(
#'   data,
#'   file = filename,
#'   sep = '\t',
#'   quote = F,
#'   col.names = c("x", "y", "new_x", "new_y", "pixel_x", "pixel_y")
#' )
#'
#' ## Parse spot file
#' parseSpotFile(filename)
#'
#' ## Delete spot file
#' unlink(filename)
parseSpotFile <- function(file){
  spots <- read.table(file, header = T)
  xcoord <- as.numeric(spots$pixel_x)
  ycoord <- as.numeric(spots$pixel_y)
  coords <- as.data.frame(cbind(pixel_x = xcoord, pixel_y = ycoord))
  rownames(coords) <- paste(spots$x, spots$y, sep = "x")
  return(coords)
}
