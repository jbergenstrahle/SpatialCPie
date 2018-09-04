context("Utility functions")


test_that("SpatialCPie::parseSpotFile parses spot files correctly", {
    data <- rbind(
        c(7, 18, 7.00, 18.07, 563.2, 947.0),
        c(8, 11, 8.00, 11.04, 612.5, 627.7)
    )
    filename <- tempfile()
    write.table(
        data,
        file = filename,
        sep = "\t",
        quote = FALSE,
        col.names = c("x", "y", "new_x", "new_y", "pixel_x", "pixel_y")
    )
    result <- parseSpotFile(filename)
    unlink(filename)

    expect_equal(
        result,
        data.frame(
            x = c(563.2, 612.5),
            y = c(947.0, 627.7),
            row.names = c("7x18", "8x11")
        )
    )
})
