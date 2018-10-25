context("Helper functions")

test_that("SpatialCPie:::.likeness outputs expected values", {
    # Same distance to all clusters
    expect_equal(
        .likeness(matrix(rep(1, 3), nrow = 1)),
        matrix(rep(1 / 3, 3), nrow = 1)
    )

    # Zero log-multiplier
    expect_equal(
        .likeness(matrix(1:3, nrow = 1), c = 0.0),
        matrix(rep(1 / 3, 3), nrow = 1)
    )

    # Entropy is decreasing with log-multiplier
    set.seed(123)
    m <- matrix(rexp(1000), nrow = 100)
    max_likeness <- function(c)
        apply(m, 1, function(xs) max(.likeness(xs, c = c)))
    expect_true(all(max_likeness(1.0) < max_likeness(1.0 + 1e-3)))
})

test_that("SpatialCPie:::.maximizeOverlap outputs expected values", {
    expect_equal(
        .maximizeOverlap(list(
            "2" = setNames(c(rep("1", 5), rep("2", 5)), letters[1:10]),
            "3" = setNames(c(rep("3", 4), rep("1", 2), rep("2", 4)), letters[1:10]),
            "4" = setNames(c(rep("2", 3), rep("4", 3), rep("3", 1), rep("1", 3)), letters[1:10])
        )),
        list(
            "2" = setNames(c(rep("1", 5), rep("2", 5)), letters[1:10]),
            "3" = setNames(c(rep("1", 4), rep("3", 2), rep("2", 4)), letters[1:10]),
            "4" = setNames(c(rep("1", 3), rep("3", 3), rep("4", 1), rep("2", 3)), letters[1:10])
        )
    )
})

test_that("SpatialCPie:::.maximizeOverlap is idempotent", {
    maxOverlap <- .maximizeOverlap(list(
        "2" = setNames(c(rep("1", 5), rep("2", 5)), letters[1:10]),
        "3" = setNames(c(rep("3", 4), rep("1", 2), rep("2", 4)), letters[1:10]),
        "4" = setNames(c(rep("2", 3), rep("4", 3), rep("3", 1), rep("1", 3)), letters[1:10])
    ))
    expect_equal(.maximizeOverlap(maxOverlap), maxOverlap)
})
