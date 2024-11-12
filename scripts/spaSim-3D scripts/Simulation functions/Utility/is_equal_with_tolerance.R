is_equal_with_tolerance <- function(x, y, tolerance = 1e-6) {
  abs(x - y) <= tolerance
}
