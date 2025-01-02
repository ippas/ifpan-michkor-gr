standardize_median <- function (x, factor, control.samples = NA, only.mean = TRUE) 
{
  if (!is.list(factor)) {
    stop("factor is not a list")
  }
  if (length(x) != length(factor[[1]])) {
    stop("lengths of data and factor differ")
  }
  factor <- do.call(paste, factor)
  x.new <- x
  control.bool <- F
  if (!is.na(control.samples[1])) {
    control.bool <- T
    if (length(x) != length(factor)) {
      stop("lengths of data and factor differ")
    }
    control.samples <- which(control.samples)
  }
  else {
    control.samples <- 1:length(x)
  }
  x.mean <- mean(x[control.samples])
  x.sd <- sd(x[control.samples])
  for (group in factor) {
    wh <- which(factor == group)
    control.wh <- wh
    if (control.bool) {
      control.wh <- intersect(wh, control.samples)
    }
    if (only.mean) {
      group.mean <- median(x[control.wh])
      x.new[wh] <- x[wh] - group.mean
    }
    else {
      group.mean <- mean(x[control.wh])
      group.sd <- sd(x[control.wh])
      x.new[wh] <- (x[wh] - group.mean)/group.sd
    }
  }
  if (only.mean) {
    x.new <- x.new + x.mean
  }
  else {
    x.new <- x.new * x.sd + x.mean
  }
  x.new
}