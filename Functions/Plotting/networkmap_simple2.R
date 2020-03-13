networkmap_simple2 <- function (CMAT, lat, lon, weights = rep(1, dim(CMAT)[1]), thresh = NULL, title) {
  noPts <- dim(CMAT)[1]
  rbPal <- colorRampPalette(c("red", "grey", "blue"))
  COLZ <- array(rbPal(9)[as.numeric(cut(c(-1, 1, c(CMAT)), 
                                        breaks = 9))][-c(1:2)], dim = dim(CMAT))
  if (is.numeric(thresh) && (thresh > 0) && (thresh < 1)) {
    ind <- which(abs(CMAT) < thresh, arr.ind = TRUE)
    CMAT[ind] <- NA
  }
  plot(c(-180, 180), c(-60, 80), type = "n", xlab = "long", ylab= "lat", main = title)
  maps::map("world", add = TRUE, col = "grey", interior = FALSE)
  syms <- apply(CMAT, 2, function(x) {
    all(is.na(x))
  })
  if (any(syms)) {
    symbols(x = lon[syms], y = lat[syms], squares = sqrt(weights[syms]), 
            inches = 1/20, bg = rgb(0.1, 0.9, 0.1, 0.5), fg = rgb(0.1, 
                                                                  0.9, 0.1), add = TRUE)
    symbols(x = lon[!syms], y = lat[!syms], circles = sqrt(weights[!syms]), 
            inches = 1/20, bg = rgb(0.1, 0.1, 0.9, 0.5), fg = rgb(0.1, 
                                                                  0.1, 0.9), add = TRUE)
  }
  else {
    symbols(x = lon, y = lat, circles = sqrt(weights), inches = 1/20, 
            bg = rgb(0.1, 0.1, 0.9, 0.5), fg = rgb(0.1, 0.1, 
                                                   0.9), add = TRUE)
  }
  for (i in 1:(noPts - 1)) {
    for (j in (i + 1):(noPts)) {
      if (!is.na(CMAT[i, j])) {
        lines(c(lon[i], lon[j]), c(lat[i], lat[j]), col = COLZ[i, 
                                                               j], lwd = 1.5)
      }
    }
  }
}
