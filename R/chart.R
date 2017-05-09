chart <- function (samp.width = 0.01){
  file <- list.files(path = "./vis")[1]
  vis.jpeg <- readJPEG(paste("./vis/", file, sep = ""))
  vis.red <- raster(vis.jpeg[, , 1])
  vis.green <- raster(vis.jpeg[, , 2])
  vis.blue <- raster(vis.jpeg[, , 3])
  rgb <- stack(vis.red, vis.green, vis.blue)
  options(warn = -1)
  plotRGB(rgb, scale = 1, asp = nrow(vis.red)/ncol(vis.red))
  options(warn = 0)
  chart.coords <- data.frame(x = numeric(), y = numeric())
  message("Click on all 24 color chart cells in sequence. The sequence follows left to right and starts at cell 1 (brown, top left) and finishes on cell 24 (black, bottom right).")
  for (i in 1:24) {
    options(warn = -1)
    chart.coords[i, 1:2] <- click(xy = T)[1:2]
    options(warn = 0)
  }
  sp.chart <- SpatialPoints(chart.coords)
  chart.buff <- gBuffer(sp.chart, width = samp.width, byid = T)
  plot(chart.buff, add = T, col = "green")
  return(chart.buff)
}
