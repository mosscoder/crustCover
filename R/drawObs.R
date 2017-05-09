drawObs <- function(){

  file <- list.files(path="./vis")[1]

  vis.jpeg <- readJPEG(paste("./vis/",file,sep=""))

  vis.red <- raster(vis.jpeg[,,1])
  vis.green <- raster(vis.jpeg[,,2])
  vis.blue <- raster(vis.jpeg[,,3])

  rgb <- stack(vis.red, vis.green, vis.blue)

  options(warn = -1)
  plotRGB(rgb, scale = 1, asp = nrow(vis.red)/ncol(vis.red))
  message("Click at points along the boundary of the observation area in the plotted image. Press the escape key when finished.")
  poly <- drawPoly()
  options(warn = 0)
  message("Extracting cells. Please wait.")
  cells <- data.frame(extract(vis.red, poly, cellnumbers = T))[,1]
  out <- data.frame(cells, rowColFromCell(vis.red, cells), xyFromCell(vis.red,cells))
  return(out)
}
