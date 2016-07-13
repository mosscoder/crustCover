pot <- function(){
  file() <- list.files(path="./vis")[1]

  vis.jpeg <- readJPEG(paste("./vis/",file,sep=""))

  vis.red <- raster(vis.jpeg[,,1])
  vis.green <- raster(vis.jpeg[,,2])
  vis.blue <- raster(vis.jpeg[,,3])

  rgb <- stack(vis.red, vis.green, vis.blue)

  plotRGB(rgb, scale=1)

  pot.coords <- data.frame(x=numeric(),
                           y=numeric())
  for(i in 1:16){
    pot.coords[i,1:2] <- click(xy=T)[1:2]
  }

  poly <- Polygon(pot.coords)
  ps = Polygons(list(poly),1)
  pot.poly = SpatialPolygons(list(ps))

  plot(pot.poly, col="green", add=T)
  return(pot.poly)
}
