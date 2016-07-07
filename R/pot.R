pot <- function(){
  file <- list.files(path="./vis")[1]
  
  vis.jpeg <- readJPEG(paste("./vis/",file,sep=""))
  
  vis.red <- raster(vis.jpeg[,,1])
  
  plot(vis.red, col=gray(1:100/100))
  
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