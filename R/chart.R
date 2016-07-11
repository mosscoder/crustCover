chart <- function(){
  file <- list.files(path="./vis")[1]

  vis.jpeg <- readJPEG(paste("./vis/",file,sep=""))

  vis.red <- raster(vis.jpeg[,,1])
  vis.green <- raster(vis.jpeg[,,2])
  vis.blue <- raster(vis.jpeg[,,3])

  rgb <- stack(vis.red, vis.green, vis.blue)

  plotRGB(rgb, scale=1)

  chart.coords <- data.frame(x=numeric(),
                             y=numeric())
  for(i in 1:24){
    chart.coords[i,1:2] <- click(xy=T)[1:2]
  }

  sp.chart <- SpatialPoints(chart.coords)
  chart.buff <- gBuffer(sp.chart, width = .01, byid=T)

  plot(chart.buff, add=T, col="green")
  return(chart.buff)
}
