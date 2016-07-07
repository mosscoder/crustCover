chart <- function(){
  file <- list.files(path="./vis")[1]

  vis.jpeg <- readJPEG(paste("./vis/",file,sep=""))

  vis.red <- raster(vis.jpeg[,,1])

  plot(vis.red, col=gray(1:100/100))

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
