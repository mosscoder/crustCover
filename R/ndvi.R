ndvi <- function(chart, pot, threshold = 0.3, agg.fact=10){
  t <- as.numeric(list(threshold)[1])
  agg <- as.numeric(list(agg.fact)[1])

  vis.files <- list.files(path = "./vis")
  nir.files <- list.files(path = "./nir")

  df <- data.frame(unit <- character(),
                   median.ndvi <- numeric(),
                   mean.ndvi <- numeric(),
                   threshold <- numeric(),
                   percent.cover <- numeric())

  colnames(df) <- c("unit","median.ndvi","mean.ndvi","threshold","percent.cover")

  write.csv(df, "ndvi.data.csv", row.names = F)

  if(file.exists("names.csv")){
    names <- c(as.character(read.csv("names.csv")[,1]))
  }else{
    names <- c(names = paste("pot_", 1:length(vis.files), sep=""))
  }

  data(chart.vals)

  calcs <- function(x){

    vis.jpeg <- readJPEG(paste("./vis/",vis.files[x],sep=""))

    vis.red <- raster(vis.jpeg[,,1])
    vis.green <- raster(vis.jpeg[,,2])
    vis.blue <- raster(vis.jpeg[,,3])

    vis.red.df <- as.data.frame(rasterToPoints(vis.red))[,3]
    vis.green.df <- as.data.frame(rasterToPoints(vis.green))[,3]
    vis.blue.df <- as.data.frame(rasterToPoints(vis.blue))[,3]

    vis.rgb.df <- data.frame(vis.red.df,vis.green.df,vis.blue.df)

    vis.color.obj <- colorspace::RGB(vis.rgb.df[,1],vis.rgb.df[,2],vis.rgb.df[,3])
    vis.lab.unequal <- coords(as(vis.color.obj, "LAB"))
    vis.lab.equal <- histeq((vis.lab.unequal[,1]/100)*256)

    vis.template <- rasterToPoints(vis.red)

    vis.equal.ras <- rasterFromXYZ(as.matrix(data.frame(vis.template[,1:2],vis.lab.equal)))

    vis.cal <- data.frame(val = numeric())

    for(i in 1:24){
      poly <- chart[i]
      vis.cal[i,1] <- mean(unlist(extract(vis.equal.ras, poly)))/256
    }

    vis.cal.df <- data.frame(vis.cal, chart.vals[,2])
    colnames(vis.cal.df) <- c("x","y")
    vis.cal.df <- vis.cal.df[order(vis.cal.df$x),]

    vis.mod <- nls(y ~ I(a*exp(b*x)), data = vis.cal.df, start = list(a = 0.01, b = 1))

    vis.predictions <- predict(vis.mod, list(x = vis.lab.equal/256))

    vis.color.lab.obj <- colorspace::LAB(vis.lab.unequal[,1], vis.lab.unequal[,2], vis.predictions*100)
    vis.equal.rgb <- coords(as(vis.color.lab.obj, "RGB"))
    vis.equal.r <- rasterFromXYZ(as.matrix(data.frame(vis.template[,1:2],vis.equal.rgb[,1])))

    #### NIR ####

    nir.jpeg <- readJPEG(paste("./nir/",nir.files[x],sep=""))

    nir.red <- raster(nir.jpeg[,,1])
    nir.green <- raster(nir.jpeg[,,2])
    nir.blue <- raster(nir.jpeg[,,3])

    nir.red.df <- as.data.frame(rasterToPoints(nir.red))[,3]
    nir.green.df <- as.data.frame(rasterToPoints(nir.green))[,3]
    nir.blue.df <- as.data.frame(rasterToPoints(nir.blue))[,3]

    nir.rgb.df <- data.frame(nir.red.df,nir.green.df,nir.blue.df)

    nir.color.obj <- colorspace::RGB(nir.rgb.df[,1],nir.rgb.df[,2],nir.rgb.df[,3])
    nir.lab.unequal <- coords(as(nir.color.obj, "LAB"))
    nir.lab.equal <- histeq((nir.lab.unequal[,1]/100)*256)

    nir.template <- rasterToPoints(nir.red)
    nir.equal.ras <- rasterFromXYZ(as.matrix(data.frame(nir.template[,1:2],nir.lab.equal)))

    nir.cal <- data.frame(val = numeric())

    for(i in 1:24){
      poly <- chart[i]
      nir.cal[i,1] <- mean(unlist(extract(nir.equal.ras, poly)))/256
    }

    nir.cal.df <- data.frame(nir.cal, chart.vals[,3])
    colnames(nir.cal.df) <- c("x","y")
    nir.cal.df <- nir.cal.df[order(nir.cal.df$x),]

    nir.mod <- nls(y ~ I(a*exp(b*x)), data = nir.cal.df, start = list(a = 0.01, b = 1))

    nir.predictions <- predict(nir.mod, list(x = nir.lab.equal/256))

    nir.equal.refl <- rasterFromXYZ(as.matrix(data.frame(nir.template[,1:2],nir.predictions)))+(10/256)

    #### NDVI ####

    ndvi <- aggregate((nir.equal.refl - vis.equal.r)/(nir.equal.refl + vis.equal.r), fact=agg)
    pot.ndvi <- mask(crop(ndvi, pot), pot)
    writeRaster(pot.ndvi, paste("./ndvi/", names[x], "_ndvi_ras.tif", sep = ""), format = "GTiff", overwrite = T)
    pot.cut <- pot.ndvi > t
    writeRaster(pot.cut, paste("./ndvi/", names[x], "_cover_ras.tif", sep = ""), format = "GTiff", overwrite = T)

    png(filename=paste("./ndvi/", names[x],".png",sep=""), w=5, h=15, units="in", res=600)
    par(mfrow=c(3,1))
    plot(pot.ndvi, col=gray(1:100/100), main=paste(names[x],"NDVI Values"),axes=FALSE, box=FALSE)
    hist(pot.ndvi, breaks=1000, main="Distribution")
    plot(pot.cut, col=c("black","green"), legend=F, main=paste("Binary Cover at",threshold),axes=FALSE, box=FALSE)
    dev.off()

    dat <- read.csv("ndvi.data.csv")
    ndvi.pts <- rasterToPoints(pot.ndvi)[,3]
    ndvi.mean <- mean(ndvi.pts)
    ndvi.median <- median(ndvi.pts)
    percent.cover <- nrow(rasterToPoints(reclassify(pot.cut, rcl=cbind(0,NA))))/length(ndvi.pts)

    new.dat <- data.frame(names[x], ndvi.median, ndvi.mean, t, percent.cover)
    colnames(new.dat) <- colnames(dat)
    dat.bind <- rbind(dat, new.dat)
    write.csv(dat.bind, "ndvi.data.csv", row.names = F)

  }
  lapply(FUN=calcs, X=1:length(names))
}
