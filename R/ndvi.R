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
    names <- c(names = paste("obs.", 1:length(vis.files), sep=""))
  }

  data(chart.vals)

  calcs <- function(x){

    vis.jpeg <- readJPEG(paste("./vis/",vis.files[x],sep=""))

    vis.red <- raster(vis.jpeg[,,1])
    vis.red.df <- as.data.frame(rasterToPoints(vis.red))[,3]
    vis.template <- rasterToPoints(vis.red)

    vis.cal.df <- data.frame(x=numeric(),
                             y=numeric())

    for(i in 1:24){
      poly <- chart[i]
      df <- data.frame(x = extract(vis.red, poly), y = chart.vals[i,2])
      df.samp <- df[sample(x=1:nrow(df),size=10,replace=F),]
      colnames(df.samp) <- c("x","y")
      vis.cal.df <- rbind(vis.cal.df,df.samp )
    }

    vis.cal.df <- vis.cal.df[order(vis.cal.df$x),]

    vis.svm_tune <- tune(svm, train.x=vis.cal.df$x, train.y=vis.cal.df$y,
                         kernel="radial", ranges=list(cost=1:10, gamma=c(0.1,1,10,100,1000,10000)))

    vis.svm.model <- svm(y ~ x,
                         data = vis.cal.df,
                         cost = vis.svm_tune$best.parameters$cost,
                         gamma = vis.svm_tune$best.parameters$gamma)
   # plot(vis.cal.df)
    #lines(vis.cal.df$x, predict(vis.svm.model, new.data = (list(x = vis.cal.df$x))), col = "green")
    vis.predictions <- predict(vis.svm.model, data.frame(x = vis.red.df))

    vis.cal.r <- rasterFromXYZ(as.matrix(data.frame(vis.template[,1:2],vis.predictions)))

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
    nir.lab <- coords(as(nir.color.obj, "LAB"))[,1]/100

    nir.template <- rasterToPoints(nir.red)
    nir.refl.ras <- rasterFromXYZ(as.matrix(data.frame(nir.template[,1:2], nir.lab)))

    nir.cal.df <- data.frame(x=numeric(),
                             y=numeric())

    for(i in 1:24){
      poly <- chart[i]
      df <- data.frame(x = extract(nir.refl.ras, poly), y = chart.vals[i,2])
      df.samp <- df[sample(x=1:nrow(df),size=10,replace=F),]
      colnames(df.samp) <- c("x","y")
      nir.cal.df <- rbind(nir.cal.df, df.samp )
    }

    nir.cal.df <- nir.cal.df[order(nir.cal.df$x),]

    nir.svm_tune <- tune(svm, train.x=nir.cal.df$x, train.y=nir.cal.df$y,
                         kernel="radial", ranges=list(cost=1:10, gamma=c(0.1,1,10,100,1000,10000)))

    nir.svm.model <- svm(y ~ x,
                         data = nir.cal.df,
                         cost = nir.svm_tune$best.parameters$cost,
                         gamma = nir.svm_tune$best.parameters$gamma)
    #plot(nir.cal.df)
    #lines(nir.cal.df$x, predict(nir.svm.model, new.data = (list(x = nir.cal.df$x))), col = "green")
    nir.predictions <- predict(nir.svm.model, data.frame(x = nir.lab))

    nir.equal.refl <- rasterFromXYZ(as.matrix(data.frame(nir.template[,1:2], nir.predictions)))

    #### NDVI ####
    ndvi <- aggregate((nir.equal.refl - vis.cal.r)/(nir.equal.refl + vis.cal.r), fact=agg)
    pot.ndvi <- mask(crop(ndvi, pot), pot)
    writeRaster(pot.ndvi, paste("./ndvi/", names[x], "_ndvi_ras.tif", sep = ""), format = "GTiff", overwrite = T)
    pot.cut <- pot.ndvi > t
    writeRaster(pot.cut, paste("./ndvi/", names[x], "_cover_ras.tif", sep = ""), format = "GTiff", overwrite = T)

    png(filename=paste("./ndvi/", names[x],".png",sep=""), w=5, h=15, units="in", res=600)
    par(mfrow=c(3,1))
    plot(pot.ndvi, col=gray(1:100/100), main=paste(names[x],"NDVI Values"),axes=FALSE, box=FALSE)
    hist(pot.ndvi, breaks=1000, main="Distribution")
    plot(pot.cut, col=c("black","green"), legend=F, main=paste("Binary Cover at", t),axes=FALSE, box=FALSE)
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
