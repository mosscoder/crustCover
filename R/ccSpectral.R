ccSpectral <- function(chart, obs.area, rasters = F, ml = F, ml.cutoff = 0.9, thresholds = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)){

  vis.files <- list.files(path = "./vis")
  nir.files <- list.files(path = "./nir")

  out.dir <- paste("output",Sys.time())
  dir.create(out.dir)

  df <- data.frame(unit <- character(),
                   vis.file <- character(),
                   nir.file <- character(),
                   red.rsq <- numeric(),
                   green.rsq <- numeric(),
                   blue.rsq <- numeric(),
                   nir.rsq <- numeric(),
                   ndvi.median <- numeric(),
                   ndvi.mean <- numeric(),
                   ndvi.threshold <- numeric(),
                   ndvi.cover <- numeric(),
                   vi.median <- numeric(),
                   vi.mean <- numeric(),
                   vi.threshold <- numeric(),
                   vi.cover <- numeric(),
                   msavi.median <- numeric(),
                   msavi.mean <- numeric(),
                   msavi.threshold <- numeric(),
                   msavi.cover <- numeric(),
                   evi.median <- numeric(),
                   evi.mean <- numeric(),
                   evi.threshold <- numeric(),
                   evi.cover <- numeric(),
                   bsci.median <- numeric(),
                   ci.mean <- numeric(),
                   ci.threshold <- numeric(),
                   ci.cover <- numeric(),
                   bsci.median <- numeric(),
                   bsci.mean <- numeric(),
                   bsci.threshold <- numeric(),
                   bsci.cover <- numeric(),
                   bi.median <- numeric(),
                   bi.mean <- numeric(),
                   bi.threshold <- numeric(),
                   bi.cover <- numeric())

  colnames(df) <- c("unit", "vis.file","nir.file", "red.rsq", "green.rsq","blue.rsq","nir.rsq",
                    "ndvi.median","ndvi.mean","ndvi.threshold", "ndvi.cover",
                    "sr.median","sr.mean","sr.threshold", "sr.cover",
                    "msavi.median","msavi.mean","msavi.threshold", "msavi.cover",
                    "evi.median","evi.mean","evi.threshold", "evi.cover",
                    "ci.median","ci.mean","ci.threshold", "ci.cover",
                    "bsci.median","bsci.mean","bsci.threshold", "bsci.cover",
                    "bi.median","bi.mean","bi.threshold", "bi.cover")

  write.csv(df, paste(out.dir,"/summary_data.csv", sep = ""), row.names = F)

  if(file.exists("names.csv")){
    names <- c(as.character(read.csv("names.csv")[,1]))
  }else{
    names <- c(names = paste("obs_", 1:length(vis.files), sep=""))
  }

  calcs <- function(x){

    vis.jpeg <- readJPEG(paste("./vis/",vis.files[x],sep=""))

    vis.red <-raster(vis.jpeg[,,1])
    vis.green <-raster(vis.jpeg[,,2])
    vis.blue <-raster(vis.jpeg[,,3])

    nir.jpeg <- readJPEG(paste("./nir/",nir.files[x],sep=""))

    nir.blue <- raster(nir.jpeg[,,3]) + 10/256

    asp = nrow(vis.red)/ncol(vis.red)

    all.bands <- stack(vis.red, vis.green, vis.blue, nir.blue)
    names(all.bands) <- c("vis.red", "vis.green", "vis.blue", "nir.blue")

    obs.ext <- extent(
      min(obs.area$x),
      max(obs.area$x),
      min(obs.area$y),
      max(obs.area$y))

    temp.mat <- raster(matrix(data = NA,
                              nrow = nrow(all.bands),
                              ncol = ncol(all.bands),
                              byrow = T))

    bands.df <- data.frame(extract(all.bands, obs.area$cells))
    colnames(bands.df) <- c("vis.red", "vis.green", "vis.blue", "nir.blue")

    train.df <- data.frame()

    chart.vals <- data.frame(red.chart = c(0.17,0.63,0.15,0.11,0.31,0.2,0.63,0.12,0.57,0.21,0.33,0.67,0.04,0.1,0.6,0.79,0.7,0.07,0.93,0.59,0.36,0.18,0.08,0.03),
                             green.chart = c(0.1,0.32,0.19,0.14,0.22,0.47,0.27,0.11,0.13,0.06,0.48,0.4,0.06,0.27,0.07,0.62,0.13,0.22,0.95,0.62,0.38,0.2,0.09,0.03),
                             blue.chart = c(0.07,0.24,0.34,0.06,0.42,0.42,0.06,0.36,0.12,0.14,0.1,0.06,0.24,0.09,0.04,0.08,0.31,0.38,0.93,0.62,0.39,0.2,0.09,0.02),
                             nir.chart = c(0.43,0.87,0.86,0.18,0.86,0.43,0.85,0.54,0.54,0.79,0.49,0.66,0.52,0.44,0.72,0.82,0.88,0.42,0.91,0.51,0.27,0.13,0.06,0.02))

    for(i in c(1:24)[-3]){
      poly <- chart[i]
      options(warn = -1)
      df.samp <- data.frame(chart.vals[i,], extract(all.bands, poly))
      options(warn = 0)

      if(nrow(df.samp) >= 50){
        df.samp <- df.samp[sample(x=1:nrow(df.samp),size=50,replace=F),]
      }
      train.df <- rbind(train.df, df.samp )
    }


    #### Correct Visible Red Band ####
    red.nls <- nls(red.chart ~ (a*exp(b*vis.red)), trace = F, data = train.df,
                   start = c(a = 0.1, b = 0.1))
    red.preds <- predict(red.nls, bands.df)
    red.rsq <- 1-sum((train.df$red.chart - predict(red.nls, train.df))^2)/(length(train.df$red.chart)*var(train.df$red.chart))

    red.mat <- temp.mat
    values(red.mat)[obs.area$cells] <- red.preds
    red.mat <- crop(red.mat, extent(obs.ext))

    #### Correct Visible Green Band ####
    green.nls <- nls(green.chart ~ (a*exp(b*vis.green)), trace = F, data = train.df,
                   start = c(a = 0.1, b = 0.1))
    green.preds <- predict(green.nls, bands.df)
    green.rsq <- 1-sum((train.df$green.chart - predict(green.nls, train.df))^2)/(length(train.df$green.chart)*var(train.df$green.chart))

    green.mat <- temp.mat
    values(green.mat)[obs.area$cells] <- green.preds
    green.mat <- crop(green.mat, extent(obs.ext))

    #### Correct Visible Blue Band ####
    blue.nls <- nls(blue.chart ~ (a*exp(b*vis.blue)), trace = F, data = train.df,
                     start = c(a = 0.1, b = 0.1))
    blue.preds <- predict(blue.nls, bands.df)
    blue.rsq <- 1-sum((train.df$blue.chart - predict(blue.nls, train.df))^2)/(length(train.df$blue.chart)*var(train.df$blue.chart))

    blue.mat <- temp.mat
    values(blue.mat)[obs.area$cells] <- blue.preds
    blue.mat <- crop(blue.mat, extent(obs.ext))

    ##### Correct NIR Blue Band ####

    nir.nls <- nls(nir.chart ~  (a*exp(b*nir.blue)), trace = F, data = train.df,
                   start = c(a = 0.1, b = 0.1))
    nir.preds <- predict(nir.nls, bands.df)
    nir.rsq <- 1-sum((train.df$nir.chart - predict(nir.nls, train.df))^2)/(length(train.df$nir.chart)*var(train.df$nir.chart))

    nir.mat <- temp.mat
    values(nir.mat)[obs.area$cells] <- nir.preds
    nir.mat <- crop(nir.mat, extent(obs.ext))

    if(ml == T){
      rf.mod <- rfsrc(cbind(nir.chart + red.chart + green.chart + blue.chart) ~ vis.red + vis.green + vis.blue,
                             data = train.df,
                             ntree = 250,
                             tree.err = T,
                             set.seed=123)
      rf.preds <- predict(rf.mod, bands.df)

      nir.ml <- temp.mat
      red.ml <- temp.mat
      green.ml <- temp.mat
      blue.ml <- temp.mat

      nir.rsq.ml <- 1-sum((train.df$nir.chart - rf.mod$regrOutput$nir.chart$predicted.oob)^2)/(length(train.df$nir.chart)*var(train.df$nir.chart))
      red.rsq.ml <- 1-sum((train.df$red.chart - rf.mod$regrOutput$red.chart$predicted.oob)^2)/(length(train.df$red.chart)*var(train.df$red.chart))
      green.rsq.ml <- 1-sum((train.df$green.chart - rf.mod$regrOutput$green.chart$predicted.oob)^2)/(length(train.df$green.chart)*var(train.df$green.chart))
      blue.rsq.ml <- 1-sum((train.df$blue.chart - rf.mod$regrOutput$blue.chart$predicted.oob)^2)/(length(train.df$blue.chart)*var(train.df$blue.chart))

      values(nir.ml)[obs.area$cells] <- rf.preds$regrOutput$nir.chart$predicted
      values(red.ml)[obs.area$cells] <- rf.preds$regrOutput$red.chart$predicted
      values(green.ml)[obs.area$cells] <- rf.preds$regrOutput$green.chart$predicted
      values(blue.ml)[obs.area$cells] <- rf.preds$regrOutput$blue.chart$predicted

      nir.ml <- crop(nir.ml, extent(obs.ext))
      red.ml <- crop(red.ml, extent(obs.ext))
      green.ml <- crop(green.ml, extent(obs.ext))
      blue.ml <- crop(blue.ml, extent(obs.ext))

      if(nir.rsq < ml.cutoff){
        nir.mat <- nir.ml
        nir.rsq <- nir.rsq.ml
      }
      if(red.rsq < ml.cutoff){
        red.mat <- red.ml
        red.rsq <- red.rsq.ml
      }
      if(green.rsq < ml.cutoff){
        green.mat <- green.ml
        green.rsq <- green.rsq.ml
      }
      if(blue.rsq < ml.cutoff){
        blue.mat <- blue.ml
        blue.rsq <- blue.rsq.ml
      }
    }

    #### Spectral Calculations ####
    ndvi <-  (nir.mat - red.mat)/(nir.mat + red.mat)
    sr <- nir.mat/red.mat
    msavi <- (2*nir.mat + 1 - sqrt((2*nir.mat + 1)^2 - 8*(nir.mat - red.mat)))/2
    evi <- 2.5*((nir.mat - red.mat)/(nir.mat + 6*red.mat - 7.5*blue.mat + 1))
    ci <- 1 - (red.mat - blue.mat)/(red.mat + blue.mat)
    bsci <- (1 - 2*abs(red.mat - green.mat))/mean(stack(green.mat, red.mat, nir.mat))
    bi <- sqrt(green.mat^2 + red.mat^2 + nir.mat^2)

    #Apply Thresholds
    ndvi.cut <- ndvi >= thresholds[1]
    sr.cut <- sr >= thresholds[2]
    msavi.cut <- msavi >= thresholds[3]
    evi.cut <- evi >= thresholds[4]
    ci.cut <- ci >= thresholds[5]
    bsci.cut <- bsci >= thresholds[6]
    bi.cut <- bi <= thresholds[7]

    #Write data to raster format
    if(rasters == T){
      writeRaster(round(stack(ndvi,sr,msavi,evi,ci,bsci,bi)*1000),
                  paste(out.dir,"/" , names[x], "_spectral_stack.tif", sep = ""), format = "GTiff", overwrite = T)
      writeRaster(stack(ndvi.cut, sr.cut, msavi.cut, evi.cut, ci.cut, bsci.cut, bi.cut),
                  paste(out.dir,"/", names[x], "_cover_stack.tif", sep = ""), format = "GTiff", overwrite = T)
    }

    #Save Plots
    pal <- colorRampPalette(colors = rev(brewer.pal(11, "Spectral")))(100)
    pdf(file=paste(out.dir,"/", names[x],".pdf",sep=""), w=10, h=25)

    par(mfrow=c(7,3))

    hist(ndvi, breaks=1000, main="NDVI Distribution")
    plot(ndvi, col = pal, main="NDVI Values",axes=FALSE, box=FALSE, asp = asp)
    plot(ndvi.cut, col=c("black","green"), legend=F, main=paste("NDVI Binary Cover threshold", thresholds[1]),axes=FALSE, box=FALSE, asp = asp)

    hist(sr, breaks=1000, main="SR Distribution")
    plot(sr, col = pal, main="SR Values",axes=FALSE, box=FALSE, asp = asp)
    plot(sr.cut, col=c("black","green"), legend=F, main=paste("SR Binary Cover threshold", thresholds[2]),axes=FALSE, box=FALSE, asp = asp)


    hist(msavi, breaks=1000, main="MSAVI Distribution")
    plot(msavi, col = pal, main="MSAVI Values",axes=FALSE, box=FALSE, asp = asp)
    plot(msavi.cut, col=c("black","green"), legend=F, main=paste("MSAVI Binary Cover threshold", thresholds[3]),axes=FALSE, box=FALSE, asp = asp)

    hist(evi, breaks=1000, main="EVI Distribution")
    plot(evi, col = pal, main="EVI Values",axes=FALSE, box=FALSE, asp = asp)
    plot(evi.cut, col=c("black","green"), legend=F, main=paste("EVI Binary Cover threshold", thresholds[4]),axes=FALSE, box=FALSE, asp = asp)

    hist(ci, breaks=1000, main="Crust Index Distribution")
    plot(ci, col = pal, main="Crust Index Values",axes=FALSE, box=FALSE, asp = asp)
    plot(ci.cut, col=c("black","green"), legend=F, main=paste("Crust Index Binary Cover threshold", thresholds[5]),axes=FALSE, box=FALSE, asp = asp)

    hist(bsci, breaks=1000, main="BSCI Index Distribution")
    plot(bsci, col = pal, main="BSC Index Values",axes=FALSE, box=FALSE, asp = asp)
    plot(bsci.cut, col=c("black","green"), legend=F, main=paste("BSC Index Binary Cover threshold", thresholds[6]),axes=FALSE, box=FALSE, asp = asp)


    hist(bi, breaks=1000, main="Brightness Index Distribution")
    plot(bi, col = pal, main="Brightness Index Values",axes=FALSE, box=FALSE, asp = asp)
    plot(bi.cut, col=c("black","green"), legend=F, main=paste("Brightness Index Binary Cover threshold", thresholds[7]),axes=FALSE, box=FALSE, asp = asp)

    dev.off()

    #Save summary data
    dat <- read.csv(paste(out.dir,"/summary_data.csv", sep = ""))

    ndvi.mean <- cellStats(ndvi, stat = "mean")
    ndvi.median <- median(na.omit(values(ndvi)))
    ndvi.cover <- nrow(rasterToPoints(reclassify(ndvi.cut, rcl=cbind(0,NA))))/nrow(obs.area)

    sr.mean <- cellStats(sr, stat = "mean")
    sr.median <- median(na.omit(values(sr)))
    sr.cover <- nrow(rasterToPoints(reclassify(sr.cut, rcl=cbind(0,NA))))/nrow(obs.area)


    msavi.mean <- cellStats(msavi, stat = "mean")
    msavi.median <- median(na.omit(values(msavi)))
    msavi.cover <- nrow(rasterToPoints(reclassify(msavi.cut, rcl=cbind(0,NA))))/nrow(obs.area)

    evi.mean <- cellStats(evi, stat = "mean")
    evi.median <- median(na.omit(values(evi)))
    evi.cover <- nrow(rasterToPoints(reclassify(evi.cut, rcl=cbind(0,NA))))/nrow(obs.area)

    ci.mean <- cellStats(ci, stat = "mean")
    ci.median <- median(na.omit(values(ci)))
    ci.cover <- nrow(rasterToPoints(reclassify(ci.cut, rcl=cbind(0,NA))))/nrow(obs.area)

    bsci.mean <- cellStats(bsci, stat = "mean")
    bsci.median <- median(na.omit(values(bsci)))
    bsci.cover <- nrow(rasterToPoints(reclassify(bsci.cut, rcl=cbind(0,NA))))/nrow(obs.area)

    bi.mean <- cellStats(bi, stat = "mean")
    bi.median <- median(na.omit(values(bi)))
    bi.cover <- nrow(rasterToPoints(reclassify(bi.cut, rcl=cbind(0,NA))))/nrow(obs.area)


    new.dat <- data.frame(names[x], vis.files[x], nir.files[x], red.rsq, green.rsq, blue.rsq, nir.rsq,
                          ndvi.median, ndvi.mean,thresholds[1], ndvi.cover,
                          sr.median, sr.mean, thresholds[2], sr.cover,
                          msavi.mean, msavi.median, thresholds[3], msavi.cover,
                          evi.mean, evi.median, thresholds[4], evi.cover,
                          ci.mean, ci.median, thresholds[5], ci.cover,
                          bsci.mean, bsci.median, thresholds[6], bsci.cover,
                          bi.mean, bi.median, thresholds[7], bi.cover
                          )

    colnames(new.dat) <- colnames(dat)
    dat.bind <- rbind(dat, new.dat)
    write.csv(dat.bind, paste(out.dir,"/summary_data.csv", sep = ""), row.names = F)
    message(paste(names[x], "processed..."))
  }
  lapply(FUN=calcs, X=1:length(names))
  message("Processed files may be found at: ",
          paste(getwd(),"/output ",Sys.time(), sep = ""))
}



