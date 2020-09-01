
#this code uses parts of baggedETS function from the forecast package
MBB <- function(x, window_size) {

  bx = array(0, (floor(length(x)/window_size)+2)*window_size)
  for (i in 1:(floor(length(x)/window_size)+2)){
    c <- sample(1:(length(x)-window_size+1),1)
    bx[((i-1)*window_size+1):(i*window_size)] <- x[c:(c+window_size-1)]
  }
  start_from <- sample(0:(window_size-1),1) + 1
  bx[start_from:(start_from+length(x)-1)]
}



bld.mbb.bootstrap<-function (x, num, block_size = if (frequency(x) > 1) 2 * frequency(x) else 8)
{
  freq <- frequency(x)
  xs <- list()
  xs[[1]] <- x
  if (num > 1) {
    lambda <- BoxCox.lambda(x, lower = 0, upper = 1)
    x.bc <- BoxCox(x, lambda)
    if (freq > 1) {
      x.stl <- stl(ts(x.bc, frequency = freq), "per")$time.series
      seasonal <- x.stl[, 1]
      trend <- x.stl[, 2]
      remainder <- x.stl[, 3]
    }
    else {
      trend <- 1:length(x)
      suppressWarnings(x.loess <- loess(x.bc ~ trend, span = 6/length(x),
                                        degree = 1))
      seasonal <- rep(0, length(x))
      trend <- x.loess$fitted
      remainder <- x.loess$residuals
    }
    for (i in 2:num) {
      xs[[i]] <- InvBoxCox(trend + seasonal + MBB(remainder,
                                                  block_size), lambda)
    }
  }
  xs
}




forecast.baggedClusterETS<-function (object, cores=detectCores()-1,h = ifelse(frequency(object$x) > 1, 2 * frequency(object$x),
                             10), ...)
{

  if(Sys.info()[1]=="Windows"){
    out <- list(model = object, series = object$series, x = object$y,
                method = object$method)
    tspx <- tsp(out$x)

    cl <- makeCluster(getOption("cl.cores", cores))
    clusterExport(cl=cl, varlist=c("forecast"))
    forecasts_boot <- parLapply(cl,out$model$models,function(mod) {
      forecast(mod, PI = FALSE, h = h)$mean
    })

    forecasts_boot <- as.matrix(as.data.frame(forecasts_boot))
    colnames(forecasts_boot) <- NULL
    if (!is.null(tspx))
      start.f <- tspx[2] + 1/frequency(out$x)
    else start.f <- length(out$x) + 1
    out$forecasts_boot <- forecasts_boot
    out$mean <- ts(apply(forecasts_boot, 1, mean), frequency = frequency(out$x),
                   start = start.f)
    out$median <- ts(apply(forecasts_boot, 1, median))
    out$lower <- ts(apply(forecasts_boot, 1, min))
    out$upper <- ts(apply(forecasts_boot, 1, max))
    out$level <- 100
    tsp(out$median) <- tsp(out$lower) <- tsp(out$upper) <- tsp(out$mean)
    class(out) <- "forecast"
    out
  }else{
    out <- list(model = object, series = object$series, x = object$y,
                method = object$method)
    tspx <- tsp(out$x)

    forecasts_boot <- mclapply(out$model$models,mc.cores=cores ,function(mod) {
      forecast(mod, PI = FALSE, h = h)$mean
    })
    forecasts_boot <- as.matrix(as.data.frame(forecasts_boot))
    colnames(forecasts_boot) <- NULL
    if (!is.null(tspx))
      start.f <- tspx[2] + 1/frequency(out$x)
    else start.f <- length(out$x) + 1
    out$forecasts_boot <- forecasts_boot
    out$mean <- ts(apply(forecasts_boot, 1, mean), frequency = frequency(out$x),
                   start = start.f)
    out$median <- ts(apply(forecasts_boot, 1, median))
    out$lower <- ts(apply(forecasts_boot, 1, min))
    out$upper <- ts(apply(forecasts_boot, 1, max))
    out$level <- 100
    tsp(out$median) <- tsp(out$lower) <- tsp(out$upper) <- tsp(out$mean)
    class(out) <- "forecast"
    out
  }
}






