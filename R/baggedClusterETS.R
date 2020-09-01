#' Perform BaggedClusterETS method
#'
#'
#' BaggedClusterETS method
#' @param y a time series
#' @param cores (optional) number of cores used to do the calculations
#' @param nclusters (optional) number of clusters used. If not selected it uses the Silhouette method to define it
#' @param distance distance function used to create the clusters. Default is the Euclidian distance but one can choose among many from the TSclust package
#' @param silhouette If TRUE the method performs the Silhouette method and ignores the nclusters selection
#' @param h_pseudo number of observations to the validation set that will be used to select the best series
#' @param boot_samples number of bootstrap series generated to perform the cluster phase
#' @param ... other parameters that affect the ETS estimation
#'
#' @return a list with several elements such as the generated bootstrapped series, fitted method to each of the bootstrapped series and so on
#' @import forecast
#' @import TSclust
#' @import parallel
#' @export
#' @examples \dontrun{
#' model <- baggedCLusterETS(gas)
#' plot(forecast(model))
#' }
baggedClusterETS<-function (y,cores=detectCores()-1,nclusters=5,distance="EUCL",silhouette=F,
                            h_pseudo = NULL,boot_samples=1000,
                            ...)
{
  #Windows OS
  if(Sys.info()[1]=="Windows"){
    cl <- makeCluster(getOption("cl.cores", cores))

    bootstrapped_series <- bld.mbb.bootstrap(y, boot_samples)

    if (is.null(h_pseudo)==T){
      h_pseudo <- ifelse(frequency(y) > 1, 2 * frequency(y), 10)
    }


    start.i <- tsp(y)[1]
    start.f <- tsp(y)[2] + 1/frequency(y)


    real.pre<-y[((length(y)-h_pseudo+1):length(y))]


    if(length(y)-h_pseudo>h_pseudo){
      pseudo<-parLapply(cl,bootstrapped_series, function(x) {
        ts(x[(1:(length(y)-h_pseudo))],frequency=frequency(y),start=start.i)
      })


      clusterExport(cl=cl, varlist=c("forecast","ets"))
      forecasts_pseudo <- parLapply(cl,pseudo, function(x) {
        mod <- forecast(ets(x),h=h_pseudo)$mean
      })



      resultado.pre.lista<-list()

      #clusterExport(cl=cl, varlist=c("real.pre"))
      resultado.pre.lista<-parLapply(cl, forecasts_pseudo, function(x) {
        100*sum(abs((real.pre-as.numeric(x))/real.pre))/length(real.pre)})

      resultado.pre<-unlist(resultado.pre.lista)

      selec=which(rank(resultado.pre)<300)


      matseries<-ts(matrix(unlist(bootstrapped_series), ncol = length(bootstrapped_series), byrow = F),frequency=frequency(y),start=start.i)

      selecao<-NULL
      selecClus<-NULL
      eucl<-diss(matseries[,selec],distance)

      #silhouette
      clusterExport(cl=cl, varlist=c("pam"))
      if (silhouette==T){
        k<-(2:100)
        teste<-parSapply(cl,k,function(x){pam(eucl, k=x) $ silinfo $ avg.width})
        k.best<-which.max(teste)+1
      }else{k.best=nclusters}

      eucl.pamclus <- pam(eucl, k = k.best)$clustering


      Nh<-NULL
      n<-100
      nh<-NULL
      Sh<-NULL
      selecClus<-NULL
      selecao3<-list()

      t<-(1:k.best)
      Nh<-as.numeric(table(eucl.pamclus))

      nh<-round(Nh/299*100)
      nh<-ifelse(nh==0,1,nh)


      for (t in 1:k.best){
        selecao2<-NULL
        teste2<-names(eucl.pamclus[eucl.pamclus==t])

        for (i in (1:length(teste2))){
          pre_sele2<-as.numeric(strsplit(teste2," ")[[i]][2])
          selecao2<-c(selecao2,pre_sele2)
        }
        selecao3[[t]]<-selecao2

      }


      for (t in 1:k.best){
        selecClus<-c(selecClus,selecao3[[t]][which(rank(resultado.pre[selecao3[[t]]],ties.method ="first")<(nh[t]+1))])
      }
      bootstrapped_series_ori<-bootstrapped_series
      bootstrapped_series<-bootstrapped_series[selecClus]
    }else{
      bootstrapped_series_ori<-bootstrapped_series
      bootstrapped_series<-bootstrapped_series[sample(1:1000,100)]
      k.best<-NA
    }
    ###########################################


    mod_boot <- parLapply(cl,bootstrapped_series, function(x) {
      mod <- ets(x,...)
    })
    out <- list()
    out$y <- as.ts(y)
    out$selec<-selec
    out$resultado.pre<-resultado.pre
    out$clusters<-selecao3
    out$bootstrapped_series <- bootstrapped_series
    out$bootstrapped_series_ori <- bootstrapped_series_ori
    out$models <- mod_boot
    out$etsargs <- list(...)
    fitted_boot <- lapply(out$models, fitted)
    fitted_boot <- as.matrix(as.data.frame(fitted_boot))
    out$fitted <- ts(apply(fitted_boot, 1, mean))
    tsp(out$fitted) <- tsp(out$y)
    out$residuals <- out$y - out$fitted
    out$series <- deparse(substitute(y))
    out$k<-k.best
    out$method <- "baggedClusterETS"
    out$call <- match.call()
    return(structure(out, class = c("baggedClusterETS")))
  }else{

    #UNIX systems (OSX, Linux)
    start.i <- tsp(y)[1]
    start.f <- tsp(y)[2] + 1/frequency(y)

    bootstrapped_series <- bld.mbb.bootstrap(y, boot_samples)

    if (is.null(h_pseudo)==T){
      h_pseudo <- ifelse(frequency(y) > 1, 2 * frequency(y), 10)
    }



    real.pre<-y[((length(y)-h_pseudo+1):length(y))]


    if(length(y)-h_pseudo>h_pseudo){
      #real.pre<-y[((length(y)-h_pseudo+1):length(y))]
      pseudo <- mclapply(bootstrapped_series,mc.cores=cores, function(x) {
        ts(x[(1:(length(y)-h_pseudo))],frequency=frequency(y),start=start.i)
      })


      forecasts_pseudo <- mclapply(pseudo,mc.cores=cores, function(x) {
        mod <- forecast(ets(x),h=h_pseudo)$mean
      })

      #  forecasts_pseudo <- mclapply(pseudo,mc.cores=cores, function(x) {
      #      mod <- hw(x,h=h_pseudo)$mean
      #    })


      resultado.pre.lista<-list()

      resultado.pre.lista<-mclapply(forecasts_pseudo,mc.cores=cores, function(x) {
        100*sum(abs((real.pre-as.numeric(x))/real.pre))/length(real.pre)})

      resultado.pre<-unlist(resultado.pre.lista)

      selec=which(rank(resultado.pre)<300)


      matseries<-ts(matrix(unlist(bootstrapped_series), ncol = length(bootstrapped_series), byrow = F),frequency=frequency(y),start=start.i)

      selecao<-NULL
      selecClus<-NULL
      eucl<-diss(matseries[,selec],distance)

      #silhouette
      if (silhouette==T){
        k<-(2:100)
        teste<-sapply(k,function(x){pam(eucl, k=x) $ silinfo $ avg.width})
        k.best<-which.max(teste)+1
      }else{k.best=nclusters}

      eucl.pamclus <- pam(eucl, k = k.best)$clustering


      Nh<-NULL
      n<-100
      nh<-NULL
      Sh<-NULL
      selecClus<-NULL
      selecao3<-list()

      t<-(1:k.best)
      Nh<-as.numeric(table(eucl.pamclus))

      nh<-round(Nh/299*100)
      nh<-ifelse(nh==0,1,nh)


      for (t in 1:k.best){
        selecao2<-NULL
        teste2<-names(eucl.pamclus[eucl.pamclus==t])

        for (i in (1:length(teste2))){
          pre_sele2<-as.numeric(strsplit(teste2," ")[[i]][2])
          selecao2<-c(selecao2,pre_sele2)
        }
        selecao3[[t]]<-selecao2

      }


      for (t in 1:k.best){
        selecClus<-c(selecClus,selecao3[[t]][which(rank(resultado.pre[selecao3[[t]]],ties.method ="first")<(nh[t]+1))])
      }
      bootstrapped_series_ori<-bootstrapped_series
      bootstrapped_series<-bootstrapped_series[selecClus]
    }else{
      bootstrapped_series_ori<-bootstrapped_series
      bootstrapped_series<-bootstrapped_series[sample(1:1000,100)]
      k.best<-NA
    }
    ###########################################


    mod_boot <- mclapply(bootstrapped_series,mc.cores=cores, function(x) {
      mod <- ets(x)
    })
    out <- list()
    out$y <- as.ts(y)
    out$selec<-selec
    out$resultado.pre<-resultado.pre
    out$clusters<-selecao3
    out$bootstrapped_series <- bootstrapped_series
    out$bootstrapped_series_ori <- bootstrapped_series_ori
    out$models <- mod_boot
    out$etsargs <- list(...)
    fitted_boot <- lapply(out$models, fitted)
    fitted_boot <- as.matrix(as.data.frame(fitted_boot))
    out$fitted <- ts(apply(fitted_boot, 1, mean))
    tsp(out$fitted) <- tsp(out$y)
    out$residuals <- out$y - out$fitted
    out$series <- deparse(substitute(y))
    out$k<-k.best
    out$method <- "baggedClusterETS"
    out$call <- match.call()
    return(structure(out, class = c("baggedClusterETS")))
  }
}


