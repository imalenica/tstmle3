#' Plot estimates for the optimal treatment regime parameter
#'
#' Function to plot the estimated Phi for the optimal treatment regime,
#' and the relevant confidence intervals for each of the rules considered.
#'
#' @param res Result of \code{tstmleOPT} function.
#'
#' @export
#'

plot_OPT <- function(res) {

  #Create one dataframe
  all<-cbind.data.frame(Psi=res$tmlePsi,sd=res$tmleSD,lower=res$tmleCI$lower,
                        upper=res$tmleCI$upper)

  ggplot2::ggplot(all, aes(y = row.names(all), x = Psi, xmin = lower,
                  xmax = upper)) + ggplot2::geom_point() + ggplot2::geom_errorbarh() + ylab("Rule")
}

#' Plot for a single time-series
#'
#' Function to plot the simulated or actual time-series.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param plot create a plot of type "xts" or "pa".
#'
#' @importFrom xts plot.xts as.xts
#' @importFrom PerformanceAnalytics chart.TimeSeries
#'
#' @export
#'

plot_ts<-function(data, plot="xts"){

  orig<-data

  if(nrow(data)>100){
    ts<-data.frame(data=data[1:100,])
    row.names(ts)<-row.names(data)[1:100]
    names<-row.names(ts)
    data<-ts(ts)

  }else{
    names<-row.names(data)
    data<-ts(data)
  }

  if(plot=="xts"){
    xts::plot.xts(xts::as.xts(data), main = "Time-Series", col="black")
  }else if(plot=="pa"){

    par(mfrow=c(1, 2))
    PerformanceAnalytics::chart.TimeSeries(xts::as.xts(data), auto.grid = FALSE, main = "Time-Series", xaxis.labels=names)

    #Plot processes separately:
    Y <- orig[grep("Y", row.names(orig), value = TRUE), ]
    A <- orig[grep("A", row.names(orig), value = TRUE), ]
    W <- orig[grep("W", row.names(orig), value = TRUE), ]

    all<-ts(cbind.data.frame(Y=Y,A=A))

    PerformanceAnalytics::chart.TimeSeries(xts::as.xts(all), auto.grid = FALSE, main = "A and Y Process")
  }

}



