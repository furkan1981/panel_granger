# Traditional Granger Causality Test
# written by Furkan Emirmahmutoglu
# emirfurkan@gmail.com; f.emirmahmutoglu@hbv.edu.tr

GC <- function(series,lag.max = NULL, ic = c("AIC", "BIC"), type = c("constant","trend")) {
  TYGC(series, dmax = 0, lag.max, ic, type)
}
