require(survival)
require(cobs)
require(latex2exp)

#' fit baseline cumulative hazard
fit_smooth_baseline = function(model, 
                               nknots = 20,
                               start_time = 0,
                               lambda = 0,
                               toler.kn = 0,
                               verbose = FALSE){
  # get baseline hazard Breslow
  bh = basehaz(model, centered = FALSE)
  t <- c(start_time,bh$time)
  Lambda0 <- c(0,bh$hazard)
  
  # Smooth version of Lambda0
  Lambda0S <- cobs(t,
                   Lambda0,
                   constraint=c("increase"),
                   pointwise=matrix(c(0,start_time,0),nrow=1),
                   nknots=nknots,
                   lambda=lambda,
                   toler.kn=toler.kn,
                   print.mesg = verbose)
  return(Lambda0S)
}


#' plot baseline hazard (breslow and smooth estimate)
plot_baseline_hazard = function(model,Lambda0S, start_time = 0,xlegend = NA, ylegend = NA){
    if (is.na(xlegend)){xlegend = 0.75*(max(Lambda0S$x) - start_time)}
    if (is.na(ylegend)){ylegend = 0.25*max(Lambda0S$y)}
    
    # get baseline hazard Breslow
    bh = basehaz(model, centered = FALSE)
    t <- c(0,bh$time)
    Lambda0 <- c(start_time,bh$hazard)
    
    # Plot comparison between basic and smooth estimate
    plot(Lambda0S$x,Lambda0S$fitted,type="l",
         main="",ylab="Baseline cumulative hazard",xlab="time", col = 'red')
    points(t,Lambda0,type="s",lty=2, col = 'blue')
    legend(xlegend, ylegend,
           legend=c(TeX("$\\tilde{\\Lambda}_0 (Smoothed)"), TeX("$\\hat{\\Lambda}_0 (Breslow)")),
           col=c("red", "blue"), lty=1:2, cex=0.7,y.intersp=3,box.lty=0)
}
