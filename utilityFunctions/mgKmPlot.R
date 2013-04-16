mgKmPlot <- function(dat, mg, subset=NULL){
  
  if( !is.null(subset) ){
    dat <- dat[ subset, ]
  }
  dat$mg <- dat[, mg]
  dat$mgGrp <- ifelse(dat$mg > median(dat$mg), "high", "low")
  
  coxFit <- coxph(Surv(osYears, osStatus) ~ mg, data=dat)
  logtestPval <- summary(coxFit)$logtest["pvalue"]
  fit <- survfit(Surv(osYears, osStatus) ~ factor(mgGrp), data=dat)
  
  df <- data.frame(
    time    = fit$time,
    surv    = fit$surv,
    strata  = gsub("factor(mgGrp)=", "", summary(fit, censored = T)$strata, fixed=T),
    upper   = fit$upper,
    lower   = fit$lower
  )
  zeros <- data.frame(time = 0, surv = 1, strata = gsub("factor(mgGrp)=", "", levels(summary(fit)$strata), fixed=T), upper = 1, lower = 1)
  df <- rbind(zeros, df)
  
  kmPlot <- ggplot(df, aes(time, surv, colour = strata)) +
    geom_step() +
    xlim(0, max(fit$time)) +
    ylim(0, 1) +
    labs(data=NULL, colour=mg) +
    xlab("years from biopsy") +
    ylab("survival probability") +
    geom_text(data=NULL, x=5, y=0.15, colour="black",
              label=paste("log likelihood p-value = ", round(logtestPval, 6), sep="")) +
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=12))
  
  return(kmPlot)
  
}

