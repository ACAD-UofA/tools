# data exploration

require(ggplot2)
require(cowplot)

# Use a different logistic regression method, from the arm package.
# bayesglm() uses a weakly informative ~Cauchy prior for the coefficients.
# This is more robust to 'complete separation', which is problematic in
# some of our data due to low sample numbers in some categories.
# NOTE: some results do change c.f glm(). I.e. for bears, I(lat)^4 is significant for
# glm() but not bayesglm(). I consider that this is not biologically useful,
# and so is of little concern.
# Interestingly, the VIFs drop with bayesglm() compared with glm().
require(arm)
glm <- bayesglm

theme_set(
  theme_bw() +
    theme(plot.title=element_text(hjust=0.5))
    #theme(text=element_text(size=15), plot.title=element_text(size=16,face="bold",hjust=0.5))
    #theme(legend.position="none")
)

plt.hist <- function(d, x, y, binsize=10) {
  brk = seq(round(min(d[[x]],na.rm=T)), round(max(d[[x]],na.rm=T))+1, binsize)
  ggplot(data=d, aes_string(x, colour=y)) +
    #geom_histogram(aes(y=..density..), breaks=brk) +
    geom_density() +
    #theme(legend.position=c(.5,.85), legend.title=element_blank()) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(ncol = 1)) +
    geom_rug(data=subset(d,d[[y]]==levels(d[[y]])[1]), aes_string(x=x)) +
    geom_rug(data=subset(d,d[[y]]==levels(d[[y]])[2]), aes_string(x=x), sides="t")

}

plt.resid <- function (model) {
  ggplot(data=model, aes(x=.fitted,y=.resid)) +
    geom_point() +
    geom_smooth()
}

fit.logistic <- function(d, x, y, poly=1) {
  terms = lapply(1:poly, function(a){if (a==1) {x} else {paste("I(",x,"^",a,")",sep="")}})
  f = reformulate(termlabels=as.character(terms), response=y)
  m = glm(f, data=d, family=binomial)
  return (m)
}

plt.fit <- function(d, x, y, m) {
  # see https://www.fromthebottomoftheheap.net/2017/05/01/glm-prediction-intervals-i/

  # 100 equally spaced points for fitting to the model
  d2 = data.frame(seq(min(d[[x]],na.rm=T), max(d[[x]],na.rm=T), length=100))
  colnames(d2) = x

  # Use model prediction to get fit, StdErr, and calculate Wald CI
  d2 = cbind(d2, predict(m, d2, type="link", se.fit=T))
  link = m$family$linkinv
  d2 = transform(d2, wald.lwr=link(fit-2*se.fit), wald.upr=link(fit+2*se.fit), fit=link(fit))

  l1 = levels(d[[y]])[2]
  p = ggplot(data=d, aes_string(x=x,y=as.integer(d[[y]])-1)) +
          geom_point() +
          stat_smooth(colour="black", fill=NA) +
          geom_line(data=d2, aes_string(x,"fit"), inherit.aes=F, colour="red") +
          geom_ribbon(data=d2, aes_string(x=x, ymin="wald.lwr", ymax="wald.upr"), inherit.aes=F, fill="red",alpha=0.2) +
          ggtitle(deparse(m$formula)) +
          ylim(0,1) +
          ylab(paste("proportion",l1))

  return(p)
}

explore.continuous <- function(d, x, y, poly=1) {
  p1 = plt.hist(d, x, y)
  m = fit.logistic(d,x,y,poly)
  p2 = plt.fit(d,x,y,m)
  p = plot_grid(p1, p2, ncol=1,align='h',axis='b')
  print(p)

  # likelihood ratio test, comparing to intercept only model
  m.alt1 = update(m, as.formula(paste(y,"~ 1")), data=m$model)
  print(anova(m.alt1, m, test="LRT"))

  print(summary(m))
  return(m)
}

explore.categorical <- function(d, x, y) {
  print(table(d[[x]], d[[y]]))

  #print(chisq.test(d[[x]], d[[y]], simulate.p.value=T, B=100000))
  m = fit.logistic(d, x, y)

  d2 = data.frame(factor(levels(d[[x]])))
  colnames(d2) = x

  # Use model prediction to get fit, StdErr, and calculate Wald CI
  d2 = cbind(d2, predict(m, d2, type="link", se.fit=T))
  link = m$family$linkinv
  critval = qnorm(0.975) # 1.96
  d2 = transform(d2, wald.lwr=link(fit-critval*se.fit), wald.upr=link(fit+critval*se.fit), fit=link(fit))

  # 95% CI using a Beta() distribution with non-informative (Jeffereys) prior
  d2 = cbind(d2, unstack(data.frame(table(d[[y]], d[[x]])), Freq~Var1))
  l1 = levels(d[[y]])[2]
  l2 = levels(d[[y]])[1]
  d2 = transform(d2, beta.lwr=qbeta(0.025, get(l1)+0.5, get(l2)+0.5), beta.upr=qbeta(0.975, get(l1)+0.5, get(l2)+0.5))

  p = ggplot(d2, aes_string(x, "fit")) +
    geom_errorbar(aes(ymin=wald.lwr, ymax=wald.upr), colour="blue", width=0.4, size=1) +
    geom_errorbar(aes(ymin=beta.lwr, ymax=beta.upr), colour="red", width=0.2, size=1) +
    geom_point(size=3) +
    ggtitle(deparse(m$formula)) +
    ylim(0,1) +
    ylab(paste("proportion",l1))

  print(p)

  # likelihood ratio test, comparing to intercept only model
  print(anova(m, test="LRT"))
  print(summary(m))
  return(m)
}
