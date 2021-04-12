rm(list=ls())

# inference for infectiousness (results for Fig 1c)

#--- data ---
# package: readxl
data = data.frame(readxl::read_xlsx("temporal dynamics in viral shedding/Fig1c_data.xlsx"))
data <- read_excel("~/Google Drive/RESEARCH/COVID Colombia/replications/temporal dynamics in viral shedding/Fig1c_data.xlsx",
                   col_types = c("numeric", "date", "date", "date")) %>% 
  mutate(across(c(x.lb, x.ub, y), as.Date))

ref.date = as.Date("2020-01-01")

data$x.lb <- as.numeric(as.Date(data$x.lb)-ref.date)
data$x.ub <- as.numeric(as.Date(data$x.ub)-ref.date)
data$y <- as.numeric(as.Date(data$y)-ref.date)

# data: (x.lb, x.ub): lower and upper bounds of infectors symtpom onset dates
# y: symptom onset dates of infectee

#--- functions ---

#--- CDF of serial interval ---
p.Z  = function(lb, ub, gpar, lnpar) {
  
  length(lb)
  
  #--- infectiousness, gamma distribution ---
  # gpar[1:2]: hyper-parameters (gamma)
  # x        : infection time of infectee w.r.t onset time of infector
  f.Xc = function(x, gpar) { dgamma(x, gpar[1], gpar[2]) }
  
  #--- incubation, from Li et al NEJM 2020 ---
  # lnpar[1:2]: hyper-parameter (logNormal)
  # y         : length of incubation period of infectee
  f.Y  = function(y, lnpar) { dlnorm(y, lnpar[1], lnpar[2]) }
  
  #--- convolution between incubation and infectiousness profile ---
  # gpar[3]: shift c days before symptom onset of infector
  # z      : length of serial interval
  f.Z = function(z, gpar, lnpar) {
    integrate(
      f = function(x, z, gpar, lnpar) { f.Y(z+gpar[3]-x, lnpar)*f.Xc(x, gpar) },
      lower = -Inf, 
      upper = Inf,
      z     = z,
      gpar  = gpar,
      lnpar = lnpar
    )$value
  } 
  f.Z2 = Vectorize(f.Z, vectorize.args = "z")
  
  # print(lb)
  
  #--- p.Z ---
  integrate(
    f = function(x, gpar, lnpar) { f.Z2(x, gpar, lnpar) },
    lower = lb,
    upper = ub,
    gpar  = gpar,
    lnpar = lnpar
  )$value
}
p.Z2 = Vectorize(p.Z, vectorize.args = c("lb", "ub"))

#--- logLikelihood for the observed serial intervals ---
# x.lb: lower bound of infectors symtpom onset dates
# x.ub: upper bound of infectors symtpom onset dates
# y   : symptom onset dates of infectee
# 0.5 : continuity correction
data %>%
  mutate(lli = pmap_dbl(list(x.lb, x.ub, y),
                        .f = ~ log(p.Z(
                          lb = ..3 - (..2 + 0.5),
                          ub = ..3 - (..1 - 0.5),
                          gpar = c(10, 0.5, 20),
                          lnpar = c(ln.par1, ln.par2)
                        )))) %>%
  print_all


lli.fx = function(gpar, x.lb, x.ub, y, lnpar) {
  # print(x.ub)
  
  lli = log(p.Z2(lb = y-(x.ub+0.5), ub = y-(x.lb-0.5), gpar, lnpar))
  out <- -sum(lli)
  print(out)
  return(out)
  
}


# OTHER VERSION which I am 65% sure is the wrong one
# lli.fx = function(gpar, x.lb, x.ub, y, lnpar) {
#   lli = log(p.Z2(y-(x.lb-0.5), gpar, lnpar) - p.Z2(y-(x.ub+0.5), gpar, lnpar))
#   return(-sum(lli[!is.infinite(lli)]))
# }




# median = lambda (ln 2)^(1/k)
# >> lambda = median / (ln 2)^(1/k)

# MY VERSION
lli.fx_map <- function(gpar, x.lb, x.ub, y, lnpar) {
  
  # print(x.lb)
  
  lli <- pmap_dbl(list(x.lb, x.ub, y),
           .f = ~ log(p.Z(
             lb = ..3 - (..2 + 0.5),
             ub = ..3 - (..1 - 0.5),
             gpar = gpar,
             lnpar = c(lnpar[1], lnpar[2])
           )))
  
  return(-sum(lli))
  
}

# TEST
lli.fx_map(gpar = c(10, 0.5, 20), x.lb = data$x.lb, x.ub = data$x.ub, y = data$y, lnpar = c(ln.par1, ln.par2))


lli.fx(gpar = c(10, 0.5, 20), x.lb = data$x.lb, x.ub = data$x.ub, y = data$y, lnpar = c(ln.par1, ln.par2))
gpar_default <- c(10, 0.5, 20)

p.Z2(lb = data$y-(data$x.ub+0.5), ub = data$y-(data$x.lb-0.5), gpar = gpar_default, lnpar = c(ln.par1, ln.par2))

#--- incubation period ---
# from Li et al NEJM 2020
# lognormal mean = 5.2; 95% CI = c(4.1, 7.0)
ln.par1 = 1.434065
ln.par2 = 0.6612

#
#--- estimation ---
#


?optim

fit = optim(
  par = c(10, 0.5, 20), 
  fn = lli.fx, 
  method = "L-BFGS-B",
  lower  = 1e-3 + c(1,0,4),
  x.lb = data$x.lb, 
  x.ub = data$x.ub,
  y    = data$y,
  lnpar = c(ln.par1, ln.par2)
)


#--- result ---
inf.par = fit$par       #20.5, 1.59, 12.272
inf.par[3]                                 # shift c
(inf.par[1]-1)/inf.par[2] - inf.par[3]     # mode
pgamma(inf.par[3], inf.par[1], inf.par[2]) # proportion of pre-symptomatic transmission


rgamma(n = 1000, inf.par[1], inf.par[2])


# Estimate serial interval
#--- functions ---

#--- logLikelihood for the observed serial intervals ---
# x.lb: lower bound of infectors symtpom onset dates
# x.ub: upper bound of infectors symtpom onset dates
# y   : symptom onset dates of infectee
# 0.5 : continuity correction

# fit lognormal 
# allow for negative serial intervals for (shifted) gamma distribution 
min.serial <- min((data$y-data$x.ub))

lli.g = function(par, x.lb, x.ub, y){
  lli = log(pgamma(y-(x.lb-0.5)-min.serial, par[1], par[2])-
              pgamma(pmax(y-(x.ub+0.5)-min.serial,0.1), par[1], par[2]))
  return(-sum(lli))
}

ser.fit = optim(c(5, 1), lli.g, x.lb=data$x.lb, x.ub=data$x.ub, y=data$y)
ser.par <- ser.fit$par

# mean and median of the serial intervals
ser.par[1]/ser.par[2]+min.serial 				# mean
median(rgamma(1000000,ser.par[1], ser.par[2]))+min.serial	# median

df_randomly_drawn <- tibble(
  serial = rgamma(1000000,ser.par[1], ser.par[2]) + min.serial - 0.5
)

df_in_plot <- tibble(
  x = seq(-4, 35, 0.1),
  y = dgamma.shift(x+0.5, min.serial, ser.par[1], ser.par[2])
)


ggplot() + 
  geom_density(data = df_randomly_drawn, aes(x = serial)) + 
  geom_line(data = df_in_plot, aes(x = x, y = y), linetype = "dashed", colour = "red")
  


# % of negative serial intervals
pgamma(-min.serial, ser.par[1], ser.par[2])


# Figure 1c
dgamma.shift <- function(x,min.serial,gpar1,gpar2) dgamma(x-min.serial,gpar1,gpar2)

windows(height=12, width=6)
par(mfrow=c(3,1), mar=c(5,5,1,1))
# Estimated serial interval distribution
plot(NA, axes=F, ann=F, xlim=c(-5,20), ylim=c(0,0.15))
axis(1, at=0:13*2-5, label=0:13*2-4, pos=-0.005)
axis(2, las=1, at=0:3*0.05, lab=paste0(0:3*5,'%'))
curve(dgamma.shift(x+0.5, min.serial, ser.par[1], ser.par[2]), from=-4-0.5, to=20-0.5, add=T)
mtext('Serial interval (days)', 1, line=2.5)
mtext('Density',2, line=3)

# Inferred distribution of infectiousness
plot(NA, axes=F, ann=F, ylim=c(0,0.3), xlim=c(-10,8))
axis(2, las=1, at=0:3*0.1, lab=paste0(0:3*10,'%'))
abline(v=0, col=gray(0.8))
curve(dgamma(x+inf.par[3], inf.par[1], inf.par[2]), from=-10, to=8, add=T)
axis(1, at=(-5:4)*2, lab=(-5:4)*2, cex.axis=1)
mtext('Density', 2, line=3)
mtext('Days after symptom onset', 1, line=2.5)

# Incubation period
curve(dlnorm(x, ln.par1, ln.par2), from=0, to=14, axes=F, ann=F, ylim=c(0,0.3), xlim=c(0,14))
axis(2, las=1, at=0:3*0.1, lab=paste0(0:3*10,'%'))
axis(1, at=0:5*3, lab=0:5*3, cex.axis=1)
mtext('Density', 2, line=3)
mtext('Days from infection to symptom onset', 1, line=2.5)


#
#--- END ---
#


#--- bootstrap 95%CI ---

#n = 1100
#fit.list = lapply(
#  1:n, function(i) {
#    set.seed(i)
#    ind = sample(1:nrow(data), replace = T)
#    d = data[ind, ]
#    fit = optim(
#      c(10, 0.5, 20), lli.fx, 
#      method = "L-BFGS-B",
#      lower  = 1e-3 + c(1,0,4),
#      x.lb = d[, "x.lb"], 
#      x.ub = d[, "x.ub"],
#      y    = d[, "y"],
#      lnpar = c(ln.par1, ln.par2)
#    )
#    return(fit)
#  }
#)

