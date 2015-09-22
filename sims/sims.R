library(xtable)

## sample size
n <- 1000000

## rate of block discovery
lambda <- 1/10

## means and variances of lognormal distribution
mus <- c(0.035, 0.04, 0.045, 0.05)
sigs <- c(0,10^(1:7-4))

out <- matrix(0,nrow=length(mus),ncol=length(sigs))
colnames(out) <- paste("$\\sigma = ", sigs, "$", sep="")
rownames(out) <- paste("$\\mu = ", mus, "$", sep="")
for(i in 1:length(mus)){
  for(j in 1:length(sigs)){
    mu <- mus[i]
    sig <- sigs[j]
    mbs <- mu*arrive + sig*sqrt(arrive)*rnorm(n)
    out[i,j] <- mean(mbs > 1)
  }
}


xout <- xtable(out, digits=4)
print.xtable(xout, digits=4, sanitize.text.function=function(x){x})
