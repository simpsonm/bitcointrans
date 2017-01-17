# install.packages("devtools")
# devtools::install_github("RcppCore/Rcpp")
# devtools::install_github("rstats-db/DBI")
# devtools::install_github("rstats-db/RPostgres")


### COMMENTED OUT TO PREVENT OVERWRITING DATA
# library(DBI)
# con <- dbConnect(RPostgres::Postgres(), dbname = "toshi", host = "toshi.cn6zzwcfsto5.us-east-1.rds.amazonaws.com", 
	# port = 5432, user = "readonly", password = "password")

# # LOAD ALL THE BLOCK HEADERS FROM THE MAIN CHAIN
# res <- dbSendQuery(con, "select * from blocks where branch=0")
# while (!dbHasCompleted(res)) {
	# chunk <- dbFetch(res, n = 5)
	# print(nrow(chunk))
# }
# blocks <- dbFetch(res)
# dbClearResult(res)
# dbDisconnect(con)


# # library(bit64)
# # unsigned <- function(x) {
	# # if (x >= 0) {
		# # return(as.integer64(x))
	# # } else {
		# # return(as.integer64(-x + 2147483648))
	# # }
# # }

# # CLEANING UP AND SAVING THE DATA
# blocks$created_at <- NULL
# blocks$id <- NULL
# blocks <- blocks[order(blocks$height), ]
# # foo <- lapply(blocks$fees, unsigned)
# save(blocks, file = "~/Documents/Research/Block Size/bitcointrans/data/blocks.Rda")
# blocks$work <- vapply(blocks$work, paste, collapse = ", ", character(1L))
# write.table(blocks, file = "~/Documents/Research/Block Size/bitcointrans/data/blocks.csv", sep = ",", col.names = NA, 
	# qmethod = "double")


load("~/Documents/Research/Block Size/bitcointrans/data/blocks.Rda")


# PLOT BLOCK SIZES VS BLOCK HEIGHT
plot(blocks$height, blocks$size, main = "Bitcoin block sizes over time", xlab = "Block Height ", ylab = "Size in Bytes ", 
	pch = ".")

# LET'S VISUALIZE THE SOFT LIMITS IN USE
d <- density(blocks$size, from = 0, to = 1e+06)
plot(d)

# LET'S SEE HOW DENSITY CHANGES OVER TIME
d1 <- density(blocks$size[blocks$height <= 1e+05], from = 0, to = 1e+06)
d2 <- density(blocks$size[blocks$height <= 2e+05 & blocks$height > 1e+05], from = 0, to = 1e+06)
d3 <- density(blocks$size[blocks$height <= 3e+05 & blocks$height > 2e+05], from = 0, to = 1e+06)
d4 <- density(blocks$size[blocks$height > 3e+05], from = 0, to = 1e+06)
png(file = "~/Documents/Research/Block Size/bitcointrans/charts/density100k.png")
plot(d1)
dev.off()
png(file = "~/Documents/Research/Block Size/bitcointrans/charts/density200k.png")
plot(d2)
dev.off()
png(file = "~/Documents/Research/Block Size/bitcointrans/charts/density300k.png")
plot(d3)
dev.off()
png(file = "~/Documents/Research/Block Size/bitcointrans/charts/density400k.png")
plot(d4)
dev.off()

# OUT OF CURIOSITY, LET'S EXPLORE FEES
# plot(blocks$height, foo, main = "Bitcoin transaction fees per block over time", xlab = "Block Height ", 
	# ylab = "Fees in Satoshis ", pch = ".")
# we have a problem because fees should be an unsigned integer

# LET'S PLOT MEAN vs. VARIANCE FOR BLOCK SIZES
maxheight <- length(blocks$height)
# overwriting
maxheight <- 358000

id <- integer()
mean <- numeric()
var <- numeric()
sd <- numeric()
summary <- data.frame(c(id, mean, var, sd))

sequence <- seq(1,maxheight,by=500)
for (i in 1:as.integer(length(blocks$height)/500)) {
	mean <- mean(blocks$size[blocks$height >= sequence[i] & blocks$height < sequence[i] + 500])
	var <- var(blocks$size[blocks$height >= sequence[i] & blocks$height < sequence[i] + 500])
	sd <- sd(blocks$size[blocks$height >= sequence[i] & blocks$height < sequence[i] + 500])
	summary <- rbind(summary,c(i,mean,var,sd))
}
colnames(summary) <- c("id", "mean", "var", "sd")
plot(summary$mean, summary$sd)

varreg <- lm(mean ~ var + I(var^2),summary)
sdreg <- lm(mean ~ sd,summary)
bothreg <-lm(mean ~ sd + var + I(var^2),summary)
summary(varreg)
summary(sdreg)
summary(bothreg)

plot(summary$mean, summary$var)
meanreg <- lm(var ~ mean + I(mean^2) + I(mean^3) - 1, data = summary)



weights <- diff(blocks$time)
weights[weights < 0] <- 1
sizes <- blocks$size[-1]
heights <- blocks$height[-1]

library(Hmisc)

sequence <- seq(1,maxheight,by=500)
weightedsummary <- NULL
for(i in 1:length(sequence)){
  tempsizes <- sizes[heights >= sequence[i] & heights < sequence[i] + 500]
  tempweights <- weights[heights >= sequence[i] & heights < sequence[i] + 500]
  tempmean <- wtd.mean(tempsizes, tempweights, TRUE)
  tempvar <- wtd.var(tempsizes, tempweights, TRUE)
  tempsd <- sqrt(tempvar)
  totalweight <- sum(tempweights)
  weightedsummary <- rbind(weightedsummary,c(i,tempmean,tempvar,tempsd,totalweight))
}
colnames(weightedsummary) <- c("id", "mean", "var", "sd", "wt")
weightedsummary <- data.frame(weightedsummary)

plot(weightedsummary$mean, weightedsummary$var)

weightedreg1 <- lm(var ~ mean - 1, data = weightedsummary, weights = wt)
weightedreg2 <- lm(var ~ mean + I(mean^2) - 1, data = weightedsummary, weights = wt)
weightedreg3 <- lm(var ~ mean + I(mean^2) + I(mean^3) - 1, data = weightedsummary, weights = wt)
weightedreg4 <- lm(var ~ mean + I(mean^2) + I(mean^3) + I(mean^4) - 1, data = weightedsummary, weights = wt)

summary(weightedreg1)
summary(weightedreg2)
summary(weightedreg3)
summary(weightedreg4)


### added for checking basic distributional assumptions about transaction sizes
### upshot: gamma distribution seems reasonable (variance of trans size is linear in mean)
weights <- blocks$transactions_count
sizes <- blocks$size
means <- weights/sizes
heights <- blocks$height
maxheight <- max(heights)

sequence <- seq(1,maxheight,by=500)
weightedsummary <- NULL
for(i in 1:length(sequence)){
  tempsizes <- means[heights >= sequence[i] & heights < sequence[i] + 500]
  tempweights <- weights[heights >= sequence[i] & heights < sequence[i] + 500]
  tempmean <- wtd.mean(tempsizes, tempweights, TRUE)
  tempvar <- wtd.var(tempsizes, tempweights, TRUE)
  tempsd <- sqrt(tempvar)
  totalweight <- sum(tempweights)
  weightedsummary <- rbind(weightedsummary,c(i,tempmean,tempvar,tempsd,totalweight))
}
colnames(weightedsummary) <- c("id", "mean", "var", "sd", "wt")
weightedsummary <- data.frame(weightedsummary)

plot(weightedsummary$mean, weightedsummary$var)

weightedreg1 <- lm(var ~ mean, data = weightedsummary, weights = wt)
weightedreg2 <- lm(var ~ mean + I(mean^2), data = weightedsummary, weights = wt)
weightedreg3 <- lm(var ~ mean + I(mean^2) + I(mean^3), data = weightedsummary, weights = wt)
weightedreg4 <- lm(var ~ mean + I(mean^2) + I(mean^3) + I(mean^4), data = weightedsummary, weights = wt)

summary(weightedreg1)
summary(weightedreg2)
summary(weightedreg3)
summary(weightedreg4)
