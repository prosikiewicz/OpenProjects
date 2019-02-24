### file that used to test my GitHub account
###	2018.02.24

# generate artificial data for galaxy plot
data_x <- sort(c(rnorm(n = 1000, mean = 10, sd = 2)))
data_y <- sapply(data_x, function(x){x+rnorm(1, mean=0, sd=2)})

# two plots on one page
par(mfrow=c(1,2))

# plot 1
plot(data_x, data_y, pch=16, cex=0.4, ylim=c(0,20), xlim=c(0,20))
points(data_x, data_x, pch=16, cex=0.5, col="red")

#plot 2
# plot(data_x, c(1:1000), pch=16, cex=0.4)

# plot 3
plot(data_x, data_y, pch=16, cex=0.4, ylim=c(0,20), xlim=c(0,20), col="red")
points(data_x, data_x, pch=16, cex=0.5, col="blue")

# enjoy
