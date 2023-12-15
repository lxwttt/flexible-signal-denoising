source("./main.R")

# experiment 1: motor data

data <- read.csv("./datasets/motor.dat", header = TRUE, sep = "\t", skip = 36) # Replace read_csv with read.table or read.delim if appropriate

times <- data$times
accel <- data$accel
J <- log2(length(accel))
n <- length(accel)
accel.est <- smash.gaus(accel)

var.est <- smash.gaus((accel - accel.est)^2)
sigma <- sqrt(var.est)

library(ggplot2)
p <- ggplot(data, aes(x = times, y = accel)) +
    geom_point() +
    geom_line(aes(y = accel.est), color = "red") +
    theme_minimal() +
    labs(x = "Time", y = "Acceleration") +
    geom_line(aes(y = accel.est + 2 * sigma), color = "blue") +
    geom_line(aes(y = accel.est - 2 * sigma), color = "blue")
print(p)

# experiment 2:
library(smashr)
data(treas)
y <- ar(treas, FALSE, 5)$res
y <- y[!is.na(y)]
x <- sort(treas[5:(length(treas) - 1)])
y <- y[order(treas[5:(length(treas) - 1)])]
detach("package:smashr", unload = TRUE)

x.mod <- unique(x)
y.mod <- 0
for (i in 1:length(x.mod)) {
    y.mod[i] <- median(y[x == x.mod[i]])
}
y.exp <- c(y.mod, y.mod[length(y.mod):(2 * length(y.mod) - 2^10 + 1)])
y.final <- c(y.exp, y.exp[length(y.exp):1])
y.est <- smash.gaus(y.final)

par(mfrow = c(2, 2))
plot(treas, xlab = "year", ylab = "interest rate", type = "l")
plot(x, y)
lines(x.mod, y.est[1:(length(y.mod))], xlab = "X", ylab = "Y", col = 2)
plot(x.mod, y.est[1:(length(y.mod))],
    xlab = "X", ylab = "Y", type = "l",
    ylim = c(-0.3, 0.3)
)
