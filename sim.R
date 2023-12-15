library(ggplot2)
source("./main.R")

# sim1: spike data

n <- 2^10 + 1
t <- 1:n / n
spike.f <- function(x) {
    (0.75 * exp(-500 * (x - 0.23)^2) +
        1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) +
        2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
}
mu.s <- spike.f(t)
mu.t <- (1 + mu.s) / 5
var.fn <- (0.0001 + 4 * (exp(-550 * (t - 0.2)^2) + exp(-200 * (t - 0.5)^2) +
    exp(-950 * (t - 0.8)^2))) / 1.35
rsnr <- sqrt(5)
sigma.t <- sqrt(var.fn) / mean(sqrt(var.fn)) * sd(mu.t) / rsnr^2
X.s <- rnorm(n, mu.t, sigma.t)
mu.est <- smash.gaus(X.s)

data <- data.frame(x = t, y = X.s, mu.est = mu.est)
p <- ggplot(data, aes(x = t, y = X.s)) +
    geom_point(shape = 21, fill = "white", color = "black", size = 1) +
    geom_line(aes(y = mu.est), color = "blue", linewidth = 1) +
    geom_line(aes(y = mu.t), color = "red", linewidth = 1) +
    theme_minimal()
print(p)

# sim2: "Corner" mean signal with "Doppler" variance

n <- 2^10
t <- 1:n / n

corner_mean <- ifelse(t < 0.5, 2 * t, 2 * (1 - t))
doppler_variance <- 0.1 + (sin(2 * pi / ((t + 0.4)^2)))^2

# Create a data frame to store the signals
X.s <- rnorm(n, corner_mean, doppler_variance)
mu.est <- smash.gaus(X.s)

data <- data.frame(x = t, y = X.s, mu.est = mu.est, mu = corner_mean)
p <- ggplot(data, aes(x = t, y = X.s)) +
    geom_point(shape = 21, color = "black", size = 1) +
    geom_line(aes(x = t, y = mu), color = "blue") +
    geom_line(aes(x = t, y = mu.est), color = "red") +
    theme_minimal()
print(p)

var.est <- smash.gaus((X.s - mu.est)^2)
sigma <- sqrt(var.est)
# plot the doppler variance and estimated variance
data <- data.frame(x = t, y = doppler_variance, sigma = sigma)
p <- ggplot(data, aes(x = t, y = doppler_variance)) +
    geom_line(color = "blue") +
    geom_line(aes(y = sigma), color = "red") +
    theme_minimal()
print(p)
