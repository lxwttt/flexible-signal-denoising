library(ggplot2)
source("./main.R")

# self exp 1: spike data

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
x <- rnorm(n, mu.t, sigma.t)

if (!ispowerof2(length(x))) {
    reflect.res <- reflect(x)
    idx <- reflect.res$idx
    x <- reflect.res$x
} else {
    idx <- 1:length(x)
}
J <- log2(length(x))
n <- length(x)

ashparam <- list(
    pointmass = TRUE, prior = "nullbiased",
    gridmult = 2, mixcompdist = "normal",
    nullweight = 10, outputlevel = 2, fixg = FALSE
)


tsum <- sum(x)
x.w.d <- sdtable(x)$difftable


var.est1.ini <- (x - lshift(x))^2 / 2
var.est2.ini <- (rshift(x) - x)^2 / 2
var.est.ini <- (var.est1.ini + var.est2.ini) / 2
mu.est <- smash.mu(
    x.w.d, var.est.ini, tsum, ashparam, J, n
)

data <- data.frame(x = t, y = x[idx], mu.est = mu.est[idx])
p <- ggplot(data, aes(x = t, y = x)) +
    # geom_point(shape = 21, fill = "white", color = "black", size = 1) +
    geom_line(aes(y = mu.est), color = "blue", linewidth = 1) +
    geom_line(aes(y = mu.t), color = "red", linewidth = 1) +
    theme_minimal()
# print(p)

var.est <- (x - mu.est)^2
var.var.est <- 2 / 3 * var.est^2
var.est <- smash.var(
    var.est, var.var.est, var.est.ini,
    ashparam, 0.5, J, n
)
var.est[var.est <= 0] <- 1e-08
sigma <- sqrt(var.est)


ashparam.mean <- ashparam
ashparam.mean$gridmult <- 64
mu.est <- smash.mu(
    x.w.d, sigma^2, tsum, ashparam.mean, J, n
)

data <- data.frame(x = t, y = x[idx], mu.est = mu.est[idx])
p <- ggplot(data, aes(x = t, y = x)) +
    geom_line(aes(y = mu.est), color = "blue", linewidth = 1) +
    geom_line(aes(y = mu.t), color = "red", linewidth = 1) +
    theme_minimal()
# print(p)

var.est <- (x - mu.est)^2
var.var.est <- 2 / 3 * var.est^2
var.est2 <- smash.var(
    var.est, var.var.est, var.est.ini,
    ashparam, 0.5, J, n
)
var.est[var.est <= 0] <- 1e-08
sigma <- sqrt(var.est)


ashparam.mean <- ashparam
ashparam.mean$gridmult <- 64
mu.est <- smash.mu(
    x.w.d, sigma^2, tsum, ashparam.mean, J, n
)

data <- data.frame(x = t, y = x[idx], mu.est = mu.est[idx])
p <- ggplot(data, aes(x = t, y = x)) +
    geom_line(aes(y = mu.est), color = "blue", linewidth = 1) +
    geom_line(aes(y = mu.t), color = "red", linewidth = 1) +
    theme_minimal()
print(p)

# self exp 2: sin(1/x^2) data
n <- 2^10 + 1
t <- 1:n / n
mu.t <- sin(1 / t^2)
x <- rnorm(n, mu.t, 1)
mu.est <- smash.gaus(x)

data <- data.frame(x = t, y = x, mu.est = mu.est)
p <- ggplot(data, aes(x = t, y = x)) +
    geom_line(aes(y = mu.t), color = "red", linewidth = 1) +
    geom_line(aes(y = mu.est), color = "blue", linewidth = 1) +
    theme_minimal()
print(p)
