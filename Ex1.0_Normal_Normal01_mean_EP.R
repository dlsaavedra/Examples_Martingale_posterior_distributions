# Martingalas Bayes Prediction Ejemplo 1 con empirical predictive
library(parallel)
f_y_yn_EP <- function(theta_n, Y, N, R){
  
  n = length(Y)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  theta_matrix[,1] = theta_n
  Y_matrix = matrix(nrow = R, ncol = N + n)
  Y_matrix[,1:n] = matrix(rep(Y, R), nrow = R, byrow = TRUE)
  
  for(j in 1:N){
    # Sample row indices
    sampled_rows <- sample(1:((j-1)+n), R, replace = T)
    # Extract sampled rows
    Y_matrix[, j + n] = Y_matrix[cbind(1:R, sampled_rows)]
    theta_matrix[,j + 1] = (theta_matrix[,j] * (n+j) + Y_matrix[, j + n])/(n+j +1)
    
  }
  rm(Y_matrix)
  return(theta_matrix)
}

f_y_yn_EP_prob <- function(theta_n, Y, N, R){
  
  n = length(Y)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  theta_matrix[,1] = theta_n
  prob_matrix = matrix(1,nrow = R, ncol = n)
  
  for(j in 1:N){
    # Sample row indices
    
    sampled_rows <- sapply(1:R, function(i) {
      sample(1:n, 1, prob = prob_matrix[i, ])
    })
    
    prob_matrix[cbind(1:R, sampled_rows)] = prob_matrix[ cbind(1:R, sampled_rows)]  + 1
    # Extract sampled rows
    Y_new = Y[sampled_rows]
    theta_matrix[,j + 1] = (theta_matrix[,j] * (n+j) + Y_new)/(n+j +1)
    
  }
  return(theta_matrix)
}

f_y_yn_EP_recur <- function(theta_n, Y, N, R){
  
  n = length(Y)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  theta_matrix[,1] = theta_n
  sampled_rows = matrix(nrow = R, ncol = N + n)
  sampled_rows[,1:n] = matrix(rep(1:n, R), nrow = R, byrow = TRUE)
  
  for(j in (n+1):(N + n)){
    # Sample row indices
    sampled_rows[,j] <- sample(1:((j-1)), R, replace = T)
    sampled_rows[, j] = sampled_rows[cbind(1:R, sampled_rows[,j])]
    # Extract sampled rows
    theta_matrix[, j - n + 1] = (theta_matrix[,j - n] * (j) + Y[sampled_rows[, j]])/(j +1)
  }
  
  
  
    
  return(theta_matrix)
}


# Tiempo ----
system.time(g())
system.time(f())

# Datos simulados-------
n = 10
theta = rnorm(1, 0, 1)
Y = rnorm(n, theta, 1)
S_y = sum(Y)
theta_n = S_y/(n+1)
var_n = 1/(n+1)

# Correr las trazas-----
N = 5000
R = 1000
theta_N = f_y_yn_EP(theta_n, Y, N,R)[,N]
plot(density(theta_N))
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "red")

mean(theta_N)
theta_n
var(theta_N)
var_n

acf(theta_N)
hist(Y)
# Traceplot ----
matplot(t(f_y_yn_EP(theta_n, Y, N=10000,R = 20)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot de las filas de M")
# Tiempo----
system.time(f_y_yn_EP(theta_n, Y, N=2000,R = 10000))

# Comparar el tiempo de la predictiva con empirical predictiva-----
#source("Ex1.0_Normal_Normal01_mean.R")
n = 1e03
theta = rnorm(1, 0, 1)
Y = rnorm(n, theta, 1)
S_y = sum(Y)
theta_n = S_y/(n+1)
var_n = 1/(n+1)


N = 5000
R = 1000

theta_N = f_y_yn(theta_n, n, N,R)[,N]
theta_N_EP = f_y_yn_EP(theta_n, Y, N, R)[,N]
theta_N_EP_recur = f_y_yn_EP_recur(theta_n, Y, N,R)[,N]

plot(density(theta_N), ylim = c(0, 20), col = "red", lwd = 2, xlab = "Theta", main = "Normal - Normal n = 1e04, N = 6000")
lines(density(theta_N_EP), col = "blue" , lwd = 2)
lines(density(theta_N_EP_recur), col = "purple" , lwd = 2)
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "black", lwd = 3)
legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical_predictive"), 
       col = c("black", "blue", "red"), cex = .5, lty = 1)



# Traceplot ----
matplot(t(f_y_yn_EP_recur(theta_n, Y, N=20000,R = 5)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot EP n = 1e04")
matplot(t(f_y_yn_EP(theta_n, Y, N=20000,R = 5)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot EP n = 1e04")
matplot(t(f_y_yn(theta_n, n, N=20000,R = 5)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot Predictive n = 1e04")

# Tiempo ----
system.time(f_y_yn(theta_n, n, N=50000,R = 1000))
system.time(f_y_yn_EP(theta_n, Y, N=500000,R = 100))
system.time(f_y_yn_EP_recur(theta_n, Y, N=500000,R = 100))

# Graficos informe -----

# Datos simulados
n = 100
theta = rnorm(1, 0, 1)
Y = rnorm(n, theta, 1)
S_y = sum(Y)
theta_n = S_y/(n+1)
var_n = 1/(n+1)
# Correr las trazas
N1 = 500
N2 = 1000
N3 = 5000
R = 10000
theta_N1 = f_y_yn_EP(theta_n, Y, N1,R)[,N1]
theta_N2 = f_y_yn_EP(theta_n, Y, N2,R)[,N2]
theta_N3 = f_y_yn_EP(theta_n, Y, N3,R)[,N3]
# comparación distintos N ----
plot(density(theta_N1), xlab = expression(theta), col = "red", main ="", lwd = 2)
lines(density(theta_N2), col = "blue", lwd = 2)
lines(density(theta_N3), col = "green", lwd = 2)
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "black", lwd = 3,lty = 2)
legend("topright",c( "Posterior", "Martingale Empirical_predictive N = 500",
                     "Martingale Empirical_predictive N = 1000", "Martingale Empirical_predictive N = 5000"), 
       col = c("black", "red", "blue","green"), cex = .5, lty = c(2,1,1,1))

# Traceplot de la matriz M
matplot(t(f_y_yn_EP(theta_n, Y, 6000,R = 50)),
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")


## Comparación predictiva con empirical predictiva
N = 5000
theta_N = f_y_yn(theta_n, n, N,R)[,N]
theta_N_EP = f_y_yn_EP(theta_n, Y, N, R)[,N]
#theta_N_EP_recur = f_y_yn_EP_recur(theta_n, Y, N,R)[,N]

plot(density(theta_N),ylim = c(0,4.5) ,
     col = "red", lwd = 2, xlab = expression(theta), main = "")
lines(density(theta_N_EP), col = "blue" , lwd = 2)
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "black", lwd = 3, lty = 2)
legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical_predictive"), 
       col = c("black", "red","blue"), lty = c(2,1,1),cex = .6)


# Datos simulados n = 1e04
n = 10000
theta = rnorm(1, 0, 1)
Y = rnorm(n, theta, 1)
S_y = sum(Y)
theta_n = S_y/(n+1)
var_n = 1/(n+1)

N_max = 200000
sec = seq(1,N_max, by = 20)
# Traceplot de la matriz M
matplot(sec, t(f_y_yn(theta_n, n, N_max,R = 50))[sec, ],
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")
# Traceplot de la matriz 

matplot(sec,t(f_y_yn_EP(theta_n, Y, N_max,R = 50))[sec,],
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")

## Comparación predictiva con empirical predictiva en paralelo
# Cargar el paquete parallel


R = 100
N = 15000

theta_N <- unlist(mclapply(1:100, function(x)  f_y_yn(theta_n, n, N,R)[,N], mc.cores = detectCores() - 1))
theta_N_EP <- unlist(mclapply(1:100, function(x)  f_y_yn_EP(theta_n, Y, N,R)[,N], mc.cores = detectCores() - 1))
#theta_N = as.vector(replicate(10, f_y_yn(theta_n, n, N,R)[,N]))
#theta_N_EP = as.vector(replicate(100,f_y_yn_EP(theta_n, Y, N, R)[,N]))


plot(density(theta_N),ylim = c(0,55) ,
     col = "red", lwd = 2, xlab = expression(theta), main = "")
lines(density(theta_N_EP), col = "blue" , lwd = 2)
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "black", lwd = 3, lty = 2)
legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical_predictive"), 
       col = c("black", "red","blue"), lty = c(2,1,1),cex = .58)


