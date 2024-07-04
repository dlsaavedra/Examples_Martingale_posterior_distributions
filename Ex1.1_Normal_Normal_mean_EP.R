# Martingalas Bayes Prediction  Ejemplo 1 variando mu0 y sigma02 con empirical predictive
library(parallel)
f_y_yn_EP <- function(mu_new, sigma2_new, sigma2, Y, N, R){
  
  mu_aux = rep(mu_new, R)
  sigma2_aux = rep(sigma2_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  Y_matrix = matrix(nrow = R, ncol = N + n)
  Y_matrix[,1:n] = matrix(rep(Y, R), nrow = R, byrow = TRUE)
  
  for(j in 1:N){
    
    theta_matrix[,j] = mu_aux
    # Sample row indices
    
    sampled_rows <- sample(1:((j-1)+n), R, replace = T)
    # Extract sampled rows
    Y_matrix[, j + n] = Y_matrix[cbind(1:R, sampled_rows)]
    
    mu_aux = mu_aux/sigma2_aux
    sigma2_aux = 1/(1/sigma2_aux + 1/sigma2)
    mu_aux = (mu_aux + Y_matrix[, j + n]/sigma2)*sigma2_aux
    
  }
  theta_matrix[,N + 1] = mu_aux
  return(theta_matrix)
}


# Datos simulados
n = 100
mu0 = 10
sigma02 = 2
theta = rnorm(1, mu0, sqrt(sigma02))
sigma2 = 1
Y = rnorm(n, theta, sqrt(sigma2))
theta_n = (mu0/sigma02 + sum(Y)/sigma2)/(1/sigma02+ n/sigma2)
var_n = 1/(1/sigma02+ n/sigma2)

# Correr las trazas
N = 5000
R = 100
theta_N <- unlist(mclapply(1:100, function(x)  f_y_yn(theta_n, var_n, sigma2, N,R)[,N+1], mc.cores = detectCores() - 1))
theta_N_EP <- unlist(mclapply(1:100, function(x)  f_y_yn_EP(theta_n, var_n, sigma2,Y, N,R)[,N+1], mc.cores = detectCores() - 1))
#theta_N_EP = f_y_yn_EP(theta_n, var_n, sigma2,Y, N,R)[,N+1]
#theta_N = f_y_yn(theta_n, var_n, sigma2, N,R)[,N+1]

plot(density(theta_N), col = "red", xlab = expression(theta), main = "", 
     lwd = 2, ylim = c(0,5))
lines(density(theta_N_EP), col = "blue", lwd = 2)
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "black", lwd = 3, lty = 2)

legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical predictive"), 
       col = c("black", "red", "blue"), cex = .6, lty = c(2,1,1))

# Traceplot

N_max = 5000
sec = seq(1,N_max, by = 10)
# Traceplot de la matriz M
matplot(sec, t(f_y_yn(theta_n, var_n, sigma2,N_max,R = 50))[sec, ],
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")
# Traceplot de la matriz 

matplot(sec,t(f_y_yn_EP(theta_n, var_n, sigma2,Y, N_max,R = 50))[sec,],
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")

# Datos simulados n grande 1e4 -----
n = 10000
mu0 = 10
sigma02 = 2
theta = rnorm(1, mu0, sqrt(sigma02))
sigma2 = 1
Y = rnorm(n, theta, sqrt(sigma2))
theta_n = (mu0/sigma02 + sum(Y)/sigma2)/(1/sigma02+ n/sigma2)
var_n = 1/(1/sigma02+ n/sigma2)

# Traceplot

N_max = 500000
sec = seq(1,N_max, by = 100)
# Traceplot de la matriz M
matplot(sec, t(f_y_yn(theta_n, var_n, sigma2,N_max,R = 50))[sec, ],
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")
# Traceplot de la matriz 

matplot(sec,t(f_y_yn_EP(theta_n, var_n, sigma2,Y, N_max,R = 50))[sec,],
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")

# Correr las trazas
N = 200000
R = 100
theta_N <- unlist(mclapply(1:100, function(x)  f_y_yn(theta_n, var_n, sigma2, N,R)[,N+1], mc.cores = detectCores() - 1))
theta_N_EP <- unlist(mclapply(1:100, function(x)  f_y_yn_EP(theta_n, var_n, sigma2,Y, N,R)[,N+1], mc.cores = detectCores() - 1))
#theta_N_EP = f_y_yn_EP(theta_n, var_n, sigma2,Y, N,R)[,N+1]
#theta_N = f_y_yn(theta_n, var_n, sigma2, N,R)[,N+1]

plot(density(theta_N), col = "red", xlab = expression(theta), main = "", 
     lwd = 2, ylim = c(0,45))
lines(density(theta_N_EP), col = "blue", lwd = 2)
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "black", lwd = 3, lty = 2)

legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical predictive"), 
       col = c("black", "red", "blue"), cex = .6, lty = c(2,1,1))




mean(theta_N)
theta
var(theta_N)
var_n

acf(theta_N)


# Traceplot 
matplot(t(f_y_yn_EP(theta_n, var_n, sigma2,Y, N=20000,R = 10)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot de las filas de M")

# Traceplot 
matplot(t(f_y_yn(theta_n, var_n, sigma2, N=20000,R = 10)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot de las filas de M")
