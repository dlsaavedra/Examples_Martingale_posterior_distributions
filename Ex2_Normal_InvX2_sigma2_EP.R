# Martingalas Bayes Prediction Ejemplo 1 con empirical predictive
library(parallel)
f_y_yn_EP <- function(v_new, sigma2_new, mu, Y, N, R){
  
  n = length(Y)
  v_aux = rep(v_new, R)
  sigma2_aux = rep(sigma2_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  Y_matrix = matrix(nrow = R, ncol = N + n)
  Y_matrix[,1:n] = matrix(rep(Y, R), nrow = R, byrow = TRUE)
  
  for(j in 1:N){
    
    theta_matrix[,j]  = v_aux*sigma2_aux/(v_aux - 2)
    # Sample row indices
    sampled_rows <- sample(1:((j-1)+n), R, replace = T)
    # Extract sampled rows
    Y_matrix[, j + n] = Y_matrix[cbind(1:R, sampled_rows)]
    sigma2_aux = sigma2_aux*v_aux
    v_aux = v_aux + 1 
    sigma2_aux = (sigma2_aux + (Y_matrix[, j + n]-mu)^2)/v_aux
    
  }
  theta_matrix[,N+1]  = v_aux*sigma2_aux/(v_aux - 2)
  
  rm(Y_matrix)
  return(theta_matrix)
}


# Datos simulados-------
n = 100
v0 = 5
sigma02 = 3
mu = 0

sigma2 = 1/rgamma(1, shape = v0/2, rate = v0*sigma02/2)
Y = rnorm(n, mu, sqrt(sigma2))
v = sum(Y-mu)^2/n

v_new = v0 + n
sigma2_new = (v0*sigma02 + sum((Y-mu)^2))/(v0 + n)


N = 5000
R = 100

theta_N <- unlist(mclapply(1:100, function(x)  f_y_yn(alpha_new = v_new,
                                                      beta2_new = sigma2_new,
                                                      mu = mu,
                                                      N,R)[,N +1], 
                           mc.cores = detectCores() - 1))

theta_N_EP <- unlist(mclapply(1:100, function(x) f_y_yn_EP(v_new = v_new,
                                                           sigma2_new = sigma2_new,
                                                           mu = mu,
                                                           Y = Y,
                                                           N,R)[,N +1], 
                              mc.cores = detectCores() - 1))

#theta_N = f_y_yn(alpha_new = v_new,
#                 beta2_new = sigma2_new,
#                 mu = mu,
#                 N,R)[,N +1]
#theta_N_EP = f_y_yn_EP(v_new = v_new,
#                    sigma2_new = sigma2_new,
#                    mu = mu,
#                    Y = Y,
#                    N,R)[,N +1]

plot(density(theta_N), col = "red", xlab = "sigma2", main = "", 
     lwd = 2, ylim = c(0,.5))
lines(density(theta_N_EP), col = "blue", lwd = 2)
lines(density(1/(rgamma(100000, shape = v_new/2,
                        rate =v_new* sigma2_new/2))),col = "black", lwd = 3, lty = 2)
legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical predictive"), 
       col = c("black", "red", "blue"), cex = .5, lty = c(2,1,1))




mean(theta_N)
var(theta_N)
#2*alpha_new^2*beta2_new^2/(alpha_new-2)^2/(alpha_new-4)

acf(theta_N)
hist(theta_N)
# Traceplot ----

N_max = 10000
sec = seq(1,N_max, by = 10)
matplot(sec, t(f_y_yn(alpha_new = v_new, beta2_new = sigma2_new,
                    mu = mu, N = N_max, R=50))[sec, ], 
        type = "l" ,lty = 1,xlab = "Step i", 
        ylab = expression("Posteriori Mean " * sigma^2), main = "")

matplot(sec, t(f_y_yn_EP(v_new = v_new,sigma2_new = sigma2_new,
                    mu = mu, Y = Y, N = N_max,R=50))[sec, ], 
        type = "l" ,lty = 1, xlab = "Step i", 
        ylab = expression("Posteriori Mean " * sigma^2), main = "")



#n grande n = 1e4-------
n = 10000
v0 = 5
sigma02 = 3
mu = 0

sigma2 = 1/rgamma(1, shape = v0/2, rate = v0*sigma02/2)
Y = rnorm(n, mu, sqrt(sigma2))
v = sum(Y-mu)^2/n

v_new = v0 + n
sigma2_new = (v0*sigma02 + sum((Y-mu)^2))/(v0 + n)

N_max = 500000
sec = seq(1,N_max, by = 100)
matplot(sec, t(f_y_yn(alpha_new = v_new, beta2_new = sigma2_new,
                      mu = mu, N = N_max, R=50))[sec, ], 
        type = "l" ,lty = 1,xlab = "Step i", 
        ylab = expression("Posteriori Mean " * sigma^2), main = "")

matplot(sec, t(f_y_yn_EP(v_new = v_new,sigma2_new = sigma2_new,
                         mu = mu, Y = Y, N = N_max,R=50))[sec, ], 
        type = "l" ,lty = 1, xlab = "Step i", 
        ylab = expression("Posteriori Mean " * sigma^2), main = "")

N = 200000
R = 1000
theta_N <- unlist(mclapply(1:10, function(x)  f_y_yn(alpha_new = v_new,
                                                      beta2_new = sigma2_new,
                                                      mu = mu,
                                                      N,R)[,N +1], 
                           mc.cores = detectCores() - 1))

theta_N_EP <- unlist(mclapply(1:10, function(x) f_y_yn_EP(v_new = v_new,
                                                           sigma2_new = sigma2_new,
                                                           mu = mu,
                                                           Y = Y,
                                                           N,R)[,N +1], 
                              mc.cores = detectCores() - 1))




plot(density(theta_N), col = "red", xlab = "sigma2", main = "", 
     lwd = 2, ylim = c(0,20))
lines(density(theta_N_EP), col = "blue", lwd = 2)
lines(density(1/(rgamma(100000, shape = v_new/2,
                        rate =v_new* sigma2_new/2))),col = "black", lwd = 3, lty = 2)
legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical predictive"), 
       col = c("black", "red", "blue"), cex = .5, lty = c(2,1,1))





# Tiempo----
system.time(f_y_yn_EP(v_new = v_new,sigma2_new = sigma2_new,
                      mu = mu, Y = Y, N=N_max,R = 50))











# Comparar el tiempo de la predictiva con empirical predictiva-----
source("Ex1_Normal_Normal_sigma2.R")
# Datos simulados-------
n = 1e3
v0 = 1/2
sigma02 = 1
mu = 1

sigma2 = 1/rgamma(1, shape = v0/2, rate = v0*sigma02/2)
Y = rnorm(n, mu, sqrt(sigma2))
plot(density(Y))
alpha_new = v0 + n
beta2_new = (v0*sigma02 + sum((Y-mu)^2))/(v0 + n)
theta_n = alpha_new*beta2_new/(alpha_new - 2)


N = n+10000
R = 1000

theta_N = f_y_yn(alpha_new = alpha_new, beta2_new = beta2_new,
                 mu = mu, N,R)[,N +1]
theta_N_EP = f_y_yn_EP(alpha_new = alpha_new,beta2_new = beta2_new,
                       mu = mu, Y = Y, N,R)[,N+1]

plot(density(theta_N))
lines(density(theta_N_EP), col = "blue")
lines(density(1/(rgamma(10000, shape = alpha_new/2,
                        rate =alpha_new* beta2_new/2))),col = "red")
legend("topright",c("predictive", "Empirical_predictive", "real"), 
       col = c("black", "blue", "red"), cex = .4, lty = 1)


# Traceplot ----
matplot(t(f_y_yn(alpha_new = alpha_new,beta2_new = beta2_new,
                    mu = mu, N = 100000,R=20)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot de las filas de M")
matplot(t(f_y_yn_EP(alpha_new = alpha_new,beta2_new = beta2_new,
                    mu = mu, Y = Y, N = 100000,R=20)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot de las filas de M")

# Tiempo ----
system.time(f_y_yn(alpha_new = alpha_new,beta2_new = beta2_new,
                   mu = mu, N=10000,R = 1000))
system.time(f_y_yn_EP(v_new = alpha_new,beta2_new = beta2_new,
                      mu = mu, Y =Y, N=10000,R = 1000))
