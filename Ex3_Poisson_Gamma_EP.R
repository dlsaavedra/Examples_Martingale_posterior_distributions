# Ejemplo 3 Modelo Poisson- Gamma con empirical predictive----
library(parallel)

f_y_yn_EP <- function(alpha_new, beta_new, Y ,N, R){
  
  n = length(Y)
  alpha_aux = rep(alpha_new, R)
  beta_aux = rep(beta_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  Y_matrix = matrix(nrow = R, ncol = N + n)
  Y_matrix[,1:n] = matrix(rep(Y, R), nrow = R, byrow = TRUE)
  
  for(j in 1:N){
    
    theta_matrix[,j]  = alpha_aux/beta_aux
    # Sample row indices
    sampled_rows <- sample(1:((j-1)+n), R, replace = T)
    # Extract sampled rows
    Y_matrix[, j + n] = Y_matrix[cbind(1:R, sampled_rows)]

    alpha_aux = alpha_aux + Y_matrix[, j + n] 
    beta_aux = beta_aux + 1
  }
  
  theta_matrix[,N+1]  = alpha_aux/beta_aux
  rm(Y_matrix)
  return(theta_matrix)
}


n = 100
alpha = 2
beta = 5
theta = rgamma(1, alpha, rate = beta)
Y = rpois(n, theta)
alpha_new = alpha + sum(Y)
beta_new = beta + n
theta_n = alpha_new/beta_new
#curve(dgamma(x, alpha_new, shape = beta_new ), xlim = c(0,2))



N = 5000
R = 100

theta_N <- unlist(mclapply(1:100, function(x) f_y_yn(alpha_new = alpha_new,
                                                        beta_new = beta_new, 
                                                        N,R)[,N], 
                              mc.cores = detectCores() - 1))

theta_N_EP <- unlist(mclapply(1:100, function(x)  f_y_yn_EP(alpha_new = alpha_new, beta_new = beta_new, 
                                                         Y=Y, N,R)[,N], 
                           mc.cores = detectCores() - 1))

#theta_N_EP = f_y_yn_EP(alpha_new = alpha_new, beta_new = beta_new, 
#                 Y=Y, N,R)[,N+1]
#theta_N = f_y_yn(alpha_new = alpha_new,beta_new = beta_new,
#                 N,R)[,N]



plot(density(theta_N), col = "red", xlab = expression(theta), main = "", 
     lwd = 2, ylim = c(0,6.5))
lines(density(theta_N_EP), col = "blue", lwd = 2)
curve(dgamma(x,alpha_new, rate = beta_new), add = T, col = "black", lwd = 3, lty = 2)
legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical predictive"), 
       col = c("black", "red", "blue"), cex = .6, lty = c(2,1,1))



# Traceplot ----

N_max = 5000
sec = seq(1,N_max, by = 10)

matplot(sec, t(f_y_yn(alpha_new = alpha_new,
                      beta_new = beta_new, 
                      N_max,R=50))[sec, ], 
        type = "l" ,lty = 1, xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")

matplot(sec, t(f_y_yn_EP(alpha_new = alpha_new, beta_new = beta_new, 
                          Y=Y, N_max,R = 50))[sec, ], 
        type = "l" ,lty = 1,xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")




# Muestra grande n = 1e04 ------

n = 10000
alpha = 2
beta = 5
theta = rgamma(1, alpha, rate = beta)
Y = rpois(n, theta)
alpha_new = alpha + sum(Y)
beta_new = beta + n
theta_n = alpha_new/beta_new


N_max = 500000
sec = seq(1,N_max, by = 100)

matplot(sec, t(f_y_yn(alpha_new = alpha_new,
                      beta_new = beta_new, 
                      N_max,R=50))[sec, ], 
        type = "l" ,lty = 1, xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")

matplot(sec, t(f_y_yn_EP(alpha_new = alpha_new, beta_new = beta_new, 
                         Y=Y, N_max,R = 50))[sec, ], 
        type = "l" ,lty = 1,xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")


N = 200000
R = 100

theta_N <- unlist(mclapply(1:100, function(x) f_y_yn(alpha_new = alpha_new,
                                                     beta_new = beta_new, 
                                                     N,R)[,N], 
                           mc.cores = detectCores() - 1))

theta_N_EP <- unlist(mclapply(1:100, function(x)  f_y_yn_EP(alpha_new = alpha_new, beta_new = beta_new, 
                                                            Y=Y, N,R)[,N], 
                              mc.cores = detectCores() - 1))

#theta_N_EP = f_y_yn_EP(alpha_new = alpha_new, beta_new = beta_new, 
#                 Y=Y, N,R)[,N+1]
#theta_N = f_y_yn(alpha_new = alpha_new,beta_new = beta_new,
#                 N,R)[,N]



plot(density(theta_N), col = "red", xlab = expression(theta), main = "", 
     lwd = 2, ylim = c(0,6.5))
lines(density(theta_N_EP), col = "blue", lwd = 2)
curve(dgamma(x,alpha_new, rate = beta_new), add = T, col = "black", lwd = 3, lty = 2)
legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical predictive"), 
       col = c("black", "red", "blue"), cex = .6, lty = c(2,1,1))







mean(theta_N)
theta_n
var(theta_N)
alpha_new/beta_new^2

acf(theta_N)

# Traceplot 
matplot(t(f_y_yn_EP(alpha_new = alpha_new, beta_new = beta_new,
                 Y= Y, N=50000,R = 10)), 
        type = "l" ,lty = 1, xlab = "Índicet", 
        ylab = "Valor", main = "Traceplot de las filas de M")

