# Ejemplo 4 Modelo Beta-Binomial con empirical predictive  ----
library(parallel)


f_y_yn_EP <- function(m,alpha_new, beta_new,Y , N, R){

  n = length(Y)
  alpha_aux = rep(alpha_new, R)
  beta_aux = rep(beta_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  Y_matrix = matrix(nrow = R, ncol = N + n)
  Y_matrix[,1:n] = matrix(rep(Y, R), nrow = R, byrow = TRUE)
  
  for(j in 1:N){
    
    theta_matrix[,j]  = alpha_aux/(alpha_aux + beta_aux)
    # Sample row indices
    sampled_rows <- sample(1:((j-1)+n), R, replace = T)
    # Extract sampled rows
    Y_matrix[, j + n] = Y_matrix[cbind(1:R, sampled_rows)]

    alpha_aux = alpha_aux + Y_matrix[, j + n] 
    beta_aux = beta_aux + m - Y_matrix[, j + n]
    
  }
  theta_matrix[,N + 1]  = alpha_aux/(alpha_aux + beta_aux)
  
  rm(Y_matrix)
  return(theta_matrix)
}

# Crear Datos Simulados
n = 100

m = 10
alpha0 = 2
beta0 = 3
theta = rbeta(1, alpha0,beta0)
Y = rbinom(n,m, theta)

alpha_new = alpha0 + sum(Y)
beta_new = n*m - sum(Y) + beta0
theta_n = alpha_new/(alpha_new + beta_new)
#curve(dgamma(x, alpha_new, shape = beta_new ), xlim = c(0,2))

N = 5000
R = 100

theta_N <- unlist(mclapply(1:100, function(x) f_y_yn(alpha_new = alpha_new,
                                                     beta_new = beta_new,
                                                     m = m,
                                                     N,R)[,N], 
                           mc.cores = detectCores() - 1))

theta_N_EP <- unlist(mclapply(1:100, function(x)  f_y_yn_EP(alpha_new = alpha_new,
                                                            beta_new = beta_new,
                                                            m = m, Y = Y,
                                                            N,R)[,N], 
                              mc.cores = detectCores() - 1))

#theta_N = f_y_yn(alpha_new = alpha_new,
#                 beta_new = beta_new,
#                 m = m,
#                 N,R)[,N]
#theta_N_EP = f_y_yn_EP(alpha_new = alpha_new,
#                 beta_new = beta_new,
#                 m = m, Y = Y,
#                 N,R)[,N]

#plot(density(theta_N))
#curve(dbeta(x,alpha_new,beta_new), add = T, col = "red")

plot(density(theta_N), col = "red", xlab = expression(theta), main = "", 
     lwd = 2, ylim = c(0,30))
lines(density(theta_N_EP), col = "blue", lwd = 2)
curve(dbeta(x,alpha_new,beta_new), add = T, col = "black", lwd = 3, lty = 2)
legend("topright",c( "Posterior","Martingale Predictive", "Martingale Empirical predictive"), 
       col = c("black", "red", "blue"), cex = .7, lty = c(2,1,1))

# Traceplot ----

N_max = 5000
sec = seq(1,N_max, by = 10)

matplot(sec, t(f_y_yn(alpha_new = alpha_new,
                      beta_new = beta_new,
                      m = m,
                      N,R))[sec, ], 
        type = "l" ,lty = 1,xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")
matplot(sec, t(f_y_yn_EP(alpha_new = alpha_new,
                         beta_new = beta_new,
                         m = m, Y = Y,
                         N_max,R=50))[sec, ], 
        type = "l" ,lty = 1, xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")

# n grande n = 10000
n = 10000
m = 10
alpha0 = 2
beta0 = 3
theta = rbeta(1, alpha0,beta0)
Y = rbinom(n,m, theta)

alpha_new = alpha0 + sum(Y)
beta_new = n*m - sum(Y) + beta0
theta_n = alpha_new/(alpha_new + beta_new)


N_max = 500000
sec = seq(1,N_max, by = 100)

matplot(sec, t(f_y_yn(alpha_new = alpha_new,
                      beta_new = beta_new,
                      m = m,
                      N_max,R=50))[sec, ], 
        type = "l" ,lty = 1,xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")

matplot(sec, t(f_y_yn_EP(alpha_new = alpha_new,
                         beta_new = beta_new,
                         m = m, Y = Y,
                         N_max,R=50))[sec, ], 
        type = "l" ,lty = 1, xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")





mean(theta_N)
mean(theta_N_EP)
theta_n
var(theta_N)
var(theta_N_EP)
alpha_new*beta_new/((beta_new+ alpha_new)^2*(alpha_new + beta_new + 1))

acf(theta_N)

# Traceplot 
matplot(t(f_y_yn_EP(alpha_new = alpha_new,
                 beta_new = beta_new,
                 m = m, Y=Y,
                 N=50000,R = 10)), type = "l" ,lty = 1, xlab = "Índicet", ylab = "Valor", main = "Traceplot de las filas de M")

system.time(f_y_yn_EP(alpha_new = alpha_new,
                   beta_new = beta_new,
                   m = m,Y=Y,
                   N=20000,R = 1000))
