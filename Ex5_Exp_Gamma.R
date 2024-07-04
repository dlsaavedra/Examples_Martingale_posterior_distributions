# Ejemplo 5 Modelo Exp-Gamma con empirical predictive  ----
library(VGAM)

f_y_yn <- function(alpha_new, beta_new,Y , N, R){
  
  n = length(Y)
  alpha_aux = rep(alpha_new, R)
  beta_aux = rep(beta_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  Y_matrix = matrix(nrow = R, ncol = N + n)
  Y_matrix[,1:n] = matrix(rep(Y, R), nrow = R, byrow = TRUE)
  
  for(j in 1:N){
    
    theta_matrix[,j]  = alpha_aux/beta_aux
    # Extract sampled rows
    Y_matrix[, j + n] = VGAM::rlomax(R, shape3.q = alpha_aux, scale = beta_aux) # Y_matrix[cbind(1:R, sampled_rows)]
    
    alpha_aux = alpha_aux + 1
    beta_aux = beta_aux + Y_matrix[, j + n] 
  }
  
  theta_matrix[,N+1]  = alpha_aux/beta_aux
  return(theta_matrix)
}

# Crear Datos Simulados
n = 1000
alpha0 = 2
beta0 = 3
theta = rgamma(1, alpha0,rate = beta0)
Y = rexp(n, rate =  theta)

alpha_new = alpha0 + n
beta_new = beta0 + sum(Y)
theta_n = alpha_new/beta_new
#curve(dgamma(x, alpha_new, shape = beta_new ), xlim = c(0,2))



N = 5000
R = 10000
theta_N = f_y_yn(alpha_new = alpha_new,
                    beta_new = beta_new,
                    Y = Y,
                    N,R)[,N]
plot(density(theta_N))
curve(dgamma(x,alpha_new,rate = beta_new), add = T, col = "red")

mean(theta_N)
theta_n
var(theta_N)
alpha_new/beta_new^2


acf(theta_N)

# Traceplot 
matplot(t(f_y_yn_EP(alpha_new = alpha_new,
                    beta_new = beta_new,
                    Y=Y,
                    N=10000,R = 20)), type = "l" ,lty = 1, xlab = "Índicet", ylab = "Valor", main = "Traceplot de las filas de M")

system.time(f_y_yn_EP(alpha_new = alpha_new,
                      beta_new = beta_new,
                      Y=Y,
                      N=20000,R = 1000))