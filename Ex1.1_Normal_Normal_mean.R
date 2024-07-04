# Martingalas Bayes Prediction Ejemplo 1 variando mu0 y sigma02


f_y_yn <- function(mu_new, sigma2_new, sigma2, N, R){
  
  mu_aux = rep(mu_new, R)
  sigma2_aux = rep(sigma2_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  
  for(j in 1:N){
    
    theta_matrix[,j] = mu_aux
    
    Y_n_j = rnorm(R, theta_matrix[, j], sqrt(sigma2_aux + sigma2))
    
    mu_aux = mu_aux/sigma2_aux
    sigma2_aux = 1/(1/sigma2_aux + 1/sigma2)
    mu_aux = (mu_aux + Y_n_j/sigma2)*sigma2_aux
    
  }
  theta_matrix[,N + 1] = mu_aux
  return(theta_matrix)
}


# Datos simulados
n = 100
mu0 = 2000
sigma02 = 1000
theta = rnorm(1, mu0, sqrt(sigma02))
sigma2 = 10
Y = rnorm(n, theta, sqrt(sigma2))
theta_n = (mu0/sigma02 + sum(Y)/sigma2)/(1/sigma02+ n/sigma2)
var_n = 1/(1/sigma02+ n/sigma2)

# Correr las trazas
N = 5000
R = 1000
theta_N = f_y_yn(theta_n, var_n, sigma2, N,R)[,N+1]
plot(density(theta_N))
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "red")

mean(theta_N)
theta
var(theta_N)
var_n

acf(theta_N)


# Traceplot 
matplot(t(f_y_yn(theta_n, var_n, sigma2, N=5000,R = 20)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot de las filas de M")
