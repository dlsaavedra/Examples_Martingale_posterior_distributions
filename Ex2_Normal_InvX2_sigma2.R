# Ejemplo 1 Modelo Normal-Inv_X^2 ----
# Mu conocido, sigma^2 desconocido
# Posteriori de sigma2| y es inv_X^2
# Predictiva es t-student
f_y_yn <- function(alpha_new, beta2_new,mu, N, R){
  
  alpha_aux = rep(alpha_new, R)
  beta2_aux = rep(beta2_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  for(j in 1:N){
    theta_matrix[,j]  = alpha_aux*beta2_aux/(alpha_aux - 2)
    Y_n_j = rt(R, alpha_aux)*sqrt(beta2_aux) + mu
    beta2_aux = beta2_aux*alpha_aux
    alpha_aux = alpha_aux + 1 
    beta2_aux = (beta2_aux + (Y_n_j-mu)^2)/alpha_aux
  }
  theta_matrix[, N+1]  = alpha_aux*beta2_aux/(alpha_aux - 2)
  return(theta_matrix)
  
}


n = 100
v0 = 10/2
sigma02 = 5
mu = 500

sigma2 = 1/rgamma(1, shape = v0/2, rate = v0*sigma02/2)
Y = rnorm(n, mu, sqrt(sigma2))
v = sum(Y-mu)^2/n
alpha_new = v0 + n
beta2_new = (v0*sigma02 + n*v)/(v0 + n)
theta_n = alpha_new*beta2_new/(alpha_new - 2)



N = n + 5000
R = 1000
theta_N = f_y_yn(alpha_new = alpha_new,
                  beta2_new = beta2_new,
                  mu = mu,
                  N,R)[,N +1]
plot(density(theta_N))
lines(density(1/(rgamma(10000, shape = alpha_new/2,
                       rate =alpha_new* beta2_new/2))),col = "red")

# Traceplot  
matplot(t(f_y_yn(alpha_new = alpha_new,
                 beta2_new = beta2_new,
                 mu = mu,N=5000,R = 20)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot de las filas de M")


mean(theta_N)
alpha_new*beta2_new/(alpha_new - 2)

acf(theta_N)
