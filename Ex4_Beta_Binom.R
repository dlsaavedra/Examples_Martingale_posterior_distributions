# Ejemplo 4 Modelo Beta-Binomial  ----

prob_sample <- function(m, alpha_new, beta_new){
  
  c = 1/ choose(alpha_new + beta_new + m - 1, m)
  p = rep(0,m+1)
  for(x0 in 0:m){
    p[x0+1] = choose(alpha_new + x0 - 1, x0) * choose(beta_new + m - x0 - 1, m - x0) * c
  } 
  return(p)
}

f_y_yn <- function(m,alpha_new, beta_new, N, R){
  
  sampling <- Vectorize(function(alpha, beta){
    sample(0:m, size = 1, 
           replace = T, 
           prob = prob_sample(m, alpha, beta))})
  
  alpha_aux = rep(alpha_new, R)
  beta_aux = rep(beta_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  for(j in 1:N){
    theta_matrix[,j]  = alpha_aux/(alpha_aux + beta_aux)
    Y_n_j= sampling(alpha_aux, beta_aux)
    alpha_aux = alpha_aux + Y_n_j 
    beta_aux = beta_aux + m - Y_n_j
    
  }
  return(theta_matrix)
}

# Crear Datos Simulados
n = 1000

m = 20
alpha0 = 2
beta0 = 3
theta = rbeta(1, alpha0,beta0)
Y = rbinom(n,m, theta)

alpha_new = alpha0 + sum(Y)
beta_new = n*m - sum(Y) + beta0
theta_n = alpha_new/(alpha_new + beta_new)
#curve(dgamma(x, alpha_new, shape = beta_new ), xlim = c(0,2))



N = 20000
R = 1000
theta_N = f_y_yn(alpha_new = alpha_new,
                 beta_new = beta_new,
                 m = m,
                 N,R)[,N]
plot(density(theta_N))
curve(dbeta(x,alpha_new,beta_new), add = T, col = "red")

mean(theta_N)
sd(theta_N)
alpha_new/(beta_new + alpha_new)

acf(theta_N)

# Traceplot 
matplot(t(f_y_yn(alpha_new = alpha_new,
                 beta_new = beta_new,
                 m = m,
                 N=50000,R = 20)), type = "l" ,lty = 1, xlab = "Índicet", ylab = "Valor", main = "Traceplot de las filas de M")

system.time(f_y_yn(alpha_new = alpha_new,
                   beta_new = beta_new,
                   m = m,
                   N=1000,R = 1000))
