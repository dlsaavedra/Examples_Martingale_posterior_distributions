# Ejemplo 2 Modelo Poisson- Gamma ----

n = 1000
alpha = 1
beta = 20
theta = rgamma(1, alpha, rate = beta)
Y = rpois(n, theta)
alpha_new = alpha + sum(Y)
beta_new = beta + n
theta_n = alpha_new/beta_new
#curve(dgamma(x, alpha_new, shape = beta_new ), xlim = c(0,2))



f_y_yn <- function(alpha_new, beta_new, N, R){
  
  alpha_aux = rep(alpha_new, R)
  beta_aux = rep(beta_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  for(j in 1:N){
    theta_matrix[,j]  = alpha_aux/beta_aux
    Y_n_j = rnbinom(R, alpha_aux, beta_aux/(beta_aux+1))
    alpha_aux = alpha_aux + Y_n_j 
    beta_aux = beta_aux + 1
  }
  return(theta_matrix)
}

N = 30000
R = 10000
theta_N = f_y_yn(alpha_new = alpha_new,
                  beta_new = beta_new, 
                  N,R)[,N]
plot(density(theta_N))
curve(dgamma(x,alpha_new, rate = beta_new), add = T, col = "red")

mean(theta_N)
theta_n
var(theta_N)
alpha_new/beta_new^2

acf(theta_N)

# Traceplot 
matplot(t(f_y_yn(alpha_new = alpha_new,
                 beta_new = beta_new,
                 N=50000,R = 20)), 
        type = "l" ,lty = 1, xlab = "Índicet", 
        ylab = "Valor", main = "Traceplot de las filas de M")



system.time(f_y_yn(alpha_new = alpha_new,
                   beta_new = beta_new,
                   N=2000,R = 100000))
