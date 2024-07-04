# Ejemplo NO funciona bine Modelo Uniforme Pareto con empirical predictive  ----

library(EnvStats)

f_y_yn_EP <- function(alpha_new, xm_new, Y , N, R){
  
  n = length(Y)
  alpha_aux = rep(alpha_new, R)
  xm_aux = rep(xm_new, R)
  theta_matrix = matrix(nrow = R, ncol = N +1)
  Y_matrix = matrix(nrow = R, ncol = N + n)
  Y_matrix[,1:n] = matrix(rep(Y, R), nrow = R, byrow = TRUE)
  
  for(j in 1:N){
    
    theta_matrix[,j]  = alpha_aux*xm_aux/(alpha_aux - 1)
    # Sample row indices
    sampled_rows <- sample(1:((j-1)+n), R, replace = T)
    # Extract sampled rows
    Y_matrix[, j + n] = Y_matrix[cbind(1:R, sampled_rows)]
    
    xm_aux = max(xm_aux, Y_matrix[, j + n])
    alpha_aux = alpha_aux + 1 
    
  }
  
  theta_matrix[,N + 1]  = alpha_aux/(alpha_aux - 1)* xm_aux 
  rm(Y_matrix)
  return(theta_matrix)
}


# Crear Datos Simulados
n = 100

alpha0 = 2
xm0 = 3
theta = rpareto(1, xm0, alpha0)
Y = runif(n, 0, theta)
plot(density(Y))
xm_new = max(c(Y,xm0))
alpha_new = alpha0 + n
theta_n = alpha_new* xm_new/(alpha_new - 1)
#curve(dgamma(x, alpha_new, shape = beta_new ), xlim = c(0,2))

N = 100
R = 100
theta_N <- unlist(mclapply(1:100, function(x) f_y_yn_EP(alpha_new = alpha_new,
                                                        xm_new = xm_new,
                                                        Y = Y,
                                                        N,R)[,N], 
                           mc.cores = detectCores() - 1))

#theta_N = f_y_yn_EP(alpha_new = alpha_new, xm_new = xm_new,
#                    Y = Y, N,R)[,N]

#plot(density(theta_N), xlim = c(0,30))
#curve(dpareto(x, location = xm_new, shape = 2), add = T, col = "red")

plot(density(theta_N), col = "blue", xlab = expression(theta), main = "", 
     lwd = 2, ylim = c(0,1))
curve(dpareto(x, location = xm_new, shape = alpha_new), add = T, col = "black", lwd = 3)
legend("topright",c( "Posterior", "Martingale Empirical predictive"), 
       col = c("black", "blue"), cex = .6, lty = 1)

mean(theta_N)
theta_n
var(theta_N)
xm_new^2*alpha_new/((alpha_new-1)^2*(alpha_new-2))

acf(theta_N)

# Traceplot 
matplot(t(f_y_yn_EP(alpha_new = alpha_new,
                    xm_new = xm_new,
                    Y=Y,
                    N=10000,R = 50)), type = "l" ,lty = 1,xlab = "Step i", 
        ylab = expression("Posteriori Mean " * theta), main = "")

system.time(f_y_yn_EP(alpha_new = alpha_new,
                      xm_new = xm_new,
                      Y=Y,
                      N=20000,R = 1000))