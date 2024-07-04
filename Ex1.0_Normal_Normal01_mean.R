# Martingalas Bayes Prediction Ejemplo 1


f_y_yn <- function(theta_n, n, N, R){
  theta_matrix = matrix(nrow = R, ncol = N +1)
  theta_matrix[,1] = theta_n
  for(j in 1:N){
    Y_n_j = rnorm(R, theta_matrix[, j], sqrt(1/(n + j) + 1))
    theta_matrix[,j + 1] = (theta_matrix[,j] * (n+j) + Y_n_j)/(n+j +1)
  }
  
  return(theta_matrix)
}


# Datos simulados
n = 100
theta = rnorm(1, 0, 1)
Y = rnorm(n, theta, 1)
S_y = sum(Y)
theta_n = S_y/(n+1)
var_n = 1/(n+1)

# Correr las trazas
N1 = 500
N2 = 1000
N3 = 5000
R = 10000
theta_N1 = f_y_yn(theta_n, n, N1,R)[,N1]
theta_N2 = f_y_yn(theta_n, n, N2,R)[,N2]
theta_N3 = f_y_yn(theta_n, n, N3,R)[,N3]
# comparación distintos N ----
plot(density(theta_N1), xlab = expression(theta), ylim = c(0,4.5), col = "red", main ="", lwd = 2)
lines(density(theta_N2), col = "blue", lwd = 2)
lines(density(theta_N3), col = "green", lwd = 2)
curve(dnorm(x,theta_n, sqrt(var_n)), add = T, col = "black", lwd = 3,lty = 2)
legend("topright",c( "Posterior", "Martingale Predictive N = 500",
                     "Martingale Predictive N = 1000", "Martingale Predictive N = 5000"), 
       col = c("black", "red", "blue","green"), cex = .55, lty = c(2,1,1,1))

# Traceplot de la matriz M
matplot(t(f_y_yn(theta_n, n, 6000,R = 50)),
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")


#. Graficos inveros ----
n = 10
theta = rnorm(1, 0, 1)
Y = rnorm(n, theta, 1)
S_y = sum(Y)
theta_n = S_y/(n+1)
var_n = 1/(n+1)
N = 1000
R = 10000
M = f_y_yn(theta_n, n, N,R)
theta_N = M[,N]
D = density(theta_N)
# Range of values for y
y <- seq(-3, 3, length.out = 100)

plot(D$y, D$x, ylim = c(0,2), type = "l", 
     xlim =c(-.01, 1.5),lwd= 3, main = "", xlab = "Density", ylab = expression(theta), col = "red")
lines(dnorm(y,theta_n, sqrt(var_n)), y, col = "black", lwd= 3, lty = 2)
#abline(h = theta, col = "blue", lty = 2)
legend("topright", c("Posteriori theta","Martingale posteriori"),
       col = c("black", "red"), lty =c (2,1,2) ,cex = 0.8)


# Traceplot de la matriz M
matplot(t(M[1:100,]), ylim = c(0,2),
        type = "l",lty = 1, xlab = "Step i", ylab = expression("Posteriori Mean " * theta), 
        main = "")


mean(theta_N)
theta
var(theta_N)
var_n

acf(theta_N)

# Traceplot 
matplot(t(f_y_yn(theta_n, n, N=10000,R = 20)), 
        type = "l" ,lty = 1, xlab = "Índicet",
        ylab = "Valor", main = "Traceplot de las filas de M", )



## Recordar que density puede no ajustar tan bien.

n = 1000
plot(density(rnorm(n,0,1)))
curve(dnorm(x,0,1), add = T, col = "red")

