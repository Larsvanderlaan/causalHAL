
n <- 5000
W <- replicate(2, runif(n, -1,1))
colnames(W) <- c("W1", "W2")
Y <- rnorm(n, 2*sin(4*W[,1])*sin(4*W[,2]), 0.4)
form <- Y~te(W1, W2, k = 20, bs='bs', m=c(1,1))
out <- gam(form , data= data.frame(Y = Y, W1 = W[,1], W2 = W[,2]))






library(isotone)
n <- 1000
W <- runif(n, -3,3)
pi1 <- plogis(W)
pi0 <- 1-pi1
A <- rbinom(n,1,pi1)

calibrator_pi1 <- as.stepfun(isoreg(pi1, A))
pi1_star <- 1/gpava(1/pi1, A, solver = function(...) {1/weighted.mean(...)})$x
#pi1_star[pi1_star==0] <- min(pi1_star[A==1]) /2


pi0_star <- 1/gpava(1/pi0, 1-A, solver = function(...) {1/weighted.mean(...)})$x
#pi0_star[pi0_star==1] <- min(pi0_star[A==0]) /2

mean((1 )*(A/pi1_star - 1))
mean((1/pi1_star)*(A/pi1_star - 1))

 sort(unique(pi1_star))
mean((pi1_star==0.02500000)*(A/pi1_star - 1))


alpha <- ifelse(A==1, 1/pi1_star, -1/pi0_star)

g <- 1/alpha^2
mg <- pi1_star^2 - pi0_star^2
mean( (alpha *g - mg))

mean(alpha)
mean((alpha  - !is.infinite(mg)*1))



alpha^2 - 2m(alpha)
alpha + f(alpha) is monotone
f(alpha) alpha - m(f(alpha))


