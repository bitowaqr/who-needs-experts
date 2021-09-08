

# Survival analysis coefficient
coeff <- c(0.37409680, -5.49093500, -0.03670220, 0.76853600, -1.34447400)
names(coeff) <- c("lngamma", "cons", "age", "male", "NP1")


#Choilesky matrix
what <- c(0.0474501,0,0,0,0,-0.119936522789204,0.169806696467586,0,0,0,5.90093593058813E-07,-0.00461070877953144,0.00242857358178443,0,0,0.000107481333021427,-0.0426020246146981,-0.067292831300903,0.0745125704696859,0,0.00545836573579402,7.45400901505983E-05,-0.0455656518412987,-0.0386465373290773,0.377847881235649)

mat <- t(matrix(what, ncol = 5))

fun <- function(x){
  rnd <- c(qnorm(runif(n = 1,min = .0000000001,max = .99999999999),0,1),
           qnorm(runif(n = 1,min = .0000000001,max = .99999999999),0,1),
           qnorm(runif(n = 1,min = .0000000001,max = .99999999999),0,1),
           qnorm(runif(n = 1,min = .0000000001,max = .99999999999),0,1),
           qnorm(runif(n = 1,min = .0000000001,max = .99999999999),0,1))
  
  err <- t(mat %*% rnd)
  return(coeff + err)
  
}



set.seed(070921)
prob_coeff <- t(sapply(1:1000, fun))
colnames(prob_coeff) <- c("lngamma", "cons", "age", "male", "NP1")


head(prob_coeff)




