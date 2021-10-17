
# REPLICATING Ex66btemp ....
# version 0.1

# rm(list=ls())

# meta vals
set.seed(1234)
CYCLES = 60     # how many cycles to run the model for?
cDR = 0.06      # cost discount rate
oDR = 0.015     # outcome discount rate
cRatio = 2200   # willingness to pay ?
INIT_AGE = 60   # age at start
SET_MALE = F    # use males = T,  use females = F

PSA_IT = 100000  # NUMBER OF PSA ITERATIONS

n_steps = 20
np1C_means_pre = seq(0,2,length.out = n_steps)
cNP1_means_pre = seq(0, 1000, length.out = n_steps)
np1_means = expand.grid(np1C_means_pre, cNP1_means_pre)
np1_means = np1_means[rep(1:n_steps^2, length.out = PSA_IT),]
colnames(np1_means) = c("np1C_means","cNP1_means")

# sex specific mortality rates
mr_male = rep(c(0.00151, 0.00393, 0.0109, 0.0316, 0.0801, 0.1879), each = 10)
mr_female = rep(c(0.00099, 0.0026, 0.0067, 0.0193, 0.0535, 0.1548), each = 10)

# function to return sex and age specific mortality rate
getMR = function(age = 35, male = T, mr_m = mr_male, mr_f = mr_female){
  
  if(male){
    mr = mr_m
  } else {
    mr = mr_f
  }
  
  index = age - 34 +1
  
  if(index > 60){
    return(mr[length(mr)]) 
  } else {
    return(mr[index])
  }
}

# transisition prob drawing functions
# Operative mortality rate following primary THR
draw_omrPTHR <- function(n = 1){
  rbeta(n = n, shape1 = 2, shape2 = 	98)
}
# Operative mortality rate following revision THR 
draw_omrRTHR <- function(n = 1){
  rbeta(n = n, shape1 = 2, shape2 = 	98)
}
# Re-revision risk (assumed to be constant)
draw_rrr <- function(n = 1){
  rbeta(n = n, shape1 = 4, shape2 = 	96)
}

# cholesky decomposition
choleskyDecomp = function(n){
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
  
  # set.seed(070921)
  prob_coeff <- t(sapply(1:n, fun))
  colnames(prob_coeff) <- c("lngamma", "cons", "age", "male", "NP1")
  
  return(prob_coeff)
  
}

# Revision risk standard
draw_RR <- function(n = NA, t_cycle = 1, age = 60, male = 1, cholesky_res, np1C = NULL){
  
  # cholesky_res = choleskyDecomp(n)
  gammaC = exp(cholesky_res["lngamma"])
  ageC = cholesky_res["age"]
  maleC = cholesky_res["male"]
  cons = cholesky_res["cons"]
  lambda = exp(cons+ageC*age+maleC*male)
  if(is.null(np1C)){
    np1C = exp(cholesky_res["NP1"]) 
  }
  res_standard = 1-exp(lambda*((t_cycle-1)^gammaC-t_cycle^gammaC))
  res_np1 = 1-exp(lambda*np1C*((t_cycle-1)^gammaC-t_cycle^gammaC))
  cbind(res_standard, res_np1)
}


# counting QALYs
draw_uSuccessP = function(n){
  rbeta(n, shape1 = 119.57, shape2 = 21.10)
}
draw_uSuccessR = function(n){
  rbeta(n, shape1 = 87.14, shape2 = 29.05)
}
draw_uRevision = function(n){
  rbeta(n, shape1 = 69.70, shape2 = 162.63)
}

uPrimary = 0
uSuccessP = draw_uSuccessP(PSA_IT)
uSuccessR = draw_uSuccessR(PSA_IT)
uRevision = draw_uRevision(PSA_IT)

countQALYs = function(mat, uPrimary = 0, uSuccessP = 0.85, uSuccessR = 0.7, uRevision = 0.32, oDR = 0.015){
  
  # utils by year
  u_by_year = mat %*% c(uPrimary, uSuccessP, uRevision, uSuccessR, 0)
  
  # discount rate by year
  v_disc <- 1/(1+oDR)^(0:(nrow(mat)-1))
  
  # discounted costs by year
  u_by_year_disc  = u_by_year * v_disc
  
  # discounted total
  disc_u_total = sum(u_by_year_disc)
  
  return(disc_u_total)  
}

# counting costs
cPrimary = 0
cSuccess = 0
draw_cRevision = function(n = 1){
  rgamma(n = n, shape = 12.67, scale = 417.67)
}
cRevision = draw_cRevision(PSA_IT)
cStandard = 394
cNP1 = 579

# comp cost of markov mat
countCosts = function(mat, cPrimary = 0, cSuccess = 0, cRevision = 5200, c_trt = 394, cDR = 0.06){
  
  # cost by year
  c_by_year = mat %*% c(c_trt, cSuccess, cRevision, cSuccess, 0)
  
  # discount rate by year
  v_disc <- 1/(1+cDR)^(0:(nrow(mat)-1))
  
  # discounted costs by year
  c_by_year_disc  = c_by_year * v_disc
  
  # discounted total
  disc_c_total = sum(c_by_year_disc)
  
  return(disc_c_total)
}

# TRANS MAT STANDARD FUN
# function to return a time/... dependant trans mat
getTransMat = function(omrPTHR = 0.02, omrRTHR = 0.02, mr = 0.05, RR = 0.003, rrr = 0.04){
  
  # 1 primary THR 
  # 2 successful primary
  # 3 revision THR
  # 4 Successful Revision
  # 5 Death
  
  tmat = matrix(
    ncol = 5, nrow = 5, byrow = T,
    data = c(
      0, 1-omrPTHR,     0,                 0,       omrPTHR,
      0, 1-(RR+mr),    RR,                 0,            mr, 
      0,         0,     0, 1- (omrRTHR + mr),  omrRTHR + mr,
      0,         0,   rrr,        1-(rrr+mr),            mr,
      0,         0,     0,                 0,              1  
    )
  )
  
  return(tmat)
}



# -----------------------------
#    SIMULATION
# -----------------------------

start_timer <- Sys.time()

# draw random params
omrPTHR = draw_omrPTHR(PSA_IT)
omrRTHR = draw_omrRTHR(PSA_IT)
rrr = draw_rrr(PSA_IT)
cholesky_res = choleskyDecomp(PSA_IT)

res_mat = matrix(
  data = NA,
  # what do we want to store?
  ncol = 6, 
  nrow = PSA_IT
)

colnames(res_mat) = c(
  "QALYs Standard","Costs Standard", "NB Standard",
  "QALYs np1","Costs np1", "NB np1"
)


for(j in 1:PSA_IT){
  
  age = INIT_AGE        # average patient age
  male = SET_MALE        # sex indicator
  
  RR_j = draw_RR(n = 1,age = INIT_AGE, male = SET_MALE, t_cycle = 1:(CYCLES+1),cholesky_res = cholesky_res[j,], np1C = rep(np1_means[j,"np1C_means"], length= CYCLES+1))
  
  # INIT MARKOV TRACE STANDARD
  mat_standard = matrix(data = NA, nrow = CYCLES+1, ncol = 5)
  mat_standard[1,] = c(1,0,0,0,0)
  
  # draw trans mats for standard
  trans_mat_standard = lapply(age:(age+CYCLES), function(x){
    getTransMat(
      omrPTHR = omrPTHR[j],
      omrRTHR = omrRTHR[j],
      rrr = rrr[j],
      RR = RR_j[x-age+1,1],
      mr = getMR(
        age = x, 
        male = male, 
        mr_m = mr_male, 
        mr_f = mr_female
      )
    )
  })
  
  # INIT MARKOV TRACE NP1
  mat_np1 = matrix(data = NA, nrow = CYCLES+1, ncol = 5)
  mat_np1[1,] = c(1,0,0,0,0)
  
  # draw trans mats for np1
  trans_mat_np1 = lapply(age:(age+CYCLES), function(x){
    getTransMat(
      omrPTHR = omrPTHR[j],
      omrRTHR = omrRTHR[j],
      rrr = rrr[j],
      RR = RR_j[x-age+1,2],
      mr = getMR(
        age = x, 
        male = male, 
        mr_m = mr_male, 
        mr_f = mr_female
      )
    )
  })
  
  
  # INNER LOOP
  for(i in 1:CYCLES){
    
    trans_mat_i_standard = trans_mat_standard[[i]]
    mat_standard[i+1,] = mat_standard[i,] %*% trans_mat_i_standard
    
    trans_mat_i_np1 = trans_mat_np1[[i]]
    mat_np1[i+1,] = mat_np1[i,] %*% trans_mat_i_np1
    
    age = age + 1
    
  }
  
  
  c_res_j_standard = countCosts(
    mat = mat_standard, 
    cPrimary = cPrimary, 
    cSuccess = cSuccess, 
    cRevision = cRevision[j], 
    c_trt = cStandard, 
    cDR = cDR
  )
  
  c_res_j_np1 = countCosts(
    mat = mat_np1, 
    cPrimary = cPrimary, 
    cSuccess = cSuccess, 
    cRevision = cRevision[j], 
    c_trt = np1_means[j,"cNP1_means"], 
    cDR = cDR
  )
  
  u_res_j_standard = countQALYs(
    mat = mat_standard, 
    uPrimary = uPrimary,
    uSuccessP = uSuccessP[j], 
    uSuccessR = uSuccessR[j], 
    uRevision = uRevision[j]
  )
  
  u_res_j_np1 = countQALYs(
    mat = mat_np1, 
    uPrimary = uPrimary,
    uSuccessP = uSuccessP[j], 
    uSuccessR = uSuccessR[j], 
    uRevision = uRevision[j]
  )
  
  nb_j_standard = u_res_j_standard * cRatio - c_res_j_standard
  nb_j_np1 = u_res_j_np1 * cRatio - c_res_j_np1
  
  res_mat[j,] <- c(
    u_res_j_standard, c_res_j_standard, nb_j_standard,
    u_res_j_np1, c_res_j_np1, nb_j_np1
  )
  
  
}

# timer
end_timer <- Sys.time() 
elapsed_time = end_timer - start_timer
elapsed_time

# mean INB
mean(res_mat[,6] - res_mat[,3])

apply(res_mat, 2, mean)


write.csv(res_mat, "mean_heatmap_dat.csv", row.names = F)


# heatmap
inb = res_mat[,6] - res_mat[,3]
heat_df = data.frame(np1_means,inb)

heat_agg_inb = aggregate(inb ~ cNP1_means + np1C_means, heat_df, mean)

heat_agg_prop = aggregate(inb ~ cNP1_means + np1C_means, heat_df, function(x){sum(x>=0)/length(x)})

# heatmap(heat_df)
library(ggplot2)

# mean INB heatmap
ggplot(data = heat_agg_inb, aes(x=cNP1_means,y=np1C_means, fill = inb)) +
  geom_tile(alpha = 1) +
  scale_fill_gradient2(low="red", high="darkgreen", guide="colorbar") +
  xlab(label = "NP1 costs") +
  ylab(label = "NP1 HR") +
  theme_minimal()


# prop INB>0 heatmap
ggplot(data = heat_agg_prop, aes(x=cNP1_means,y=np1C_means, fill = inb)) +
  geom_tile(alpha = 0.8) +
  scale_fill_gradientn(colours=c("red", "white", "darkgreen"),
                       values=c(0,0.5,1),
                       guide="colorbar", 
                       name = "% inb > 0") +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  xlab(label = "NP1 costs") +
  ylab(label = "NP1 HR") +
  theme_minimal()


