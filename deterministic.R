
# REPLICATING Ex66btemp ....
# version 0.1


# DETERMINISTIC ----- !

rm(list=ls())

# meta vals
CYCLES = 60     # how many cycles to run the model for?
cDR = 0.06      # cost discount rate
oDR = 0.015     # outcome discount rate
cRatio = 2200   # willingness to pay ?
INIT_AGE = 60   # age at start
SET_MALE = T    # use males = T,  use females = F

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

# Revision risk standard
draw_RR <- function(n = NA, t_cycle = 1, init_age = 60, male = 1, np1 = F){
  
  gammaC = exp(0.3740968)
  ageC = -0.0367022
  maleC = 0.768536
  cons = -5.490935
  lambda = exp(cons+ageC*init_age+maleC*male)
  np1C = 1
  if(np1){
    np1C = exp(-1.344474)
  }
  res = 1-exp(lambda*((t_cycle-1)^gammaC-t_cycle^gammaC))
  res = 1-exp(lambda*np1C*((t_cycle-1)^gammaC-t_cycle^gammaC))
  return(res)
}

uPrimary = 0
uSuccessP = 0.85
uSuccessR = 0.75
uRevision = 0.30

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
cRevision = 5294
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
omrPTHR = 0.02
omrRTHR = 0.02
rrr = 0.04

res_mat = matrix(
  data = NA,
  ncol = 6, 
  nrow = 1
)

colnames(res_mat) = c(
  "QALYs Standard","Costs Standard", "NB Standard",
  "QALYs np1","Costs np1", "NB np1"
)

RR_standard = lapply(INIT_AGE:(INIT_AGE+CYCLES), function(x){
  draw_RR(n = PSA_IT,init_age = INIT_AGE, male = SET_MALE, np1 = F, t_cycle = x-INIT_AGE+1)
})

RR_np1 = lapply(INIT_AGE:(INIT_AGE+CYCLES), function(x){
  draw_RR(n = PSA_IT,init_age = INIT_AGE, male = SET_MALE, np1 = T, t_cycle = x-INIT_AGE+1)
})



  age = INIT_AGE        # average patient age
  male = SET_MALE        # sex indicator
  
  # INIT MARKOV TRACE STANDARD
  mat_standard = matrix(data = NA, nrow = CYCLES+1, ncol = 5)
  mat_standard[1,] = c(1,0,0,0,0)
  
  # draw trans mats for standard
  trans_mat_standard = lapply(age:(age+CYCLES), function(x){
    getTransMat(
      omrPTHR = omrPTHR,
      omrRTHR = omrRTHR,
      rrr = rrr,
      RR = RR_standard[[x-age+1]],
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
      omrPTHR = omrPTHR,
      omrRTHR = omrRTHR,
      rrr = rrr,
      RR = RR_np1[[x-age+1]],
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
    cRevision = cRevision, 
    c_trt = cStandard, 
    cDR = cDR
  )
  
  c_res_j_np1 = countCosts(
    mat = mat_np1, 
    cPrimary = cPrimary, 
    cSuccess = cSuccess, 
    cRevision = cRevision, 
    c_trt = cNP1, 
    cDR = cDR
  )
  
  u_res_j_standard = countQALYs(
    mat = mat_standard, 
    uPrimary = uPrimary,
    uSuccessP = uSuccessP, 
    uSuccessR = uSuccessR, 
    uRevision = uRevision
  )
  
  u_res_j_np1 = countQALYs(
    mat = mat_np1, 
    uPrimary = uPrimary,
    uSuccessP = uSuccessP, 
    uSuccessR = uSuccessR, 
    uRevision = uRevision
  )
  
  nb_j_standard = u_res_j_standard * cRatio - c_res_j_standard
  nb_j_np1 = u_res_j_np1 * cRatio - c_res_j_np1
  
  res_mat[1,] <- c(
    u_res_j_standard, c_res_j_standard, nb_j_standard,
    u_res_j_np1, c_res_j_np1, nb_j_np1
  )
  
  


# timer
end_timer <- Sys.time() 
elapsed_time = end_timer - start_timer
elapsed_time

# approx 800 it per second, constant between 10k-40k
# ==> 1 million iterations in 21 mins ?!

# mean INB
mean(res_mat[,6] - res_mat[,3])

apply(res_mat, 2, mean)





