# evppi

library(mgcv)
library(nlme)

# load res_mat
res_mat = read.csv("res_mat.csv")
res_mat = data.frame(res_mat)
colnames(res_mat) = c("QALYs_Standard", "Costs_Standard", "NB_Standard", "QALYs_np1", "Costs_np1", "NB_np1")

res_mat$np1_winner = res_mat$NB_np1 > res_mat$NB_Standard
max_nb = apply(cbind(res_mat$NB_Standard,res_mat$NB_np1),1,max)
nb_diff = max_nb - res_mat$NB_np1
nb_diff = mean(nb_diff)

EVPI = mean(max_nb - max(mean(res_mat$NB_Standard), mean(res_mat$NB_np1))) 
EVPI

# cDR = 0.06
pop_disc = cDR 
pop_per_year = 40000
rollout_years = 1:10
pop_by_year = pop_per_year / (1+pop_disc)^(rollout_years-1)
pop_n = sum(pop_by_year)
EVPI*pop_n


#################



##### VOI gam models
get_voppi = function(x,y=inb,show.plot=T){
  
  require(mgcv)
  evppi=NA
  tryCatch({
    
    model <- gam(y ~ s(x))    
    fittedValues <- fitted(model)
    evppi <- mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))
  }, error=function(e){cat("-> error (probably not enough unique values)")
    # cat("unique values:",length(unique(x)))
  })
  
  if(show.plot){
    require(ggplot2)
    p1 = ggplot() +
      geom_point(aes(x=x,y=y,col="observed")) +
      geom_line(aes(x=x,y=fittedValues,col="fitted")) +
      geom_hline(yintercept = 0) +
      ylab("INB") +
      ggtitle("EVPPI=",round(evppi,2))+
      theme_minimal()
    print(evppi)
    return(p1)
  } else {
    return(evppi)  
  }
  
  
}


# # cholesky missing
# uSuccessP
# uSuccessR
# uRevision
# cRevision
# omrPTHR
# omrRTHR
# rrr


inb <- res_mat$NB_np1 - res_mat$NB_Standard 
get_voppi(x = uSuccessP, y = inb, show.plot = T)
get_voppi(x = uSuccessR, y = inb, show.plot = T)
get_voppi(x = uRevision, y = inb, show.plot = T)
get_voppi(x = cRevision, y = inb, show.plot = T)
get_voppi(x = omrPTHR, y = inb, show.plot = T)
get_voppi(x = omrRTHR, y = inb, show.plot = T)
get_voppi(x = rrr, y = inb, show.plot = T)
get_voppi(x = cholesky_res[,1], y = inb, show.plot = T)
get_voppi(x = cholesky_res[,2], y = inb, show.plot = T)
get_voppi(x = cholesky_res[,3], y = inb, show.plot = T)
get_voppi(x = cholesky_res[,4], y = inb, show.plot = T)
get_voppi(x = cholesky_res[,5], y = inb, show.plot = T)


get_comb_voppi = function(...,y=inb){
  require(mgcv)
  evppi=NA
  tryCatch({
    
    model <- gam(y ~ te(...))    
    fittedValues <- fitted(model)
    evppi <- mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))
  }, error=function(e){cat("-> error (probably not enough unique values)")
    # cat("unique values:",length(unique(x)))
  })
  
  return(evppi)
  
}


get_comb_voppi(y = inb, uSuccessP, uSuccessR, uRevision)
get_comb_voppi(y = inb, cRevision, omrRTHR,rrr)
# with more vars, R might crash...





