############################################################
# DATA AND PACKAGES
############################################################

library(quantreg)
library(parallel)

# data sets overall and exclusion af obese and smokers
df     <- readRDS("./170525_n10090_bmj")
df     <- df[,c("Geschlecht","ALTER", "REL_VO2_MAX", "VO2_MAX", "obese", "ZIGARETTEN")]
dfm    <- subset(df, Geschlecht == "Male")
dff    <- subset(df, Geschlecht == "Female")
dfm_ex <- subset(dfm, obese == "no" & ZIGARETTEN == "no")
dff_ex <- subset(dff, obese == "no" & ZIGARETTEN == "no")

# quantiles
taus       <- seq(from = 0.01, to = 0.99, by = 0.01)

# range of REL vo2max (0.1 to 99.9 percentiles of data)
rvo <- seq(from = 11, to = 61, by = 1)
avo <- seq(from = 0.9, to = 5, by = 0.1)

# range of age
ag <- 25:69

rfrm <- as.formula("REL_VO2_MAX ~ ALTER + I(ALTER^2)")
afrm <- as.formula("VO2_MAX ~ ALTER + I(ALTER^2)")

nboot <- 10000

############################################################
# FUNCTIONS FOR MATRIX "QUANTILE"
############################################################

# equivalent to predict.rq
# predict.rq cannot be applied
# as 10,000*100 quantile regressions (coeficients, residuals etc) are stored in RAM (too large)
# this function is able to predict with only coefficients and not whole rq-object
predict_rq <- function(fit){
  
  Intercept = fit[1,]
  Slope1    = fit[2,]
  Slope2    = fit[3,]
  
  newdata = matrix(rep(25:69, 99),nrow = 45)
  colnames(newdata) <- as.character(taus)
  
  prd <- function(i) Intercept[i] + Slope1[i] * newdata[,i] + Slope2[i] * newdata[,i] * newdata[,i]
  
  for(i in 1:99) newdata[,i] <- prd(i)
  
  newdata <- as.data.frame(newdata)
  newdata$ALTER <- 25:69 # adding ALTER
  
  return(newdata)
  
}

# which predicted value from which tau is the closest for given age and vo2max?
# CAVE: data set has to be called pred!
q.check.fun <- function(age, vo2){
  
  sset <- subset(pred, ALTER == age)
  sset <- sset[,!(names(sset)) == "ALTER"] # drop ALTER
  
  row <- which.min(abs(vo2 - sset))
  
  return(names(row))
  
}

# INNER FUNCTION
# for a given age: taus are calculated for each vo2max
q.inner <- function(x) {
  
  col <- lapply(vo, function(y) q.check.fun(age = x, vo2 = y))
  col <- do.call(data.frame, col)
  col <- t(col)
  row.names(col) <- vo
  
  return(col)
  }

# OUTER FUNCTION
q.outer <- function() {
  
  matrix <- lapply(ag, q.inner)
  
  matrix <- do.call(data.frame, matrix)
  colnames(matrix)  <- as.character(ag)
  row.names(matrix) <- as.character(vo)
  
  return(matrix)
}


## FUNCTION OVERALL
#fit  <- rq(frm, data = in_data, tau = taus)$coefficients
#pred <- predict_rq(fit)
#res  <- q.outer()

## FUNCTION CALLS
t1 <- Sys.time()

#####
## relvo2max
#####
frm       <- rfrm 
vo        <- rvo

# MALES
fit       <- rq(frm, data = dfm, tau = taus)$coefficients
pred      <- predict_rq(fit)
q_m_rel   <- q.outer()
write.csv(q_m_rel, "q_m_rel.csv")

fit       <- rq(frm, data = dfm_ex, tau = taus)$coefficients
pred      <- predict_rq(fit)
q_m_ex_rel<- q.outer()
write.csv(q_m_ex_rel, "q_m_ex_rel.csv")

# FEMALES
fit       <- rq(frm, data = dff, tau = taus)$coefficients
pred      <- predict_rq(fit)
q_f_rel   <- q.outer()
write.csv(q_f_rel , "q_f_rel.csv")

fit       <- rq(frm, data = dff_ex, tau = taus)$coefficients
pred      <- predict_rq(fit)
q_f_ex_rel<- q.outer()
write.csv(q_f_ex_rel , "q_f_ex_rel.csv")

#####
## absvo2max
#####
frm       <- afrm
vo        <- avo

# MALES
fit       <- rq(frm, data = dfm, tau = taus)$coefficients
pred      <- predict_rq(fit)
q_m_abs   <- q.outer()
write.csv(q_m_abs  , "q_m_abs.csv")

fit       <- rq(frm, data = dfm_ex, tau = taus)$coefficients
pred      <- predict_rq(fit)
q_m_ex_abs<- q.outer()
write.csv(q_m_ex_abs, "q_m_ex_abs.csv")

# FEMALES
fit       <- rq(frm, data = dff, tau = taus)$coefficients
pred      <- predict_rq(fit)
q_f_abs   <- q.outer()
write.csv(q_f_abs, "q_f_abs.csv")

fit       <- rq(frm, data = dff_ex, tau = taus)$coefficients
pred      <- predict_rq(fit)
q_f_ex_abs<- q.outer()
write.csv(q_f_ex_abs, "q_f_ex_abs.csv")

t2 <- Sys.time()
t2-t1 # time spent

############################################################
# BOOTSTRAP FUNCTIONS
############################################################

# function to calculate 10,000 bootstrap samples
# output are 10,000 rq coefficients!
boot_all <- function(dat){
  
  cl <- makeCluster(mc <- getOption("cl.cores", parallel::detectCores()-1))
  clusterEvalQ(cl,library(quantreg))
  clusterExport(cl,varlist = c("taus", "dfm","dff","dfm_ex", "dff_ex", "frm"))
  
  # creating bootstrap sample (replace = T)
  gen_sample <- function() dat[sample(nrow(dat), replace = T),]
  
  l1  <- parLapply(cl,1:nboot, function(z) rq(frm, data = gen_sample(), tau = taus)$coefficients )
  
  stopCluster(cl)
  return(l1)
}

# predicted values for all quantile regressions and all ages
# input must be called l
pred_all <- function() {
  
  l1 <- lapply(1:nboot, function(z) {
    
    d        <- predict_rq(l[[z]])
    return(d)
    
  })
  
  return(l1)
}


# function for checking out nearest tau of quantile regression
# list must be called l1
check.fun <- function(age, vo2){
  
  l2    <- lapply(l1, function(x) subset(as.data.frame(x), ALTER == age)[,as.character(taus)])
  qtl   <- lapply(l2, function(x) {
    
    col  <- which.min(abs(x - vo2))
    qtl  <- as.numeric(names(col))
    
    return(qtl)
  })# End of lapply
  
  qtl1 <- t(do.call(data.frame,qtl))
  lo   <- round(quantile(qtl1, probs = c(0.025))*100, digits = 0)
  hi   <- round(quantile(qtl1, probs = c(0.975))*100, digits = 0)
  ci   <- paste0(lo, "% to ", hi, "%")
  
  return(ci)
}

# apply check.fun to all achieved vo2_max (inner function) and all ages (outer function)
# INNER FUNCTION
inner <- function(x) {

  col <- lapply(vo, function(y) check.fun(age = x, vo2 = y))
  col <- do.call(data.frame, col)
  col <- t(col)
  row.names(col) <- vo
  
  return(col)
  
  }

# OUTER FUNCTION
outer <- function() {
  cl <- makeCluster(mc <- getOption("cl.cores", parallel::detectCores()-1))
  clusterExport(cl=cl, varlist=c("vo", "ag", "l1", "inner", "check.fun", "taus"))
  clusterCall(cl, function() library("quantreg"))
  ci_matrix <- parLapply(cl, ag, inner)
  stopCluster(cl)
  
  ci_matrix <- do.call(data.frame, ci_matrix)
  colnames(ci_matrix) <- as.character(ag)
  row.names(ci_matrix) <- as.character(vo)
  
  return(ci_matrix)
}


## FUNCTION OVERALL
#l   <- boot_all(dat_input)
#l1  <- pred_all()
#res <- outer()

## FUNCTION CALLS

#####
## relvo2max
#####
frm       <- rfrm 
vo        <- rvo

# MALES 
l           <- boot_all(dfm)
l1          <- pred_all()
ci_m_rel    <- outer()
write.csv(ci_m_rel, "ci_m_rel.csv")

l           <- boot_all(dfm_ex)
l1          <- pred_all()
ci_m_ex_rel <- outer()
write.csv(ci_m_ex_rel, "ci_m_ex_rel.csv")

# FEMALES
l           <- boot_all(dff)
l1          <- pred_all()
ci_f_rel    <- outer()
write.csv(ci_f_rel, "ci_f_rel.csv")

l           <- boot_all(dff_ex)
l1          <- pred_all()
ci_f_ex_rel <- outer()
write.csv(ci_f_ex_rel, "ci_f_ex_rel.csv")

#####
## absvo2max
#####
frm       <- afrm
vo        <- avo

# MALES
l           <- boot_all(dfm)
l1          <- pred_all()
ci_m_abs    <- outer()
write.csv(ci_m_abs, "ci_m_abs.csv")

l           <- boot_all(dfm_ex)
l1          <- pred_all()
ci_m_ex_abs <- outer()
write.csv(ci_m_ex_abs, "ci_m_ex_abs.csv")

# FEMALES
l           <- boot_all(dff)
l1          <- pred_all()
ci_f_abs    <- outer()
write.csv(ci_f_abs, "ci_f_abs.csv")

l           <- boot_all(dff_ex)
l1          <- pred_all()
ci_f_ex_abs <- outer()
write.csv(ci_f_ex_abs, "ci_f_ex_abs.csv")






# END
