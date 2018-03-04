#_________________________________________________
#_________________________________________________
# TABLE OF CONTENTS
#_________________________________________________
#_________________________________________________
# line
# 25   Specifications
# 77   data
# 246  Random sample
# 443  Descriptive statistics
# 737  Representativity
# 885  Reference Values
# 1327 Validation
# 1790 Multiple Quantile regression
# 2019 comparison of reference values
# 2182 web application
#_________________________________________________
#_________________________________________________
# TABLE OF CONTENTS
#_________________________________________________
#_________________________________________________

##################################################
##################################################
# Specifications
##################################################
##################################################

# loading packages
library(Hmisc)
library(compareGroups)
library(boot)
library(quantreg)
library(grid)
library(gridExtra)
library(xtable)
library(psych)
library(e1071) 
library(reshape2)
library(plyr)
library(dplyr)
library(epitools)
library(directlabels)
library(splines)

# quantiles
TAUS <- c(0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95)

# defining numer of bootstrap samples
NBOOT <- 10000

# labels for plots
REL_LAB <- expression("Relative"~dot("V")*O["2peak"]~"[mL"*O["2"]*"/min/kg"*"]")
ABS_LAB <- expression("Absolute"~dot("V")*O["2peak"]~"[L"*O["2"]*"/min"*"]")

# function for saving plots
SAVE.PLOT <- function(PLOT, FILE.NAME, HEIGHT, WIDTH){
  
  # path
  P1 <- paste0("./results/fig/", FILE.NAME)
  
  # format
  png(filename    = P1, 
      units       = "in", 
      width       = WIDTH, 
      height      = HEIGHT, 
      res         = 300)
  
  plot(PLOT)
  
  dev.off()
  
}

##################################################
##################################################
# Data
##################################################
##################################################

##################################################
# reading data
##################################################

# reading full data
df <- read.table(
  "./data/160301_Originaldaten.csv",
  header = T, sep = ";", dec = ",")

# adding missing data
df2 <- read.table(
  "./data/missings161020.csv",
  header = T, dec = ",", sep = ";")

##################################################
# Calculating new variables
##################################################

# renaming id
names(df)[which(names(df) == "ï..ID")]    <- "id"
names(df2)[which(names(df2) == "ï..ID")]  <- "id"

# adding missing data
df_miss             <- subset(df, id %in% df2$id)
df_miss             <- df_miss[order(df_miss$id),]
df2                 <- df2[order(df2$id),]
df_miss$ALTER       <- df2$ALTER # adding age from df2
df_miss$REL_VO2_MAX <- df2$REL_VO2_MAX # adding rel VO2 from df2
df_miss$VO2_MAX     <- df2$VO2_MAX # adding abs VO2 from df2
df_nomiss           <- subset(df, !(df$id %in% df2$id))
df                  <- rbind(df_miss, df_nomiss)

# define date variable
df$datum  <- as.Date(df$DATUM, "%d.%m.%Y")
df        <- df[order(df$datum),]
df$DATUM  <- NULL

# recode binary variables with -1 / 0
# defining relevant variables
vars <- c("ZIGARETTEN","EXRAUCHER")

# function for changing -1 to 1
f <- function(x){ 
  ifelse(x == -1, 1 , x)}

# apply function to all relevant variables
df[, vars] <- sapply(df[, vars], f)

# computing vars 
# diabetes, hypertension, overweight, obesity
df$MANIFEST_DM        <- ifelse(df$BLUTZUCKER >= 126 | df$HBA1C >= 6.5, 1, 0)
df$BLUTHOCHDRUCK_WHO  <- ifelse(df$SYST >= 140 | df$DIAST >= 90, 1, 0)
df$BMI_GE_25          <- ifelse(df$BMI >= 25, 1, 0)
df$obese              <- ifelse(df$BMI >= 30, 1, 0)

##################################################
# adding study center
##################################################

sto <- read.csv("./data/160810_standort.csv")
df <- merge(df, sto, by = "id", all.x = T)
df$standort <- factor(df$standort, 
                      levels = c("F", "M"), 
                      labels = c("Frankfurt/Ruedesheim", "Munich"))

##################################################
# Implausible cases
##################################################

# Function for assinging NA to implausible cases
f   <- function(var, mi = min(df[var]), ma = max(df[var])){
  d <- ifelse(df[var] > ma | df[var] <= mi, NA, df[[var]])
  d <- as.vector(d)
  return(d)
}

df$ALTER                <- f("ALTER", 15, 100)
df$KOERPERFETT_TANITA   <- f("KOERPERFETT_TANITA", 0, 100)
df$KOERPERFETT_CALIPER  <- f("KOERPERFETT_CALIPER", 0, 100)
df$SYST                 <- f("SYST", 50, 400)
df$DIAST                <- f("DIAST", 30, 200)
df$VO2_MAX              <- f("VO2_MAX", 0, 7)
df$REL_VO2_MAX          <- f("REL_VO2_MAX", 0, 100)
df$GESAMTCHOLESTERIN    <- f("GESAMTCHOLESTERIN", 0, 600)
df$HDL_CHOLESTERIN      <- f("HDL_CHOLESTERIN", 10, 150)
df$TRIGLYCERIDE         <- f("TRIGLYCERIDE", 0, 2000)
df$HBA1C                <- f("HBA1C", 0, 100)

# implausible dates
df <- df[order(df$datum), ]
df$datum[df$datum < as.Date("2000-01-01")] <- NA

# defining factor variables in df
f <- function(var){
  out <- factor(df[,var], levels = c(0,1), labels = c("no", "yes"))
}

df$Geschlecht         <- factor(df$Geschlecht, 
                                levels = c(0,1), 
                                labels = c("Male", "Female"))
df$BMI_GE_25          <- f("BMI_GE_25")
df$obese              <- f("obese")
df$ZIGARETTEN         <- f("ZIGARETTEN")
df$EXRAUCHER          <- f("EXRAUCHER")
df$MANIFEST_DM        <- f("MANIFEST_DM")
df$BLUTHOCHDRUCK_WHO  <- f("BLUTHOCHDRUCK_WHO")

# create variable ex smoker
df$EXRAUCHER[which(df$ZIGARETTEN == "yes")] <- "no"
df$EXRAUCHER[which(df$ZIGARETTEN == "no" & is.na(df$EXRAUCHER))] <- "no"

# add variable age group
df$ageclass <-cut(df$ALTER, 
                  breaks = c(25,30,35,40,45,50,55,60,65,69),
                  include.lowest = TRUE, include.highest = F, right = F)
# creating numerical age class used as linear predictor in regression (validation)
df$ageclass_double <- as.double(df$ageclass)

##################################################
# cleaning data
##################################################

keeps <- c("id", "Geschlecht", "ALTER", 
           "GEWICHT", "GROESSE", 
           "BMI", "BAUCHUMFANG", "KOERPERFETT_CALIPER", "KOERPERFETT_TANITA",
           "SYST", "DIAST",
           "VO2_MAX", "REL_VO2_MAX", "MAX_LEISTUNG",
           "GESAMTCHOLESTERIN","HDL_CHOLESTERIN","LDL_CHOLESTERIN", "TRIGLYCERIDE",
           "ZIGARETTEN", "EXRAUCHER",
           "MANIFEST_DM", "BLUTHOCHDRUCK_WHO", "BMI_GE_25", "obese",
           "ageclass", "ageclass_double",
           "BLUTZUCKER", "HBA1C", "datum", "standort")

df <- df[,which(names(df) %in% keeps)]

##################################################
# Adding Hmisc variable labels
##################################################
label(df$Geschlecht) <- "Sex"
label(df$ALTER) <- "Age [years]"
label(df$ageclass) <- "Age class"
label(df$REL_VO2_MAX) <- "Relative $\\dot{V}O_{2peak}$ [$\\unitfrac{mL}{min*kg}$]"
label(df$VO2_MAX) <- "Absolute $\\dot{V}O_{2peak}$ [$\\unitfrac{L}{min}$]"
label(df$GEWICHT) <- "Weight [kg]"
label(df$GROESSE) <- "Height [cm]"
label(df$BMI) <- "BMI [$\\unitfrac{kg}{m^2}$]"
label(df$KOERPERFETT_CALIPER) <- "Body fat Caliper [\\%]"
label(df$BMI_GE_25) <- "Overweight"
label(df$obese) <- "Obese"
label(df$SYST) <- "Systolic [mmHg]"
label(df$DIAST) <- "Diastolic [mmHg]"
label(df$BLUTHOCHDRUCK_WHO) <- "Hypertension"
label(df$BLUTZUCKER) <- "Blood glucose [$\\unitfrac{mg}{dl}$]"
label(df$HBA1C) <- "$HbA1_c$ [\\%]"
label(df$MANIFEST_DM) <- "Diabetes mellitus"
label(df$GESAMTCHOLESTERIN) <- "Total cholesterol [$\\unitfrac{mg}{dl}$]"
label(df$HDL_CHOLESTERIN) <- "HDL Cholesterol [$\\unitfrac{mg}{dl}$]"
label(df$LDL_CHOLESTERIN) <- "LDL Cholesterol [$\\unitfrac{mg}{dl}$]"
label(df$TRIGLYCERIDE) <- "Triglycerides [$\\unitfrac{mg}{dl}$]"
label(df$ZIGARETTEN) <- "Smoker"
label(df$EXRAUCHER) <- "Ex smoker"
label(df$standort) <- "Study center"

##################################################
##################################################
# Random sample
##################################################
##################################################

# draw random sample
set.seed(1)
spl<-df[sample(nrow(df),252),]

# comparing random sample and full data set

# df for comparison of df and spl
spl$spl       <- 1
df_no_spl     <- subset(df, !(df$id %in% spl$id))
df_no_spl$spl <- 0
test_spl      <- rbind(df_no_spl, spl)
test_spl$spl  <- factor(test_spl$spl, 
                        levels = c(0,1), 
                        labels = c("Full data", "Random sample"))

# defining quatitative and categorical variables
qvars<-c("ALTER","REL_VO2_MAX","VO2_MAX", "GEWICHT","GROESSE","BMI")
cvars<-c("Geschlecht","ZIGARETTEN","BMI_GE_25", "standort")

# creating formula
frm <- as.formula(paste0("spl ~", paste(c(qvars, cvars), collapse = " + ")))

# creating bivariate table for sample vs. full data
tab <- compareGroups(frm , data = test_spl, method = 2)
tab <- createTable(tab)
export2latex(tab)

##################################################
# reading file with manually added information
##################################################

# read sample Ruedesheim and create two data sets
# sr: only sample characteristics
# sf: merged with original data 
# (without excluded cases as random IDs were drawn prior to exclusion)
sr <- read.csv("./data/sample_Rued_ergaenzt_160622.csv", 
               sep=";", dec=",", header=T)
sr <- merge(sr, df, by.x = "ID", by.y = "id", all.x = T)

##################################################
# measures of exertion
##################################################

# calculating age-predicted maximal heart rate (APMHR, pmid: 11153730)
# the criterion for maximal effort is 90% of APMHR (pmid: 18027991)
sr$th_hr_max <- (208 - 0.7 * sr$ALTER) * 0.9

# measure of exertion: lac>8 OR RER >1.1 OR HR >90% of APMHR
sr$lac_ge8 <- ifelse(sr$Laktat >= 8,1,0)
sr$rer_ge1 <- ifelse(sr$RER    >= 1.1,1,0)
sr$hr_ge   <- ifelse(sr$HF.max >= sr$th_hr_max,1,0)

# specifying rows where lac, rer and hr are all NA
sr$na       <- rowSums(is.na(cbind(sr$rer_ge1,sr$lac_ge8,sr$hr_ge)))
sr$Ex       <- rowSums(cbind(sr$lac_ge8,sr$rer_ge1,sr$hr_ge),na.rm=T)
sr$Ex       <- ifelse(sr$na == 3, NA, sr$Ex)
sr$Exertion <- ifelse(sr$Ex >=1,1,sr$Ex)

# Bootstrapped 95% confidence intervals for proportions
SAMPLE.PROPORTION <- function(x, d) {
  z = sum(x[d]==1, na.rm = T)
  n = sum(x[d]==0, na.rm = T)
  p = z/(z+n)
  return(p)
}

# function to list N(%) [95% CI] of exertion
EXERTION <- function(SUBGROUP = ""){
  
  x       <- sr[SUBGROUP, "Exertion"]
  YES     <- sum(x == 1, na.rm = T)
  NO      <- sum(x == 0, na.rm = T)
  N       <- paste0(YES, "/", YES + NO)
  PERCENT <- round((YES / (YES + NO))*100, digits = 0)
  
  set.seed(1)
  b       <- boot(x, SAMPLE.PROPORTION, R = NBOOT)
  CI      <- paste(round(quantile(b$t, probs = c(0.025, 0.975))*100, digits = 0))
  return(paste0(N, " (", PERCENT, "%)", " [", CI[1], "% to ", CI[2], "%]"))
  
}

# overall, males, females
EXERTION(sr$sex == "m" | sr$sex == "w") 
EXERTION(sr$sex == "m") 
EXERTION(sr$sex == "w")

##################################################
# does inclusion of cases with exertion == "no" change results?
##################################################

sr$Exertion <- factor(sr$Exertion, levels = c(0, 1), labels = c("no", "yes"))

# checking visually using scatter plots
# producing scatter plot with two regression lines:
# 1) all cases of the sample 2) only cases with Exertion == "yes"
# for males and females / absolute and relative vo2max (fig 4.13)

SC.PLOT <- function(SEX = "m", OUTCOME, XLAB, YLAB, TITLE = ""){
  
  DAT       <- subset(sr, sex == SEX)
  DAT.NOEX  <- subset(DAT, Exertion == "yes")
  
  FRM      <- as.formula(paste(OUTCOME, "~ ALTER"))
  
  FIT      <- rq(FRM, data = DAT,      tau = 0.5)$coefficients
  FIT.NOEX <- rq(FRM, data = DAT.NOEX, tau = 0.5)$coefficients
  
  gp <- ggplot(data = DAT, aes_string(x = "ALTER", y = OUTCOME, colour = "Exertion")) + 
    geom_point(na.rm = T) +
    geom_abline(intercept = FIT[1],      slope = FIT[2]) +
    geom_abline(intercept = FIT.NOEX[1], slope = FIT.NOEX[2], lty = 2) +
    scale_x_continuous(limits = c(25,69), breaks = seq(from = 25, to = 69, by = 5)) +
    scale_y_continuous(limits = c(min(DAT[,OUTCOME]), max(DAT[,OUTCOME]))) +
    labs(x = XLAB, y = YLAB, title = TITLE) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="none")
  
  return(gp)
}

SAVE.PLOT(
  PLOT = grid.arrange(
    SC.PLOT("m", "REL_VO2_MAX", "",            REL_LAB, "Male"),
    SC.PLOT("w", "REL_VO2_MAX", "",            REL_LAB, "Female"),
    SC.PLOT("m", "VO2_MAX",     "Age [years]", ABS_LAB, ""),
    SC.PLOT("w", "VO2_MAX",     "Age [years]", ABS_LAB, ""),
    ncol = 2),
  "sample.png", 8, 10
)

# function to create an xtable with ANCOVA- results 
# Question: does Exertion significantly alter Results?
# Median regression: 
# DV = VO2peak 
# IV = ALTER + Exertion (yes/no) + interaction term ALTER * Exertion
ANCOVA <- function(dat = srf, outcome = "REL_VO2_MAX"){
  
  # 1) data without missings (avoid error by boot.rq)
  dat      <- dat[!is.na(dat[,outcome   ]),]
  dat      <- dat[!is.na(dat[,"ALTER"   ]),]
  dat      <- dat[!is.na(dat[,"Exertion"]),]
  
  # 2) how many cases with "exertion = no"?
  print(table(dat[,"Exertion"]))
  
  # 2) extract variables
  A        <- dat[,"ALTER"]
  E        <- as.numeric(dat[, "Exertion"])
  DV       <- dat[, outcome]
  
  # 3) create formula
  frm      <- paste(outcome, "~ A + E + I(A * E)")
  frm      <- as.formula(frm)
  
  # 4) estimate for coefficient and bootstrapped P values
  fit      <- rq(frm, tau = 0.5, data = dat)
  estpv    <- summary(fit, se = "boot")$coefficients
  est      <- estpv[,1]
  pval     <- estpv[,4]
  
  # 5) bootstrapped confidence intervals
  set.seed(1)
  fit.b    <- boot.rq(cbind(1, A, E, I(A * E)), DV, tau = 0.5, R = NBOOT)
  ci       <- t(apply(fit.b$B, 2, quantile, c(0.025,0.975)))
  
  # 6) cbind to create a table
  tab      <- cbind(est, ci, pval)
  tab      <- as.data.frame(tab)
  row.names(tab) <-  c("Intercept", "Age", "Exertion(yes)", "Age * Exertion(yes)")
  
  print(xtable(tab,digits=c(2,2,2,2,3)), include.rownames = T)
  
} 

ANCOVA(dat = subset(sr, sex == "m"), outcome = "REL_VO2_MAX")
ANCOVA(dat = subset(sr, sex == "m"), outcome = "VO2_MAX")
ANCOVA(dat = subset(sr, sex == "w"), outcome = "REL_VO2_MAX")
ANCOVA(dat = subset(sr, sex == "w"), outcome = "VO2_MAX")

## PLAUSIBILITY CHECK:
## visualise the regression lines for exertion = yes/no seperately
#pdat                  <- data.frame(age  = seq(from = 20, to = 80))
#pdat$pred_exertion    <- 18.27 + (0.11*pdat$age) + (5.85*2) + (-0.07*2*pdat$age) 
#pdat$pred_no_exertion <- 18.27 + (0.11*pdat$age) + (5.85*1) + (-0.07*1*pdat$age) 

#ggplot() + 
#  geom_point(aes(x = ALTER, y = REL_VO2_MAX), data = srf, na.rm = T) + 
#  geom_line(aes(x = age, y = pred_exertion),    data = pdat, colour = "red") + 
#  geom_line(aes(x = age, y = pred_no_exertion), data = pdat, colour = "blue")


##################################################
##################################################
# Descriptive statistics
##################################################
##################################################

##################################################
# Numbers of cases for flow chart
##################################################

# function for counting non-missings
n   <- function(x) sum(!is.na(x))
ex  <- function(x) sum(is.na(x))

# crude n
nrow(df)
table(df$Geschlecht)

# excluded cases due to missing values
ex(df$REL_VO2_MAX | df$VO2_MAX)
ex(df[!is.na(df$REL_VO2_MAX), "standort"])
ex(df[!is.na(df$REL_VO2_MAX) &
        !is.na(df$standort),"ALTER"])

# n after exclusion
n(df[!is.na(df$REL_VO2_MAX) &
       !is.na(df$standort) & 
       !is.na(df$ALTER),"REL_VO2_MAX"])
tapply(df[!is.na(df$REL_VO2_MAX) &
            !is.na(df$standort) &
            !is.na(df$ALTER), "REL_VO2_MAX"],
       df[!is.na(df$REL_VO2_MAX) &
            !is.na(df$standort) &
            !is.na(df$ALTER), "Geschlecht"], n)

# exclusion of age>=70 or <25
ex(df[!is.na(df$REL_VO2_MAX) &
        !is.na(df$standort) & 
        !is.na(df$ALTER), "ageclass"])

n(df[!is.na(df$REL_VO2_MAX) &
       !is.na(df$standort) & 
       !is.na(df$ALTER) &
       !is.na(df$ageclass),"REL_VO2_MAX"])

# exclude cases
df <- df[! is.na(df$REL_VO2_MAX),]
df <- df[! is.na(df$standort),]
df <- df[! is.na(df$ALTER),]
df <- df[! is.na(df$ageclass),]
tapply(df$VO2_MAX,df$Geschlecht,n)

# split data by location
df_FR <- subset(df, standort == "Frankfurt/Ruedesheim")
df_M  <- subset(df, standort == "Munich")

# split data by sex
dfm <- subset(df, Geschlecht == "Male")
dff <- subset(df, Geschlecht == "Female")

# split data by sex and location
dfm_FR <- subset(df_FR, Geschlecht == "Male")
dff_FR <- subset(df_FR, Geschlecht == "Female")
dfm_M  <- subset(df_M, Geschlecht == "Male")
dff_M  <- subset(df_M, Geschlecht == "Female")

##################################################
# Date and location
##################################################

SAVE.PLOT(
  
  # density function of date of examination
  ggplot(df, aes(x = datum, colour = Geschlecht)) + 
    stat_bin(data=subset(df,Geschlecht=="Male"),aes(y=cumsum(..count..)),geom="step", na.rm = T)+
    stat_bin(data=subset(df,Geschlecht=="Female"),aes(y=cumsum(..count..)),geom="step", na.rm = T)+
    scale_x_date(breaks = seq(as.Date("2001-01-01"), as.Date("2015-12-31"), by="2 years"), date_labels = "%Y") +
    scale_y_continuous(breaks = seq(from = 0, to = 6000, by = 1000)) +
    scale_colour_discrete(name  ="Sex",
                          breaks=c("Female", "Male")) +
    geom_vline(xintercept = as.Date("2001-01-01", "%Y-%m-%d"), colour = "grey60", size = 0.3 ) +
    geom_vline(xintercept = as.Date("2015-12-31", "%Y-%m-%d"), colour = "grey60", size = 0.3 ) +
    labs(y = "Cumulative number of participants", x = "Date")
  
  , "time.png", 5, 5
)

##################################################
# Bivariate descriptive Table
##################################################

# bivariate desctiptive table

# exclude missing values from
# age, relative and absolute VO2peak
# use only subset of variables
vars <- c(
  # id, independent and dependent variables
  "Geschlecht", "REL_VO2_MAX", "VO2_MAX", "ALTER", "ageclass",
  # descriptive statistics and multivariable regression
  "GEWICHT", "GROESSE", "BMI", "KOERPERFETT_CALIPER", "BMI_GE_25", "obese",
  "SYST", "DIAST", "BLUTHOCHDRUCK_WHO", "BLUTZUCKER", "HBA1C", "MANIFEST_DM",
  "GESAMTCHOLESTERIN", "HDL_CHOLESTERIN", "LDL_CHOLESTERIN", "TRIGLYCERIDE",
  "ZIGARETTEN", "EXRAUCHER", "standort")

df1 <- df[, vars]

# description of full data
tab <- compareGroups(
  Geschlecht ~ ., 
  data = df1, 
  method = 2)
export2latex(createTable(tab),
             header=c(p.overall='P value'),
             caption = "Descriptive statistics by sex.")

# description of data for reference values
df1m <- subset(df1, Geschlecht == "Male")
df1f <- subset(df1, Geschlecht == "Female")

males <- compareGroups(
  standort ~ . -Geschlecht, 
  data = df1m, 
  method = 2)

females <- compareGroups(
  standort ~ . -Geschlecht, 
  data = df1f, 
  method = 2)

males   <- createTable(males)
females <- createTable(females)

export2latex(cbind(Males = males, Females = females),
             header=c(p.overall='P value'),
             caption = "Descriptive statistics by sex and study center.")


##################################################
# Description of all variables in data set
##################################################

# description of all variables in original data set

# quantitative
qvars <- c(
  # target variables
  "ALTER" ,"REL_VO2_MAX", "VO2_MAX",
  # other quantitative variables
  "GEWICHT", "GROESSE", "BMI", "KOERPERFETT_CALIPER", "KOERPERFETT_TANITA", "BAUCHUMFANG",
  "SYST", "DIAST", "BLUTZUCKER", "HBA1C",
  "GESAMTCHOLESTERIN", "HDL_CHOLESTERIN", "LDL_CHOLESTERIN", "TRIGLYCERIDE")

qvars_labels <- c(
  # target variables
  "Age" , "Relative \n VO2peak", "Absolute \n VO2peak" ,
  # other quantitative variables
  "Weight", "Height", "BMI", "Body fat \n (Caliper)", "Body fat \n (Tanita)", "Waist \n circumference",
  "Systolic", "Diastolic", "Blood \n glucose", "HbA1c",
  "Total \n cholesterol", "HDL \n cholesterol", "LDL \n cholesterol", "Tri- \n glycerides")

# equations as labels in histograms
qvars_labels2     <- qvars_labels 
qvars_labels2[2]  <- expression("Relative"~dot("V")*O["2peak"])
qvars_labels2[3]  <- expression("Absolute"~dot("V")*O["2peak"])
qvars_labels2[13] <- expression("HbA1"["c"])

# histograms for all qvars using labels
hst <- function(dat){
  
  # empty list to save in
  lst <- list()
  
  for(i in qvars){
    # calculate means, sd and the index of label
    m   <- mean(dat[, i], na.rm = T)
    s   <- sd(dat[, i], na.rm = T)
    num <- qvars_labels2[which(qvars %in% i)]
    
    # plot
    g <- ggplot(data = dat, aes_string(x = i)) +
      geom_histogram(aes(y=..density..), color="black", bins = 15,  na.rm=  T) +
      scale_x_continuous(limits = c(m - 3*s, m + 3*s)) + 
      stat_function(fun = dnorm, args = list(mean = m, sd = s), color = "red", na.rm = TRUE) +
      labs(x="", y="") +
      ggtitle(num)
    
    # save plot in list
    lst[[i]] <- g }
  
  return(lst)
}

# Arranging all plots to one

# males
lst_m <- hst(dfm_FR)
gm <- do.call(grid.arrange, c(lst_m, list(ncol = 6, top = "Males")))
SAVE.PLOT(gm, "histograms_male.png", 10,15)

# females
lst_f <- hst(dff_FR)
gf <- do.call(grid.arrange, c(lst_f, list(ncol = 6, top = "Females")))
SAVE.PLOT(gf, "histograms_female.png", 10,15)

# correlation matrices (SPLOM)
png(height = 4500, width = 6000, pointsize = 20, res = 300,
    file="./results/fig/corplot_males.png")

pairs.panels(
  setNames(dfm_FR[, qvars], qvars_labels),
  pch = ".", 
  gap = 0,
  method = "spearman",
  density = F,
  hist.col = "white",
  rug = F)
dev.off()

png(height = 4000, width = 5000, pointsize = 20, res = 300,
    file="./results/fig/corplot_females.png")

pairs.panels(
  setNames(dff_FR[, qvars], qvars_labels),
  pch = ".", 
  gap = 0,
  method = "spearman",
  density = F,
  hist.col = "white",
  rug = F)
dev.off()

# defining quantile-based skewness (Hao, 2007: 14)
qsk <- function(x){
  qsk_numerator   <- 
    quantile(x, probs = 0.90, na.rm = T) - 
    quantile(x, probs = 0.5, na.rm = T)
  qsk_denominator <- 
    quantile(x, probs = 0.5, na.rm = T) - 
    quantile(x, probs = 0.1, na.rm = T) 
  qsk <- (qsk_numerator / qsk_denominator) - 1
  return(qsk)}

# Skewness and qsk for all quantitative variables
melted <- melt(df_FR[, c(qvars, "Geschlecht")], id.vars = "Geschlecht")

grouped <- group_by(melted, Geschlecht, variable)
tab     <- summarise(grouped, skewness = skewness(value, na.rm = T), qsk=qsk(value))
tab     <- cbind(tab[tab$Geschlecht == "Male", 2:ncol(tab)],
                 tab[tab$Geschlecht == "Female", 3:ncol(tab)])

print(xtable(tab), include.rownames = F)

##################################################
# Box plots Vo2 by age class
##################################################

# Function because BPLs are also needed for validation
PLT.BPL <- function(DAT, REL = T){
  
  # ageclass is needed as "factor"
  DF1 <- ddply(DAT, c("ageclass", "Geschlecht"), transform, N = length(ageclass))
  DF1$label <- paste0(DF1$ageclass, "\n" , "(n=",DF1$N,")")
  
  # labels and axes
  if(REL == T){
    OUTCOME = "REL_VO2_MAX"
    LAB = REL_LAB
    LMT = c(0,60)
    BRK = seq(from=0,to=60,by=5)
  }else{
    OUTCOME = "VO2_MAX"
    LAB = ABS_LAB
    LMT = c(0,6) 
    BRK = seq(from = 0,to = 6,by = 0.5)
  }
  
  # absolute VO2peak
  BPL <- ggplot(DF1, aes_string(x = "label", y = OUTCOME)) + 
    geom_boxplot(na.rm = T) + 
    scale_y_continuous(limits = LMT, breaks = BRK) +
    facet_grid(.~Geschlecht, scales = "free")+
    labs(x = "", 
         y = LAB)
  
  return(BPL)
  
}

SAVE.PLOT(grid.arrange(PLT.BPL(df_FR, REL = F),
                       PLT.BPL(df_FR, REL = T),
                       ncol=1),
          "bpl.png", 8, 10)

##################################################
##################################################
# Representativity
##################################################
##################################################

# read standard population
sb <- read.table("./data/Altersverteilung_Zensus_2011.csv", 
                 header = T, sep = ";", dec = ",")

# same age groups for standard population and PF sample
df_FR$a.standard <- cut(df_FR$ALTER, 
                        breaks = c(25,30,35,40,45,50,55,60,65,69),
                        include.lowest = TRUE, include.highest = F, right = F)

sb$a.standard <- cut(sb$Alter, 
                     breaks = c(25,30,35,40,45,50,55,60,65,69),
                     include.lowest = TRUE, include.highest = F, right = F)

# restricting census age groups to the margins of the present sample
sb <- subset(sb, !is.na(a.standard))

# aggregating standard population by age group
sb <- aggregate(cbind(Male, Female) ~ a.standard, data = sb, sum)

# bootstrap function
ias <- function(dat, sex, var, b, standard = sb){
  
  set.seed(1)
  
  bt <- dat[b,]
  
  # variables
  repvars <- c("ZIGARETTEN",
               "EXRAUCHER", 
               "BMI_GE_25", 
               "obese", 
               "BLUTHOCHDRUCK_WHO")
  
  # aggregating PF sample by age group
  events <- aggregate(
    x = bt[,repvars],
    by = bt[,c("a.standard", "Geschlecht")],
    function(x) sum(x=="yes", na.rm =T))
  
  n_per_AG <- aggregate(
    x = bt[,repvars],
    by = bt[,c("a.standard","Geschlecht")],
    function(x) sum(!is.na(x)))
  
  e <- subset(events, Geschlecht == sex)
  n <- subset(n_per_AG, Geschlecht == sex)
  
  # event rate
  er <- epitools::ageadjust.direct(
    e[,var], n[,var], stdpop = standard[,sex])
  
  return(er*100)
}

# function to execute bootstrap for sexes and variables
bt <- function(sx, vr){
  
  b = boot::boot(
    data = df_FR,
    sex = sx, var = vr, statistic = ias, 
    R = NBOOT
  )
  
  return(
    c(b$t0[2],
      quantile(b$t[,2], probs = c(0.025,0.975))
    )
  )
}

# Creating overall table

PF_MALE <- data.frame(
  
  Smoker        = bt("Male", "ZIGARETTEN"),
  Former_Smoker = bt("Male", "EXRAUCHER"),
  Overweight    = bt("Male", "BMI_GE_25"),
  Obesity       = bt("Male", "obese"),
  Hypertension  = bt("Male", "BLUTHOCHDRUCK_WHO")
  
)

PF_FEMALE <- data.frame(
  
  Smoker        = bt("Female", "ZIGARETTEN"),
  Former_Smoker = bt("Female", "EXRAUCHER"),
  Overweight    = bt("Female", "BMI_GE_25"),
  Obesity       = bt("Female", "obese"),
  Hypertension  = bt("Female", "BLUTHOCHDRUCK_WHO")
  
)

DEGS_MALE <- data.frame(
  
  Smoker        = c(26.1, 24.0, 28.2),
  Former_Smoker = c(33.7, 31.9, 35.5),
  Overweight    = c(67.1, 65.0, 69.2),
  Obesity       = c(23.3, 21.2, 25.4),
  Hypertension  = c(33.3, 31.1, 35.6)
  
)

DEGS_FEMALE <- data.frame(
  
  Smoker        = c(21.4, 19.7, 23.1),
  Former_Smoker = c(22.8, 21.4, 24.2),
  Overweight    = c(53.0, 50.8, 55.1),
  Obesity       = c(23.9, 22.0, 25.9),
  Hypertension  = c(29.9, 28.1, 31.9)
  
)

TAB <- round(as.data.frame(cbind(
  
  t(PF_MALE),
  t(DEGS_MALE),
  t(PF_FEMALE),
  t(DEGS_FEMALE)
  
)), digits = 1)

LTAB <- cbind(
  paste0(TAB[,1], " (", TAB[,2]," to ", TAB[,3], ")"),
  paste0(TAB[,4], " (", TAB[,5]," to ", TAB[,6], ")"),
  paste0(TAB[,7], " (", TAB[,8]," to ", TAB[,9], ")"),
  paste0(TAB[,10], " (", TAB[,11]," to ", TAB[,12], ")")
)

LTAB <- as.data.frame(LTAB)
names(LTAB) <- c("PF_M", "DEGS_M", "PF_F", "DEGS_F")
row.names(LTAB) <- row.names(TAB)

# saving results as a list on hard drive
TAB_REPRESENTATIVITY <- list(
  
  DATA = TAB,
  LATEX = xtable(LTAB)
  
)

saveRDS(TAB_REPRESENTATIVITY, "./results/dat/tab_representativity.rds")

##################################################
##################################################
# Reference values
##################################################
##################################################

##################################################
# Quantiles
##################################################

# relative VO2peak

# aggregating by TAUS and age class
tab1 <- aggregate(REL_VO2_MAX ~ ageclass + Geschlecht, 
                  data = df_FR, 
                  FUN = function(x) quantile(x, probs = TAUS))

tab2 <- aggregate(REL_VO2_MAX ~ ageclass + Geschlecht, 
                  data = df_FR, 
                  FUN = function(x) length(x))

tab  <- do.call(data.frame, merge(tab1, tab2, by = c("Geschlecht", "ageclass")))
names(tab) <- c("Sex", "Age class", TAUS, "N")

print(xtable(tab, 
             digits = c(0, 0, rep(1, 15))), 
      include.rownames = F)

## PLAUSIBILITY CHECK:
# sum(df_FR$Geschlecht == "Female" & df_FR$ageclass == "[25,30)")
# quantile(subset(dfm_FR, ageclass == "[65,69]")$REL_VO2_MAX, probs = TAUS)
# correct!

# absolute VO2peak
tab1 <- aggregate(VO2_MAX ~ ageclass + Geschlecht, 
                  data = df_FR, 
                  FUN = function(x) quantile(x, probs = TAUS))

tab2 <- aggregate(VO2_MAX ~ ageclass + Geschlecht, 
                  data = df_FR, 
                  FUN = function(x) length(x))

tab  <- do.call(data.frame, merge(tab1, tab2, by = c("Geschlecht", "ageclass")))
names(tab) <- c("Sex", "Age class", TAUS, "N")

print(xtable(tab,
             digits = c(0, 0, rep(2, 15))), 
      include.rownames = F)

## PLAUSIBILITY CHECK:
# quantile(
#  subset(df, 
#          Geschlecht == "Male" & 
#          ALTER >= 25 & 
#          ALTER < 30)$REL_VO2_MAX,
#   probs = c(0.05, 0.5, 0.95)
# )
## --> correct!

##################################################
# Nomograms: calculation data
##################################################

PLT.DAT <- function(dat, var){
  
  # empty list to save results
  lst <- list()
  
  # loop over all taus and calculate predicted values
  for(i in TAUS){
    
    # fit quantile regression
    fit   <- rq(dat[,var] ~ ALTER + I(ALTER^2), data = dat, tau = i)
    
    # create empty prediction data set
    pdat  <- data.frame(ALTER = seq(from = 25, to = 70, length.out = 1000))
    
    # calculate predicted values and confidence intervals, and adding tau
    set.seed(1)
    pdat  <- 
      cbind(pdat, tau = i,
            predict.rq(fit, newdata = pdat, interval = "confidence", level = .95,
                       se = "boot", type = "percentile", R = NBOOT)
      )
    
    # creating list with all results
    nam <- paste("pdat", i, sep = "")
    lst[[nam]] <- pdat
  }
  
  ############################################
  ## >pdat
  ## head(lst$pdat0.95)
  
  #  ALTER    tau  fit      lower    higher
  #1 25.00000 0.95 52.41516 49.79706 54.98178
  #2 25.04505 0.95 52.40571 49.79676 54.96343
  #3 25.09009 0.95 52.39624 49.79642 54.94508
  
  ## --> all entries (taus) should be row-binded
  ############################################
  
  # row-binding for all prediction data frames (pdat)
  ds_overall <- do.call(rbind, lst)
  
  # express tau as perventile
  ds_overall$tau        <- ds_overall$tau*100
  ds_overall$percentile <- paste(ds_overall$tau, "%", sep = "")
  
  return(ds_overall)
  
}# plt.dat function end


# Apply PLT.DAT to desired data sets

# list of all desired data sets
# --> 8 data sets --> perfect for parallelisation with 4 cpus
DSETS <- list(
  
  # males and females with all observations
  dfm = dfm,
  dff = dff,
  
  # males and female, smokers and obese excluded
  # exm = exclusion of smokers&obese, male
  exm = subset(dfm, ZIGARETTEN == "no" & obese == "no"),
  exf = subset(dff, ZIGARETTEN == "no" & obese == "no"),
  
  # including only participants from Frankfurt/Ruedesheim
  # dfrm = df Ruedesheim, male
  dfm_FR = dfm_FR,
  dff_FR = dfm_FR,
  
  # including only participants from F/R and no smokers/obese
  # exrm = exclusion of smoker&obese, Ruedesheim, male
  exm_FR = subset(dfm_FR, ZIGARETTEN == "no" & obese == "no"),
  exf_FR = subset(dff_FR, ZIGARETTEN == "no" & obese == "no")
)

# initiate clusters for parallel computing
cl <- makeCluster(mc <- getOption("cl.cores", 4))
clusterExport(cl=cl, varlist=ls())
clusterCall(cl, function() library("quantreg"))

ABS <- parLapply(cl, DSETS, function(x) PLT.DAT(x, "VO2_MAX"))

REL <- parLapply(cl, DSETS, function(x) PLT.DAT(x, "REL_VO2_MAX"))

stopCluster(cl)

# saving data sets (long computing time)
saveRDS(ABS, "./results/dat/nomodata_abs.rds")
saveRDS(REL, "./results/dat/nomodata_rel.rds")


##################################################
# Nomograms: plotting
##################################################

# defining colors
COL  <- c("dodgerblue","dodgerblue1","tomato3","dodgerblue2","dodgerblue3", "dodgerblue4",
          "tomato3",
          "dodgerblue4", "dodgerblue3", "tomato3", "dodgerblue2", "dodgerblue1", "dodgerblue")

# function for fitting of ggplots
PLT <- function(data_male, data_female, rel = T, ex = F){
  
  # y-label
  if(rel == T){
    ylb <- REL_LAB
    yby <- 1
    dig <- 0
  } else {
    ylb <- ABS_LAB
    yby <- 0.1
    dig <- 1
  }
  
  # breaks are not dependent on a "from =" and "to =" argument
  brk_fun <- function(k) {
    step <- k
    function(y) seq(floor(min(y)), ceiling(max(y)), by = step)       
  }
  
  dm     <- data_male
  df     <- data_female
  dm$sex <- "Male"
  df$sex <- "Female"
  d      <- rbind(dm, df)
  
  ggp <- ggplot(data = d, aes(x = ALTER, y = fit)) +
    
    # defining background of ggplot
    theme(panel.grid.major = element_line(colour = "grey60", size = 0.3), 
          panel.background = element_rect(fill = "white"),
          axis.text.x      = element_text(angle = 90, vjust=0.5),
          plot.title       = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(25,69),
                       breaks = seq(from = 25, to = 69, by = 1))+  
    scale_y_continuous(breaks = brk_fun(yby)) +
    labs(x     = "Age [years]", 
         y     = ylb) +
    geom_vline(xintercept = c(30, 40, 50, 60), 
               colour     = "gray25", 
               size       = 0.6) +
    facet_wrap(~sex, scales = "free_y", ncol = 1)
  
  # adding confidence bands
  ggp <- ggp + 
    geom_ribbon(aes(x = ALTER, ymin = lower, ymax = higher, fill = percentile),
                alpha = 0.3, inherit.aes = F, show.legend = F, na.rm = T) +
    scale_fill_manual(values = COL)
  
  # adding quantile curves
  ggp <- ggp + 
    geom_line(aes(x = ALTER, y = fit, colour = percentile), 
              show.legend = F, size = 1, na.rm = T) +
    scale_color_manual(values = COL)
  
  # adding percentile labels to curves
  ggp <- direct.label(ggp, list("last.polygons", colour = "white", cex=0.65))
  
  return(ggp)
  
}

# saving plots

# Only Frankfurt and Ruedesheim (Munich is for external validation)
SAVE.PLOT(PLT(REL$dfm_FR, REL$dff_FR, rel = T, ex = F), "nomo_rel_include_FR.png", 10.5, 7)
SAVE.PLOT(PLT(ABS$dfm_FR, ABS$dff_FR, rel = F, ex = F), "nomo_abs_include_FR.png", 10.5, 7)
SAVE.PLOT(PLT(REL$exm_FR, REL$exf_FR, rel = T, ex = T), "nomo_rel_exclude_FR.png", 10.5, 7)
SAVE.PLOT(PLT(ABS$exm_FR, ABS$exf_FR, rel = F, ex = T), "nomo_rel_exclude_FR.png", 10.5, 7)

# Frankfurt, Ruedesheim and Munich included
SAVE.PLOT(PLT(REL$dfm, REL$dff, rel = T, ex = F), "nomo_rel_include_OVERALL.png", 10.5, 7)
SAVE.PLOT(PLT(ABS$dfm, ABS$dff, rel = F, ex = F), "nomo_abs_include_OVERALL.png", 10.5, 7)
SAVE.PLOT(PLT(REL$exm, REL$exf, rel = T, ex = T), "nomo_rel_exclude_OVERALL.png", 10.5, 7)
SAVE.PLOT(PLT(ABS$exm, ABS$exf, rel = F, ex = T), "nomo_rel_exclude_OVERALL.png", 10.5, 7)
dev.off()

##################################################
# Nomograms: Coefficient tables
##################################################

# calculating all coefficients and all 95% CIs for coefficients

COEF <- function(DAT, OUTCOME){
  
  # two empty lists to save coefficients and ci (lower and higher CI)
  COEF <- list()
  CIL  <- list()
  CIH  <- list()
  
  # loop over all taus and calculate predicted values
  for(i in TAUS){
    
    # name
    NAM   <- paste0("fit", i)
    
    # fit quantile regression
    FIT   <- rq(DAT[,OUTCOME] ~ ALTER + I(ALTER^2), data = DAT, tau = i)
    
    # bootstrapped confidence intervals for coefficients
    BT    <- boot.rq(cbind(1, DAT[,"ALTER"], I(DAT[,"ALTER"]^2)), 
                     DAT[, OUTCOME], 
                     tau = i, R = NBOOT)
    
    COEF[[NAM]]  <- FIT$coefficients
    CIL[[NAM]]   <- apply(BT$B, 2, quantile, c(0.025))
    CIH[[NAM]]   <- apply(BT$B, 2, quantile, c(0.975))
    
  }# End of loop
  
  return(list(COEF = COEF, CIL = CIL, CIH = CIH))
  
}

# calculating coefficients for all data sets

# initiate clusters for parallel computing
cl <- makeCluster(mc <- getOption("cl.cores", 4))
clusterExport(cl=cl, varlist=ls())
clusterCall(cl, function() library("quantreg"))

COEF.ABS <- parLapply(cl, DSETS, function(x) COEF(x, "VO2_MAX"))

COEF.REL <- parLapply(cl, DSETS, function(x) COEF(x, "REL_VO2_MAX"))

stopCluster(cl)

# saving data sets (long computing time)
saveRDS(COEF.ABS, "./results/dat/coef_abs.rds")
saveRDS(COEF.REL, "./results/dat/coef_rel.rds")


# creating tables for all data sets
COEF.TAB <- function(DAT){
  
  C <- round(t(do.call(data.frame, DAT$COEF)), digits = 7)
  L <- round(t(do.call(data.frame, DAT$CIL)), digits = 7)
  H <- round(t(do.call(data.frame, DAT$CIH)), digits = 7)
  
  # for coefficient plots
  TAB <- cbind(C[,1], L[,1], H[,1],
               C[,2], L[,2], H[,2],
               C[,3], L[,3], H[,3])
  
  TAB <- as.data.frame(TAB)
  
  # Coefficients + CI, I = Intercept, a = age, a2 = age^2
  # names are needed in this format for reshape using varying (coef plots)
  names(TAB) <- c(
    "c.1", "cil.1", "cih.1", "c.2", "cil.2", "cih.2", "c.3", "cil.3", "cih.3")
  
  # for Latex
  LTAB <- cbind(
    paste0(C[,1], " (", L[,1]," to ", H[,1], ")"),
    paste0(C[,2], " (", L[,2]," to ", H[,2], ")"),
    paste0(C[,3], " (", L[,3]," to ", H[,3], ")")
  )
  
  LTAB <- as.data.frame(LTAB)
  
  names(LTAB) <- c(
    "Intercept", "age", "age2")
  
  TAB$FIT        <- row.names(C)
  row.names(TAB) <- NULL
  
  return(list(LTAB = xtable(TAB), DATA = TAB))
  
}

# Create tables

# Absolute VO2max

# All cases
COEF.TAB(COEF.ABS$dfm)
COEF.TAB(COEF.ABS$dff)

# exclusion of smokers and obese
COEF.TAB(COEF.ABS$exm)
COEF.TAB(COEF.ABS$exf)

# Frankfurt/Ruedesheim (FR)
COEF.TAB(COEF.ABS$dfrm)
COEF.TAB(COEF.ABS$dfrf)

# Frankfurt/Ruedesheim (FR)
COEF.TAB(COEF.ABS$exrm)
COEF.TAB(COEF.ABS$exrf)

# Relative VO2max

# All cases
COEF.TAB(COEF.ABS$dfm)
COEF.TAB(COEF.ABS$dff)

# exclusion of smokers and obese
COEF.TAB(COEF.ABS$exm)
COEF.TAB(COEF.ABS$exf)

# Frankfurt/Ruedesheim (FR)
COEF.TAB(COEF.ABS$dfrm)
COEF.TAB(COEF.ABS$dfrf)

# Frankfurt/Ruedesheim (FR)
COEF.TAB(COEF.ABS$exrm)
COEF.TAB(COEF.ABS$exrf)

##################################################
# coefficient plots
##################################################

PLT.COEF <- function(DAT_MALE = COEF.ABS$dfm, DAT_FEMALE = COEF.ABS$dff){
  
  DM <- COEF.TAB(DAT_MALE)$DATA
  DF <- COEF.TAB(DAT_FEMALE)$DATA
  
  # adding sex-variable
  DM$Sex <- "Male"
  DF$Sex <- "Female"
  
  # Rbindung both sexes
  COEF.DATA <- rbind(DM, DF)
  
  # >COEF.DATA 
  # colnames are: c1 = Intercept, c2 = age, c3 = age^2
  # 
  #             c.1     cil.1    cih.1    c.2       cil.2        
  # fit0.05 2.148064 1.2871575 2.833846 0.0220968 -0.0071967 ...
  # fit0.1  2.283736 1.3166107 2.540710 0.0250659  0.0004518 ...
  # fit0.2  2.407500 1.2070498 2.863562 0.0309868 -0.0089181 ...
  # ...
  
  LONG      <- reshape(COEF.DATA, dir="long", varying = 1:9, idvar = c("FIT", "Sex"))
  LONG$time <- factor(LONG$time, levels = c(1,2,3), 
                      labels = c("Age",expression(Age^"2"), "Intercept"))
  
  # >LONG
  # FIT       Sex        time c          cil        cih
  # fit0.05   Male       Age  2.1480645  1.1823432  3.1150844
  # fit0.1    Male       Age  2.2837363  1.4854133  2.9921140
  # fit0.2    Male       Age  2.4075000  1.9149587  2.6529048
  # ...
  
  
  PLT <- ggplot(LONG) +
    geom_point(    aes(x = FIT, y = c)) +
    geom_linerange(aes(x = FIT, ymin = cil, ymax = cih)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(Sex ~ time, scales = "free", labeller = label_parsed) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = "Coefficient (95% CI)", x = "Quantile regression")
  
  return(PLT)
}

# Create plots

# Absolute VO2max

# All cases
SAVE.PLOT(PLT.COEF(COEF.ABS$dfm,COEF.ABS$dff), "coefplot_abs_OVERALL.png", 4, 8)
# exclusion of smokers and obese
SAVE.PLOT(PLT.COEF(COEF.ABS$exm,COEF.ABS$exf), "coefplot_abs_EX.png", 4, 8)
# Frankfurt/Ruedesheim (FR)
SAVE.PLOT(PLT.COEF(COEF.ABS$dfm_FR,COEF.ABS$dff_FR), "coefplot_abs_FR.png", 4, 8)

# Relative VO2max

# All cases
SAVE.PLOT(PLT.COEF(COEF.REL$dfm,COEF.REL$dff), "coefplot_rel_OVERALL.png", 4, 8)
# exclusion of smokers and obese
SAVE.PLOT(PLT.COEF(COEF.REL$exm,COEF.REL$exf), "coefplot_rel_EX.png", 4, 8)
# Frankfurt/Ruedesheim (FR)
SAVE.PLOT(PLT.COEF(COEF.REL$dfm_FR,COEF.REL$dff_FR), "coefplot_rel_FR.png", 4, 8)

##################################################
##################################################
# Validation
##################################################
##################################################

##################################################
# calculating validation data
##################################################

VALI.DAT <- function(TRN.DAT = dfm_FR, VAL.DAT = dfm_FR, OUTCOME = "REL_VO2_MAX",
                     FORMULA = "ageclass + I(ageclass^2)"){
  
  # GOALS: producing all data for validation
  
  # G1) empirical quantiles (EMP.QTLS) from validation data (VAL.DAT)
  # G2) quantile regression in trainaing data (TRN.DAT)
  # G3) predicted values from quantile regressions
  # G4) 95% CI for predicted values
  # G5) AIC table for each TAUS2
  # G6) comparing empirical and estimated values in separate function
  
  # quantiles for validation
  TAUS2 <- c(0.25,0.5,0.75)
  
  # GOAL 1
  
  # empirical quantiles: tau = {0.25,0.5,0.75}
  FRM1     <- as.formula(paste(OUTCOME, "~ ageclass + Geschlecht"))
  EMP.QTLS <- aggregate(FRM1, data = VAL.DAT,
                        FUN = function(x) quantile(x, probs = TAUS2))
  EMP.QTLS <- do.call(data.frame, EMP.QTLS)
  names(EMP.QTLS) <- c("ageclass", "Geschlecht", "tau0.25", "tau0.5", "tau0.75")
  EMP.QTLS <- melt(EMP.QTLS, id.vars = c("ageclass", "Geschlecht"))
  
  # > EMP.QTLS
  # ageclass Geschlecht variable  value
  # [25,30)       Male  tau0.25 33.525
  # [30,35)       Male  tau0.25 33.400
  # [35,40)       Male  tau0.25 32.700       
  # ...
  
  # GOAL 2
  
  # ageclass variable (from trainaing data) for regression
  TRN.DAT$ageclass <- as.numeric(TRN.DAT$ageclass)
  
  # formula for quantile regressions
  FRM2 <- as.formula(paste(OUTCOME, "~", FORMULA))
  
  PDAT     <- list()
  AIC.DAT  <- list()
  
  # loop over all taus and calculate predicted values
  for(i in TAUS2){
    
    FIT   <- rq(FRM2, data = TRN.DAT, tau = i)
    
    # create empty prediction data set
    pdat  <- data.frame(ageclass = 1:9)
    
    # function to predict value and 95% CI
    PRED <- function(){
      
      set.seed(1)
      PV <- predict.rq(FIT, 
                       newdata = pdat, 
                       interval = "confidence", 
                       se = "boot", type = "percentile",
                       level = .95, R = NBOOT)
      return(PV)
    }
    
    # calculate predicted values and confidence intervals, and adding tau
    pdat  <- cbind(pdat, 
                   tau = paste0("tau",i),
                   PRED()
    )
    # > pdat
    # ageclass   tau      fit    lower   higher
    #       1   tau0.50 39.30000 38.31400 41.09600
    #       2   tau0.50 38.83333 38.24494 39.88267
    #       3   tau0.50 38.00000 37.70000 38.49750
    # ...
    
    # column must be the same name as outcome
    # in order to add the dots in the box plot
    pdat[,OUTCOME] <- pdat$fit
    
    # saving prediction data
    nam1 <- paste0("pdat", i)
    PDAT[[nam1]] <- pdat
    
    # saving aic data
    nam2 <- paste0("aic", i)
    AIC.DAT[[nam2]] <- AIC.rq(FIT)
    
  }
  
  PDAT <- do.call(rbind, PDAT)
  row.names(PDAT) <- NULL
  
  AIC.DAT <- do.call(rbind, AIC.DAT)
  # > test$AIC.DAT
  #         [,1]
  # aic0.25 32672.15
  # aic0.5  32406.82
  # aic0.75 33242.66
  
  return(list(EMP.QTLS = EMP.QTLS, PDAT = PDAT, AIC.DAT = AIC.DAT))
}

# saving all data
# Nomenclature: 
# A = a = apparent, e = external
# M = male, F = female
# R = REL_VO2_MAX, A = ABS_VO2_MAX

# __________________________________
# Main models: Quadratic
# __________________________________

# apparent validation
AMR <- VALI.DAT(dfm_FR, dfm_FR, "REL_VO2_MAX", "ageclass + I(ageclass^2)")
AFR <- VALI.DAT(dff_FR, dff_FR, "REL_VO2_MAX", "ageclass + I(ageclass^2)")
AMA <- VALI.DAT(dfm_FR, dfm_FR, "VO2_MAX",     "ageclass + I(ageclass^2)")
AFA <- VALI.DAT(dff_FR, dff_FR, "VO2_MAX",     "ageclass + I(ageclass^2)")

# external validation
EMR <- VALI.DAT(dfm_FR, dfm_M, "REL_VO2_MAX", "ageclass + I(ageclass^2)")
EFR <- VALI.DAT(dff_FR, dff_M, "REL_VO2_MAX", "ageclass + I(ageclass^2)")
EMA <- VALI.DAT(dfm_FR, dfm_M, "VO2_MAX",     "ageclass + I(ageclass^2)")
EFA <- VALI.DAT(dff_FR, dff_M, "VO2_MAX",     "ageclass + I(ageclass^2)")

# __________________________________
# Other models: Linear
# __________________________________

# apparent validation
AMR_lin <- VALI.DAT(dfm_FR, dfm_FR, "REL_VO2_MAX", "ageclass")
AFR_lin <- VALI.DAT(dff_FR, dff_FR, "REL_VO2_MAX", "ageclass")
AMA_lin <- VALI.DAT(dfm_FR, dfm_FR, "VO2_MAX",     "ageclass")
AFA_lin <- VALI.DAT(dff_FR, dff_FR, "VO2_MAX",     "ageclass")

# __________________________________
# Other models: spline
# __________________________________

# apparent validation
AMR_spl <- VALI.DAT(dfm_FR, dfm_FR, "REL_VO2_MAX", "bs(ageclass)")
AFR_spl <- VALI.DAT(dff_FR, dff_FR, "REL_VO2_MAX", "bs(ageclass)")
AMA_spl <- VALI.DAT(dfm_FR, dfm_FR, "VO2_MAX",     "bs(ageclass)")
AFA_spl <- VALI.DAT(dff_FR, dff_FR, "VO2_MAX",     "bs(ageclass)")

##################################################
# validation box plots
##################################################

# function for adding the estimates or CI to the box plot
ADD.LINES <- function(DAT = DAT, COLOUR = "blue", CI = F, ALPHA = 1){
  
  # adding confidence bands  
  if(CI == T){
    BPL <-  BPL + 
      geom_linerange(data = DAT, 
                     aes(x = ageclass, ymin = lower, ymax = higher, group = tau),
                     size = 7, colour = "red", alpha = 0.5)
  }
  
  BPL <- BPL +
    geom_line(data = DAT,
              stat = "smooth", method = "loess",
              aes(x = ageclass, y = fit, group = tau),
              size = 0.5, colour = COLOUR, alpha = ALPHA,
              na.rm = T) +
    geom_point(data = DAT, aes(x = ageclass, y = fit, shape = factor(tau)),
               size = 2, colour = COLOUR, alpha = ALPHA) +
    theme(legend.position="none")
  
  return(BPL)
} 

# __________________________________
# Apparent validation
# __________________________________

# Apparent validation, rel, males
BPL <- PLT.BPL(dfm_FR, REL = T)
BPL1 <- ADD.LINES(DAT = AMR$PDAT, COLOUR = "blue", CI = T)

# Apparent validation, rel, females
BPL <- PLT.BPL(dff_FR, REL = T)
BPL2 <- ADD.LINES(DAT = AFR$PDAT, COLOUR = "blue", CI = T)

SAVE.PLOT(grid.arrange(BPL1,BPL2),
          "vali_apparent_rel.png", 9, 6)


# Apparent validation, abs, males
BPL <- PLT.BPL(dfm_FR, REL = F)
BPL3 <- ADD.LINES(DAT = AMA$PDAT, COLOUR = "blue", CI = T)

# Apparent validation, abs, females
BPL <- PLT.BPL(dff_FR, REL = F)
BPL4 <- ADD.LINES(DAT = AFA$PDAT, COLOUR = "blue", CI = T)

SAVE.PLOT(grid.arrange(BPL3,BPL4),
          "vali_apparent_abs.png", 9, 6)

# __________________________________
# External validation
# __________________________________

# Apparent validation, rel, males
BPL <- PLT.BPL(dfm_M, REL = T)
BPL1 <- ADD.LINES(DAT = EMR$PDAT, COLOUR = "blue", CI = T)

# Apparent validation, rel, females
BPL <- PLT.BPL(dff_M, REL = T)
BPL2 <- ADD.LINES(DAT = EFR$PDAT, COLOUR = "blue", CI = T)

SAVE.PLOT(grid.arrange(BPL1,BPL2),
          "vali_external_rel.png", 9, 6)


# Apparent validation, abs, males
BPL <- PLT.BPL(dfm_M, REL = F)
BPL3 <- ADD.LINES(DAT = EMA$PDAT, COLOUR = "blue", CI = T)

# Apparent validation, abs, females
BPL <- PLT.BPL(dff_M, REL = F)
BPL4 <- ADD.LINES(DAT = EFA$PDAT, COLOUR = "blue", CI = T)

SAVE.PLOT(grid.arrange(BPL3,BPL4),
          "vali_external_abs.png", 9, 6)


# __________________________________
# Apparent validation including linear and spline models
# __________________________________

# Apparent validation, rel, males
BPL <- PLT.BPL(dfm_FR, REL = T)
BPL <- ADD.LINES(DAT = AMR_lin$PDAT, COLOUR = "red", CI = F, ALPHA = 0.5)
BPL <- ADD.LINES(DAT = AMR$PDAT,     COLOUR = "blue", CI = F, ALPHA = 0.5)
BPL1 <- ADD.LINES(DAT = AMR_spl$PDAT, COLOUR = "green", CI = F, ALPHA = 0.5)

# Apparent validation, rel, females
BPL <- PLT.BPL(dff_FR, REL = T)
BPL <- ADD.LINES(DAT = AFR_lin$PDAT, COLOUR = "red", CI = F, ALPHA = 0.5)
BPL <- ADD.LINES(DAT = AFR$PDAT,     COLOUR = "blue", CI = F, ALPHA = 0.5)
BPL2 <- ADD.LINES(DAT = AFR_spl$PDAT, COLOUR = "green", CI = F, ALPHA = 0.5)

SAVE.PLOT(grid.arrange(BPL1,BPL2),
          "vali_apparent_rel_comparison.png", 9, 6)


# Apparent validation, abs, males
BPL <- PLT.BPL(dfm_FR, REL = F)
BPL <- ADD.LINES(DAT = AMA_lin$PDAT, COLOUR = "red", CI = F, ALPHA = 0.5)
BPL <- ADD.LINES(DAT = AMA$PDAT,     COLOUR = "blue", CI = F, ALPHA = 0.5)
BPL3 <- ADD.LINES(DAT = AMA_spl$PDAT, COLOUR = "green", CI = F, ALPHA = 0.5)

# Apparent validation, abs, females
BPL <- PLT.BPL(dff_FR, REL = F)
BPL <- ADD.LINES(DAT = AFA_lin$PDAT, COLOUR = "red", CI = F, ALPHA = 0.5)
BPL <- ADD.LINES(DAT = AFA$PDAT,     COLOUR = "blue", CI = F, ALPHA = 0.5)
BPL4 <- ADD.LINES(DAT = AFA_spl$PDAT, COLOUR = "green", CI = F, ALPHA = 0.5)

SAVE.PLOT(grid.arrange(BPL3,BPL4),
          "vali_apparent_abs_comparison.png", 9, 6)

##################################################
# validation calibration plots
##################################################

# build basic calibration plot without lines
PLT.CALI <- function(REL = T){
  
  if(REL == T){ MIN = 20; MAX = 50; BY = 5  }else{
    MIN = 1 ; MAX = 4 ; BY = 0.5}
  
  P <- ggplot(data = D, aes(x = fit, y = value)) +
    geom_abline(intercept = 0, slope = 1, lty = 3) +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(MIN,MAX), breaks = seq(MIN, MAX, by = BY)) +
    scale_y_continuous(limits = c(MIN,MAX), breaks = seq(MIN, MAX, by = BY)) +
    labs(x = "Predicted", y = "Observed") 
  
  return(P)
}

# add regression lines and dots
LINE.CALI <- function(X = AMR$PDAT, Y = AMR$EMP.QTLS, PLOT = CP, COLOUR = "blue", ALPHA = 1){
  
  Y$ageclass <- as.numeric(Y$ageclass)
  
  D <- merge(X, Y, 
             by.x = c("ageclass", "tau"),
             by.y = c("ageclass", "variable"))
  
  CP <- CP + 
    geom_line(stat = "smooth", method = "lm", se = F,
              data = D, 
              aes(x = fit, y = value), 
              colour = COLOUR, alpha = ALPHA, na.rm = T) +
    geom_point(data = D, aes(x = fit, y = value, shape = tau),
               colour = COLOUR, alpha = ALPHA, na.rm = T) +
    facet_grid(~Geschlecht)
  
  return(CP)
  
}

# __________________________________
# Apparent validation 
# __________________________________

CP  <- PLT.CALI(REL = T)
CP1 <- LINE.CALI(AMR$PDAT, AMR$EMP.QTLS)
CP2 <- LINE.CALI(AFR$PDAT, AFR$EMP.QTLS)
CP  <- PLT.CALI(REL = F)
CP3 <- LINE.CALI(AMA$PDAT, AMA$EMP.QTLS)
CP4 <- LINE.CALI(AFA$PDAT, AFA$EMP.QTLS)
SAVE.PLOT(grid.arrange(
  arrangeGrob(CP1,CP2, top = textGrob(REL_LAB)), 
  arrangeGrob(CP3,CP4, top = textGrob(ABS_LAB)),
  ncol = 2), "cali_apparent.png", 7,6)

# __________________________________
# External validation 
# __________________________________
CP  <- PLT.CALI(REL = T)
CP1 <- LINE.CALI(EMR$PDAT, EMR$EMP.QTLS)
CP2 <- LINE.CALI(EFR$PDAT, EFR$EMP.QTLS)
CP  <- PLT.CALI(REL = F)
CP3 <- LINE.CALI(EMA$PDAT, EMA$EMP.QTLS)
CP4 <- LINE.CALI(EFA$PDAT, EFA$EMP.QTLS)
SAVE.PLOT(grid.arrange(
  arrangeGrob(CP1,CP2, top = textGrob(REL_LAB)), 
  arrangeGrob(CP3,CP4, top = textGrob(ABS_LAB)),
  ncol = 2), "cali_external.png", 7,6)

# __________________________________
# Apparent validation comparison of regression models 
# __________________________________
CP  <- PLT.CALI(REL = T) # Males, REL_VO2_MAX
CP  <- LINE.CALI(AMR$PDAT,     AMR$EMP.QTLS    , ALPHA = 0.5, COLOUR = "blue")
CP  <- LINE.CALI(AMR_lin$PDAT, AMR_lin$EMP.QTLS, ALPHA = 0.5, COLOUR = "green")
CP1 <- LINE.CALI(AMR_spl$PDAT, AMR_spl$EMP.QTLS, ALPHA = 0.5, COLOUR = "red")
CP  <- PLT.CALI(REL = T) # Females REL_VO2_MAX
CP  <- LINE.CALI(AFR$PDAT,     AFR$EMP.QTLS    , ALPHA = 0.5, COLOUR = "blue")
CP  <- LINE.CALI(AFR_lin$PDAT, AFR_lin$EMP.QTLS, ALPHA = 0.5, COLOUR = "green")
CP2 <- LINE.CALI(AFR_spl$PDAT, AFR_spl$EMP.QTLS, ALPHA = 0.5, COLOUR = "red")
CP  <- PLT.CALI(REL = F) # Males VO2_MAX
CP  <- LINE.CALI(AMA$PDAT,     AMA$EMP.QTLS    , ALPHA = 0.5, COLOUR = "blue")
CP  <- LINE.CALI(AMA_lin$PDAT, AMA_lin$EMP.QTLS, ALPHA = 0.5, COLOUR = "green")
CP3 <- LINE.CALI(AMA_spl$PDAT, AMA_spl$EMP.QTLS, ALPHA = 0.5, COLOUR = "red")
CP  <- PLT.CALI(REL = F) # Females VO2_MAX
CP  <- LINE.CALI(AFA$PDAT,     AFA$EMP.QTLS    , ALPHA = 0.5, COLOUR = "blue")
CP  <- LINE.CALI(AFA_lin$PDAT, AFA_lin$EMP.QTLS, ALPHA = 0.5, COLOUR = "green")
CP4 <- LINE.CALI(AFA_spl$PDAT, AFA_spl$EMP.QTLS, ALPHA = 0.5, COLOUR = "red")
SAVE.PLOT(grid.arrange(
  arrangeGrob(CP1,CP2, top = textGrob(REL_LAB)), 
  arrangeGrob(CP3,CP4, top = textGrob(ABS_LAB)),
  ncol = 2), "cali_apparent_comparison.png", 7,6)

##################################################
# validation AIC tables
##################################################

# Males
AIC.TAB <- rbind(
  
  # females
  cbind(AFR_lin$AIC.DAT, AFR_spl$AIC.DAT, AFR$AIC.DAT,
        AFA_lin$AIC.DAT, AFA_spl$AIC.DAT, AFA$AIC.DAT),
  
  # males
  cbind(AMR_lin$AIC.DAT, AMR_spl$AIC.DAT, AMR$AIC.DAT,
        AMA_lin$AIC.DAT, AMA_spl$AIC.DAT, AMA$AIC.DAT)
)

AIC.TAB <- as.data.frame(AIC.TAB)
names(AIC.TAB)     <- c("AL","AS","AP","RL","RS","RP")
row.names(AIC.TAB) <- c("F.25","F.5","F.75","M.25","M.5","M.75")
xtable(AIC.TAB)

##################################################
# validation coefficient tables
##################################################

COEF.TAB2 <- function(DAT){
  
  X <- DAT["PDAT"]; X <- X$PDAT$fit
  Y <- DAT["EMP.QTLS"]; Y <- Y$EMP.QTLS$value
  
  D <- data.frame(X,Y)
  
  COEF.BT <- function(data, INDEX){
    D1  <- data[INDEX,]  
    FIT <- lm(Y ~ X, data = D1)
    ROW <- c(FIT$coefficients, summary(FIT)$r.squared)
  }
  
  BT  <- boot(data = D, statistic = COEF.BT, R = NBOOT)
  TAB <- rbind(BT$t0, apply(BT$t, 2, quantile, c(0.025,0.975)))
  TAB <- round(t(as.data.frame(TAB)), digits = 2)
  COL <- paste0(TAB[,1], " [", TAB[,2], " to ", TAB[,3], "]")
  COL <- data.frame(COL = COL)
  row.names(COL) <- c("Intercept", "Slope", "Rsquared")
  
  #> COL
  #Intercept 5.94 [3.11 to 8.69]
  #Slope     0.83 [0.76 to 0.90]
  #Rsquared  0.92 [0.88 to 0.96]
  
  return(COL)
}

# __________________________________
# Apparent validation
# __________________________________

# CAVE: females on top, absolute on the left column

TAB <- rbind(
  cbind(COEF.TAB2(AFA), COEF.TAB2(AFR)),
  cbind(COEF.TAB2(AMA), COEF.TAB2(AMR))
)
names(TAB) <- c("Absolute", "Relative")
xtable(TAB)

# __________________________________
# External validation
# __________________________________

TAB <- rbind(
  cbind(COEF.TAB2(EFA), COEF.TAB2(EFR)),
  cbind(COEF.TAB2(EMA), COEF.TAB2(EMR))
)
names(TAB) <- c("Absolute", "Relative")
xtable(TAB)

# __________________________________
# Apparent validation comparison of regression models
# __________________________________

TAB <- rbind(
  # female
  cbind(# absolute
    COEF.TAB2(AFA_lin), COEF.TAB2(AFA_spl), COEF.TAB2(AFA),
    # relative
    COEF.TAB2(AFR_lin), COEF.TAB2(AFR_spl), COEF.TAB2(AFR)),
  # male
  cbind(# absolute
    COEF.TAB2(AMA_lin), COEF.TAB2(AMA_spl), COEF.TAB2(AMA),
    # relative
    COEF.TAB2(AMR_lin), COEF.TAB2(AMR_spl), COEF.TAB2(AMR))
)
names(TAB) <- c("AL", "AS", "AP", "RL", "RS", "RP")
xtable(TAB)

##################################################
##################################################
# Multiple quantile regression
##################################################
##################################################

# histogram rubust for outliers

h <- function(var){
  
  # set margins of histogram
  m <- mean(df[,var], na.rm = T)
  s <- sd(df[,var], na.rm = T)
  
  # log transformed variable
  d   <- df
  d$l <- log(df[,var])
  ml  <- mean(d$l, na.rm = T)
  sl  <- sd(d$l, na.rm = T)
  
  # histogram with crude variable
  c <- ggplot(df, aes_string(x = var)) +
    geom_histogram(aes(y = ..density..), bins = 20, na.rm = T) +
    scale_x_continuous(limits = c(m - 3*s, m + 3*s)) +
    facet_wrap(~Geschlecht) +
    labs(title = "crude")
  
  # histogram with log-transformed variable
  l <- ggplot(d, aes(x = l)) +
    geom_histogram(aes(y = ..density..), bins = 20, na.rm = T) +
    scale_x_continuous(limits = c(ml - 3*sl, ml + 3*sl)) +
    facet_wrap(~Geschlecht) +
    labs(title = "log-transformed")
  
  gg <- grid.arrange(c,l)
  return(gg)
}

h("BMI") # --> transform
h("GEWICHT") # --> transform
h("TRIGLYCERIDE") # --> transform
h("BAUCHUMFANG") # transform
h("BLUTZUCKER") # --> transform
h("HBA1C") # --> transform

# checking correlation of variables
# remove when r>=0.75 to avoid collinearity
# based on SPLOM
# --> remove: 
# GEWICHT, BAUCHUMFANG, KOERPERFETT_CALIPER & TANITA,
# GESAMTCHOLESTERIN

# calculating new variables
df$Age_squared        = df$ALTER^2
df$log_BMI            = log(df$BMI)
df$log_triglycerides  = log(df$TRIGLYCERIDE)
df$log_blood_glucose  = log(df$BLUTZUCKER)
df$log_HbA1c          = log(df$HBA1C)

# renaming variables (for tick labels)
df$Height    <- df$GROESSE
df$Age       <- df$ALTER
df$Systolic  <- df$SYST
df$Diastolic <- df$DIAST 
df$HDL_cholesterol <- df$HDL_CHOLESTERIN
df$LDL_cholesterol <- df$LDL_CHOLESTERIN

# recoding factor variables 
# --> necessary for following function
# --> coding as no = 1, yes = 2
BIN <- function(x){
  var <- as.numeric(df[,x])
  var <- ifelse(var == 1, 0, ifelse(var == 2, 1, NA))
  return(var)
}
df$Smoker            <- BIN("ZIGARETTEN")
df$Exsmoker          <- BIN("EXRAUCHER")
df$Overweight        <- BIN("BMI_GE_25")
df$Obese             <- BIN("obese")
df$Diabetes_mellitus <- BIN("MANIFEST_DM")

# Defining final independent variables
IV <- c(
  "Age", "Height", 
  "Overweight", "Obese",
  "Systolic", "Diastolic", "HDL_cholesterol", "LDL_cholesterol", "log_triglycerides",
  "log_blood_glucose", "log_HbA1c", "Diabetes_mellitus", "Smoker", "Exsmoker"
)

# removing all cases with NAs in data set
df$exclusion <- rowSums(is.na(df[,c(IV,"VO2_MAX", "REL_VO2_MAX")]))
df1 <- subset(df, exclusion == 0)

# subsetting data including new IVs
dfm <- subset(df1, Geschlecht == "Male")
dff <- subset(df1, Geschlecht == "Female")

# Regression fitting

QR <- function(dat, outcome, qtile) {
  
  # creating formula
  # IV have to be defined!
  frm <- as.formula(paste(outcome ,"~", paste(IV, collapse = "+")))
  
  # fitting quantile regression
  fit <- step(rq(frm, data = dat, tau = qtile), direction = "both")
  
  # bootstrapped P values
  # NBOOT has to be defined!
  set.seed(1)
  tab <- summary(fit, se = "boot", R = NBOOT)
  
  # extract all variables that are in the formula of final model and remove dependent var 
  var <- setdiff(all.vars(tab$call$formula), outcome)
  
  # bottstrapping 95% confidence intervals
  set.seed(1)
  bt <- boot.rq(cbind(1, dat[, var]), dat[, outcome], tau = qtile, R = NBOOT)
  
  # add bootstrapped 95% CI to table
  tab <- cbind(
    tab$coefficients[, c(1,4)],
    t(apply(bt$B, 2, quantile, c(0.025,0.975))))
  
  # are coefficients within ci bounds?
  tab <- as.data.frame(tab)
  tab <- tab[, c("Value", "2.5%", "97.5%", "Pr(>|t|)")]
  
  coef_in_ci <- data.frame(
    hi = tab$Value<tab$`97.5%`,
    lo = tab$Value>tab$`2.5%`)
  
  # calculation of R^2 according to
  # Hao, 2007, p. 52
  # https://stat.ethz.ch/pipermail/r-help/2006-August/110386.html
  # http://stats.stackexchange.com/questions/129200/r-squared-in-quantile-regression
  fit0 <- rq(REL_VO2_MAX ~ 1, tau = qtile, data = dat)
  rho <- function(u,tau=qtile)u*(tau - (u < 0))
  R1 <- round(1 - fit$rho/fit0$rho, digits = 3)
  
  # creating name for coefficient plot
  # sex
  if(dat[3,"Geschlecht"] == "Male") sex <- "Males"
  if(dat[3,"Geschlecht"] == "Female") sex <- "Females"
  
  # return the name VO2_MAX
  out <- colnames(dat)[which(colnames(dat) %in% outcome)]
  
  # build name
  if(out == "REL_VO2_MAX") {
    lab <- bquote(atop(.(paste0(sex))*";"
                       ~Dependent~variable~"="
                       ~Relative~dot('V')*O['2peak']~ "[mL" * O["2"]*"/min/kg"*"]",
                       Quantile~"="~.(qtile)*";"
                       ~Pseudo*"-"*R^2~"="~ .(R1)))}
  if(out == "VO2_MAX") {
    lab <- bquote(atop(.(paste0(sex))*";"
                       ~Dependent~variable~"="
                       ~Absolute~dot('V')*O['2peak']~ "[L" * O["2"]*"/min"*"]",
                       Quantile~"="~.(qtile)*";"
                       ~Pseudo*"-"*R^2~"="~ .(R1)))}
  
  # add variable names
  name <- rownames(tab)
  tab2 <- tab
  tab2$name <- name
  
  # coefficient plot
  plt <-ggplot(data = tab2[2:nrow(tab2),]) +
    geom_point(aes(x = name, y = Value)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(linetype = "dotted", colour = "grey60", size = 0.3),
          plot.title = element_text(size = 10)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_linerange(aes(x = name, ymin = `2.5%`, ymax = `97.5%`)) +
    labs(x = "", y = "") +
    ggtitle(lab) +
    scale_y_continuous(expand = c(0.1,0.1), breaks = ) +
    coord_flip()
  
  l <- list(tab, plt, R1, coef_in_ci)
  
  return(l)}

# relative VO2peak

# males
rm25 <- QR(dfm, "REL_VO2_MAX", 0.25)
rm50 <- QR(dfm, "REL_VO2_MAX", 0.50)
rm75 <- QR(dfm, "REL_VO2_MAX", 0.75)

# females
rf25 <- QR(dff, "REL_VO2_MAX", 0.25)
rf50 <- QR(dff, "REL_VO2_MAX", 0.50)
rf75 <- QR(dff, "REL_VO2_MAX", 0.75)

# arranging coefficient plots
rm <- grid.arrange(rm25[[2]], rm50[[2]], rm75[[2]], ncol = 1)
rf <- grid.arrange(rf25[[2]], rf50[[2]], rf75[[2]], ncol = 1)

# relative VO2peak

# males
am25 <- QR(dfm, "VO2_MAX", 0.25)
am50 <- QR(dfm, "VO2_MAX", 0.50)
am75 <- QR(dfm, "VO2_MAX", 0.75)

# females
af25 <- QR(dff, "VO2_MAX", 0.25)
af50 <- QR(dff, "VO2_MAX", 0.50)
af75 <- QR(dff, "VO2_MAX", 0.75)

# arranging coefficient plots
am <- grid.arrange(am25[[2]], am50[[2]], am75[[2]], ncol = 1)
af <- grid.arrange(af25[[2]], af50[[2]], af75[[2]], ncol = 1)

# saving plots
SAVE.PLOT(rm, "multi_relative_m.png", 8,6)
SAVE.PLOT(rf, "multi_relative_f.png", 8,6)
SAVE.PLOT(am, "multi_absolute_m.png", 8,6)
SAVE.PLOT(af, "multi_absolute_f.png", 8,6)

# tables
print(xtable(rbind(rm25[[1]], rm50[[1]], rm75[[1]]), digits=c(3,3,3,3,3)))
print(xtable(rbind(rf25[[1]], rf50[[1]], rf75[[1]]), digits=c(3,3,3,3,3)))
print(xtable(rbind(am25[[1]], am50[[1]], am75[[1]]), digits=c(3,3,3,3,3)))
print(xtable(rbind(af25[[1]], af50[[1]], af75[[1]]), digits=c(3,3,3,3,3)))

##################################################
##################################################
# Comparison of reference values
##################################################
##################################################  

# SHIP Study (Koch, 2009)

# CODING
# m=1, f=2
# age groups: 25-34, 35-44, 45-54, 55-64, >=64
# bmi: <25 = 0, >=25 = 1

# VO2peak 5%
v.05 <- function(sex, Age, bmi){
  p <- 30.9643 + 
    (-2.5661*Age) +
    (-0.0263*I(Age^2)) +
    (-3.7224*sex) +
    (1.8765*bmi) + 
    (0.1082*Age*bmi) +
    (-2.9703*sex*bmi) +
    (0.7361*Age*sex) +
    (0.2799*Age*sex*bmi)
  return(p)
}

v.95 <- function(sex, Age, bmi){
  p <- 61.3721 + 
    (-1.9479*Age) +
    (-0.3053*I(Age^2)) +
    (-9.1229*sex) +
    (3.8892*bmi) + 
    (-1.9492*Age*bmi) +
    (-6.7455*sex*bmi) +
    (0.0716*Age*sex) +
    (1.6900*Age*sex*bmi)
  return(p)
}

v.50 <- function(sex, Age, bmi){
  p <- 47.7565 + 
    (-0.9880*Age) +
    (-0.2356*I(Age^2)) +
    (-8.8697*sex) +
    (2.3597*bmi) + 
    (-2.0308*Age*bmi) +
    (-3.7405*sex*bmi) +
    (0.2512*Age*sex) +
    (1.3797*Age*sex*bmi)
  return(p)
}

# prediction data set
ship_m <- data.frame(Age = 1:5)
ship_f <- data.frame(Age = 1:5)
own_m  <- data.frame(Age = 25:70)
own_f  <- data.frame(Age = 25:70)

# males
ship_m$m.05 <- v.05(1, ship_m$Age, 0)
ship_m$m.50 <- v.50(1, ship_m$Age, 0)
ship_m$m.95 <- v.95(1, ship_m$Age, 0)

# females
ship_f$f.05 <- v.05(2, ship_f$Age, 0)
ship_f$f.50 <- v.50(2, ship_f$Age, 0)
ship_f$f.95 <- v.95(2, ship_f$Age, 0)

ship_m$Age <- ship_f$Age <- factor(
  ship_m$Age,
  levels = 1:5, 
  labels = c("25-34", "35-44", "45-54", "55-64", ">=64"))

ship_m_long <- melt(ship_m, id.vars = "Age")
ship_f_long <- melt(ship_f, id.vars = "Age")

# FRIEND Study (Kaminsky, 2015)

# Age groups (6): 20-29, 30-39, 40-49, 50-59, 60-69, 70-79
kaminsky_m <- data.frame(
  Age = rep(1:6, 3),
  variable = c(rep("m.05",6), rep("m.50",6), rep("m.95",6)),
  value = c(
    # m.50
    41.9, 30.1, 27.1, 24.8, 22.4, 19.5))

kaminsky_f <- data.frame(
  Age = rep(1:6, 3),
  variable = c(rep("f.05",6), rep("f.50",6), rep("f.95",6)),
  value = c(
    # f.50
    31.0, 21.6, 19.4, 17.3, 16.0, 14.8))

# recoding age
kaminsky_m$Age <- kaminsky_f$Age <- factor(
  kaminsky_m$Age,
  levels = 1:6,
  labels = c("20-29", "30-30", "40-49", "50-59", "60-69", "70-79"))

# own

# males
fitm.05 <- rq(REL_VO2_MAX ~ Age + I(Age^2), data = dfm, tau = 0.05)
fitm.50 <- rq(REL_VO2_MAX ~ Age + I(Age^2), data = dfm, tau = 0.50)
fitm.95 <- rq(REL_VO2_MAX ~ Age + I(Age^2), data = dfm, tau = 0.95)

own_m$m.05 <- predict(fitm.05, newdata = own_m) 
own_m$m.50 <- predict(fitm.50, newdata = own_m) 
own_m$m.95 <- predict(fitm.95, newdata = own_m)

own_m_long <- melt(own_m, id.vars = "Age")

# females
fitf.05 <- rq(REL_VO2_MAX ~ Age + I(Age^2), data = dff, tau = 0.05)
fitf.50 <- rq(REL_VO2_MAX ~ Age + I(Age^2), data = dff, tau = 0.50)
fitf.95 <- rq(REL_VO2_MAX ~ Age + I(Age^2), data = dff, tau = 0.95)

own_f$f.05 <- predict(fitf.05, newdata = own_f) 
own_f$f.50 <- predict(fitf.50, newdata = own_f) 
own_f$f.95 <- predict(fitf.95, newdata = own_f)

own_f_long <- melt(own_f, id.vars = "Age")

g <- function(dat, tit, xtit, ytit, class = T){
  # plot
  gp <- ggplot(data = dat, aes(x = Age, y = value, group = factor(variable))) +
    geom_line() +
    theme(panel.grid.major = element_line(colour = "grey60", size = 0.3),
          panel.background = element_rect(fill = "white")) +
    scale_y_continuous(limits = c(10,55), breaks = seq(from=10, to = 55, by = 5)) +
    labs(title = tit, x = xtit, y = ytit)
  
  # modifying x axis
  if(class == F){
    gp <- gp + scale_x_continuous(limits = c(25,70), breaks = seq(from= 25, to = 70, by = 5))
  }else{
    gp <- gp + geom_point()
  }
  
  return(gp)
}

# plotting males
g1 <- g(own_m_long, "Prevention First", "Age [years]", REL_LAB, class = F)
g2 <- g(ship_m_long, "SHIP Study", "Age class [years]", "")
g3 <- g(kaminsky_m, "FRIEND Study", "Age class [years]", "")

# agganging plots
gm <- grid.arrange(g1, g2, g3, ncol = 3, top = "Males")

# plotting females
g1 <- g(own_f_long, "Prevention First", "Age [years]", REL_LAB, class = F)
g2 <- g(ship_f_long, "SHIP Study", "Age class [years]", "")
g3 <- g(kaminsky_f, "FRIEND Study", "Age class [years]", "")

# agganging plots
gf <- grid.arrange(g1, g2, g3, ncol = 3, top = "Females")

# saving plots
SAVE.PLOT(gm,"Comp_m.png", 5, 9)
SAVE.PLOT(gf,"Comp_f.png", 5, 9)

##################################################
##################################################
# Web application
##################################################
##################################################  

# see separate code
