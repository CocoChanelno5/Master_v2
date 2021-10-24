#Panel Unit Root Testing
library(dplyr)
library(plm)
library(tseries)
library(CADFtest)

setwd("~/Desktop/Magisterka/Master_git")
source("unit_root_package.R")

##################### INTERCEPT
#Levin, Lin and Chu Unit Root Test
PL_GDP_test <- PL_GDP %>% select(1,6,7)
test <- data.frame(split(PL_GDP_test$Value,PL_GDP_test$ID))
LLC1 <- purtest(test, test = "levinlin", exo = "intercept", lags = "AIC", pmax = 1)
summary(LLC1) #cannot reject the null, all of the series are I(1)
Choi_PL_GDP <- pCADFtest(Y=test, type = "drift", max.lag.y = 1, criterion = "AIC")
print(paste("Statystyka testu levinlin dla zestawu PL_GDP wynosi",LLC1$statistic,"a p-value:",LLC1$p.value,
            "natomiast p-value w tescie Choia wynosi",Choi_PL_GDP$p.value))

PL_UE_test <- PL_UE %>% select(2,8,9)
test <- data.frame(split(PL_UE_test$Value,PL_UE_test$ID))
LLC2 <- purtest(test, test = "levinlin", exo = "intercept", lags = "AIC", pmax = 4)
summary(LLC2) #cannot reject the null, all of the series are I(1)
Choi_PL_UE <- pCADFtest(Y=test, type = "drift", max.lag.y = 12, criterion = "AIC")
print(paste("Statystyka testu levinlin dla zestawu PL_UE wynosi",LLC2$statistic,"a p-value:",LLC2$p.value,
            "natomiast p-value w tescie Choia wynosi",Choi_PL_UE$p.value))

USA_GDP_test <- USA_GDP %>% select(1,7,8)
test <- data.frame(split(USA_GDP_test$Value,USA_GDP_test$ID))
LLC3 <- purtest(test, test = "levinlin", exo = "intercept", lags = "AIC", pmax = 4)
summary(LLC3) #cannot reject the null, all of the series are I(1)
Choi_USA_GDP <- pCADFtest(Y=test, type = "drift", max.lag.y = 4, criterion = "AIC")
print(paste("Statystyka testu levinlin dla zestawu USA_GDP wynosi",LLC3$statistic,"a p-value:",LLC3$p.value,
            "natomiast p-value w tescie Choia wynosi",Choi_USA_GDP$p.value))

USA_UE_test <- USA_UE %>% select(1,6,7)
test <- data.frame(split(USA_UE_test$Value,USA_UE_test$ID))
LLC4 <- purtest(test, test = "levinlin", #exo = "intercept", 
               lags = "AIC", pmax = 12)
summary(LLC4) #reject the null, at least one series is I(0)
Choi_USA_UE <- pCADFtest(Y=test, type = "drift", max.lag.y = 12, criterion = "AIC")
print(paste("Statystyka testu levinlin dla zestawu USA_UE wynosi",LLC4$statistic,"a p-value:",LLC4$p.value,
            "natomiast p-value w tescie Choia wynosi",Choi_USA_UE$p.value))


#Im Pesaran Shin Unit Root Test
PL_GDP_test <- PL_GDP %>% select(1,6,7)
test <- data.frame(split(PL_GDP_test$Value,PL_GDP_test$ID))
LLC1 <- purtest(test, test = "ips", exo = "intercept",
               lags = "AIC", pmax = 2)
summary(LLC1)
print(LLC1)#NA
print(paste("Statystyka testu ips dla zestawu PL_GDP wynosi",LLC1$statistic,"a p-value:",LLC1$p.value))

PL_UE_test <- PL_UE %>% select(2,8,9)
test <- data.frame(split(PL_UE_test$Value,PL_UE_test$ID))
LLC2 <- purtest(test, test = "ips", exo = "intercept", lags = "AIC", pmax = 4)
summary(LLC2) #cannot reject the null, all of the series are I(1)
print(LLC2)
print(paste("Statystyka testu ips dla zestawu PL_UE wynosi",LLC2$statistic,"a p-value:",LLC2$p.value))

USA_GDP_test <- USA_GDP %>% select(1,7,8)
test <- data.frame(split(USA_GDP_test$Value,USA_GDP_test$ID))
LLC3 <- purtest(test, test = "ips", exo = "intercept", lags = "AIC", pmax = 1)
summary(LLC3) #cannot reject the null, all of the series are I(1)
print(LLC3)
print(paste("Statystyka testu ips dla zestawu USA_GDP wynosi",LLC3$statistic,"a p-value:",LLC3$p.value))

USA_UE_test <- USA_UE %>% select(1,6,7)
test <- data.frame(split(USA_UE_test$Value,USA_UE_test$ID))
LLC4 <- purtest(test, test = "ips", exo = "intercept", lags = "AIC", pmax =12)
summary(LLC4) #reject the null, at least one series is I(0)
print(LLC4)
print(paste("Statystyka testu ips dla zestawu USA_UE wynosi",LLC4$statistic,"a p-value:",LLC4$p.value))


write.csv(PL_GDP_test,"~/Desktop/Magisterka/Master_git/PL_GDP.csv")
write.csv(PL_UE_test,"~/Desktop/Magisterka/Master_git/PL_UE.csv")
write.csv(USA_GDP_test,"~/Desktop/Magisterka/Master_git/USA_GDP.csv")
write.csv(USA_UE_test,"~/Desktop/Magisterka/Master_git/USA_UE.csv")

################################### SECOND APPROACH
p_data<-PL_GDP
panel_PL_GDP<-pdata.frame(p_data, index=c("ID","Date"))
test1<-purtest(panel_PL_GDP$Value, exo = "intercept",data =panel_PL_GDP,test = "levinlin", lags = "AIC", pmax = 2)
#summary(test1)
print(test1)

p_data<-PL_UE
panel_PL_UE<-pdata.frame(p_data, index=c("ID","Date"))
test2<-purtest(panel_PL_UE$Value, exo = "intercept",data = panel_PL_UE,test = "levinlin", lags = "AIC", pmax = 4)
#summary(test2)
print(test2)

p_data<-USA_GDP
panel_USA_GDP<-pdata.frame(p_data, index=c("ID","Date"))
test3<-purtest(panel_USA_GDP$Value,exo = "intercept", data = panel_USA_GDP,test = "levinlin", lags = "AIC", pmax = 1)
#summary(test3)
print(test3)

p_data<-USA_UE
panel_USA_UE<-pdata.frame(p_data, index=c("ID","Date"))
test4<-purtest(panel_USA_UE$Value, exo = "intercept",data = panel_USA_UE,test = "levinlin", lags = "AIC", pmax = 12)
#summary(test4)
print(test4)





################################### THIRD APPROACH

purtest(
  object,
  data = NULL,
  index = NULL,
  test = c("levinlin", "ips", "madwu", "Pm", "invnormal", "logit", "hadri"),
  exo = c("none", "intercept", "trend"),
  lags = c("SIC", "AIC", "Hall"),
  pmax = 10,
  Hcons = TRUE,
  q = NULL,
  dfcor = FALSE,
  fixedT = TRUE,
  ips.stat = NULL
)
## S3 method for class 'purtest'
print(x, ...)
## S3 method for class 'purtest'
summary(object, ...)
## S3 method for class 'summary.purtest'
print(x, ...)


