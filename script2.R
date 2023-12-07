library("dplyr")
library("ggplot2")
# load the adpa (Efficacy Analysis Dataset)
## One record per study identifier per subject identifier per parameter code (per analysis visit)
adpa <- read.csv("data_phase 3/ADPA.csv")
adsl <- read.csv("data_phase 3/ADSL.csv")
########################
######## Step 1 ########
########################
set.seed(1000)
subjid_sample<-sample(1:1831,1831*0.1) # randomly sample 10% of 1831 subjects
# set their visit6 PASI score (AVAL) and PASI 75 (PCHG) to NA
# adpa_mcar <- adpa %>%
#   mutate(AVAL = ifelse(SUBJID %in% subjid_sample & PARAMN == 10 & AVISITN == 6, NA, AVAL)) %>%
#   mutate(PCHG = ifelse(SUBJID %in% subjid_sample & PARAMN == 10 & AVISITN == 6, NA, PCHG))
########################
######## Step 2 ########
########################
# 2.1
# completers = patients who are not selected to have their PASI at visit6 missing
#           = patients id not in subjid_sample
# exclue SPGA data and visit2 data of PASI
completers10PData <- adpa %>%
  filter(!(SUBJID %in% subjid_sample) & PARAMN == 10 & !(AVISITN == 2)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
write.csv(completers10PData,"2a.csv",row.names = FALSE)
# # logistic regression: PCHGCA1N ~sex
# model2.1SEX <- glm(PCHGCA1N~SEX, family = binomial(), data = completersData10P)
# summary(model2.1SEX)
# # logistic regression: PCHGCA1N ~trtp
# model2.1TRTP <-glm(PCHGCA1N~TRTP, family = binomial(), data = completersData)
# summary(model2.1TRTP)
# # logistic regression: PCHGCA1N ~sex+trtp 
# model2.1SEXTRTP<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = completersData)
# summary(model2.1SEXTRTP)
# 2.2
imputed10PData <- adpa %>%
  filter(PARAMN == 10 & !(AVISITN == 2)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjid_sample,0,PCHGCA1N))%>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID")  %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
write.csv(imputed10PData,"2b.csv",row.names = FALSE)
# # logistic regression: PCHGCA1N ~sex
# model2.2SEX <- glm(PCHGCA1N~SEX, family = binomial(), data = imputedData)
# summary(model2.2SEX)
# # logistic regression: PCHGCA1N ~ trtp
# model2.2TRTP <-glm(PCHGCA1N~TRTP, family = binomial(), data = imputedData)
# summary(model2.2TRTP)
# # logistic regression: PCHGCA1N ~sex + trtp
# model2.2SEXTRTP<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = imputedData)
# summary(model2.2SEXTRTP)
# 2.3
completeData<-adpa %>%
  filter(PARAMN == 10 & !(AVISITN == 2)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
write.csv(completeData,"2c.csv",row.names = FALSE)
# # logistic regression: PCHGCA1N ~sex
# model2.3SEX <- glm(PCHGCA1N~SEX, family = binomial(), data = completeData)
# summary(model2.3SEX)
# # logistic regression: PCHGCA1N ~ trtp
# model2.3TRTP <-glm(PCHGCA1N~TRTP, family = binomial(), data = completeData)
# summary(model2.3TRTP)
# # logistic regression: PCHGCA1N ~ sex + trtp
# model2.3SEXTRTP<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = completeData)
# summary(model2.3SEXTRTP)
########################
######## Step 3 ########
########################
# 20% completers
subjid20Psample<-sample(1:1831,1831*0.2) 
completers20PData <- adpa %>%
  filter(!(SUBJID %in% subjid20Psample) & PARAMN == 10 & !(AVISITN == 2)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
write.csv(completers20PData,"3-20-a.csv",row.names = FALSE)
# 20% imputed
imputed20PData <- adpa %>%
  filter(PARAMN == 10 & !(AVISITN == 2)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjid20Psample,0,PCHGCA1N)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID")  %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
write.csv(imputed20PData,"3-20-b.csv",row.names = FALSE)
# 30% completers
subjid30Psample<-sample(1:1831,1831*0.3)
completers30PData <- adpa %>%
  filter(!(SUBJID %in% subjid30Psample) & PARAMN == 10 & !(AVISITN == 2)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
write.csv(completers30PData,"3-30-a.csv",row.names = FALSE)
# 30% imputed
imputed30PData <- adpa %>%
  filter(PARAMN == 10 & !(AVISITN == 2)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjid30Psample,0,PCHGCA1N)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID")  %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
write.csv(imputed30PData,"3-30-b.csv",row.names = FALSE)
########################
######## Step 4 ########
########################
# subjects with 30% missing visit4
subjidMARvisit4high<- adpa %>%
  filter(AVISITN == 3 & PARAMN == 10 & PCHG < 10) %>%
  pull(SUBJID)
subjidMARvisit4high_sample <- sample(subjidMARvisit4high,0.3*length(subjidMARvisit4high))
# subjects with 5% missing visit4
allSubjid <- 1:1831
subjidMARvisit4low <- setdiff(allSubjid,subjidMARvisit4high) 
subjidMARvisit4low_sample <- sample(subjidMARvisit4low,0.05*length(subjidMARvisit4low))
# set the selected subjects' PASI score to NA on visit 4 and on subsequent visits
adpa_mar <- adpa %>% 
  left_join(select(adsl,SUBJID,SEX),by="SUBJID")  %>%
  mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit4high_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,AVAL)) %>%
  mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit4high_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,PCHG)) %>%
  mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit4low_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,AVAL)) %>%
  mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit4low_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,PCHG)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit4low_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,PCHGCA1N)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit4high_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,PCHGCA1N)) 
# subjects with 30% missing visit5
subjidMARvisit5high<- adpa_mar %>%
  filter(AVISITN == 4 & PARAMN == 10 & PCHG < 10) %>%
  pull(SUBJID)
subjidMARvisit5high_sample <- sample(subjidMARvisit5high,0.3*length(subjidMARvisit5high))
# subjects with 5% missing visit5
subjidMARvisit5low <- setdiff(allSubjid,subjidMARvisit5high) 
subjidMARvisit5low_sample <- sample(subjidMARvisit5low,0.05*length(subjidMARvisit5low))
# set the selected subjects' PASI score to NA on visit 5 and 6
adpa_mar <- adpa_mar %>% 
  mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,AVAL)) %>%
  mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHG)) %>%
  mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,AVAL)) %>%
  mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHG)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHGCA1N)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHGCA1N)) 
# subjects with 30% missing visit6
subjidMARvisit6high<- adpa_mar %>%
  filter(AVISITN == 5 & PARAMN == 10 & PCHG < 10) %>%
  pull(SUBJID)
subjidMARvisit6high_sample <- sample(subjidMARvisit6high,0.3*length(subjidMARvisit6high))
# subjects with 5% missing visit5
subjidMARvisit6low <- setdiff(allSubjid,subjidMARvisit6high) 
subjidMARvisit6low_sample <- sample(subjidMARvisit6low,0.05*length(subjidMARvisit6low))
adpa_mar <- adpa_mar %>% 
  mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN == 6,NA,AVAL)) %>%
  mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN ==6,NA,PCHG)) %>%
  mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,AVAL)) %>%
  mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,PCHG)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,PCHGCA1N)) %>%
  mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN == 6,NA,PCHGCA1N)) %>%
  select(c(1:20,29))%>%
  filter(PARAMN == 10 & !(AVISITN == 2)) 
write.csv(adpa_mar,"4.csv",row.names = FALSE)
