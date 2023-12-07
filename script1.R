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
adpa_mcar <- adpa %>%
  mutate(AVAL = ifelse(SUBJID %in% subjid_sample & PARAMN == 10 & AVISITN == 6, NA, AVAL)) %>%
  mutate(PCHG = ifelse(SUBJID %in% subjid_sample & PARAMN == 10 & AVISITN == 6, NA, PCHG))
########################
######## Step 2 ########
########################
# 2.1
# completers = patients who are not selected to have their PASI at visit6 missing
#           = patients id not in subjid_sample
# exclue SPGA data and visit2 data of PASI
completersData <- adpa %>%
  filter(!(SUBJID %in% subjid_sample) & PARAMN == 10 & !(AVISITN == 2)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  mutate(SEXN = ifelse(SEX == "M",1,0)) %>% # add a column SEXN, male = 1, female = 0
  mutate(SEXNTRTPN = case_when((SEXN == 1 & TRTPN == 1) ~ 1,# add a column, sex and treatment, 8 labels
                               (SEXN == 1 & TRTPN == 2) ~ 2,
                               (SEXN == 1 & TRTPN == 3) ~ 3,
                               (SEXN == 1 & TRTPN == 4) ~ 4,
                               (SEXN == 0 & TRTPN == 1) ~ 5,
                               (SEXN == 0 & TRTPN == 2) ~ 6,
                               (SEXN == 0 & TRTPN == 3) ~ 7,
                               (SEXN == 0 & TRTPN == 4) ~ 8,
    ))
completersData$SEXNTRTPN <-factor(completersData$SEXNTRTPN,level = unique(completersData$SEXNTRTPN))
# logistic regression: regress PASI75 on sex (male =1 , female =0) only
model2.1SEX <- glm(PCHGCA1N~SEX, family = binomial(), data = completersData)
summary(model2.1SEX)
# logistic regression: regress on TRTP only
model2.1TRTP <-glm(PCHGCA1N~TRTP, family = binomial(), data = completersData)
summary(model2.1TRTP)
# logistic regression on Treatment and Sex (8 categories)
model2.1SEXNTRTPN<-glm(PCHGCA1N~SEXNTRTPN, family = binomial(), data = completersData)
summary(model2.1SEXNTRTPN)
# logistic regression: PCHGCA1N ~sex+trtp 
model2.1SEXTRTPN<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = completersData)
summary(model2.1SEXTRTPN)
imputedData <- adpa_mcar %>%
  filter(PARAMN == 10 & !(AVISITN == 2)) %>%
  mutate(PCHGCA1N = ifelse(is.na(AVAL),0,PCHGCA1N)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  mutate(SEXN = ifelse(SEX == "M",1,0)) %>% # add a column SEXN, male = 1, female = 0
  mutate(SEXNTRTPN = case_when((SEXN == 1 & TRTPN == 1) ~ 1,# add a column, sex and treatment, 8 labels
                               (SEXN == 1 & TRTPN == 2) ~ 2,
                               (SEXN == 1 & TRTPN == 3) ~ 3,
                               (SEXN == 1 & TRTPN == 4) ~ 4,
                               (SEXN == 0 & TRTPN == 1) ~ 5,
                               (SEXN == 0 & TRTPN == 2) ~ 6,
                               (SEXN == 0 & TRTPN == 3) ~ 7,
                               (SEXN == 0 & TRTPN == 4) ~ 8,
  ))
imputedData$SEXNTRTPN <-factor(imputedData$SEXNTRTPN,level = unique(imputedData$SEXNTRTPN))
model2.2SEX <- glm(PCHGCA1N~SEX, family = binomial(), data = imputedData)
summary(model2.2SEX)
# logistic regression: regress on TRTP only
model2.2TRTP <-glm(PCHGCA1N~TRTP, family = binomial(), data = imputedData)
summary(model2.2TRTP)
# logistic regression on Treatment and Sex (8 categories)
model2.2SEXNTRTPN<-glm(PCHGCA1N~SEXNTRTPN, family = binomial(), data = imputedData)
summary(model2.2SEXNTRTPN)
#2.3
completeData<-adpa %>%
  filter(PARAMN == 10 & !(AVISITN == 2)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  mutate(SEXN = ifelse(SEX == "M",1,0)) %>% # add a column SEXN, male = 1, female = 0
  mutate(SEXNTRTPN = case_when((SEXN == 1 & TRTPN == 1) ~ 1,# add a column, sex and treatment, 8 labels
                               (SEXN == 1 & TRTPN == 2) ~ 2,
                               (SEXN == 1 & TRTPN == 3) ~ 3,
                               (SEXN == 1 & TRTPN == 4) ~ 4,
                               (SEXN == 0 & TRTPN == 1) ~ 5,
                               (SEXN == 0 & TRTPN == 2) ~ 6,
                               (SEXN == 0 & TRTPN == 3) ~ 7,
                               (SEXN == 0 & TRTPN == 4) ~ 8,
  ))
completeData$SEXNTRTPN <-factor(completeData$SEXNTRTPN,level = unique(completeData$SEXNTRTPN))
model2.3SEX <- glm(PCHGCA1N~SEX, family = binomial(), data = completeData)
summary(model2.3SEX)
# logistic regression: regress on TRTP only
model2.3TRTP <-glm(PCHGCA1N~TRTP, family = binomial(), data = completeData)
summary(model2.3TRTP)
# logistic regression on Treatment and Sex (8 categories)
model2.3SEXNTRTPN<-glm(PCHGCA1N~SEXNTRTPN, family = binomial(), data = completeData)
summary(model2.3SEXNTRTPN)

model2.1.1<-glm(PCHGCA1N~SEXN + TRTPN + SEXN * TRTPN, family = binomial(), data = completersData)
