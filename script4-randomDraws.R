library("dplyr")
library("ggplot2")
library("epitools")
library("descr")
library("ggstance")
adpa <- read.csv("data_phase 3/ADPA.csv")
adsl <- read.csv("data_phase 3/ADSL.csv")
############ calculate OR via epitools
MCAR_get_OR <- function(SEED,percentage){
  set.seed(SEED)
  subjid_sample<-sample(1:1831,1831*percentage)
  ##### completers ######
  completersData <- adpa %>%
    filter(!(SUBJID %in% subjid_sample) & PARAMN == 10 & !(AVISITN == 2)) %>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
    select(c(1:20,29)) %>%
    filter(AVISITN == 6)
  #threeWayTableCompleters<-ftable(table(completersData$TRTPN,completersData$SEX,completersData$PCHGCA1N))
  ctCompletersSex =CrossTable(completersData$SEX,completersData$PCHGCA1N)
  ORcompletersSex<- oddsratio.wald(ctCompletersSex$t)
  ctCompletersTrtp =CrossTable(completersData$TRTPN,completersData$PCHGCA1N)
  ORcompletersTrtp<-oddsratio.wald(ctCompletersTrtp$t)
  ##### imputed ######
  imputedData <- adpa %>%
    filter(PARAMN == 10 & !(AVISITN == 2)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjid_sample,0,PCHGCA1N))%>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID")  %>%
    select(c(1:20,29)) %>%
    filter(AVISITN == 6)
  ctImputedSex =CrossTable(imputedData$SEX,imputedData$PCHGCA1N)
  ORimputedSex<- oddsratio.wald(ctImputedSex$t)
  ctImputedTrtp =CrossTable(imputedData$TRTPN,imputedData$PCHGCA1N)
  ORimputedTrtp<-oddsratio.wald(ctImputedTrtp$t)
  ##### complete ######
  completeData<-adpa %>%
    filter(PARAMN == 10 & !(AVISITN == 2)) %>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
    select(c(1:20,29)) %>%
    filter(AVISITN == 6)
  ctCompleteSex =CrossTable(completeData$SEX,completeData$PCHGCA1N)
  ORcompleteSex<- oddsratio.wald(ctCompleteSex$t)
  ctCompleteTrtp =CrossTable(completeData$TRTPN,completeData$PCHGCA1N)
  ORcompleteTrtp<-oddsratio.wald(ctCompleteTrtp$t)
  ####
  # return
  c(list(ORcompletersSex = ORcompletersSex,ORcompletersTrtp=ORcompletersTrtp,ORimputedSex=ORimputedSex,ORimputedTrtp = ORimputedTrtp,ORcompleteSex = ORcompleteSex,ORcompleteTrtp=ORcompleteTrtp))
}
# test functions
tmp<-lapply(c(999,1000),
       FUN = MCAR_get_OR,
       0.1)
#####################################################
############### MCAR ################################
#####################################################
############### calculate OR via glm ################
MCAR_get_OR <- function(SEED,percentage){
  set.seed(SEED)
  subjid_sample<-sample(1:1831,1831*percentage)
  ##### completers ######
  completersData <- adpa %>%
    filter(!(SUBJID %in% subjid_sample) & PARAMN == 10 & !(AVISITN == 2)) %>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
    select(c(1:20,29)) %>%
    filter(AVISITN == 6)
  # threeWayTableCompleters<-ftable(table(completersData$TRTPN,completersData$SEX,completersData$PCHGCA1N))
  model1 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = completersData)
  coef.model1 <- exp(coefficients(model1))
  ORcompletersSexM2F <- as.numeric(round(coef.model1[2],4)) # OR male to female
  ORcompletersAvsP <- as.numeric(round(1/coef.model1[3],4)) # OR active control vs placebo
  ORcompleters140vsP <- as.numeric(round(coef.model1[4]/coef.model1[3],4)) # OR 140mg vs placebo
  ORcompleters210vsP <- as.numeric(round(coef.model1[5]/coef.model1[3],4)) # OR 140mg vs placebo
  ORcompleters <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
  ##### imputed ######
  imputedData <- adpa %>%
    filter(PARAMN == 10 & !(AVISITN == 2)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjid_sample,0,PCHGCA1N))%>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID")  %>%
    select(c(1:20,29)) %>%
    filter(AVISITN == 6)
  # threeWayTableImputed<-ftable(table(imputedData$TRTPN,imputedData$SEX,imputedData$PCHGCA1N))
  model2 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = imputedData)
  coef.model2 <- exp(coefficients(model2))
  ORimputedSexM2F <- as.numeric(round(coef.model2[2],4)) # OR male to female
  ORimputedAvsP <- as.numeric(round(1/coef.model2[3],4)) # OR active control vs placebo
  ORimputed140vsP <- as.numeric(round(coef.model2[4]/coef.model2[3],4)) # OR 140mg vs placebo
  ORimputed210vsP <- as.numeric(round(coef.model2[5]/coef.model2[3],4)) # OR 140mg vs placebo
  ORimputed <- list(SexM2F = ORimputedSexM2F,AvsP = ORimputedAvsP,T140vsP = ORimputed140vsP, T210vsP =ORimputed210vsP)
  #######
  # return
  c(list(ORcompleters = ORcompleters,ORimputed=ORimputed))
}
n = 1000 # take seeds 1 to 1000
############# 10% missing ##################
MCAR_10P_result<-lapply(1:1000,
            FUN = MCAR_get_OR,
            0.1)
tmp.1<-unlist(MCAR_10P_result)
MCAR_10P_meanOR<- c(by(tmp.1, names(tmp.1), mean))
MCAR_10P_varOR<- c(by(tmp.1, names(tmp.1), var))
MCAR_10P_result <- matrix(,nrow=8,ncol = 4)
MCAR_10P_result[,1]<-MCAR_10P_meanOR
MCAR_10P_result[,2]<-MCAR_10P_varOR
MCAR_10P_result[,3]<-MCAR_10P_meanOR-qnorm(0.975) * sqrt(MCAR_10P_varOR)/sqrt(n)
MCAR_10P_result[,4]<-MCAR_10P_meanOR+qnorm(0.975) * sqrt(MCAR_10P_varOR)/sqrt(n)
MCAR_10P_result = data.frame(MCAR_10P_result)
names(MCAR_10P_result)<- c("mean","var","low","high")
########### 20% missing ###################
MCAR_20P_result<-lapply(1:1000,
                        FUN = MCAR_get_OR,
                        0.2)
tmp.2<-unlist(MCAR_20P_result)
MCAR_20P_meanOR<- c(by(tmp.2, names(tmp.2), mean))
MCAR_20P_varOR<- c(by(tmp.2, names(tmp.2), var))
MCAR_20P_result <- matrix(,nrow=8,ncol = 4)
MCAR_20P_result[,1]<-MCAR_20P_meanOR
MCAR_20P_result[,2]<-MCAR_20P_varOR
MCAR_20P_result[,3]<-MCAR_20P_meanOR-qnorm(0.975) * sqrt(MCAR_20P_result[,2]-MCAR_20P_varOR)/sqrt(n)
MCAR_20P_result[,4]<-MCAR_20P_meanOR+qnorm(0.975) * sqrt(MCAR_20P_result[,2]+MCAR_20P_varOR)/sqrt(n)
MCAR_20P_result = data.frame(MCAR_20P_result)
names(MCAR_20P_result)<- c("mean","var","low","high")
########### 30% missing ###################
MCAR_30P_result<-lapply(1:1000,
                        FUN = MCAR_get_OR,
                        0.3)
tmp.3<-unlist(MCAR_30P_result)
MCAR_30P_meanOR<- c(by(tmp.3, names(tmp.3), mean))
MCAR_30P_varOR<- c(by(tmp.3, names(tmp.3), var))
MCAR_30P_result <- matrix(,nrow=8,ncol = 4)
MCAR_30P_result[,1]<-MCAR_30P_meanOR
MCAR_30P_result[,2]<-MCAR_30P_varOR
MCAR_30P_result[,3]<-MCAR_30P_meanOR-qnorm(0.975) * sqrt(MCAR_30P_result[,2]-MCAR_30P_varOR)/sqrt(n)
MCAR_30P_result[,4]<-MCAR_30P_meanOR+qnorm(0.975) * sqrt(MCAR_30P_result[,2]+MCAR_30P_varOR)/sqrt(n)
MCAR_30P_result = data.frame(MCAR_30P_result)
names(MCAR_30P_result)<- c("mean","var","low","high")
#####################################################
############### MAR #################################
#####################################################
MAR_get_OR <- function(SEED,percentage)
###### complete ######
completeData<-adpa %>%
  filter(PARAMN == 10 & !(AVISITN == 2)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
model3 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = completeData)
coef.model3 <- exp(coefficients(model3))
ORcompleteSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
ORcompleteAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
ORcomplete140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
ORcomplete210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
ORcomplete <- list(SexM2F = ORcompleteSexM2F,AvsP = ORcompleteAvsP,T140vsP = ORcomplete140vsP, T210vsP =ORcomplete210vsP)
############## plot ###########################
boxLabels=factor(c(1,2,3,4,1,2,3,4,1,2,3,4))
boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
df = data.frame(
  yAxis = boxLabels,
  boxOdds = c(MCAR_10P_result$mean[1],MCAR_10P_result$mean[3],MCAR_10P_result$mean[4],MCAR_10P_result$mean[2],
              MCAR_10P_result$mean[5],MCAR_10P_result$mean[7],MCAR_10P_result$mean[8],MCAR_10P_result$mean[6],
              ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCILow = c(MCAR_10P_result$low[1],MCAR_10P_result$low[3],MCAR_10P_result$low[4],MCAR_10P_result$low[2],
               MCAR_10P_result$low[5],MCAR_10P_result$low[7],MCAR_10P_result$low[8],MCAR_10P_result$low[6],
               ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCIHigh = c(MCAR_10P_result$high[1],MCAR_10P_result$high[3],MCAR_10P_result$high[4],MCAR_10P_result$high[2],
                MCAR_10P_result$high[5],MCAR_10P_result$high[7],MCAR_10P_result$high[8],MCAR_10P_result$high[6],
                ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  group=factor(c('Completers Data','Completers Data','Completers Data','Completers Data',
                 'Imputed Data','Imputed Data','Imputed Data','Imputed Data',
                 'Complete Data','Complete Data','Complete Data','Complete Data')),
  VAL=c(round(MCAR_10P_result$mean[1],4),round(MCAR_10P_result$mean[3],4),round(MCAR_10P_result$mean[4],4),round(MCAR_10P_result$mean[2],4),
        round(MCAR_10P_result$mean[5],4),round(MCAR_10P_result$mean[7],4),round(MCAR_10P_result$mean[8],4),round(MCAR_10P_result$mean[6],4),
        round(ORcomplete$AvsP,4),round(ORcomplete$T140vsP,4),round(ORcomplete$T210vsP,4),round(ORcomplete$SexM2F,4))
)

p1 = ggplot(df, aes(x = boxOdds, y = yAxis, colour = group))
p1 + geom_errorbarh(aes(y = yAxis, xmax = boxCIHigh, xmin = boxCILow), size = 3.5, height = .2, position = position_dodgev(height=0.9)) +
  geom_point(aes(x = boxOdds, y = yAxis),size = 1,  position = position_dodgev(height=0.9))+
  geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
  theme(panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = c(1,seq(0,100,20)))+
  coord_trans(x = 'log10')+
  ylab("") +
  xlab("Odds ratio (log scale)") +
  ggtitle("Odds Ratios when 10% Patients' Data MCAR, 95% Walds CI")
