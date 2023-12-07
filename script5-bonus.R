library("dplyr")
library("ggplot2")
library("epitools")
library("descr")
library("ggstance")
adpa <- read.csv("data_phase 3/ADPA.csv")
adsl <- read.csv("data_phase 3/ADSL.csv")
# Six b
# a concave function to represent prob missing
# 0 - 0.3
corrected_logistic<-function(x){
  0.6*(1/(1+exp(-0.06*(100-x)))-0.5)
}
x=0:100
plot(x,corrected_logistic(x),type="l",
     xlab='PCHG', ylab='Missing Probability')

improvement = adpa%>%
  filter(PARAMN == 10 & AVISITN == 6)%>%
  mutate(Possiblity = corrected_logistic(PCHG))%>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID")
Missdata=function(x,p){
  r= runif(1)
  if (r<= p) {x=NA}
  x
}
improvement$PCHG <-mapply(Missdata,improvement$PCHG,improvement$Possiblity)
improvement_completers <- improvement %>%
  filter(!is.na(PCHG))
improvement_imputed <- improvement %>%
  mutate(PCHGCA1N = ifelse(is.na(PCHG),0,PCHGCA1N))
######### completers ##########
###### MNAR large sample ######
MNAR6b_get_OR <- function(SEED){
  set.seed(SEED)
  mnar_6b = adpa%>%
    filter(PARAMN == 10 & AVISITN == 6)%>%
    mutate(Possibility = corrected_logistic(PCHG))%>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID")%>%
    select(c(1:20,29,30))
  mnar_6b$PCHG <-mapply(Missdata,mnar_6b$PCHG,mnar_6b$Possibility)
  mnar_6b_completers <- mnar_6b %>%
    filter(!is.na(PCHG))
  mnar_6b_imputed <- mnar_6b %>%
    mutate(PCHGCA1N = ifelse(is.na(PCHG),0,PCHGCA1N))
  # completers
  model1 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = mnar_6b_completers)
  coef.model1 <- exp(coefficients(model1))
  ORcompletersSexM2F <- as.numeric(round(coef.model1[2],4)) # OR male to female
  ORcompletersAvsP <- as.numeric(round(1/coef.model1[3],4)) # OR active control vs placebo
  ORcompleters140vsP <- as.numeric(round(coef.model1[4]/coef.model1[3],4)) # OR 140mg vs placebo
  ORcompleters210vsP <- as.numeric(round(coef.model1[5]/coef.model1[3],4)) # OR 140mg vs placebo
  ORcompleters <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
  # imputed
  model2 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = mnar_6b_imputed)
  coef.model2 <- exp(coefficients(model2))
  ORimputedSexM2F <- as.numeric(round(coef.model2[2],4)) # OR male to female
  ORimputedAvsP <- as.numeric(round(1/coef.model2[3],4)) # OR active control vs placebo
  ORimputed140vsP <- as.numeric(round(coef.model2[4]/coef.model2[3],4)) # OR 140mg vs placebo
  ORimputed210vsP <- as.numeric(round(coef.model2[5]/coef.model2[3],4)) # OR 140mg vs placebo
  ORimputed <- list(SexM2F = ORimputedSexM2F,AvsP = ORimputedAvsP,T140vsP = ORimputed140vsP, T210vsP =ORimputed210vsP)
  # complete #
  model3 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = completeData)
  coef.model3 <- exp(coefficients(model3))
  ORcompleteSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
  ORcompleteAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
  ORcomplete140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
  ORcomplete210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
  ORcomplete <- list(SexM2F = ORcompleteSexM2F,AvsP = ORcompleteAvsP,T140vsP = ORcomplete140vsP, T210vsP =ORcomplete210vsP)
  # return #
  c(list(ORcompleters = ORcompleters,ORimputed=ORimputed, ORcomplete=ORcomplete))
}
# 1000 seeds
MNAR_6b_result1 <- lapply(1:1000,MNAR6b_get_OR)
n=1000
tmp.mnar<-unlist(MNAR_6b_result1)
MNAR_meanOR<- c(by(tmp.mnar, names(tmp.mnar), mean))
MNAR_varOR<- c(by(tmp.mnar, names(tmp.mnar), var))
MNAR_result <- matrix(,nrow=12,ncol = 4)
MNAR_result[,1]<-MNAR_meanOR
MNAR_result[,2]<-MNAR_varOR
MNAR_result[,3]<-MNAR_meanOR-qnorm(0.975) * sqrt(MNAR_result[,2]<-MNAR_varOR)/sqrt(n)
MNAR_result[,4]<-MNAR_meanOR+qnorm(0.975) * sqrt(MNAR_result[,2]<-MNAR_varOR)/sqrt(n)
MNAR_result = data.frame(MNAR_result)
names(MNAR_result)<- c("mean","var","low","high")
row.names(MNAR_result)<-names(MNAR_meanOR)
# plot MNAR #
boxLabels=factor(c(1,2,3,4,1,2,3,4,1,2,3,4))
boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
# plot 1: MCAR
df.mnar = data.frame(
  yAxis = boxLabels,
  boxOdds = c(MNAR_result$mean[5],MNAR_result$mean[7],MNAR_result$mean[8],MNAR_result$mean[6],
              MNAR_result$mean[9],MNAR_result$mean[11],MNAR_result$mean[12],MNAR_result$mean[10],
              MNAR_result$mean[1],MNAR_result$mean[3],MNAR_result$mean[4],MNAR_result$mean[2]),
  boxCILow = c(MNAR_result$low[5],MNAR_result$low[7],MNAR_result$low[8],MNAR_result$low[6],
               MNAR_result$low[9],MNAR_result$low[11],MNAR_result$low[12],MNAR_result$low[10],
               MNAR_result$low[1],MNAR_result$low[3],MNAR_result$low[4],MNAR_result$low[2]),
  boxCIHigh = c(MNAR_result$high[5],MNAR_result$high[7],MNAR_result$high[8],MNAR_result$high[6],
                MNAR_result$high[9],MNAR_result$high[11],MNAR_result$high[12],MNAR_result$high[10],
                MNAR_result$high[1],MNAR_result$high[3],MNAR_result$high[4],MNAR_result$high[2]),
  group=factor(c('Completers Data','Completers Data','Completers Data','Completers Data',
                 'Imputed Data','Imputed Data','Imputed Data','Imputed Data',
                 'Complete Data','Complete Data','Complete Data','Complete Data')),
  VAL=c(round(MNAR_result$mean[5],3),round(MNAR_result$mean[7],3),round(MNAR_result$mean[8],3),round(MNAR_result$mean[6],3),
        round(MNAR_result$mean[9],3),round(MNAR_result$mean[11],3),round(MNAR_result$mean[12],3),round(MNAR_result$mean[10],3),
        round(MNAR_result$mean[1],3),round(MNAR_result$mean[3],3),round(MNAR_result$mean[4],3),round(MNAR_result$mean[2],3))
)

p = ggplot(df.mnar, aes(x = boxOdds, y = yAxis, colour = group))
p.mnar <-p + geom_errorbarh(aes(y = yAxis, xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, position = position_dodgev(height=0.9)) +
  geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
  geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  scale_x_continuous(breaks = c(1,seq(0,100,20)))+
  coord_trans(x = 'log10')+
  ylab("") +
  xlab("Odds ratio (log scale)") +
  ggtitle("MNAR")
p.mnar  
