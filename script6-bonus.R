library("dplyr")
library("ggplot2")
library("epitools")
library("descr")
library("ggstance")
adpa = read.csv("data_phase 3/ADPA.csv")
adsl = read.csv("data_phase 3/ADSL.csv")
adae = read.csv("data_phase 3/ADAE.csv")
adae_1 = adae%>%
  mutate(SUBJID=as.numeric(substr(USUBJID,nchar(USUBJID)-3,nchar(USUBJID))))
mapping=c('VISIT 2'=0,'VISIT 3'=21,'VISIT 4'=42,'VISIT 5'=63,'VISIT 6'=84)
mapping2=c('MILD'=0.1,"MODERATE"=0.15,'SEVERE'=0.2)
mnar_6a=adpa%>%
  filter(PARAMN == 10 )%>%
  full_join(select(adae_1,SUBJID,AESTDY,AEENDY,AESER,AESEV),by='SUBJID',suffix = c("_1", "_2"))%>%
  mutate(NVISIT=mapping[AVISIT])%>%
  mutate(whethermissing=ifelse(NVISIT>=AESTDY & NVISIT<=AEENDY,1,0),
         AESERN=ifelse(AESER=='Y',1,0.3),
         AESEVN=mapping2[AESEV])%>%
  mutate(Possibility=ifelse(is.na(AESEV),0,whethermissing*AESERN*AESEVN))
Missdata=function(x,p){
  r= runif(1)
  if (r<= p) {x=NA}
  x
}
mnar_6a=mnar_6a%>%
  filter(AVISITN != 2)%>%
  mutate(PCHG = mapply(Missdata,PCHG,Possibility))
mark_missing <- function(df, SUBJID) {

    df$PCHG[df$SUBJID %in% SUBJID] <- NA

  return(df)
}
subjid=unique(mnar_6a$SUBJID[is.na(mnar_6a$PCHG)])
mnar_6a_final <- mark_missing(mnar_6a, subjid)
mnar_6a_final = mnar_6a_final%>%
  filter(AVISITN == 6)%>%
  group_by(SUBJID)%>%
  slice(1)

