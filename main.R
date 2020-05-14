library(Hmisc);library(glmnet);library(survivalROC);library(hdnom);library(boot);library(e1071);library(reshape2);library(caret)
library(rms);library(survival);library(survcomp);library(survminer);library(compareC);library(sva);library(MASS);library(ggplot2);library(TSHRC)
source('HLtest.r')
source('dca.r')
source('val.prob.ci.dec08.r')

### Read in the radiomics data and clinical data
data_pic <- read.csv('radiomics_data.csv');data_pic_ori <- data_pic
surv <- Surv(time=data_pic$time,event=data_pic$status)
data_cli <- read.csv('clinical_data.csv')

### All data is divided into four cohort according to the center
##  The variable named 'zh' is the primary cohort, 'dzl' is the validation cohort 1, 'fd' is the validation cohort 2, 'zs' is the validation cohort 3.
dzl<-which(data_pic$his=="dzl");fd<-which(data_pic$his=="fd");hn<-which(data_pic$his=="hn");yn<-which(data_pic$his=="yn");zs<-which(data_pic$his=="zs");zh<- c(hn,yn)  #×éºÏc/d
dzly <- surv[dzl,];fdy <- surv[fd,];zsy <- surv[zs,];zhy <- surv[zh,]

### Normalization
data_pic <- data_pic[,-c(1:3)]
pic_vec_num <- ncol(data_pic)
for (v in 1:pic_vec_num) {
  data_pic[zh,v] <- (data_pic[zh,v]-min(data_pic[zh,v]))/(max(data_pic[zh,v])-min(data_pic[zh,v]))
  data_pic[fd,v] <- (data_pic[fd,v]-min(data_pic[fd,v]))/(max(data_pic[fd,v])-min(data_pic[fd,v]))
  data_pic[dzl,v] <- (data_pic[dzl,v]-min(data_pic[dzl,v]))/(max(data_pic[dzl,v])-min(data_pic[dzl,v]))
  data_pic[zs,v] <- (data_pic[zs,v]-min(data_pic[zs,v]))/(max(data_pic[zs,v])-min(data_pic[zs,v]))
}
batch <- rep(1,nrow(data_pic));batch[dzl]<-2;batch[fd]<-3;batch[zs]<-4
combat_pic = ComBat(dat = data_pic[,1:pic_vec_num],batch = batch,mod = NULL,par.prior = T,mean.only = F)
data_pic[,1:pic_vec_num] <- t(combat_pic)  
zhx <- data_pic[zh,];dzlx <- data_pic[dzl,];fdx <- data_pic[fd,];zsx <- data_pic[zs,]

### Feature selection
pic_p <- c()
for (v in 1:pic_vec_num) {
  pictrain <- data.frame(zhx[,v],response=zhy)
  huozhe <- coxph(response~.,data=pictrain)
  pic_p[v] <- summary(huozhe)$coefficients[,5]
}
vecskm <- pic_p

trainx <- zhx;train_y <-zhy;testx1<-dzlx;testy1<-dzly;testx2<-fdx;testy2<-fdy;testx3<-zsx;testy3<-zsy
vectkmsig <- which(vecskm < 0.1)
vecneed <- vecskm[vectkmsig]
VECHAO1 <- order(vecneed)
VECHAO2 <- VECHAO1[1:floor(length(VECHAO1)*0.2)]
vechao <- vectkmsig[VECHAO2]
trainx0 <- trainx[,vechao];train_y <- train_y

highcorr <- findCorrelation(cor(trainx0),cutoff = 0.6,names = FALSE)
trainx1 <- trainx[,vechao];train_y <- train_y

### Lasso-cox for feature selection and radiomics signature building
jc_cox_glm <- cv.glmnet(as.matrix(trainx1),train_y,family = "cox",type.measure = "deviance",alpha=1,nfolds = nrow(trainx1))
vecp1 <- which(coef(jc_cox_glm,s="lambda.1se")!=0)

y_train_rs_pre <- predict(jc_cox_glm,as.matrix(trainx1),s=c("lambda.1se"),type="link")
y_test1_rs_pre <- predict(jc_cox_glm,as.matrix(testx1[,vechao]),s=c("lambda.1se"),type="link")
y_test2_rs_pre <- predict(jc_cox_glm,as.matrix(testx2[,vechao]),s=c("lambda.1se"),type="link")
y_test3_rs_pre <- predict(jc_cox_glm,as.matrix(testx3[,vechao]),s=c("lambda.1se"),type="link")
cdex <- rcorr.cens(y_train_rs_pre,train_y);c_train_rs <- 1-cdex[1];rs_train_upper <- cdex[3]/2*1.96+c_train_rs; rs_train_lower <- c_train_rs-cdex[3]/2*1.96
cdex <- rcorr.cens(y_test1_rs_pre,testy1);c_test1_rs <- 1-cdex[1];rs_test1_upper <- cdex[3]/2*1.96+c_test1_rs; rs_test1_lower <- c_test1_rs-cdex[3]/2*1.96
cdex <- rcorr.cens(y_test2_rs_pre,testy2);c_test2_rs <- 1-cdex[1];rs_test2_upper <- cdex[3]/2*1.96+c_test2_rs; rs_test2_lower <- c_test2_rs-cdex[3]/2*1.96
cdex <- rcorr.cens(y_test3_rs_pre,testy3);c_test3_rs <- 1-cdex[1];rs_test3_upper <- cdex[3]/2*1.96+c_test3_rs; rs_test3_lower <- c_test3_rs-cdex[3]/2*1.96

### Clinical signature building
cli_p <- c()
for (cli_i in 1:nrow(data_cli)) {
  names_single <- names(data_cli)
  clisig <- data.frame(pt=data_cli[,cli_i],response=train_y)
  huozhe <- coxph(response~.,data=clisig)
  print(names_single[cli_i])
  p_show <- summary(huozhe)$coefficients[,5]
  cli_p[cli_i] <- p_show
  print(p_show)
  hr_show <- summary(huozhe)$conf.int
  print(hr_show)
}
vect_ncol <- which(cli_p<0.05); 
clitrain <- data.frame(cli=data_cli[zh,vect_ncol],response=train_y)
clitestx1 <- data.frame(cli=data_cli[dzl,vect_ncol],response=testy1)
clitestx2 <- data.frame(cli=data_cli[fd,vect_ncol],response=testy2)
clitestx3 <- data.frame(cli=data_cli[zs,vect_ncol],response=testy3)
clicox1 <- coxph(response~.,data=clitrain)
clicox <- step(clicox1,direction = "backward",trace = 0)
y_train_cli_pre <- predict(clicox,newdata = clitrain,type = "risk")
y_test1_cli_pre <- predict(clicox,newdata = clitestx1,type = "risk")
y_test2_cli_pre <- predict(clicox,newdata = clitestx2,type = "risk")
y_test3_cli_pre <- predict(clicox,newdata = clitestx3,type = "risk")
cdex0 <- rcorr.cens(y_train_cli_pre,train_y);c_train_cli <- 1-cdex0; cli_train_upper <- cdex0[3]/2*1.96+c_train_cli; cli_train_lower <- c_train_cli-cdex0[3]/2*1.96
cdex1 <- rcorr.cens(y_test1_cli_pre,testy1);c_test1_cli <- 1-cdex1; cli_test1_upper <- cdex1[3]/2*1.96+c_test1_cli; cli_test1_lower <- c_test1_cli-cdex1[3]/2*1.96
cdex2 <- rcorr.cens(y_test2_cli_pre,testy2);c_test2_cli <- 1-cdex2; cli_test2_upper <- cdex2[3]/2*1.96+c_test2_cli; cli_test2_lower <- c_test2_cli-cdex2[3]/2*1.96
cdex3 <- rcorr.cens(y_test3_cli_pre,testy3);c_test3_cli <- 1-cdex3; cli_test3_upper <- cdex3[3]/2*1.96+c_test3_cli; cli_test3_lower <- c_test3_cli-cdex3[3]/2*1.96
  
### radiomics nomogram building
clipictrain <- data.frame(pic=y_train_rs_pre,cli=data_cli[zh,vect_ncol],response=train_y)
clipictestx1 <- data.frame(pic=y_test1_rs_pre,cli=data_cli[dzl,vect_ncol],response=testy1)
clipictestx2 <- data.frame(pic=y_test2_rs_pre,cli=data_cli[fd,vect_ncol],response=testy2)
clipictestx3 <- data.frame(pic=y_test3_rs_pre,cli=data_cli[zs,vect_ncol],response=testy3)
rncox1 <- coxph(response~.,data=clipictrain)
rncox <- step(rncox1,direction = "backward",trace = 0)
y_train_rn_pre <- predict(rncox,newdata = clipictrain,type = "risk")
y_test1_rn_pre <- predict(rncox,newdata = clipictestx1,type = "risk")
y_test2_rn_pre <- predict(rncox,newdata = clipictestx2,type = "risk")
y_test3_rn_pre <- predict(rncox,newdata = clipictestx3,type = "risk")
cdex0 <- rcorr.cens(y_train_rn_pre,train_y);c_train_rn <- 1-cdex0; rn_train_upper <- cdex0[3]/2*1.96+c_train_rn; rn_train_lower <- c_train_rn-cdex0[3]/2*1.96
cdex1 <- rcorr.cens(y_test1_rn_pre,testy1);c_test1_rn <- 1-cdex1; rn_test1_upper <- cdex1[3]/2*1.96+c_test1_rn; rn_test1_lower <- c_test1_rn-cdex1[3]/2*1.96
cdex2 <- rcorr.cens(y_test2_rn_pre,testy2);c_test2_rn <- 1-cdex2; rn_test2_upper <- cdex2[3]/2*1.96+c_test2_rn; rn_test2_lower <- c_test2_rn-cdex2[3]/2*1.96
cdex3 <- rcorr.cens(y_test3_rn_pre,testy3);c_test3_rn <- 1-cdex3; rn_test3_upper <- cdex3[3]/2*1.96+c_test3_rn; rn_test3_lower <- c_test3_rn-cdex3[3]/2*1.96

### Vincenzo's nomogram building
vntrain <- data.frame(pN=data_cli$pN[zh],pT=data_cli$pT[zh],chemo=data_cli$adjchemo[zh],surgtype=data_cli$surgerytype[zh],response=train_y)
vntestx1 <- data.frame(pN=data_cli$pN[dzl],pT=data_cli$pT[dzl],chemo=data_cli$adjchemo[dzl],surgtype=data_cli$surgerytype[dzl],response=testy1)
vntestx2 <- data.frame(pN=data_cli$pN[fd],pT=data_cli$pT[fd],chemo=data_cli$adjchemo[fd],surgtype=data_cli$surgerytype[fd],response=testy2)
vntestx3 <- data.frame(pN=data_cli$pN[zs],pT=data_cli$pT[zs],chemo=data_cli$adjchemo[zs],surgtype=data_cli$surgerytype[zs],response=testy3)
vncox <- coxph(response~.,data=vntrain)
y_train_vn_pre <- predict(rncox,newdata = vntrain,type = "risk")
y_test1_vn_pre <- predict(rncox,newdata = vntestx1,type = "risk")
y_test2_vn_pre <- predict(rncox,newdata = vntestx2,type = "risk")
y_test3_vn_pre <- predict(rncox,newdata = vntestx3,type = "risk")
cdex0 <- rcorr.cens(y_train_vn_pre,train_y);c_train_vn <- 1-cdex0; vn_train_upper <- cdex0[3]/2*1.96+c_train_vn; vn_train_lower <- c_train_vn-cdex0[3]/2*1.96
cdex1 <- rcorr.cens(y_test1_vn_pre,testy1);c_test1_vn <- 1-cdex1; vn_test1_upper <- cdex1[3]/2*1.96+c_test1_vn; vn_test1_lower <- c_test1_vn-cdex1[3]/2*1.96
cdex2 <- rcorr.cens(y_test2_vn_pre,testy2);c_test2_vn <- 1-cdex2; vn_test2_upper <- cdex2[3]/2*1.96+c_test2_vn; vn_test2_lower <- c_test2_vn-cdex2[3]/2*1.96
cdex3 <- rcorr.cens(y_test3_vn_pre,testy3);c_test3_vn <- 1-cdex3; vn_test3_upper <- cdex3[3]/2*1.96+c_test3_vn; vn_test3_lower <- c_test3_vn-cdex3[3]/2*1.96


### plot nomogram
dd=datadist(clipictrain)
options(datadist="dd")

f <- cph(response~ pic + cli.surgloca + cli.pN, x = T, y = T, data = clipictrain,surv = TRUE) 
survival <- Survival(f)
survival1 <-  function(x)survival(365,x)
survival3 <-  function(x)survival(365*2,x)
survival5 <-  function(x)survival(365*3,x)

nom2 <- nomogram(f, fun=list(survival1,survival3,survival5), fun.at = c(0.01,seq(.1,.9, by= .1),.95,.99),
                 funlabel=c('1 year Survival Probability','2 years Survival Probability','3 years Survival Probability'))
plot(nom2)

### plot KM
plottrain <- data.frame(rn=y_train_rn_pre,rs=y_train_rs_pre,cli=data_cli[zh,],time=data_pic_ori$time[zh],event=data_pic_ori$time$status[zh])
plottestx1 <- data.frame(rn=y_test1_rn_pre,rs=y_test1_rs_pre,cli=data_cli[dzl,],time=data_pic_ori$time[dzl],event=data_pic_ori$time$status[dzl])
plottestx2 <- data.frame(rn=y_test2_rn_pre,rs=y_test2_rs_pre,cli=data_cli[fd,],time=data_pic_ori$time[fd],event=data_pic_ori$time$status[fd])
plottestx3 <- data.frame(rn=y_test3_rn_pre,rs=y_test3_rs_pre,cli=data_cli[zs,],time=data_pic_ori$time[zs],event=data_pic_ori$time$status[zs])
data_plot <- rbind(plottrain,plottestx1,plottestx2,plottestx3);survy <- Surv(time=data_plot$time,event = data_plot$event)

his <- zh
veclie <- data_plot$rn[his]
hrmatrix <- hazard.ratio(x=veclie,surv.time = data_plot$time[his],surv.event = data_plot$event[his])

cutoff <- median(data_plot$rn[zh])
if(cutoff>=1){
  veclie[which(veclie <= cutoff)] <- 0
  veclie[which(veclie > cutoff)] <- 1
}
if(cutoff<1){
  veclie[which(veclie > cutoff)] <- 1
  veclie[which(veclie <= cutoff)] <- 0
}
veclie <-  as.factor(veclie)
veclie1 <- veclie
survkm <- survy[his]
kmdata <- data.frame(surv = survkm,vect1 = veclie1)
kmmodel <- survfit(surv~vect1,data=kmdata)
ggsurv <- ggsurvplot(
  fit=kmmodel,data=kmdata,  size = 1.5, risk.table = TRUE,     
  pval = TRUE,  conf.int = F, xlab = "t (days)",   break.time.by = 182.5,  xlim=c(0,365*3)
)
ggsurv$plot <- ggsurv$plot 
ggsurv <- ggpar(
  ggsurv, font.title = c(12, "bold", "black"),         
  font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"),    
  font.xtickslab = c(14, "plain", "black"),  legend = "top",
  font.ytickslab = c(14, "plain", "black"))
survp <- ggsurv
print(survp)
ggsave(paste("KM-rn-zh",".tiff"), print(survp))

### subgroup analysis
cutoff <- median(data_plot$rs[zh])
data_yazu3 <- data_plot[his];survy3 <- survy[his2]
veclie <- data_yazu3$rs
if(cutoff>=1){
  veclie[which(veclie <= cutoff)] <- 0
  veclie[which(veclie > cutoff)] <- 1
}
if(cutoff<1){
  veclie[which(veclie > cutoff)] <- 1
  veclie[which(veclie <= cutoff)] <- 0
}
veclie <-  as.factor(veclie);veclie3 <- veclie
cli_order <- c(4,13,14)
hrdata <- data.frame();hrx1 <- data.frame();siglength <- data.frame()
name_vec <- colnames(data_yazu3)
for (st in 1:3) {
  i=cli_order[st]+2
  vec_values <- data_yazu3[,i]
  vec_values <- unique(vec_values)
  for (j in 1:length(vec_values)) {
    vec_yazu <- data_yazu3[,i]
    cli_vec <- which(vec_yazu== vec_values[j]);data_yazu1 <- data_yazu3[cli_vec,];survy1 <- survy3[cli_vec,]
    veclie2 <- veclie3[cli_vec]
    # high_risk <- which(veclie2 == "1");data_yazu2 <- data_yazu1[high_risk,];survy2 <- survy1[high_risk,]
    # veclie1 <- data_yazu2$cli.afthua

    hrmatrix <- hazard.ratio(x=data_yazu1$rs,  surv.time = data_yazu1$time,  surv.event = data_yazu1$event)
    hrdata1 <- data.frame(vectname = paste(name_vec[i],vec_values[j]), hr = hrmatrix$hazard.ratio, low = hrmatrix$lower, up = hrmatrix$upper)
    hrx1 <- rbind(hrx1,hrdata1)
    
    kmdata <- data.frame(surv = survy1,vect1 = veclie2)
    kmmodel <- survfit(surv~vect1,data = kmdata)
    ggsurv <- ggsurvplot(
      fit=kmmodel,data=kmdata,  size = 1.5, risk.table = TRUE,
      pval = TRUE,  conf.int = F, xlab = "t (days)",   break.time.by = 365/2,  xlim=c(0,365*3)
    )
    ggsurv$plot <- ggsurv$plot
    ggsurv <- ggpar(
      ggsurv, font.title = c(12, "bold", "black"),
      font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"),
      font.xtickslab = c(14, "plain", "black"),  legend = "top",
      font.ytickslab = c(14, "plain", "black"))
    survp <- ggsurv
    print(survp)
    summary(kmmodel)
    ggsave(paste("KM-subgroup",name_vec[i],vec_values[j],".tiff",sep="_"), print(survp))
    ggsave(filename = paste("KM-subgroup",name_vec[i],vec_values[j],".eps",sep="_"),device = "eps",print(survp))
  }
}

### The KM analysis of adjuvant chemotherapy in subgroups
cutoff <- median(data_plot$rs[zh])
data_yazu3 <- data_plot;survy3 <- survy
veclie <- data_yazu3$rs
if(cutoff>=1){
  veclie[which(veclie <= cutoff)] <- 0
  veclie[which(veclie > cutoff)] <- 1
}
if(cutoff<1){
  veclie[which(veclie > cutoff)] <- 1
  veclie[which(veclie <= cutoff)] <- 0
}
veclie <-  as.factor(veclie);veclie3 <- veclie
cli_order <- c(4,13,14)
hrdata <- data.frame();hrx1 <- data.frame();siglength <- data.frame()
name_vec <- colnames(data_yazu3)
for (st in 1:3) {
  i=cli_order[st]+2
  vec_values <- data_yazu3[,i]
  vec_values <- unique(vec_values)
  for (j in 1:length(vec_values)) {
    vec_yazu <- data_yazu3[,i]
    cli_vec <- which(vec_yazu== vec_values[j]);data_yazu1 <- data_yazu3[cli_vec,];survy1 <- survy3[cli_vec,]
    veclie2 <- veclie3[cli_vec]
    high_risk <- which(veclie2 == "1");data_yazu2 <- data_yazu1[high_risk,];survy2 <- survy1[high_risk,]
    veclie1 <- data_yazu2$cli.adjchemo
    
    hrmatrix <- hazard.ratio(x=veclie1,  surv.time = data_yazu2$time,  surv.event = data_yazu2$event)
    hrdata1 <- data.frame(vectname = paste(name_vec[i],vec_values[j]), hr = hrmatrix$hazard.ratio, low = hrmatrix$lower, up = hrmatrix$upper)
    hrx1 <- rbind(hrx1,hrdata1)
    
    kmdata <- data.frame(surv = survy2,vect1 = veclie1)
    kmmodel <- survfit(surv~vect1,data = kmdata)
    ggsurv <- ggsurvplot(
      fit=kmmodel,data=kmdata,  size = 1.5, risk.table = TRUE,
      pval = TRUE,  conf.int = F, xlab = "t (days)",   break.time.by = 365/2,  xlim=c(0,365*3)
    )
    ggsurv$plot <- ggsurv$plot
    ggsurv <- ggpar(
      ggsurv, font.title = c(12, "bold", "black"),
      font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"),
      font.xtickslab = c(14, "plain", "black"),  legend = "top",
      font.ytickslab = c(14, "plain", "black"))
    survp <- ggsurv
    print(survp)
    summary(kmmodel)
    ggsave(paste("KM-subgroup-adjchemo",name_vec[i],vec_values[j],".tiff",sep="_"), print(survp))
    ggsave(filename = paste("KM-subgroup-adjchemo",name_vec[i],vec_values[j],".eps",sep="_"),device = "eps",print(survp))
  }
}

### Hierarchical analysis
hrdata <- data.frame()
cli_order <- c(8,13,14)
high_RS <- which(veclie== "1");data_yazu1 <- data_yazu3[high_RS,];survy1 <- survy3[high_RS,]
for (st in 1:3) {
  i=cli_order[st]+2
  veclie1 <- data_yazu1[,i]
  kmdata <- data.frame(surv = survy1,vect1 = veclie1)

  hrmatrix <- hazard.ratio(x=veclie1,  surv.time = data_yazu1$time,  surv.event = data_yazu1$event)
  hrdata1 <- data.frame(vectname = name_vec[i], hr = hrmatrix$hazard.ratio, low = hrmatrix$lower, up = hrmatrix$upper)
  hrdata <- rbind(hrdata,hrdata1)
  
  kmmodel <- survfit(survy1~veclie1)
  ggsurv <- ggsurvplot(
    fit=kmmodel,data=kmdata,  size = 1.5, risk.table = TRUE,
    pval = TRUE,  conf.int = F, xlab = "t (days)",   break.time.by = 365/2,  xlim=c(0,365*3)
  )
  ggsurv$plot <- ggsurv$plot
  ggsurv <- ggpar(
    ggsurv, font.title = c(12, "bold", "black"),
    font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"),
    font.xtickslab = c(14, "plain", "black"),  legend = "top",
    font.ytickslab = c(14, "plain", "black"))
  survp <- ggsurv
  print(survp)
  ggsave(filename = paste("KM-high-RS",name_vec[i],".eps",sep="_"),device = "eps",print(survp))
}

### Time-dependent ROC
library(timeROC)
timerocdata <- data.frame(time=data_plot$time[his],status=data_plot$event[his],vec1=as.vector(data_plot$rs[his]))
ROC_rt<- timeROC(T=timerocdata$time, delta=timerocdata$status,
                 marker=timerocdata$vec1, cause=1,
                 weighting='marginal',
                 times=c(365,365*2,365*3),ROC=TRUE)

plot(ROC_rt,time=365,col='red',title=FALSE,lwd=2)
plot(ROC_rt,time=365*2,col='green',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=365*3,col='blue',add=TRUE,title=FALSE,lwd=2)

legend('bottom',
       c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],2)),
         paste0('AUC at 2 years: ',round(ROC_rt$AUC[2],2)),
         paste0('AUC at 3 years: ',round(ROC_rt$AUC[3],2))),
       col=c('red','green','blue'),lwd=5,bty = 'n')

### Desicion curve
dcatrain <- data.frame(sig=y_train_rs_pre,cli.surgery_location=data_cli$surgloca[zh],cli.pN=data_cli$pN[zh],time=data_plot$time[zh],event=data_plot$event[zh])
dcatestx1 <- data.frame(sig=y_test1_rs_pre,cli.surgery_location=data_cli$surgloca[dzl],cli.pN=data_cli$pN[dzl],time=data_plot$time[dzl],event=data_plot$event[dzl])
dcatestx2 <- data.frame(sig=y_test2_rs_pre,cli.surgery_location=data_cli$surgloca[fd],cli.pN=data_cli$pN[fd],time=data_plot$time[fd],event=data_plot$event[fd])
dcatestx3 <- data.frame(sig=y_test3_rs_pre,cli.surgery_location=data_cli$surgloca[zs],cli.pN=data_cli$pN[zs],time=data_plot$time[zs],event=data_plot$event[zs])
data_dca <- rbind(dcatrain,dcatestx1,dcatestx2,dcatestx3)
dcs <- data_dca[,c(ncol(data_dca)-1,ncol(data_dca))]
dcs_sta<- dcs[,2]
dcs_sta[which((dcs[,1]>(365*3))&(dcs[,2]==1))] <- 0

par(bty="n",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lwd=3,lty=1,col=1)
dca.train   <- dca(yvar= dcs_sta[his], xmatrix=1-survest(f,newdata=data_dca[his,],times = 365*3)$surv, 
                   prob="Y",ymin=-0.1,ymax=0.3) 

train.xindex<-dca.train$threshold
train.yindex<-dca.train$modelp1
plot(xlim=c(0,100),ylim=c(-0.1,0.3),x=dca.train$threshold,y=dca.train$none, 
     col="black",lwd=2,type = "l",
     xlab = "Threshold Probability(%)", ylab= "Net benefit",
     cex=1,cex.axis=1,cex.lab=1)
lines(x=train.xindex, y=train.yindex,col="#FF0B97",lwd=2,type = "l")
lines(x=dca.train$threshold, y=dca.train$all,col="blue",lwd=2,type = "l")
legend("topright",bty="n",cex=0.8,legend=c("None","All","Nomogram"), col=c("black","blue","#FF0B97"),lwd=2)

### calibration curve
his <- zh;evaluaten=50
par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
traincal <- val.surv(f,S=Surv(data_dca$time[his],data_dca$event[his]),newdata=data_dca[his,], u=365*3, evaluate=100)

res = groupkm(traincal$p, Surv(data_dca$time[his],data_dca$event[his]), m=evaluaten, u=1095, pl=T, add=F,xlim=c(0,1),
              ylim=c(0,1),errbar=T, errbar.col="blue",cex.axis=2,cex.lab=2,font=1,
              xlab="Nomogram-Predicted Probability of DM",lwd = 2,
              ylab="Actual DM (proportion)",cex.subtitle=F,col="blue")
abline(0,1,lty=2)
lines(res[,c('x','KM')],type= 'o',lwd = 2,col="blue",pch = 16)

##2 years
par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
traincal <- val.surv(f,S=Surv(data_dca$time[his],data_dca$event[his]),newdata=data_dca[his,], u=365*2, evaluate=evaluaten)

res = groupkm(traincal$p, Surv(data_dca$time[his],data_dca$event[his]), m=evaluaten, u=365*2, pl=T, add=T,xlim=c(0,1),
              ylim=c(0,1),errbar=T, errbar.col="green",cex.axis=2,cex.lab=2,font=1,
              xlab="Nomogram-Predicted Probability of DM",lwd = 2,
              ylab="Actual DM (proportion)",cex.subtitle=F,col="green")
abline(0,1,lty=2)
lines(res[,c('x','KM')],type= 'o',lwd = 2,col="green",pch = 16)

##1 year
par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
traincal <- val.surv(f,S=Surv(data_dca$time[his],data_dca$event[his]),newdata=data_dca[his,], u=365*1, evaluate=evaluaten)

res = groupkm(traincal$p, Surv(data_dca$time[his],data_dca$event[his]), m=evaluaten, u=365*1, pl=T, add=T,xlim=c(0,1),
              ylim=c(0,1),errbar=T, errbar.col="red",cex.axis=2,cex.lab=2,font=1,
              xlab="Nomogram-Predicted Probability of DM",lwd = 2,
              ylab="Actual DM (proportion)",cex.subtitle=F,col="red")
abline(0,1,lty=2)
lines(res[,c('x','KM')],type= 'o',lwd = 2,col="red",pch = 16)

legend('topleft',c(paste0('1-year survival '),paste0('2-year survival '),paste0('3-year survival ')),cex=1.5,
       col=c('red','green','blue'),lwd=5,bty = 'n')




