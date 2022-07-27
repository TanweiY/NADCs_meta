########################################## Data preparation ########################################################
SIRR <- read_excel("SIRR.xlsx", col_types = c("numeric",
                                              "text", "text", "numeric", "text", "text",
                                              "text", "text", "text", "text", "text",
                                              "numeric", "text", "text", "text", "numeric",
                                              "text", "text", "text", "numeric", "numeric",
                                              "numeric", "numeric", "numeric", "text",
                                              "numeric", "text", "text", "text", "text",
                                              "text", "text", "numeric", "numeric",
                                              "numeric", "text"))


SIRR<-SIRR %>% remove_empty(which = c('rows', 'cols')) # 去掉空行

dput(names(SIRR))
cols<-c("Study", "Cohortna",  "Publicationty", "Studydesig",
        "Country2", "WHO region", "Incomecountry", "Setting",
        "Enrollment Duration", "Follow-up/linkage Duration (Totai/Frequency)",
        "cancer case identification", "HIV/AIDS identification",
        "Selectedpotime",
        "origNADCstyp", "primNADCstyp", "secNADCstyp",
        "Adjustment")

SIRR[cols]<-lapply(SIRR[cols], factor) # 统一设置类型为factor
length(unique(SIRR$Study)) ## 76


##############  SIR overall meta analysis ##############

siroverall_update<- escalc(measure="IRLN", xi=O, ti=E, data=siroverall_update)
write_xlsx(siroverall_update, path = "SIR/Overall/siroverall_update.xlsx")
overall_prim<-siroverall_update[!is.na(siroverall_update$primNADCstyp), ]
overall_prim$primNADCstyp<-as.factor(as.character(overall_prim$primNADCstyp))
summary(overall_prim$primNADCstyp)
levels(overall_prim$primNADCstyp)

## calculate robust for primary cancer
num<-length(levels(overall_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(overall_prim, overall_prim$primNADCstyp==levels(overall_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(overall_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

sirresult<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)


### 2. calculate robust for secondary cancers
overall_sec<-subset(siroverall_update, !is.na(siroverall_update$secNADCstyp))

overall_sec<-subset(overall_sec, overall_sec$secNADCstyp!='Gallbladder' &
                      overall_sec$secNADCstyp!='Nasopharynx'&
                      overall_sec$secNADCstyp!='Tonsil')

overall_sec$secNADCstyp<-as.factor(as.character(overall_sec$secNADCstyp))

summary(overall_sec$secNADCstyp)

num_2<-length(levels(overall_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(overall_sec, secNADCstyp==levels(overall_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(overall_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

sirresult_2<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

sirresult_2$Analysis<-'robu_sec'

#### 3. calculate cancers no need robust

sir_secnorob<-subset(siroverall_update, secNADCstyp =='Gallbladder' |
                       secNADCstyp =='Nasopharynx'|
                       secNADCstyp =='Tonsil'
)

sir_secnorob$secNADCstyp<-as.factor(as.character(sir_secnorob$secNADCstyp))
summary(sir_secnorob$secNADCstyp)

num_3<-length(levels(sir_secnorob$secNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(sir_secnorob, secNADCstyp==levels(sir_secnorob$secNADCstyp)[i]),method="ML")

  Name[i]<-levels(sir_secnorob$secNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

sirresult_3<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)


sirresult_final<-bind_rows(sirresult, sirresult_2, sirresult_3)

sirresult_final[,2:4] <-round(sirresult_final[,2:4],2)
sirresult_final[,7] <-round(sirresult_final[,7],0)

sirresult_final$'SIR (95%CI)' <-paste0(sirresult_final$SIR," ","(", sirresult_final$LL, "-", sirresult_final$UL, ")")
sirresult_final$lgSIR<-log10(sirresult_final$SIR)
sirresult_final$lgul<-log10(sirresult_final$UL)
sirresult_final$lgll<-log10(sirresult_final$LL)

sirresult_final<-sirresult_final[match(sir_name$Name, sirresult_final$Name),]


## add number of cancers
siroverall_update <- read_excel("SIR/Overall/siroverall_update.xlsx")

num<-data.frame(tapply(siroverall_update$O, siroverall_update$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)


num_sec<-data.frame(tapply(siroverall_update$O, siroverall_update$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:14)
num_sec<-num_sec[-4, ]

num<-bind_rows(num, num_sec)
colnames(sirresult_final)[1]<-'Name'
sirresult_final<-merge(sirresult_final, num, by='Name')

save(sirresult_final, file='sirresult_final.RData')
write_xlsx(sirresult_final, path = "SIR/Overall/sirresult_final.xlsx")

write_xlsx(sirresult_final, path = 'Data_analysis_new/Result_SIRoverall.xlsx')

###################### calculate SMR overall ######################
SMRoverall <- read_excel("SMRoverall.xlsx")
SMRoverall<- escalc(measure="IRLN", xi=O, ti=E, data=SMRoverall) # calculate effect size and variance
write_xlsx(SMRoverall, path = "SMRoverall.xlsx")
dput(names(SMRoverall))
cols<-c("study_id",  "primNADCstyp", "time_duration")
SMRoverall[cols]<-lapply(SMRoverall[cols], factor)
SMR_data<-subset(SMRoverall, select = c("study_id", "Study", "Cohortame", "Cohortab", "time_duration",
                                        "group", "origNADCstyp", "primNADCstyp",
                                        "secNADCstyp", "O", "E", "extract_SMR", "extract_ll", "extract_ul",
                                        "Adjustment", "yi", "vi"))

## robust
options(scipen = 9999999)

SMR_data_robu<-subset(SMR_data, SMR_data$primNADCstyp =='All NADCs' |
                        SMR_data$primNADCstyp =='Colon and Rectum' |
                        SMR_data$primNADCstyp =='Liver' |
                        SMR_data$primNADCstyp =='Trachea, bronchus, and lung')
SMR_data_robu$primNADCstyp<-as.factor(as.character(SMR_data_robu$primNADCstyp))
summary(SMR_data_robu$primNADCstyp)

num<-length(levels(SMR_data_robu$primNADCstyp))
SMR<-0
Name<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsmr<-robu(formula = yi ~ 1,
               data = subset(SMR_data_robu, SMR_data_robu$primNADCstyp==levels(SMR_data_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(SMR_data_robu$primNADCstyp)[i]
  SMR[i]<-exp(fitsmr$reg_table[1, 2])
  LL[i]<-exp(fitsmr$reg_table[1, 7])
  UL[i]<-exp(fitsmr$reg_table[1, 8])
  num_study[i]<-fitsmr$N
  num_outcome[i]<-fitsmr$M
  I2[i]<-fitsmr$mod_info$I.2
  Tau[i]<-fitsmr$mod_info$tau.sq
}

smrresult_robu<-data.frame(Name, SMR, LL, UL,num_study, num_outcome, I2, Tau)



# rma
"Brain and CNS"|
  "Leukemia "|
  "Melanoma of skin"|

  SMR_data_rma<-subset(SMR_data, SMR_data$primNADCstyp == "Brain and CNS"|
                         SMR_data$primNADCstyp == 'Leukemia'|
                         SMR_data$primNADCstyp == 'Melanoma of skin'|
                         SMR_data$primNADCstyp == "Pancreas")

SMR_data_rma$primNADCstyp<-as.factor(as.character(SMR_data_rma$primNADCstyp))
summary(SMR_data_rma$primNADCstyp)

num<-length(levels(SMR_data_rma$primNADCstyp))

Name<-0
SMR<-0
LL<-0
UL<-0

num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
PHeterogenity<-0
I2<-0

for (i in 1: num) {
  fit_norob<-rma(yi, vi, measure = "IRLN",
                 data =subset(SMR_data_rma, SMR_data_rma$primNADCstyp==levels(SMR_data_rma$primNADCstyp)[i]),method="ML")
  Name[i]<-levels(SMR_data_rma$primNADCstyp)[i]
  SMR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_level1[i]<-fit_norob$k
  num_level2[i]<-fit_norob$k
  num_level3[i]<-fit_norob$k


  I2[i]<-fit_norob$I2

}

smrresult_rma<-data.frame(Name, SMR, LL, UL, num_outcome, num_level1, num_level2,
                          num_level3, I2)




# mutlilevel"

SMR_data_mul<-subset(SMR_data, primNADCstyp=='Anus and anal canal'|
                       primNADCstyp== "Hodgkin lymphoma"|
                       primNADCstyp== "Stomach")

SMR_data_mul$primNADCstyp<-as.factor(as.character(SMR_data_mul$primNADCstyp))

num<-length(levels(SMR_data_mul$primNADCstyp))
SMR<-0
Name<-0
SMR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
PHeterogenity<-0
I2<-0

for (i in 1:num) {
  data = subset(SMR_data_mul, SMR_data_mul$primNADCstyp==levels(SMR_data_mul$primNADCstyp)[i])
  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | study_id/time_duration,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(SMR_data_mul$primNADCstyp)[i]
  SMR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]

  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))

}

SMR_data_mul_result<-data.frame(Name, SMR, LL, UL, num_outcome, num_level1, num_level2,
                                num_level3,I2)

smr_result<-bind_rows(smrresult_robu,
                      smrresult_rma,
                      SMR_data_mul_result)


######## summarize number of cancer cases
num<-data.frame(tapply(SMR_data$O, SMR_data$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)

smr<-merge(smr_result, num, by='Name', all.x = T)
smr<-within(smr, {
  lgsmr<-log10(SMR)
  lgll<-log10(LL)
  lgul<-log10(UL)
  final_SMR <-paste0(round(SMR, 2)," ","(", round(LL, 2), "-", round(UL,2), ")")
  I2<-round(I2, digits = 0)

})

write_xlsx(smr, path ="SMR_result.xlsx")

############## subgroup_analysis_ categorical variables ##############
######################################## gender ########################################
gender <- read_excel("SIR/subgroup/Categorical/gender.xlsx")
attach(gender)
gender$study_id<-ifelse(!is.na(Cohortab), Cohortab, Study)
detach(gender)
gender<- escalc(measure="IRLN", xi=O, ti=E, data=gender)

write_xlsx(gender, path = 'SIR/subgroup/Categorical/gender.xlsx')

######################### Men ###############
men<-subset(gender, groupsex=='Men')

# Primary robust
## rma
men_prim<-subset(men, primNADCstyp !='Bone and joints' &
                   primNADCstyp !='Breast' &
                   primNADCstyp !='Head and neck' &
                   primNADCstyp !='Mesothelial and soft tissue'&
                   primNADCstyp !='Nasal cavity, middle ear, and accessory sinuses'&
                   primNADCstyp !='Small intestine'&
                   primNADCstyp !='Thymus, heart, mediastinum, and pleura')
men_prim$primNADCstyp<-as.factor(as.character(men_prim$primNADCstyp))
summary(men_prim$primNADCstyp)

num<-length(levels(men_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(men_prim, men_prim$primNADCstyp==levels(men_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(men_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

men_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
men_sirpri$group<-'Men'

### robust for secondary
men_sec<-subset(men, !is.na(men$secNADCstyp))

men_sec<-subset(men_sec, men_sec$secNADCstyp!='Gallbladder' &
                  men_sec$secNADCstyp!='Lip' &
                  men_sec$secNADCstyp!='Testis seminoma' &
                  men_sec$secNADCstyp!='Nasopharynx'&
                  men_sec$secNADCstyp!='Tongue'&
                  men_sec$secNADCstyp!='Tonsil')

men_sec$secNADCstyp<-as.factor(as.character(men_sec$secNADCstyp))

summary(men_sec$secNADCstyp)

num_2<-length(levels(men_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(men_sec, secNADCstyp==levels(men_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(men_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

men_sirsec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
men_sirsec$group<-'Men'

###rma

men_indep<-subset(men, primNADCstyp =='Bone and joints' |
                    primNADCstyp =='Breast' |
                    primNADCstyp =='Head and neck' |
                    primNADCstyp =='Mesothelial and soft tissue'|
                    primNADCstyp =='Nasal cavity, middle ear, and accessory sinuses'|
                    primNADCstyp =='Small intestine'|
                    primNADCstyp =='Thymus, heart, mediastinum, and pleura')

### primary
men_indep$primNADCstyp<-as.factor(as.character(men_indep$ primNADCstyp))
summary(men_indep$primNADCstyp)

num_3<-length(levels(men_indep$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(men_indep, primNADCstyp==levels(men_indep$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(men_indep$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

men_sirs_indep_prim<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

men_sirs_indep_prim$group<-'Men'


### secondary

men_indepsec<-subset(men, secNADCstyp=='Nasopharynx'|
                       secNADCstyp=='Tongue'|
                       secNADCstyp=='Tonsil')

men_indepsec$secNADCstyp<-as.factor(as.character(men_indepsec$secNADCstyp))
summary(men_indepsec$secNADCstyp)

num_3<-length(levels(men_indepsec$secNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(men_indepsec, secNADCstyp==levels(men_indepsec$secNADCstyp)[i]),method="ML")

  Name[i]<-levels(men_indepsec$secNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

men_sirs_indep_sec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

men_sirs_indep_sec$group<-'Men'

mensir_final<-bind_rows(men_sirpri, men_sirsec, men_sirs_indep_prim, men_sirs_indep_sec) # bind results together

mensir_final[,2:4] <-round(mensir_final[,2:4],2)
mensir_final[,7] <-round(mensir_final[,7],0)

######## summarize number of cancer cases
num<-data.frame(tapply(men$O, men$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)

num_sec<-data.frame(tapply(men$O, men$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:10)

num<-bind_rows(num, num_sec)
mensf<-merge(mensir_final, num, by='Name', all.x = T)

save(mensf, file='mensir.RData')

####################################################### Women #############################################
gender <- read_excel("SIR/subgroup/Categorical/gender.xlsx")

women<-subset(gender, groupsex=='Women')

# Primary robust
## rma
women_prim<-subset(women, primNADCstyp !='Bone and joints' &
                     primNADCstyp !='Endometrium' &  # remove
                     primNADCstyp !='Esophagus' &
                     primNADCstyp !='Head and neck' &
                     primNADCstyp !='Kidney and other urinary' &
                     primNADCstyp !='Larynx' &
                     primNADCstyp !='Ovary'&
                     primNADCstyp !='Pancreas'&
                     primNADCstyp !='Thyroid'&
                     primNADCstyp !='Eye and adnexa'& #
                     primNADCstyp !='Nasal cavity, middle ear, and accessory sinuses'& #remove
                     primNADCstyp !='Small intestine'&
                     primNADCstyp !='Thymus, heart, mediastinum, and pleura')# remove
women_prim$primNADCstyp<-as.factor(as.character(women_prim$primNADCstyp))
summary(women_prim$primNADCstyp)

num<-length(levels(women_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(women_prim, women_prim$primNADCstyp==levels(women_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(women_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

women_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
women_sirpri$group<-'Women'

### robust for secondary
women_sec<-subset(women, !is.na(women$secNADCstyp))

women_sec<-subset(women_sec,
                  secNADCstyp!='Kidney and renal pelvis' &
                    secNADCstyp!='Women breast' & # remove
                    secNADCstyp!='Nasopharynx'&
                    secNADCstyp!='Tongue')

women_sec$secNADCstyp<-as.factor(as.character(women_sec$secNADCstyp))

summary(women_sec$secNADCstyp)

num_2<-length(levels(women_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(women_sec, secNADCstyp==levels(women_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(women_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

women_sirsec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
women_sirsec$group<-'Women'

###rma

women_indep<-subset(women, primNADCstyp =='Bone and joints' |
                      primNADCstyp =='Esophagus' |
                      primNADCstyp =='Head and neck' |
                      primNADCstyp =='Kidney and other urinary' |
                      primNADCstyp =='Larynx' |
                      primNADCstyp =='Pancreas'|
                      primNADCstyp =='Thyroid'|
                      primNADCstyp =='Eye and adnexa')

### primary
women_indep$primNADCstyp<-as.factor(as.character(women_indep$ primNADCstyp))
summary(women_indep$primNADCstyp)

num_3<-length(levels(women_indep$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(women_indep, primNADCstyp==levels(women_indep$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(women_indep$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

women_sirs_indep_prim<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

women_sirs_indep_prim$group<-'Women'


### secondary

women_indep_sec<-subset(women, secNADCstyp=='Kidney and renal pelvis' |
                          secNADCstyp=='Nasopharynx'|
                          secNADCstyp=='Tongue')


women_indep_sec$secNADCstyp<-as.factor(as.character(women_indep_sec$secNADCstyp))
summary(women_indep_sec$secNADCstyp)

num_3<-length(levels(women_indep_sec$secNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(women_indep_sec, secNADCstyp==levels(women_indep_sec$secNADCstyp)[i]),method="ML")

  Name[i]<-levels(women_indep_sec$secNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

women_sirs_indep_sec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

women_sirs_indep_sec$group<-'Women'

womensir_final<-bind_rows(women_sirpri, women_sirsec, women_sirs_indep_prim, women_sirs_indep_sec) # bind results together

womensir_final[,2:4] <-round(womensir_final[,2:4],2)
womensir_final[,7] <-round(womensir_final[,7],0)

######## summarize number of cancer cases
num<-data.frame(tapply(women$O, women$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)

num_sec<-data.frame(tapply(women$O, women$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:7)

num<-bind_rows(num, num_sec)
womensf<-merge(womensir_final, num, by='Name', all.x = T)

save(womensf, file='womensir.RData')

gendersf<-bind_rows(mensf, womensf)

gendersf<-arrange(gendersf, Name)

gendersf$'SIR (95%CI)' <-paste0(gendersf$SIR," ","(", gendersf$LL, "-", gendersf$UL, ")")
gendersf$lgSIR<-log10(gendersf$SIR)
gendersf$lgul<-log10(gendersf$UL)
gendersf$lgll<-log10(gendersf$LL)


write_xlsx(gendersf, path ="SIR/subgroup/Categorical/sirgender.xlsx" )

###################################### HIV transmission group ######################################
transmigroup <- read_excel("SIR/subgroup/Categorical/transmigroup.xlsx")
attach(transmigroup)
transmigroup$study_id<-ifelse(!is.na(Cohortab), Cohortab, Study)
detach(transmigroup)
transmigroup<- escalc(measure="IRLN", xi=O, ti=E, data=transmigroup)

write_xlsx(transmigroup, path = 'SIR/subgroup/Categorical/transmigroup.xlsx')

####################### MSM ###########
msm<-subset(transmigroup, groupmsm=='MSM')

# Primary robust
msm_prim<-subset(msm,
                 primNADCstyp !='Bladder' &
                   primNADCstyp !='Bone and joints' &
                   primNADCstyp !='Kidney and other urinary' &
                   primNADCstyp !='Leukaemia' &
                   primNADCstyp !='Lip, oral cavity and pharynx'& ## rma
                   primNADCstyp !='Liver, gallbladder, and biliary tract'& ## 直接liver
                   primNADCstyp !='Melanoma of skin'& # rma
                   primNADCstyp !='Mesothelial and soft tissue' &
                   primNADCstyp !='Multiple myeloma' &
                   primNADCstyp !='Pancreas'&
                   primNADCstyp !='Penis'&
                   primNADCstyp !='Eye and adnexa'&
                   primNADCstyp !='Brain and CNS'&
                   primNADCstyp !='Colon, rectosigmoid junction, rectum, anus and anal canal'&
                   primNADCstyp !='Skin nonmelanoma'& ## rma
                   primNADCstyp !='Thyroid')
msm_prim$primNADCstyp<-as.factor(as.character(msm_prim$primNADCstyp))
summary(msm_prim$primNADCstyp)

num<-length(levels(msm_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(msm_prim, msm_prim$primNADCstyp==levels(msm_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(msm_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

msm_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
msm_sirpri$group<-'MSM'

### robust for secondary
msm_sec<-subset(msm, !is.na(msm$secNADCstyp))

msm_sec<-subset(msm_sec,
                msm_sec$secNADCstyp!='Kidney and renal pelvis' &
                  msm_sec$secNADCstyp!='Liver' & # rma
                  msm_sec$secNADCstyp!='Salivary glands' &
                  msm_sec$secNADCstyp!='Testis seminoma')

msm_sec$secNADCstyp<-as.factor(as.character(msm_sec$secNADCstyp))

summary(msm_sec$secNADCstyp)

num_2<-length(levels(msm_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(msm_sec, secNADCstyp==levels(msm_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(msm_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

msm_sirsec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
msm_sirsec$group<-'MSM'

###rma
msm_indep<-subset(msm, primNADCstyp =='Lip, oral cavity and pharynx' |
                    primNADCstyp =='Melanoma of skin' |
                    primNADCstyp =='Skin nonmelanoma'|
                    primNADCstyp =='Brain and CNS' )


### primary
msm_indep$primNADCstyp<-as.factor(as.character(msm_indep$ primNADCstyp))
summary(msm_indep$primNADCstyp)

num_3<-length(levels(msm_indep$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(msm_indep, primNADCstyp==levels(msm_indep$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(msm_indep$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

msm_sirs_indep_prim<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

msm_sirs_indep_prim$group<-'MSM'


### secondary

msm_indepsec<-subset(msm, secNADCstyp=='Liver')

msm_indepsec$secNADCstyp<-as.factor(as.character(msm_indepsec$secNADCstyp))
summary(msm_indepsec$secNADCstyp)

num_3<-length(levels(msm_indepsec$secNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(msm_indepsec, secNADCstyp==levels(msm_indepsec$secNADCstyp)[i]),method="ML")

  Name[i]<-levels(msm_indepsec$secNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

msm_sirs_indep_sec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

msm_sirs_indep_sec$group<-'MSM'

msmsir_final<-bind_rows(msm_sirpri, msm_sirsec, msm_sirs_indep_prim, msm_sirs_indep_sec) # bind results together

msmsir_final[,2:4] <-round(msmsir_final[,2:4],2)
msmsir_final[,7] <-round(msmsir_final[,7],0)

######## summarize number of cancer cases
num<-data.frame(tapply(msm$O, msm$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)

num_sec<-data.frame(tapply(msm_sec$O, msm_sec$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:10)

num<-bind_rows(num, num_sec)
msmsf<-merge(msmsir_final, num, by='Name', all.x = T)
msmsf <- arrange(msmsf, desc(SIR))
save(msmsf, file='msmsir.RData')

#################### IVDU ############################################################
options(scipen = 9999999)
transmigroup <- read_excel("SIR/subgroup/Categorical/transmigroup.xlsx")

ivdu<-subset(transmigroup, groupivdu=='Injection drug users')

# Primary robust
ivdu_prim<-subset(ivdu, primNADCstyp !='Brain and CNS' & # rma
                    primNADCstyp !='Breast' &
                    primNADCstyp !='Colon, rectosigmoid junction, rectum, anus and anal canal' &
                    primNADCstyp !='Eye and adnexa'&
                    primNADCstyp !='Leukaemia'&
                    primNADCstyp !='Liver, gallbladder, and biliary tract'&
                    primNADCstyp !='Melanoma of skin' &
                    primNADCstyp !='Mesothelial and soft tissue' &
                    primNADCstyp !='Penis'&
                    primNADCstyp !='Skin nonmelanoma'&
                    primNADCstyp !='Testis'&
                    primNADCstyp !='Prostate'&
                    primNADCstyp !='Vulva and vagina')
ivdu_prim$primNADCstyp<-as.factor(as.character(ivdu_prim$primNADCstyp))
summary(ivdu_prim$primNADCstyp)

num<-length(levels(ivdu_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(ivdu_prim, ivdu_prim$primNADCstyp==levels(ivdu_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(ivdu_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

ivdu_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
ivdu_sirpri$group<-'ivdu'

### robust for secondary
ivdu_sec<-subset(ivdu, !is.na(ivdu$secNADCstyp))

ivdu_sec<-subset(ivdu_sec, ivdu_sec$secNADCstyp!='Colon and Rectum' &
                   secNADCstyp!='Nasopharynx'&
                   secNADCstyp!='Salivary glands' &
                   secNADCstyp!='Testis seminoma')

ivdu_sec$secNADCstyp<-as.factor(as.character(ivdu_sec$secNADCstyp))

summary(ivdu_sec$secNADCstyp)

num_2<-length(levels(ivdu_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(ivdu_sec, secNADCstyp==levels(ivdu_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(ivdu_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

ivdu_sirsec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
ivdu_sirsec$group<-'ivdu'

###rma
ivdu_indep<-subset(ivdu, primNADCstyp=='Prostate'|
                     primNADCstyp=='Brain and CNS')

ivdu_indep$primNADCstyp<-as.factor(as.character(ivdu_indep$ primNADCstyp))
summary(ivdu_indep$primNADCstyp)

num_3<-length(levels(ivdu_indep$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(ivdu_indep, primNADCstyp==levels(ivdu_indep$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(ivdu_indep$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}


ivdu_sirs_indep_prim<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

ivdu_sirs_indep_prim$group<-'ivdu'

ivdusir_final<-bind_rows(ivdu_sirpri, ivdu_sirsec, ivdu_sirs_indep_prim) # bind results together


ivdusir_final[,2:4] <-round(ivdusir_final[,2:4],2)
ivdusir_final[,7] <-round(ivdusir_final[,7],0)

######## summarize number of cancer cases
num<-data.frame(tapply(ivdu$O, ivdu$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)

num_sec<-data.frame(tapply(ivdu$O, ivdu$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:10)

num<-bind_rows(num, num_sec)
ivdusf<-merge(ivdusir_final, num, by='Name', all.x = T)

save(ivdusf, file='ivdusir.RData')

############## heterosexual contact ##############
options(scipen = 9999999)
transmigroup <- read_excel("SIR/subgroup/Categorical/transmigroup.xlsx")

het<-subset(transmigroup, groupivdu=='Heterosexual contact')

# Primary robust
het_prim<-subset(het, primNADCstyp !='Brain and CNS' &
                   primNADCstyp !='Breast' &
                   primNADCstyp !='Colon, rectosigmoid junction, rectum, anus and anal canal' &
                   primNADCstyp !='Eye and adnexa'&
                   primNADCstyp !='Lip, oral cavity and pharynx'& ## rma
                   primNADCstyp !='Liver, gallbladder, and biliary tract'&
                   primNADCstyp !='Melanoma of skin' &
                   primNADCstyp !='Mesothelial and soft tissue' &
                   primNADCstyp !='Prostate'&
                   primNADCstyp !='Penis'&
                   primNADCstyp !='Skin nonmelanoma'&
                   primNADCstyp !='Testis'&
                   primNADCstyp !='Vulva and vagina')
het_prim$primNADCstyp<-as.factor(as.character(het_prim$primNADCstyp))
summary(het_prim$primNADCstyp)

num<-length(levels(het_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(het_prim, het_prim$primNADCstyp==levels(het_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(het_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

het_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
het_sirpri$group<-'Heterosexual contact'

### robust for secondary
het_sec<-subset(het, !is.na(het$secNADCstyp))

het_sec<-subset(het_sec, het_sec$secNADCstyp!='Colon and Rectum' &
                  secNADCstyp!='Liver'&
                  secNADCstyp!='Tonsil' &
                  secNADCstyp!='Testis seminoma')

het_sec$secNADCstyp<-as.factor(as.character(het_sec$secNADCstyp))

summary(het_sec$secNADCstyp)

num_2<-length(levels(het_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(het_sec, secNADCstyp==levels(het_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(het_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

het_sirsec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
het_sirsec$group<-'Heterosexual contact'

###rma
het_indep<-subset(het, primNADCstyp =='Lip, oral cavity and pharynx'|
                    primNADCstyp =='Liver, gallbladder, and biliary tract')

het_indep$primNADCstyp<-as.factor(as.character(het_indep$ primNADCstyp))
summary(het_indep$primNADCstyp)

num_3<-length(levels(het_indep$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(het_indep, primNADCstyp==levels(het_indep$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(het_indep$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}


het_sirs_indep_prim<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

het_sirs_indep_prim$group<-'Heterosexual contact'

hetsir_final<-bind_rows(het_sirpri, het_sirsec, het_sirs_indep_prim) # bind results together


hetsir_final[,2:4] <-round(hetsir_final[,2:4],2)
hetsir_final[,7] <-round(hetsir_final[,7],0)

######## summarize number of cancer cases
num<-data.frame(tapply(het$O, het$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)

num_sec<-data.frame(tapply(het$O, het$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:10)

num<-bind_rows(num, num_sec)
hetsf<-merge(hetsir_final, num, by='Name', all.x = T)

save(hetsf, file='hetsir.RData')

load("~/Desktop/Academia/Projects/NADCs-risk and mortality/NADCs-R project/hetsir.RData")
load("~/Desktop/Academia/Projects/NADCs-risk and mortality/NADCs-R project/ivdusir.RData")
load("~/Desktop/Academia/Projects/NADCs-risk and mortality/NADCs-R project/msmsir.RData")

transmigpsf<-bind_rows(hetsf, ivdusf, msmsf)
transmigpsf<-arrange(transmigpsf, Name)
transmigpsf$'SIR (95%CI)' <-paste0(transmigpsf$SIR," ","(", transmigpsf$LL, "-", transmigpsf$UL, ")")
transmigpsf$lgSIR<-log10(transmigpsf$SIR)
transmigpsf$lgul<-log10(transmigpsf$UL)
transmigpsf$lgll<-log10(transmigpsf$LL)


write_xlsx(transmigpsf, path ="SIR/subgroup/Categorical/sirtransmission_result.xlsx" )


############################################ HIV vs AIDS ############################################
hivaids <- read_excel("SIR/subgroup/Categorical/hivaids.xlsx")
attach(hivaids)
hivaids$study_id<-ifelse(!is.na(Cohortab), Cohortab, Study)
detach(hivaids)
hivaids<- escalc(measure="IRLN", xi=O, ti=E, data=hivaids)

write_xlsx(hivaids, path = 'SIR/subgroup/Categorical/hivaids.xlsx')

####################### AIDS ###########
options(scipen = 9999999)
aid<-subset(hivaids, group2=='AIDS')

# rma
aid_rma<-subset(aid, primNADCstyp== "Anus and anal canal"|
                  primNADCstyp== "Bladder"|
                  primNADCstyp== "Bone and joints"|
                  primNADCstyp=="Brain and CNS"|
                  primNADCstyp=="Breast"|
                  primNADCstyp=="Esophagus"|
                  primNADCstyp=="Hodgkin’s lymphoma"|
                  primNADCstyp=="Kidney and renal pelvis"|
                  primNADCstyp=="Larynx"|
                  primNADCstyp=="Leukaemia"|
                  primNADCstyp=="Liver"|
                  primNADCstyp=="Melanoma of skin"|
                  primNADCstyp=="Multiple myeloma"|
                  primNADCstyp=="Nasal cavity, middle ear, and accessory sinuses"|
                  primNADCstyp=="Pancreas"|
                  primNADCstyp=="Penis"|
                  primNADCstyp=="Prostate"|
                  primNADCstyp=="Skin nonmelanoma"|
                  primNADCstyp=="Small intestine"|
                  primNADCstyp=="Stomach"|
                  primNADCstyp=="Testis"|
                  primNADCstyp=="Thymus, heart, mediastinum, and pleura"|
                  primNADCstyp=="Trachea, bronchus, and lung" )

aid_rma$primNADCstyp<-as.factor(as.character(aid_rma$primNADCstyp))
summary(aid_rma$primNADCstyp)

num_3<-length(levels(aid_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(aid_rma, primNADCstyp==levels(aid_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(aid_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

aid_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

aid_rma_result$group<-'AIDS'

# multilevel
aid_ml<-subset(aid, primNADCstyp=="Colon and Rectum"|
                 primNADCstyp=="Mesothelial and soft tissue")

aid_ml$primNADCstyp<-as.factor(as.character(aid_ml$primNADCstyp))

num<-length(levels(aid_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
PHeterogenity<-0
I2<-0

for (i in 1:num) {
  data = subset(aid_ml, aid_ml$primNADCstyp==levels(aid_ml$primNADCstyp)[i])
  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | Study/origNADCstyp,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(aid_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]

  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))

}


aid_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                          num_level3,I2)

aid_ml_result$group<-'AIDS'

# robust
aid_robu<-subset(aid, primNADCstyp =='Lip, oral cavity and pharynx')

aid_robu$primNADCstyp<-as.factor(as.character(aid_robu$primNADCstyp))
summary(aid_robu$primNADCstyp)

num<-length(levels(aid_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(aid_robu, aid_robu$primNADCstyp==levels(aid_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(aid_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

aid_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
aid_sirpri$group<-'AIDS'

############## HIV not AIDS ##############
hiv<-subset(hivaids, group2=='HIV not AIDS')

# rma
hiv_rma<-subset(hiv, primNADCstyp== "Anus and anal canal"|
                  primNADCstyp== "Bladder"|
                  #primNADCstyp== "Bone and joints"|
                  primNADCstyp=="Brain and CNS"|
                  primNADCstyp=="Breast"|
                  primNADCstyp=="Esophagus"|
                  primNADCstyp=="Hodgkin’s lymphoma"|
                  primNADCstyp=="Kidney and renal pelvis"|
                  primNADCstyp=="Larynx"|
                  primNADCstyp=="Leukaemia"|
                  primNADCstyp=="Liver"|
                  primNADCstyp=="Melanoma of skin"|
                  primNADCstyp=="Multiple myeloma"|
                  primNADCstyp=="Nasal cavity, middle ear, and accessory sinuses"|
                  primNADCstyp=="Pancreas"|
                  primNADCstyp=="Penis"|
                  primNADCstyp=="Prostate"|
                  primNADCstyp=="Skin nonmelanoma"|
                  primNADCstyp=="Small intestine"|
                  # primNADCstyp=="Stomach"|
                  primNADCstyp=="Testis"|
                  #primNADCstyp=="Thymus, heart, mediastinum, and pleura"|
                  primNADCstyp=="Trachea, bronchus, and lung" )

hiv_rma$primNADCstyp<-as.factor(as.character(hiv_rma$primNADCstyp))
summary(hiv_rma$primNADCstyp)

num_3<-length(levels(hiv_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(hiv_rma, primNADCstyp==levels(hiv_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(hiv_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

hiv_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

hiv_rma_result$group<-'HIV not AIDS'

# multilevel
hiv_ml<-subset(hiv, primNADCstyp=="Colon and Rectum"|
                 primNADCstyp=="Mesothelial and soft tissue")

hiv_ml$primNADCstyp<-as.factor(as.character(hiv_ml$primNADCstyp))

num<-length(levels(hiv_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
I2<-0

for (i in 1:num) {
  data<-subset(hiv_ml, hiv_ml$primNADCstyp==levels(hiv_ml$primNADCstyp)[i])
  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | Study/origNADCstyp,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(hiv_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]

  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))

}


hiv_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                          num_level3,PHeterogenity)

hiv_ml_result$group<-'HIV not AIDS'

# robust
hiv_robu<-subset(hiv, primNADCstyp =='Lip, oral cavity and pharynx')

hiv_robu$primNADCstyp<-as.factor(as.character(hiv_robu$primNADCstyp))
summary(hiv_robu$primNADCstyp)

num<-length(levels(hiv_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(hiv_robu, hiv_robu$primNADCstyp==levels(hiv_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(hiv_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

hiv_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
hiv_sirpri$group<-'HIV not AIDS'

hvsir_result<-bind_rows(aid_rma_result,
                        aid_ml_result,aid_sirpri,
                        hiv_rma_result, hiv_ml_result,
                        hiv_sirpri)
hvsir_result<-arrange(hvsir_result, Name)

######## summarize number of cancer cases
num<-data.frame(tapply(hivaids$O, hivaids$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:33)

hvsf<-merge(hvsir_result, num, by='Name', all.x = T)

write_xlsx(hvsf, path ="SIR/subgroup/Categorical/hivaids_result.xlsx")

############################ HAART use ####################################################################################
art <- read_excel("SIR/subgroup/Categorical/haart.xlsx")
attach(art)
art$study_id<-ifelse(!is.na(Cohortab), Cohortab, Study)
detach(art)
art<- escalc(measure="IRLN", xi=O, ti=E, data=art)

write_xlsx(art, path = 'SIR/subgroup/Categorical/haart.xlsx')

art$primNADCstyp[art$primNADCstyp=='Liver, gallbladder, and biliary tract']<-'Liver'
arty<-subset(art, group2=='HAART use')
arty$primNADCstyp<-as.factor(arty$primNADCstyp)
summary(arty$primNADCstyp)

art_meta<-subset(arty, primNADCstyp=='Hodgkin’s lymphoma'|
                   primNADCstyp=='Liver'|
                   primNADCstyp=='Trachea, bronchus, and lung')

art_meta$primNADCstyp<-as.factor(as.character(art_meta$primNADCstyp))
summary(art_meta$primNADCstyp)

num_3<-length(levels(art_meta$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(art_meta, primNADCstyp==levels(art_meta$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(art_meta$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

artuse<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

artuse$group<-'HAART use'

# add cancer cases
num<-data.frame(tapply(art_meta$O, art_meta$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)

artuse<-merge(artuse, num, by='Name', all.x = T)


## no ART use

artn<-subset(art, group2=='No HAART use')
artn$primNADCstyp<-as.factor(artn$primNADCstyp)
summary(artn$primNADCstyp)

nart_meta<-subset(artn, primNADCstyp=='Hodgkin’s lymphoma'|
                    primNADCstyp=='Liver'|
                    primNADCstyp=='Trachea, bronchus, and lung')

nart_meta$primNADCstyp<-as.factor(as.character(nart_meta$primNADCstyp))
summary(nart_meta$primNADCstyp)

num_3<-length(levels(nart_meta$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(nart_meta, primNADCstyp==levels(nart_meta$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(nart_meta$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

nartuse<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

nartuse$group<-'No HAART use'

# add cancer cases
num<-data.frame(tapply(nart_meta$O, nart_meta$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)

nartuse<-merge(nartuse, num, by='Name', all.x = T)

art<-bind_rows(artuse, nartuse)

art<-arrange(art, Name)

write_xlsx(art, path = 'SIR/subgroup/Categorical/HAART use_result.xlsx')

############################# continous variables ####################################################
#################################### HAART year ####################################
# data preparation
harrtyear <- read_excel("SIR/subgroup/Continuous/harrtyear.xlsx")

harrtyear$study_id<-ifelse(!is.na(harrtyear$Cohortab), harrtyear$Cohortab, harrtyear$Study)
harrtyear$SIR_extract<-as.numeric(harrtyear$SIR_extract)

harrtyear$logsir<-log(harrtyear$SIR_extract)
harrtyear$logul<-log(harrtyear$ul_extract)


harrtyear$O<-ifelse(is.na(harrtyear$O),
                    round((1.96/(harrtyear$logsir-harrtyear$logul))^2,digits = 0), harrtyear$O)

harrtyear$E<-ifelse(is.na(harrtyear$E),
                    harrtyear$O/harrtyear$SIR_extract, harrtyear$E)


harrtyear<- escalc(measure="IRLN", xi=O, ti=E, data=harrtyear)

write_xlsx(harrtyear, path = 'SIR/subgroup/Continuous/harrtyear.xlsx')

## overall
harrtyear <- read_excel("SIR/subgroup/Continuous/harrtyear.xlsx")
dput(names(harrtyear))
cols<-c("study_id", "group_yearfinal",
        "primNADCstyp")

harrtyear[cols]<-lapply(harrtyear[cols], factor)

################## 1980-1989 ################
cal1<-subset(harrtyear, group_yearfinal=='1980-1989')
cal1_rma<-subset(cal1, primNADCstyp== "Anus and anal canal"|
                   primNADCstyp=="Brain and CNS"|
                   primNADCstyp=="Colon and Rectum"|
                   primNADCstyp=="Hodgkin lymphoma"|
                   primNADCstyp=="Kidney and renal pelvis"|
                   primNADCstyp=="Larynx"|
                   primNADCstyp=="Melanoma of skin"|
                   primNADCstyp=="Oral cavity and pharynx"|
                   primNADCstyp=="Prostate"|
                   primNADCstyp=="Stomach"|
                   primNADCstyp=="Testis"|
                   primNADCstyp=="Trachea, bronchus, and lung")

cal1_rma$primNADCstyp<-as.factor(as.character(cal1_rma$primNADCstyp))
summary(cal1_rma$primNADCstyp)

num_3<-length(levels(cal1_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(cal1_rma, primNADCstyp==levels(cal1_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(cal1_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

cal1_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

cal1_rma_result$group<-'1980-1989'


# robust: "Leukaemia",

cal1_robu<-subset(cal1, primNADCstyp== 'Leukaemia' )
cal1_robu$primNADCstyp<-as.factor(as.character(cal1_robu$primNADCstyp))
summary(cal1_robu$primNADCstyp)

num<-length(levels(cal1_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(cal1_robu, cal1_robu$primNADCstyp==levels(cal1_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(cal1_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

cal1_roburesult<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
cal1_roburesult$group<-'1980-1989'



# single study: "Bladder", "Bones and joints", "Breast", "Esophagus", "Eye and adnexa", "Head and neck", "Lip", "Liver", "Mesothelial and soft tissue",
"Multiple myeloma", "Ovary", "Pancreas", "Thyroid",

########################################### 1990-1995 ######################################################
cal2<-subset(harrtyear, group_yearfinal=='1990-1995')
# rma: ,

cal2_rma<-subset(cal2, primNADCstyp== "Bones and joints"|
                   primNADCstyp=="Esophagus"|
                   primNADCstyp=="Head and neck"|
                   primNADCstyp=="Larynx"|
                   primNADCstyp=="Melanoma of skin"|
                   primNADCstyp=="Mesothelial and soft tissue"|
                   primNADCstyp=="Ovary"|
                   primNADCstyp=="Pancreas"|
                   primNADCstyp=="Penis"|
                   primNADCstyp=="Skin nonmelanoma"|
                   primNADCstyp=="Stomach"|
                   primNADCstyp=="Testis"|
                   primNADCstyp=="Thyroid"|
                   primNADCstyp=="Uterus"|
                   primNADCstyp=="Vulva and vagina")

cal2_rma$primNADCstyp<-as.factor(as.character(cal2_rma$primNADCstyp))
summary(cal2_rma$primNADCstyp)

num_3<-length(levels(cal2_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(cal2_rma, primNADCstyp==levels(cal2_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(cal2_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

cal2_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

cal2_rma_result$group<-'1990-1995'


# robust:

cal2_robu<-subset(cal2, primNADCstyp== "Anus and anal canal"|
                    primNADCstyp=="Brain and CNS"|
                    primNADCstyp=="Breast"|
                    primNADCstyp=="Colon and Rectum"|
                    primNADCstyp=="Hodgkin lymphoma"|
                    primNADCstyp=="Kidney and renal pelvis"|
                    primNADCstyp=="Leukaemia"|
                    primNADCstyp=="Liver"|
                    primNADCstyp=="Multiple myeloma"|
                    primNADCstyp=="Oral cavity and pharynx"|
                    primNADCstyp=="Prostate"|
                    primNADCstyp=="Trachea, bronchus, and lung")
cal2_robu$primNADCstyp<-as.factor(as.character(cal2_robu$primNADCstyp))
summary(cal2_robu$primNADCstyp)

num<-length(levels(cal2_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(cal2_robu, cal2_robu$primNADCstyp==levels(cal2_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(cal2_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

cal2_roburesult<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
cal2_roburesult$group<-'1990-1995'

# single: "Eye and adnexa",

################ 1996-2001 ################################
cal3<-subset(harrtyear, group_yearfinal=='1996-2001')
# rma: ,

cal3_rma<-subset(cal3, primNADCstyp== "Bladder"|
                   primNADCstyp=="Bones and joints"|
                   primNADCstyp=="Esophagus"|
                   primNADCstyp== "Head and neck"|
                   primNADCstyp==  "Larynx"|
                   primNADCstyp== "Leukaemia"|
                   primNADCstyp== "Lip"|
                   primNADCstyp== "Melanoma of skin"|
                   primNADCstyp=="Multiple myeloma"|
                   primNADCstyp=="Ovary"|
                   primNADCstyp=="Penis"|
                   primNADCstyp=="Skin nonmelanoma"|
                   primNADCstyp=="Thyroid"|
                   primNADCstyp=="Uterus")

cal3_rma$primNADCstyp<-as.factor(as.character(cal3_rma$primNADCstyp))
summary(cal3_rma$primNADCstyp)

num_3<-length(levels(cal3_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(cal3_rma, primNADCstyp==levels(cal3_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(cal3_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

cal3_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

cal3_rma_result$group<-'1996-2001'


############# robust:
cal3_robu<-subset(cal3, primNADCstyp== "Anus and anal canal"|
                    primNADCstyp== "Brain and CNS"|
                    primNADCstyp== "Breast"|
                    primNADCstyp== "Colon and Rectum"|
                    primNADCstyp== "Hodgkin lymphoma"|
                    primNADCstyp== "Kidney and renal pelvis"|
                    primNADCstyp== "Liver"|
                    primNADCstyp== "Mesothelial and soft tissue"|
                    primNADCstyp=="Oral cavity and pharynx"|
                    primNADCstyp=="Prostate"|
                    primNADCstyp=="Stomach"|
                    primNADCstyp=="Testis"|
                    primNADCstyp=="Trachea, bronchus, and lung"|
                    primNADCstyp=="Vulva and vagina"|
                    primNADCstyp=="Pancreas")
cal3_robu$primNADCstyp<-as.factor(as.character(cal3_robu$primNADCstyp))
summary(cal3_robu$primNADCstyp)

num<-length(levels(cal3_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(cal3_robu, cal3_robu$primNADCstyp==levels(cal3_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(cal3_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

cal3_roburesult<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
cal3_roburesult$group<-'1996-2001'

# single: "Eye and adnexa",

######################## 2002-2013 ####################

cal4<-subset(harrtyear, group_yearfinal=='2002-2013')

cal4_rma<-subset(cal4, primNADCstyp== "Bladder"|
                   primNADCstyp== "Brain and CNS"|
                   primNADCstyp== "Breast"|
                   primNADCstyp== "Colon and Rectum"|
                   primNADCstyp== "Head and neck"|
                   primNADCstyp== "Kidney and renal pelvis"|
                   primNADCstyp== "Leukaemia"|
                   primNADCstyp== "Melanoma of skin"|
                   #primNADCstyp== "Multiple myeloma"|
                   primNADCstyp== "Oral cavity and pharynx"|
                   primNADCstyp== "Testis")

cal4_rma$primNADCstyp<-as.factor(as.character(cal4_rma$primNADCstyp))
summary(cal4_rma$primNADCstyp)

num_3<-length(levels(cal4_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(cal4_rma, primNADCstyp==levels(cal4_rma$primNADCstyp)[i]), method="ML")

  Name[i]<-levels(cal4_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

cal4_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

cal4_rma_result$group<-'2002-2013'

# unique
cal4_u<-subset(cal4, primNADCstyp== "Multiple myeloma")

cal4_u$primNADCstyp<-as.factor(as.character(cal4_u$primNADCstyp))
summary(cal4_u$primNADCstyp)

fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                    data =cal4_u, method="FE")

Name<-levels(cal4_u$primNADCstyp)
SIR<-exp(fit_norob$beta)
LL<-exp(fit_norob$ci.lb)
UL<-exp(fit_norob$ci.ub)

num_study<-fit_norob$k
num_outcome<-fit_norob$k

I2<-fit_norob$I2
Tau<-fit_norob$tau2


cal4_u_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

cal4_u_result$group<-'2002-2013'


# robust:

cal4_robu<-subset(cal4, primNADCstyp== "Anus and anal canal"|
                    primNADCstyp=="Hodgkin lymphoma"|
                    primNADCstyp=="Liver"|
                    primNADCstyp=="Prostate"|
                    primNADCstyp=="Trachea, bronchus, and lung")

cal4_robu$primNADCstyp<-as.factor(as.character(cal4_robu$primNADCstyp))
summary(cal4_robu$primNADCstyp)

num<-length(levels(cal4_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(cal4_robu, cal4_robu$primNADCstyp==levels(cal4_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(cal4_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

cal4_roburesult<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
cal4_roburesult$group<-'2002-2013'

cal_final<-bind_rows(cal1_rma_result, cal1_roburesult,
                     cal2_rma_result, cal2_roburesult,
                     cal3_rma_result, cal3_roburesult,
                     cal4_rma_result, cal4_u_result,
                     cal4_roburesult)

cal_final<-arrange(cal_final, Name)

# add cancer case

cal_final<-within(cal_final, {
  lgsir<-log10(SIR)
  lgll<-log10(LL)
  lgul<-log10(UL)
  Name<-factor(Name)
  group<-factor(group)
  final_SIR <-paste0(round(SIR, 2)," ","(", round(LL, 2), "-", round(UL,2), ")")
  case_group<-paste0(Name, group)

}
)

harrtyear$case_group<-paste0(harrtyear$primNADCstyp, harrtyear$group_yearfinal)
num<-data.frame(tapply(harrtyear$O, harrtyear$case_group, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
colnames(num)[2]<-'case_group'

harrty<-merge(cal_final, num, by='case_group', all.x = T)

# single: "Eye and adnexa", "LIP" "Ovary", "Pancreas","Skin nonmelanoma",  "Stomach", "Vulva and vagina"

write_xlsx(harrty, path = 'SIR/subgroup/Continuous/harrtyear_result.xlsx')




aidsperiod <- read_excel("SIR/subgroup/Continuous/aidsperiod.xlsx")
aidsperiod$study_id<-ifelse(!is.na(aidsperiod$Cohortab), aidsperiod$Cohortab, aidsperiod$Study)
aidsperiod<- escalc(measure="IRLN", xi=O, ti=E, data=aidsperiod)
write_xlsx(aidsperiod, path = 'SIR/subgroup/Continuous/aidsperiod.xlsx')

a<-5
a*12


##################AIDS relative was defined by authors#####################
aidsperiod <- read_excel("SIR/subgroup/Continuous/aidsperiod.xlsx")

cols<-c("study_id",  "group_final", "origNADCstyp",
        "primNADCstyp")

aidsperiod[cols]<-lapply(aidsperiod[cols], factor)

# Period1
# rma
period1<-subset(aidsperiod, group_final=='period1')
period1_rma<-subset(period1, primNADCstyp== "Anus and anal canal"|
                      primNADCstyp=="Brain and CNS"|
                      primNADCstyp=="Colon and Rectum"|
                      primNADCstyp=="Hodgkin’s lymphoma"|
                      primNADCstyp=="Lip, oral cavity and pharynx"|
                      primNADCstyp=="Liver"|
                      primNADCstyp=="Melanoma of skin"|
                      primNADCstyp=="Mesothelial and soft tissue"|
                      primNADCstyp=="Prostate"|
                      primNADCstyp=="Testis")

period1_rma$primNADCstyp<-as.factor(as.character(period1_rma$primNADCstyp))
summary(period1_rma$primNADCstyp)

num_3<-length(levels(period1_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(period1_rma, primNADCstyp==levels(period1_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(period1_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

period1_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

period1_rma_result$group<-'period1'


# multilevel

period1_ml<-subset(period1, primNADCstyp =="Leukaemia"|  # 4
                     primNADCstyp =="Trachea, bronchus, and lung")

period1_ml$primNADCstyp<-as.factor(as.character(period1_ml$primNADCstyp))

num<-length(levels(period1_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
I2<-0

for (i in 1:num) {
  data = subset(period1_ml, period1_ml$primNADCstyp==levels(period1_ml$primNADCstyp)[i])
  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | Study/origNADCstyp,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(period1_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]
  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))



}


period1_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                              num_level3,I2)

period1_ml_result$group<-'period1'

# period 2
# rma
period2<-subset(aidsperiod, group_final=='period2')
period2_rma<-subset(period2, primNADCstyp== "Anus and anal canal"|
                      primNADCstyp=="Brain and CNS"|
                      primNADCstyp=="Hodgkin’s lymphoma"|
                      primNADCstyp=="Lip, oral cavity and pharynx"|
                      primNADCstyp=="Liver"|
                      primNADCstyp=="Melanoma of skin"|
                      primNADCstyp=="Multiple myeloma"|
                      primNADCstyp=="Testis")

period2_rma$primNADCstyp<-as.factor(as.character(period2_rma$primNADCstyp))
summary(period2_rma$primNADCstyp)

num_3<-length(levels(period2_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(period2_rma, primNADCstyp==levels(period2_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(period2_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

period2_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

period2_rma_result$group<-'period2'

# FE
period2_fe<-subset(period2, primNADCstyp=="Colon and Rectum"|
                     primNADCstyp=="Prostate")

period2_fe$primNADCstyp<-as.factor(as.character(period2_fe$primNADCstyp))
summary(period2_fe$primNADCstyp)

num_3<-length(levels(period2_fe$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(period2_fe, primNADCstyp==levels(period2_fe$primNADCstyp)[i]),method="FE")

  Name[i]<-levels(period2_fe$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

period2_fe_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

period2_fe_result$group<-'period2'

# multilevel

period2_ml<-subset(period2, primNADCstyp =="Leukaemia"|  # 4
                     primNADCstyp =="Trachea, bronchus, and lung")

period2_ml$primNADCstyp<-as.factor(as.character(period2_ml$primNADCstyp))

num<-length(levels(period2_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
I2<-0

for (i in 1:num) {
  data = subset(period2_ml, period2_ml$primNADCstyp==levels(period2_ml$primNADCstyp)[i])

  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | Study/origNADCstyp,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(period2_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]

  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))
}


period2_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                              num_level3,PHeterogenity)


period2_ml_result$group<-'period2'


# period 3
# rma
period3<-subset(aidsperiod, group_final=='period3')
period3_rma<-subset(period3, primNADCstyp== "Anus and anal canal"|
                      primNADCstyp=="Brain and CNS"|
                      primNADCstyp=="Colon and Rectum"|
                      primNADCstyp=="Hodgkin’s lymphoma"|
                      primNADCstyp=="Lip, oral cavity and pharynx"|
                      primNADCstyp=="Liver"|
                      primNADCstyp=="Melanoma of skin"|
                      primNADCstyp=="Mesothelial and soft tissue"|
                      primNADCstyp=="Multiple myeloma"|
                      primNADCstyp=="Prostate"|
                      primNADCstyp=="Testis")

period3_rma$primNADCstyp<-as.factor(as.character(period3_rma$primNADCstyp))
summary(period3_rma$primNADCstyp)

num_3<-length(levels(period3_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(period3_rma, primNADCstyp==levels(period3_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(period3_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

period3_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

period3_rma_result$group<-'period3'

# multilevel

period3_ml<-subset(period3, primNADCstyp =="Leukaemia"|  # 4
                     primNADCstyp =="Trachea, bronchus, and lung")

period3_ml$primNADCstyp<-as.factor(as.character(period3_ml$primNADCstyp))

num<-length(levels(period3_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
PHeterogenity<-0
I2<-0

for (i in 1:num) {
  data = subset(period3_ml, period3_ml$primNADCstyp==levels(period3_ml$primNADCstyp)[i])

  fourlevel<-rma.mv(yi,
                    vi,
                    data = data,
                    random = ~ 1 | Study/origNADCstyp,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(period3_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]
  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))
}


period3_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                              num_level3,PHeterogenity)


period3_ml_result$group<-'period3'

# period 4
# rma
period4<-subset(aidsperiod, group_final=='period4')
period4_rma<-subset(period4, primNADCstyp== "Anus and anal canal"|
                      primNADCstyp=="Colon and Rectum"|
                      primNADCstyp=="Hodgkin’s lymphoma"|
                      primNADCstyp=="Lip, oral cavity and pharynx"|
                      primNADCstyp=="Melanoma of skin"|
                      primNADCstyp=="Mesothelial and soft tissue"|
                      primNADCstyp=="Multiple myeloma"|
                      primNADCstyp=="Prostate"|
                      primNADCstyp=="Testis")

period4_rma$primNADCstyp<-as.factor(as.character(period4_rma$primNADCstyp))
summary(period4_rma$primNADCstyp)

num_3<-length(levels(period4_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(period4_rma, primNADCstyp==levels(period4_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(period4_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

period4_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

period4_rma_result$group<-'period4'

# multilevel

period4_ml<-subset(period4, primNADCstyp =="Leukaemia"|  # 4
                     primNADCstyp =="Trachea, bronchus, and lung")

period4_ml$primNADCstyp<-as.factor(as.character(period4_ml$primNADCstyp))

num<-length(levels(period4_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
PHeterogenity<-0
I2<-0

for (i in 1:num) {
  data = subset(period4_ml, period4_ml$primNADCstyp==levels(period4_ml$primNADCstyp)[i])

  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | Study/origNADCstyp,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(period4_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]
  PHeterogenity[i]<-fourlevel$QEp
  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))
}


period4_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                              num_level3,PHeterogenity)

period4_ml_result$group<-'period4'

# period 5
origNADCsty
# rma
period5<-subset(aidsperiod, group_final=='period5')

# multilevel

period5_ml<-subset(period5, primNADCstyp =="Anus and anal canal")

period5_ml$primNADCstyp<-as.factor(as.character(period5_ml$primNADCstyp))

num<-length(levels(period5_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
I2<-0

for (i in 1:num) {
  data = subset(period5_ml, period5_ml$primNADCstyp==levels(period5_ml$primNADCstyp)[i])
  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | Study/origNADCstyp,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(period5_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]
  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))
}


period5_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                              num_level3,PHeterogenity)

period5_ml_result$group<-'period5'

aidsperiod_result<-bind_rows(period1_rma_result,
                             period1_ml_result,
                             period2_rma_result,
                             period2_fe_result,
                             period2_ml_result,
                             period3_rma_result,
                             period3_ml_result,
                             period4_rma_result,
                             period4_ml_result,
                             period5_ml_result)


aidsperiod_result<-arrange(aidsperiod_result, Name)
write_xlsx(aidsperiod_result, path='SIR/subgroup/Continuous/aidsperiod_result.xlsx')

######## age ########
age <- read_excel("SIR/subgroup/Continuous/age.xlsx")
age$Adjustment<-NULL
age$study_id<-ifelse(!is.na(age$Cohortab), age$Cohortab, age$Study)
age$logsir<-log(age$SIR_extract)
age$logul<-log(age$ul_extract)

age$O<-ifelse(is.na(age$O),
              round((1.96/(age$logsir-age$logul))^2,digits = 0), age$O)

age$E<-ifelse(is.na(age$E),
              age$O/age$SIR_extract, age$E)


age<- escalc(measure="IRLN", xi=O, ti=E, data=age)

write_xlsx(age, path = 'SIR/subgroup/Continuous/age.xlsx')
age$agefinal<-ifelse(age$agemid<35, '<35',
                     ifelse(age$agemid>=50, '>=50', '35-49'))

summary(as.factor(age$ageorig))

## <35, 36-49; >=50
age <- read_excel("SIR/subgroup/Continuous/age.xlsx")
cols<-c("study_id",  "grouppop1",
        "grouppop2",  "agefinal",
        "primNADCstyp")
age[cols]<-lapply(age[cols], factor)

## overall
age_all<-subset(age, grouppop1=='Overall')

# <35
age_all1<-subset(age_all, agefinal=='<35')

# robust
age_all1_robu<-subset(age_all1, primNADCstyp== "Anus and anal canal"|
                        primNADCstyp== "Trachea, bronchus, and lung")

age_all1_robu$primNADCstyp<-as.factor(as.character(age_all1_robu$primNADCstyp))
summary(age_all1_robu$primNADCstyp)

num<-length(levels(age_all1_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(age_all1_robu, age_all1_robu$primNADCstyp==levels(age_all1_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(age_all1_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

age_all1_roburesult<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
age_all1_roburesult$group<-'<35'

# rma:
age_all1_rma<-subset(age_all1, primNADCstyp==  "Hodgkin’s lymphoma"|
                       primNADCstyp=="Liver")

age_all1_rma$primNADCstyp<-as.factor(as.character(age_all1_rma$primNADCstyp))
summary(age_all1_rma$primNADCstyp)

num_3<-length(levels(age_all1_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(age_all1_rma, primNADCstyp==levels(age_all1_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(age_all1_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

age_all1_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

age_all1_rma_result$group<-'<35'


# 36-49
age_all2<-subset(age_all, agefinal=='35-49')

# robu
age_all2_robu<-subset(age_all2, primNADCstyp== "Anus and anal canal"|
                        primNADCstyp=="Colon and Rectum")
age_all2_robu$primNADCstyp<-as.factor(as.character(age_all2_robu$primNADCstyp))
summary(age_all2_robu$primNADCstyp)

num<-length(levels(age_all2_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(age_all2_robu, age_all2_robu$primNADCstyp==levels(age_all2_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(age_all2_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

age_all2_roburesult<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
age_all2_roburesult$group<-'35-49'

# rma:
age_all2_rma<-subset(age_all2, primNADCstyp== "Hodgkin’s lymphoma"|
                       primNADCstyp=="Prostate"|
                       primNADCstyp=="Trachea, bronchus, and lung")

age_all2_rma$primNADCstyp<-as.factor(as.character(age_all2_rma$primNADCstyp))
summary(age_all2_rma$primNADCstyp)

num_3<-length(levels(age_all2_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(age_all2_rma, primNADCstyp==levels(age_all2_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(age_all2_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

age_all2_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

age_all2_rma_result$group<-'35-49'


## >=50

age_all3<-subset(age_all, agefinal=='>=50')

# robu
age_all3_robu<-subset(age_all3, primNADCstyp== "Anus and anal canal"|
                        primNADCstyp== "Colon and Rectum"|
                        primNADCstyp== "Hodgkin’s lymphoma"|
                        primNADCstyp== "Prostate"|
                        primNADCstyp== "Trachea, bronchus, and lung")

age_all3_robu$primNADCstyp<-as.factor(as.character(age_all3_robu$primNADCstyp))
summary(age_all3_robu$primNADCstyp)

num<-length(levels(age_all3_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(age_all3_robu, age_all3_robu$primNADCstyp==levels(age_all3_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(age_all3_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

age_all3_roburesult<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
age_all3_roburesult$group<-'>=50'

# rma:
age_all3_rma<-subset(age_all3, primNADCstyp== "Liver")

age_all3_rma$primNADCstyp<-as.factor(as.character(age_all3_rma$primNADCstyp))
summary(age_all3_rma$primNADCstyp)

num_3<-length(levels(age_all3_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(age_all3_rma, primNADCstyp==levels(age_all3_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(age_all3_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

age_all3_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

age_all3_rma_result$group<-'>=50'

age_result<-bind_rows(age_all1_roburesult,
                      age_all1_rma_result,
                      age_all2_roburesult,
                      age_all2_rma_result,
                      age_all3_roburesult,
                      age_all3_rma_result)


age_result<-arrange(age_result, Name)
write_xlsx(age_result, path='SIR/subgroup/Continuous/age_result.xlsx')

####################### all forestplot code ###############################################################
################################# SIR overall ##################################################################

install.packages(Cairo)
Library (Cario)

sfp <- read_excel("SIR/Overall/SIRallforcl.xlsx")

sfp<-within.data.frame(sfp,{ # 小数点统一

  `Type of Cancer`[11:16]<-paste0('  ', `Type of Cancer`[11:16])
  SIR<-as.character(SIR+0.001)
  UL<-as.character(UL+0.001)
  LL<-as.character(LL+0.001)

  SIR<-str_sub(SIR,1,nchar(SIR)-1)
  UL<-str_sub(UL,1,nchar(UL)-1)
  LL<-str_sub(LL,1,nchar(LL)-1)

  `SIR (95%CI)`<-ifelse(!is.na(SIR), paste0(SIR," ","(", LL, "-", UL, ")"),
                        NA)

})

write_xlsx(sfp, path="SIR/Overall/SIRallforcl.xlsx")

sfp_data <-
  subset(sfp, select = c('lgSIR', 'lgll', 'lgul'))
colnames(sfp_data)<-c("mean", "lower", "upper")

sfp_data<-sfp_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

dput(names(sfp))
sir_text<-subset(sfp, select = c("Type of Cancer", "SIR (95%CI)", "num_study", "num_outcome",
                                 "Observed number of cancer", "I2"))

sir_table <- cbind(
  c("Type of Cancer",  sir_text$`Type of Cancer`),
  c("SIR (95%CI)", sir_text$`SIR (95%CI)`),
  c("No. Studies", sir_text$num_study),
  c("No. Outcomes", sir_text$num_outcome),
  c("Cancer cases", sir_text$`Observed number of cancer`),
  c("Heterogenity", sir_text$I2))


f<-forestplot(sir_table,
              sfp_data, new_page = TRUE,
              clip= c(0.1, 28),
              lwd.ci=1.2,  lwd.xaxis=1,  lwd.zero=1,
              zero=0,
              graph.pos = 3,
              xlab = "SIR",
              boxsize = .4,
              hrzl_lines = list("2" = gpar(lty = 1)),
              line.margin = .5,
              align=c("l","l","l","l","l"),
              is.summary = c(TRUE, TRUE,
                             rep(FALSE,20),
                             TRUE,
                             rep(FALSE,20)),
              vertices = TRUE,
              ci.vertices.height= 0.2,
              xticks=c(-1, 0, 1, 2, 3),
              xticks.digits=0,
              graphwidth= unit(9, "cm"),
              txt_gp=fpTxtGp(label=gpar(fontsize=15, cex =0.8, fontfamily="Arial"),
                             ticks=gpar(fontsize=15, cex=0.8, fontfamily="Arial"),
                             xlab=gpar(fontsize=15,cex = 0.8, fontfamily="Arial")),
              col = fpColors(
                box = "#08519C",
                line = "#08519C",
                zero = 'black'))

## export as pdf  11.3 * 16

################################# SMR overall ##################################################################
smr <- read_excel("SMR_plot.xlsx")

smr<-within.data.frame(smr,{

  I2<-round(I2, digits = 0)
  SMR<-as.character(round(SMR, digits = 2)+0.001)
  UL<-as.character(round(UL, digits = 2) +0.001)
  LL<-as.character(round(LL, digits = 1) +0.001)

  SMR<-str_sub(SMR,1,nchar(SMR)-1)
  UL<-str_sub(UL,1,nchar(UL)-1)
  LL<-str_sub(LL,1,nchar(LL)-1)

  final_SMR<-ifelse(!is.na(SMR), paste0(SMR," ","(", LL, "-", UL, ")"),
                    NA)

})


smr_data <-
  subset(smr, select = c('lgsmr', 'lgll', 'lgul'))
colnames(smr_data)<-c("mean", "lower", "upper")

smr_data<-smr_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

dput(names(smr))
smr_text<-subset(smr, select = c("Name", "final_SMR", "num_study", "num_outcome",
                                 "Observed number of cancer", "I2"))

smr_table <- cbind(
  c("Type of Cancer",  smr_text$Name),
  c("SMR (95%CI)", smr_text$final_SMR),
  c("No. Studies", smr_text$num_study),
  c("No. Outcomes", smr_text$num_outcome),
  c("Cancer cases", smr_text$`Observed number of cancer`),
  c("Heterogenity", smr_text$I2))


forestplot(smr_table,
           smr_data, new_page = TRUE,
           clip= c(0.1, 28),
           lwd.ci=1.2,  lwd.xaxis=1,  lwd.zero=1,
           zero=0,
           graph.pos = 3,
           xlab = "SMR",
           boxsize = .25,
           hrzl_lines = list("2" = gpar(lty = 1)),
           line.margin = .5,
           align=c("l","l","l","l","l"),
           is.summary = c(TRUE, TRUE,
                          rep(FALSE,4),
                          TRUE,
                          rep(FALSE,5)),
           vertices = TRUE,
           ci.vertices.height= 0.2,
           xticks=c(-1, 0, 1, 2, 3),
           xticks.digits=0,
           graphwidth= unit(9, "cm"),
           txt_gp=fpTxtGp(label=gpar(fontsize=15, cex =0.8, fontfamily="Arial"),
                          ticks=gpar(fontsize=15, cex=0.8,fontfamily="Arial"),
                          xlab=gpar(fontsize=15,cex = 0.8,fontfamily="Arial")),
           col = fpColors(
             box = "#08519C",
             line = "#08519C",
             zero = 'black'))
## export as pdf  7 * 13

########################################## SIR_gender ###############################################################
sex <- read_excel("genderplot.xlsx")
sex<-within(sex, {
  SIR<-as.character(SIR+0.001)
  UL<-as.character(UL+0.001)
  LL<-as.character(LL+0.001)

  SIR<-str_sub(SIR,1,nchar(SIR)-1)
  UL<-str_sub(UL,1,nchar(UL)-1)
  LL<-str_sub(LL,1,nchar(LL)-1)

  final_SIR<-ifelse(!is.na(Name), paste0(SIR," ","(", LL, "-", UL, ")"),
                    final_SIR)
})


write_xlsx(sex, path = 'genderplot.xlsx')


# plot
# text
sex_table <- cbind(
  c("Cancer and SIR (95%CI)",  sex$final_SIR),
  c("Cancer cases", sex$`Observed number of cancer`),
  c("No. Studies", sex$num_study),
  c("No. Outcomes", sex$num_outcome),
  c("Heterogenity", round(sex$I2, 0)))

# data
sex_data <-
  subset(sex, select = c('lgsir', 'lgll', 'lgul'))
colnames(sex_data)<-c("mean", "lower", "upper")

sex_data<-sex_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

styles <- fpShapesGp(
  lines = list(
    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),


    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#A50026")),

  box = list(
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#A50026", fill = "#A50026"))
)


forestplot(sex_table,
           sex_data, new_psex = TRUE,
           clip= c(0.1, 2.5),
           lwd.ci=1.2,  lwd.xaxis=1,  lwd.zero=1,
           zero=0,
           graph.pos = 3,
           xlab = "SIR",
           boxsize = .5,
           hrzl_lines = list("2" = gpar(lty = 1, lwd=1)),
           line.margin = .5,
           align=c("l","l","l","l","l"),
           is.summary = c(rep(TRUE,3),
                          rep(FALSE,2), TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2), TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2), TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2), TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2), TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2), TRUE,
                          rep(FALSE,2),TRUE, TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2)),
           vertices = TRUE,
           ci.vertices.height= 0.2,
           xticks=c(-1, 0, 1, 2,3),
           xticks.digits=0,
           graphwidth= unit(9, "cm"),
           txt_gp=fpTxtGp(label=gpar(fontsize=15, cex =0.8, fontfamily="Arial"),
                          ticks=gpar(fontsize=15, cex=0.8, fontfamily="Arial"),
                          xlab=gpar(fontsize=15,cex = 0.8, fontfamily="Arial")),
           col = fpColors(zero = "black"),
           shapes_gp = styles)


###################################### 从以下开始修改代码和重新跑 ######################################
###################################### HIV transmission group ######################################
transmigroup <- read_excel("SIR/subgroup/Categorical/transmigroup.xlsx")
attach(transmigroup)
transmigroup$study_id<-ifelse(!is.na(Cohortab), Cohortab, Study)
detach(transmigroup)
transmigroup<- escalc(measure="IRLN", xi=O, ti=E, data=transmigroup)

write_xlsx(transmigroup, path = 'SIR/subgroup/Categorical/transmigroup.xlsx')

####################### MSM ###########
msm<-subset(transmigroup, groupmsm=='MSM')

# Primary robust
msm_prim<-subset(msm,
                 primNADCstyp !='Bladder' &
                   primNADCstyp !='Bone and joints' &
                   primNADCstyp !='Kidney and other urinary' &
                   primNADCstyp !='Leukaemia' &
                   primNADCstyp !='Lip, oral cavity and pharynx'& ## rma
                   primNADCstyp !='Liver, gallbladder, and biliary tract'& ## 直接liver
                   primNADCstyp !='Melanoma of skin'& # rma
                   primNADCstyp !='Mesothelial and soft tissue' &
                   primNADCstyp !='Multiple myeloma' &
                   primNADCstyp !='Pancreas'&
                   primNADCstyp !='Penis'&
                   primNADCstyp !='Eye and adnexa'&
                   primNADCstyp !='Brain and CNS'&
                   primNADCstyp !='Colon, rectosigmoid junction, rectum, anus and anal canal'&
                   primNADCstyp !='Skin nonmelanoma'& ## rma
                   primNADCstyp !='Thyroid')
msm_prim$primNADCstyp<-as.factor(as.character(msm_prim$primNADCstyp))
summary(msm_prim$primNADCstyp)

num<-length(levels(msm_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(msm_prim, msm_prim$primNADCstyp==levels(msm_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(msm_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

msm_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
msm_sirpri$group<-'MSM'

### robust for secondary
msm_sec<-subset(msm, !is.na(msm$secNADCstyp))

msm_sec<-subset(msm_sec,
                msm_sec$secNADCstyp!='Kidney and renal pelvis' &
                  msm_sec$secNADCstyp!='Liver' & # rma
                  msm_sec$secNADCstyp!='Salivary glands' &
                  msm_sec$secNADCstyp!='Testis seminoma')

msm_sec$secNADCstyp<-as.factor(as.character(msm_sec$secNADCstyp))

summary(msm_sec$secNADCstyp)

num_2<-length(levels(msm_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(msm_sec, secNADCstyp==levels(msm_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(msm_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

msm_sirsec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
msm_sirsec$group<-'MSM'

###rma
msm_indep<-subset(msm, primNADCstyp =='Lip, oral cavity and pharynx' |
                    primNADCstyp =='Melanoma of skin' |
                    primNADCstyp =='Skin nonmelanoma'|
                    primNADCstyp =='Brain and CNS' )


### primary
msm_indep$primNADCstyp<-as.factor(as.character(msm_indep$ primNADCstyp))
summary(msm_indep$primNADCstyp)

num_3<-length(levels(msm_indep$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(msm_indep, primNADCstyp==levels(msm_indep$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(msm_indep$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

msm_sirs_indep_prim<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

msm_sirs_indep_prim$group<-'MSM'


### secondary

msm_indepsec<-subset(msm, secNADCstyp=='Liver')

msm_indepsec$secNADCstyp<-as.factor(as.character(msm_indepsec$secNADCstyp))
summary(msm_indepsec$secNADCstyp)

num_3<-length(levels(msm_indepsec$secNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(msm_indepsec, secNADCstyp==levels(msm_indepsec$secNADCstyp)[i]),method="ML")

  Name[i]<-levels(msm_indepsec$secNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

msm_sirs_indep_sec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

msm_sirs_indep_sec$group<-'MSM'

msmsir_final<-bind_rows(msm_sirpri, msm_sirsec, msm_sirs_indep_prim, msm_sirs_indep_sec) # bind results together

msmsir_final[,2:4] <-round(msmsir_final[,2:4],2)
msmsir_final[,7] <-round(msmsir_final[,7],0)

######## summarize number of cancer cases
num<-data.frame(tapply(msm$O, msm$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)

num_sec<-data.frame(tapply(msm_sec$O, msm_sec$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:10)

num<-bind_rows(num, num_sec)
msmsf<-merge(msmsir_final, num, by='Name', all.x = T)
msmsf <- arrange(msmsf, desc(SIR))
save(msmsf, file='msmsir.RData')

#################### IVDU ############################################################
options(scipen = 9999999)
transmigroup <- read_excel("SIR/subgroup/Categorical/transmigroup.xlsx")

ivdu<-subset(transmigroup, groupivdu=='Injection drug users')

# Primary robust
ivdu_prim<-subset(ivdu, primNADCstyp !='Brain and CNS' & # rma
                    primNADCstyp !='Breast' &
                    primNADCstyp !='Colon, rectosigmoid junction, rectum, anus and anal canal' &
                    primNADCstyp !='Eye and adnexa'&
                    primNADCstyp !='Leukaemia'&
                    primNADCstyp !='Liver, gallbladder, and biliary tract'&
                    primNADCstyp !='Melanoma of skin' &
                    primNADCstyp !='Mesothelial and soft tissue' &
                    primNADCstyp !='Penis'&
                    primNADCstyp !='Skin nonmelanoma'&
                    primNADCstyp !='Testis'&
                    primNADCstyp !='Prostate'&
                    primNADCstyp !='Vulva and vagina')
ivdu_prim$primNADCstyp<-as.factor(as.character(ivdu_prim$primNADCstyp))
summary(ivdu_prim$primNADCstyp)

num<-length(levels(ivdu_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(ivdu_prim, ivdu_prim$primNADCstyp==levels(ivdu_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(ivdu_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

ivdu_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
ivdu_sirpri$group<-'ivdu'

### robust for secondary
ivdu_sec<-subset(ivdu, !is.na(ivdu$secNADCstyp))

ivdu_sec<-subset(ivdu_sec, ivdu_sec$secNADCstyp!='Colon and Rectum' &
                   secNADCstyp!='Nasopharynx'&
                   secNADCstyp!='Salivary glands' &
                   secNADCstyp!='Testis seminoma')

ivdu_sec$secNADCstyp<-as.factor(as.character(ivdu_sec$secNADCstyp))

summary(ivdu_sec$secNADCstyp)

num_2<-length(levels(ivdu_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(ivdu_sec, secNADCstyp==levels(ivdu_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(ivdu_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

ivdu_sirsec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
ivdu_sirsec$group<-'ivdu'

###rma
ivdu_indep<-subset(ivdu, primNADCstyp=='Prostate'|
                     primNADCstyp=='Brain and CNS')

ivdu_indep$primNADCstyp<-as.factor(as.character(ivdu_indep$ primNADCstyp))
summary(ivdu_indep$primNADCstyp)

num_3<-length(levels(ivdu_indep$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(ivdu_indep, primNADCstyp==levels(ivdu_indep$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(ivdu_indep$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}


ivdu_sirs_indep_prim<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

ivdu_sirs_indep_prim$group<-'ivdu'

ivdusir_final<-bind_rows(ivdu_sirpri, ivdu_sirsec, ivdu_sirs_indep_prim) # bind results together


ivdusir_final[,2:4] <-round(ivdusir_final[,2:4],2)
ivdusir_final[,7] <-round(ivdusir_final[,7],0)

######## summarize number of cancer cases
num<-data.frame(tapply(ivdu$O, ivdu$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)

num_sec<-data.frame(tapply(ivdu$O, ivdu$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:10)

num<-bind_rows(num, num_sec)
ivdusf<-merge(ivdusir_final, num, by='Name', all.x = T)

save(ivdusf, file='ivdusir.RData')

############## heterosexual contact ##############
options(scipen = 9999999)
transmigroup <- read_excel("SIR/subgroup/Categorical/transmigroup.xlsx")

het<-subset(transmigroup, groupivdu=='Heterosexual contact')

# Primary robust
het_prim<-subset(het, primNADCstyp !='Brain and CNS' &
                   primNADCstyp !='Breast' &
                   primNADCstyp !='Colon, rectosigmoid junction, rectum, anus and anal canal' &
                   primNADCstyp !='Eye and adnexa'&
                   primNADCstyp !='Lip, oral cavity and pharynx'& ## rma
                   primNADCstyp !='Liver, gallbladder, and biliary tract'&
                   primNADCstyp !='Melanoma of skin' &
                   primNADCstyp !='Mesothelial and soft tissue' &
                   primNADCstyp !='Prostate'&
                   primNADCstyp !='Penis'&
                   primNADCstyp !='Skin nonmelanoma'&
                   primNADCstyp !='Testis'&
                   primNADCstyp !='Vulva and vagina')
het_prim$primNADCstyp<-as.factor(as.character(het_prim$primNADCstyp))
summary(het_prim$primNADCstyp)

num<-length(levels(het_prim$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(het_prim, het_prim$primNADCstyp==levels(het_prim$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(het_prim$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

het_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
het_sirpri$group<-'Heterosexual contact'

### robust for secondary
het_sec<-subset(het, !is.na(het$secNADCstyp))

het_sec<-subset(het_sec, het_sec$secNADCstyp!='Colon and Rectum' &
                  secNADCstyp!='Liver'&
                  secNADCstyp!='Tonsil' &
                  secNADCstyp!='Testis seminoma')

het_sec$secNADCstyp<-as.factor(as.character(het_sec$secNADCstyp))

summary(het_sec$secNADCstyp)

num_2<-length(levels(het_sec$secNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num_2) {
  fitsir_2<-robu(formula = yi ~ 1, data = subset(het_sec, secNADCstyp==levels(het_sec$secNADCstyp)[i]),
                 studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(het_sec$secNADCstyp)[i]
  SIR[i]<-exp(fitsir_2$reg_table[1, 2])
  LL[i]<-exp(fitsir_2$reg_table[1, 7])
  UL[i]<-exp(fitsir_2$reg_table[1, 8])
  num_study[i]<-fitsir_2$N
  num_outcome[i]<-fitsir_2$M
  I2[i]<-fitsir_2$mod_info$I.2
  Tau[i]<-fitsir_2$mod_info$tau.sq
}

het_sirsec<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
het_sirsec$group<-'Heterosexual contact'

###rma
het_indep<-subset(het, primNADCstyp =='Lip, oral cavity and pharynx'|
                    primNADCstyp =='Liver, gallbladder, and biliary tract')

het_indep$primNADCstyp<-as.factor(as.character(het_indep$ primNADCstyp))
summary(het_indep$primNADCstyp)

num_3<-length(levels(het_indep$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(het_indep, primNADCstyp==levels(het_indep$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(het_indep$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}


het_sirs_indep_prim<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

het_sirs_indep_prim$group<-'Heterosexual contact'

hetsir_final<-bind_rows(het_sirpri, het_sirsec, het_sirs_indep_prim) # bind results together


hetsir_final[,2:4] <-round(hetsir_final[,2:4],2)
hetsir_final[,7] <-round(hetsir_final[,7],0)

######## summarize number of cancer cases
num<-data.frame(tapply(het$O, het$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:28)

num_sec<-data.frame(tapply(het$O, het$secNADCstyp, FUN=sum))
colnames(num_sec)[1]<-'Observed number of cancer'
num_sec[,1]  <-round(num_sec[,1], 0)
num_sec$Name<-rownames(num_sec)
rownames(num_sec)<-c(1:10)

num<-bind_rows(num, num_sec)
hetsf<-merge(hetsir_final, num, by='Name', all.x = T)

save(hetsf, file='hetsir.RData')

load("~/Desktop/Academia/Projects/NADCs-risk and mortality/NADCs-R project/hetsir.RData")
load("~/Desktop/Academia/Projects/NADCs-risk and mortality/NADCs-R project/ivdusir.RData")
load("~/Desktop/Academia/Projects/NADCs-risk and mortality/NADCs-R project/msmsir.RData")

transmigpsf<-bind_rows(hetsf, ivdusf, msmsf)
transmigpsf<-arrange(transmigpsf, Name)
transmigpsf$'SIR (95%CI)' <-paste0(transmigpsf$SIR," ","(", transmigpsf$LL, "-", transmigpsf$UL, ")")
transmigpsf$lgSIR<-log10(transmigpsf$SIR)
transmigpsf$lgul<-log10(transmigpsf$UL)
transmigpsf$lgll<-log10(transmigpsf$LL)


write_xlsx(transmigpsf, path ="SIR/subgroup/Categorical/sirtransmission_result.xlsx" )


############################################ HIV vs AIDS ############################################
hivaids <- read_excel("SIR/subgroup/Categorical/hivaids.xlsx")
attach(hivaids)
hivaids$study_id<-ifelse(!is.na(Cohortab), Cohortab, Study)
detach(hivaids)
hivaids<- escalc(measure="IRLN", xi=O, ti=E, data=hivaids)

write_xlsx(hivaids, path = 'SIR/subgroup/Categorical/hivaids.xlsx')

####################### AIDS ###########
options(scipen = 9999999)
aid<-subset(hivaids, group2=='AIDS')

# rma
aid_rma<-subset(aid, primNADCstyp== "Anus and anal canal"|
                  primNADCstyp== "Bladder"|
                  primNADCstyp== "Bone and joints"|
                  primNADCstyp=="Brain and CNS"|
                  primNADCstyp=="Breast"|
                  primNADCstyp=="Esophagus"|
                  primNADCstyp=="Hodgkin’s lymphoma"|
                  primNADCstyp=="Kidney and renal pelvis"|
                  primNADCstyp=="Larynx"|
                  primNADCstyp=="Leukaemia"|
                  primNADCstyp=="Liver"|
                  primNADCstyp=="Melanoma of skin"|
                  primNADCstyp=="Multiple myeloma"|
                  primNADCstyp=="Nasal cavity, middle ear, and accessory sinuses"|
                  primNADCstyp=="Pancreas"|
                  primNADCstyp=="Penis"|
                  primNADCstyp=="Prostate"|
                  primNADCstyp=="Skin nonmelanoma"|
                  primNADCstyp=="Small intestine"|
                  primNADCstyp=="Stomach"|
                  primNADCstyp=="Testis"|
                  primNADCstyp=="Thymus, heart, mediastinum, and pleura"|
                  primNADCstyp=="Trachea, bronchus, and lung" )

aid_rma$primNADCstyp<-as.factor(as.character(aid_rma$primNADCstyp))
summary(aid_rma$primNADCstyp)

num_3<-length(levels(aid_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(aid_rma, primNADCstyp==levels(aid_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(aid_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

aid_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

aid_rma_result$group<-'AIDS'

# multilevel
aid_ml<-subset(aid, primNADCstyp=="Colon and Rectum"|
                 primNADCstyp=="Mesothelial and soft tissue")

aid_ml$primNADCstyp<-as.factor(as.character(aid_ml$primNADCstyp))

num<-length(levels(aid_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
PHeterogenity<-0
I2<-0

for (i in 1:num) {
  data = subset(aid_ml, aid_ml$primNADCstyp==levels(aid_ml$primNADCstyp)[i])
  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | Study/origNADCstyp,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(aid_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]

  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))

}


aid_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                          num_level3,I2)

aid_ml_result$group<-'AIDS'

# robust
aid_robu<-subset(aid, primNADCstyp =='Lip, oral cavity and pharynx')

aid_robu$primNADCstyp<-as.factor(as.character(aid_robu$primNADCstyp))
summary(aid_robu$primNADCstyp)

num<-length(levels(aid_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(aid_robu, aid_robu$primNADCstyp==levels(aid_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(aid_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

aid_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
aid_sirpri$group<-'AIDS'

############## HIV not AIDS ##############
hiv<-subset(hivaids, group2=='HIV not AIDS')

# rma
hiv_rma<-subset(hiv, primNADCstyp== "Anus and anal canal"|
                  primNADCstyp== "Bladder"|
                  #primNADCstyp== "Bone and joints"|
                  primNADCstyp=="Brain and CNS"|
                  primNADCstyp=="Breast"|
                  primNADCstyp=="Esophagus"|
                  primNADCstyp=="Hodgkin’s lymphoma"|
                  primNADCstyp=="Kidney and renal pelvis"|
                  primNADCstyp=="Larynx"|
                  primNADCstyp=="Leukaemia"|
                  primNADCstyp=="Liver"|
                  primNADCstyp=="Melanoma of skin"|
                  primNADCstyp=="Multiple myeloma"|
                  primNADCstyp=="Nasal cavity, middle ear, and accessory sinuses"|
                  primNADCstyp=="Pancreas"|
                  primNADCstyp=="Penis"|
                  primNADCstyp=="Prostate"|
                  primNADCstyp=="Skin nonmelanoma"|
                  primNADCstyp=="Small intestine"|
                  # primNADCstyp=="Stomach"|
                  primNADCstyp=="Testis"|
                  #primNADCstyp=="Thymus, heart, mediastinum, and pleura"|
                  primNADCstyp=="Trachea, bronchus, and lung" )

hiv_rma$primNADCstyp<-as.factor(as.character(hiv_rma$primNADCstyp))
summary(hiv_rma$primNADCstyp)

num_3<-length(levels(hiv_rma$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(hiv_rma, primNADCstyp==levels(hiv_rma$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(hiv_rma$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

hiv_rma_result<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

hiv_rma_result$group<-'HIV not AIDS'

# multilevel
hiv_ml<-subset(hiv, primNADCstyp=="Colon and Rectum"|
                 primNADCstyp=="Mesothelial and soft tissue")

hiv_ml$primNADCstyp<-as.factor(as.character(hiv_ml$primNADCstyp))

num<-length(levels(hiv_ml$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
num_outcome<-0
num_level1<-0
num_level2<-0
num_level3<-0
I2<-0

for (i in 1:num) {
  data<-subset(hiv_ml, hiv_ml$primNADCstyp==levels(hiv_ml$primNADCstyp)[i])
  fourlevel<-rma.mv(yi,
                    vi,
                    random = ~ 1 | Study/origNADCstyp,
                    data = data,
                    method = "REML",control=list(optimizer="optim"),
                    digits=2)

  Name[i]<-levels(hiv_ml$primNADCstyp)[i]
  SIR[i]<-exp(fourlevel$beta) # point estimate
  LL[i]<-exp(fourlevel$ci.lb) # 95%CI, LL
  UL[i]<-exp(fourlevel$ci.ub) # 95%CI, UL
  num_outcome[i]<-fourlevel$k ## numer of effect size
  num_level1[i]<-fourlevel$s.nlevels[1]
  num_level2[i]<-fourlevel$s.nlevels[2]
  num_level3[i]<-fourlevel$s.nlevels[3]

  W <- diag(1/data$vi)

  X <- model.matrix(fourlevel)

  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  I2[i] <- 100 * sum(fourlevel$sigma2) / (sum(fourlevel$sigma2) + (fourlevel$k-fourlevel$p)/sum(diag(P)))

}


hiv_ml_result<-data.frame(Name, SIR, LL, UL, num_outcome, num_level1, num_level2,
                          num_level3,PHeterogenity)

hiv_ml_result$group<-'HIV not AIDS'

# robust
hiv_robu<-subset(hiv, primNADCstyp =='Lip, oral cavity and pharynx')

hiv_robu$primNADCstyp<-as.factor(as.character(hiv_robu$primNADCstyp))
summary(hiv_robu$primNADCstyp)

num<-length(levels(hiv_robu$primNADCstyp))
SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1:num) {
  fitsir<-robu(formula = yi ~ 1,
               data = subset(hiv_robu, hiv_robu$primNADCstyp==levels(hiv_robu$primNADCstyp)[i]),
               studynum = study_id, var.eff.size = vi)
  Name[i]<-levels(hiv_robu$primNADCstyp)[i]
  SIR[i]<-exp(fitsir$reg_table[1, 2])
  LL[i]<-exp(fitsir$reg_table[1, 7])
  UL[i]<-exp(fitsir$reg_table[1, 8])
  num_study[i]<-fitsir$N
  num_outcome[i]<-fitsir$M
  I2[i]<-fitsir$mod_info$I.2
  Tau[i]<-fitsir$mod_info$tau.sq
}

hiv_sirpri<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)
hiv_sirpri$group<-'HIV not AIDS'

hvsir_result<-bind_rows(aid_rma_result,
                        aid_ml_result,aid_sirpri,
                        hiv_rma_result, hiv_ml_result,
                        hiv_sirpri)
hvsir_result<-arrange(hvsir_result, Name)

######## summarize number of cancer cases
num<-data.frame(tapply(hivaids$O, hivaids$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
rownames(num)<-c(1:33)

hvsf<-merge(hvsir_result, num, by='Name', all.x = T)

write_xlsx(hvsf, path ="SIR/subgroup/Categorical/hivaids_result.xlsx")

############################ HAART use ####################################################################################
art <- read_excel("SIR/subgroup/Categorical/haart.xlsx")
attach(art)
art$study_id<-ifelse(!is.na(Cohortab), Cohortab, Study)
detach(art)
art<- escalc(measure="IRLN", xi=O, ti=E, data=art)

write_xlsx(art, path = 'SIR/subgroup/Categorical/haart.xlsx')

art$primNADCstyp[art$primNADCstyp=='Liver, gallbladder, and biliary tract']<-'Liver'
arty<-subset(art, group2=='HAART use')
arty$primNADCstyp<-as.factor(arty$primNADCstyp)
summary(arty$primNADCstyp)

art_meta<-subset(arty, primNADCstyp=='Hodgkin’s lymphoma'|
                   primNADCstyp=='Liver'|
                   primNADCstyp=='Trachea, bronchus, and lung')

art_meta$primNADCstyp<-as.factor(as.character(art_meta$primNADCstyp))
summary(art_meta$primNADCstyp)

num_3<-length(levels(art_meta$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(art_meta, primNADCstyp==levels(art_meta$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(art_meta$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

artuse<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

artuse$group<-'HAART use'

# add cancer cases
num<-data.frame(tapply(art_meta$O, art_meta$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)

artuse<-merge(artuse, num, by='Name', all.x = T)


## no ART use

artn<-subset(art, group2=='No HAART use')
artn$primNADCstyp<-as.factor(artn$primNADCstyp)
summary(artn$primNADCstyp)

nart_meta<-subset(artn, primNADCstyp=='Hodgkin’s lymphoma'|
                    primNADCstyp=='Liver'|
                    primNADCstyp=='Trachea, bronchus, and lung')

nart_meta$primNADCstyp<-as.factor(as.character(nart_meta$primNADCstyp))
summary(nart_meta$primNADCstyp)

num_3<-length(levels(nart_meta$primNADCstyp))

SIR<-0
Name<-0
SIR<-0
LL<-0
UL<-0
I2<-0
Tau<-0
num_study<-0
num_outcome<-0

for (i in 1: num_3) {
  fit_norob<-rma.glmm(xi=O, ti=E, measure = "IRLN",
                      data =subset(nart_meta, primNADCstyp==levels(nart_meta$primNADCstyp)[i]),method="ML")

  Name[i]<-levels(nart_meta$primNADCstyp)[i]
  SIR[i]<-exp(fit_norob$beta)
  LL[i]<-exp(fit_norob$ci.lb)
  UL[i]<-exp(fit_norob$ci.ub)

  num_study[i]<-fit_norob$k
  num_outcome[i]<-fit_norob$k

  I2[i]<-fit_norob$I2
  Tau[i]<-fit_norob$tau2

}

nartuse<-data.frame(Name, SIR, LL, UL,num_study, num_outcome, I2, Tau)

nartuse$group<-'No HAART use'

# add cancer cases
num<-data.frame(tapply(nart_meta$O, nart_meta$primNADCstyp, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)

nartuse<-merge(nartuse, num, by='Name', all.x = T)

art<-bind_rows(artuse, nartuse)

art<-arrange(art, Name)

write_xlsx(art, path = 'SIR/subgroup/Categorical/HAART use_result.xlsx')


#
########################################## HAART year ##########################
# data preparation
hp <- read_excel("haartyplot.xlsx")
hp<-within(hp, {
  lgsir<-log10(SIR)
  lgll<-log10(LL)
  lgul<-log10(UL)
  Name<-factor(Name)
  group<-factor(group)
  final_SIR <-paste0(round(SIR, 2)," ","(", round(LL, 2), "-", round(UL,2), ")")
  case_group<-paste0(Name, group)

}

# add cancer case
harrtyear <- read_excel("SIR/subgroup/Continuous/harrtyear.xlsx")
harrtyear$case_group<-paste0(harrtyear$primNADCstyp, harrtyear$group_yearfinal)
num<-data.frame(tapply(harrtyear$O, harrtyear$case_group, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
colnames(num)[2]<-'case_group'

harrty<-merge(hp, num, by='case_group', all.x = T)
write_xlsx(harrty, path = "haartyplot.xlsx")

# plot figure
hp <- read_excel("haartyplotnew.xlsx")

# text
hcal_table <- cbind(
  c("Cancer and SIR (95%CI)",  hp$final_SIR),
  c("Cancer cases", hp$`Observed number of cancer`),
  c("No. Studies", hp$num_study),
  c("No. Outcomes", hp$num_outcome),
  c("Heterogenity", round(hp$I2, 0)))

# data
hcal_data <-subset(hp, select = c('lgsir', 'lgll', 'lgul'))
colnames(hcal_data)<-c("mean", "lower", "upper")

hcal_data<-hcal_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

styles <- fpShapesGp(
  lines = list(
    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),


    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),

    gpar(col = "#FB9A99"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4"),

    gpar(col = "#E31A1C"),
    gpar(col = "#E31A1C"),
    gpar(col = "#FB9A99"),
    gpar(col = "#A6CEE3"),
    gpar(col = "#1F78B4")

  ),

  box = list(
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),
    gpar(col = "#1F78B4",fill = "#1F78B4"),

    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),
    gpar(col = "#1F78B4",fill = "#1F78B4"),

    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),
    gpar(col = "#1F78B4",fill = "#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),
    gpar(col = "#1F78B4",fill = "#1F78B4"),

    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),

    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),
    gpar(col = "#1F78B4",fill = "#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),

    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99",fill = "#FB9A99"),
    gpar(col = "#A6CEE3",fill = "#A6CEE3"),

    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),

    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),

    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"),

    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#E31A1C", fill = "#E31A1C"),
    gpar(col = "#FB9A99", fill = "#FB9A99"),
    gpar(col = "#A6CEE3", fill = "#A6CEE3"),
    gpar(col = "#1F78B4", fill="#1F78B4"))
)


forestplot(hcal_table,
           hcal_data, new_page = TRUE,
           clip= c(0.1, 5),
           lwd.ci=1.2,  lwd.xaxis=1,  lwd.zero=1,
           zero=0,
           graph.pos = 3,
           xlab = "SIR",
           boxsize = .5,
           hrzl_lines = list("2" = gpar(lty = 1,lwd=1)),
           line.margin = .5,
           align=c("l","l","l","l","l"),
           is.summary = c(rep(TRUE,3),
                          rep(FALSE,4), TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,3),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,3),TRUE,
                          rep(FALSE,3),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,3),
                          rep(TRUE,2),
                          rep(FALSE,2), TRUE,
                          rep(FALSE,3),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,4),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,3), TRUE,
                          rep(FALSE,4)),

           vertices = TRUE,
           ci.vertices.height= 0.2,
           xticks=c(-1, 0, 1, 2, 3),
           xticks.digits=0,
           graphwidth= unit(9, "cm"),
           txt_gp=fpTxtGp(label=gpar(fontsize=15, cex =0.8),
                          ticks=gpar(fontsize=15, cex=0.8),
                          xlab=gpar(fontsize=15,cex = 0.8)),
           col = fpColors(zero = "black"),
           shapes_gp = styles)


###################################################################### age ####################################################################################
# data preparation
age <- read_excel("SIR/subgroup/Continuous/age_result.xlsx")
age<-within(age, {
  lgsir<-log10(SIR)
  lgll<-log10(LL)
  lgul<-log10(UL)
  Name<-factor(Name)
  group<-factor(group)
  final_SIR <-paste0(round(SIR, 2)," ","(", round(LL, 2), "-", round(UL,2), ")")
  case_group<-paste0(Name, group)

})

# add cancer case
age_raw <- read_excel("SIR/subgroup/Continuous/age.xlsx")
age_raw$case_group<-paste0(age_raw$primNADCstyp, age_raw$agefinal)
num<-data.frame(tapply(age_raw$O, age_raw$case_group, FUN=sum))
colnames(num)[1]<-'Observed number of cancer'
num[,1]  <-round(num[,1], 0)
num$Name<-rownames(num)
colnames(num)[2]<-'case_group'

agep<-merge(age, num, by='case_group', all.x = T)
write_xlsx(agep, path = "ageplot.xlsx")


display.brewer.pal(n = 6, name = 'Blues')
brewer.pal(n = 6, name = "Blues")
"#08519C" "#9ECAE1" "#3182BD" "#08519C"
# plot figure
age <- read_excel("ageplot.xlsx")

# text
age_table <- cbind(
  c("Cancer and SIR (95%CI)",  age$final_SIR),
  c("Cancer cases", age$`Observed number of cancer`),
  c("No. Studies", age$num_study),
  c("No. Outcomes", age$num_outcome),
  c("Heterogenity", round(age$I2, 0)))

# data
age_data <-
  subset(age, select = c('lgsir', 'lgll', 'lgul'))
colnames(age_data)<-c("mean", "lower", "upper")

age_data<-age_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

styles <- fpShapesGp(
  lines = list(
    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#3182BD"),
    gpar(col = "#9ECAE1"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#3182BD"),
    gpar(col = "#9ECAE1"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#9ECAE1"),

    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#08519C"),
    gpar(col = "#3182BD"),
    gpar(col = "#9ECAE1"),

    gpar(col = "#3182BD"),
    gpar(col = "#3182BD"),
    gpar(col = "#9ECAE1"),

    gpar(col = "#3182BD"),
    gpar(col = "#3182BD"),
    gpar(col = "#9ECAE1")),

  box = list(
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#3182BD", fill = "#3182BD"),
    gpar(col = "#9ECAE1", fill = "#9ECAE1"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#3182BD", fill = "#3182BD"),
    gpar(col = "#9ECAE1", fill = "#9ECAE1"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#9ECAE1", fill = "#9ECAE1"),

    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#08519C", fill = "#08519C"),
    gpar(col = "#3182BD", fill = "#3182BD"),
    gpar(col = "#9ECAE1", fill = "#9ECAE1"),

    gpar(col = "#3182BD", fill = "#3182BD"),
    gpar(col = "#3182BD", fill = "#3182BD"),
    gpar(col = "#9ECAE1", fill = "#9ECAE1"),

    gpar(col = "#3182BD", fill = "#3182BD"),
    gpar(col = "#3182BD", fill = "#3182BD"),
    gpar(col = "#9ECAE1", fill = "#9ECAE1")
  )
)


forestplot(age_table,
           mean=age_data$mean, lower = age_data$lower, upper = age_data$upper,
           new_page = TRUE,
           clip= c(0.1, 2.5),
           lwd.ci=1.2,  lwd.xaxis=1,  lwd.zero=1,
           zero=0,
           graph.pos = 3,
           xlab = "SIR",
           boxsize = .4,
           hrzl_lines = list("2" = gpar(lty = 1)),
           line.margin = .5,
           align=c("l","l","l","l","l"),
           is.summary = c(rep(TRUE,3),
                          rep(FALSE,3), TRUE,
                          rep(FALSE,3),TRUE,
                          rep(FALSE,2),TRUE, TRUE,
                          rep(FALSE,3),TRUE,
                          rep(FALSE,2),TRUE,
                          rep(FALSE,2)),
           vertices = TRUE,
           ci.vertices.height= 0.2,
           xticks=c( -1, 0, 1, 2, 3),
           xticks.digits=0,
           graphwidth= unit(9, "cm"),
           txt_gp=fpTxtGp(label=gpar(fontsize=15, cex =0.8),
                          ticks=gpar(fontsize=15, cex=0.8),
                          xlab=gpar(fontsize=15,cex = 0.8)),
           col = fpColors(zero = "black"),
           shapes_gp = styles)





















































































































































