
library(dplyr)


d0_RA<- read.csv(file.path("/data/as/processed/clinical/database_views/","v_proj1_RA_filtered_20221020.csv"))
# d0_RA<- read.csv(file.path("/data/as/processed/clinical/database_views/","v_proj1_RA_filtered_20221005.csv"))
d0_PsA<- read.csv(file.path("/data/as/processed/clinical/database_views/","v_proj1_PsA_20220930.csv"))

PsA_study<-c("CAIN457F2306", "CAIN457F2312", "CAIN457F2318", "CAIN457F2336", "CAIN457F2342", "CAIN457F2366")
RA_study<-c("CAIN457F2201", "CAIN457F2206", "CAIN457F2208", "CAIN457F2302", "CAIN457F2309", "CAIN457F2311")
param_nms<-c("CRPSI","DAS28CRP","DAS28ESR","ESRSI","HAQDI","MCS","PADA01","PADA02","PAP01","PCS","SWJT28","TNJT28")

## for RA new Release
d1_RA<-d0_RA %>% mutate(PARAMCD2=ifelse(PARAMCD=="RADA01","PADA01",
                                        ifelse(PARAMCD=="RADA02","PADA02",
                                               ifelse(PARAMCD=="RAP01","PAP01",PARAMCD))) ) %>% 
  select(-PARAMCD) %>% dplyr::rename(PARAMCD=PARAMCD2) %>% 
  mutate(TRT01P=ifelse(STUDYID=="CAIN457F2208",TRT01A,TRT01P)) %>% 
  mutate(AVISITN=as.numeric(as.character(AVISITN))) %>%
  filter(STUDYID %in% RA_study)%>%
  filter(PARAMCD %in% param_nms)

# longitudinal data d_step1_RA
d_step1_RA <- d1_RA[!duplicated(d1_RA %>% select(STUDYID, USUBJID, AVISITN, PARAMCD)), ]

# long to wide format
# wide1 <- dcast(d_step1_RA, STUDYID+ ARM+TRT01P + USUBJID + AVISITN+TNFRES ~ PARAMCD, value.var = "AVAL") %>% tbl_df() %>%
#   as.data.frame()%>% filter(AVISITN %in% c(0,1,2,3,4,8,12,16,20,24,28,32,36,40,44,48,52))

d_step1_RA<- d_step1_RA %>% select(STUDYID,ARM,TRT01P,USUBJID,AVISITN,TNFRES,PARAMCD,AVAL)

d_step1_RA <- spread(d_step1_RA, "PARAMCD", "AVAL")%>% tbl_df() %>%
  as.data.frame()%>% filter(AVISITN %in% c(0,1,2,3,4,8,12,16,20,24,28,32,36,40,44,48,52))



## for PsA new Release
# d1_PsA<-d0_PsA %>% mutate(PARAMCD2=ifelse(PARAMCD=="RADA01","PADA01",ifelse(PARAMCD=="RADA02","PADA02",ifelse(PARAMCD=="RAP01","PAP01",PARAMCD))) ) %>% select(-PARAMCD) %>% dplyr::rename(PARAMCD=PARAMCD2) %>% 
d1_PsA<-d0_PsA %>% mutate(PARAMCD2=ifelse(PARAMCD=="RADA01","PADA01",
                                          ifelse(PARAMCD=="RADA02","PADA02",
                                                 ifelse(PARAMCD=="RAP01","PAP01",ifelse(PARAMCD=="PAP02","PAP01",PARAMCD)))) ) %>% select(-PARAMCD) %>% dplyr::rename(PARAMCD=PARAMCD2) %>% 
  mutate(TRT01P=ifelse(STUDYID=="CAIN457F2208",TRT01A,TRT01P)) %>% 
  filter(STUDYID %in% PsA_study)%>%
  mutate(AVISITN=as.numeric(as.character(AVISITN)))%>%
  filter(PARAMCD %in% c(param_nms,"PSTSCO"))

# longitudinal data d_step1_PsA
d_step1_PsA <- d1_PsA[!duplicated(d1_PsA %>% select(STUDYID, USUBJID, AVISITN, PARAMCD)), ]

sort(names(d1_PsA))

# long to wide format
# d_step1_PsA <- dcast(d_step1_PsA, STUDYID+ ARM+TRT01P + USUBJID + AVISITN+NTNFGR1+TNFRES ~ PARAMCD, value.var = "AVAL") %>% tbl_df()%>%
#   as.data.frame()%>% filter(AVISITN %in% c(0,1,2,3,4,8,12,16,20,24,28,32,36,40,44,48,52))
d_step1_PsA<- d_step1_PsA %>% select(STUDYID,ARM,TRT01P,USUBJID,AVISITN,TNFRES,PARAMCD,AVAL)

d_step1_PsA <- spread(d_step1_PsA, "PARAMCD", "AVAL")%>% tbl_df() %>%
  as.data.frame()%>% filter(AVISITN %in% c(0,1,2,3,4,8,12,16,20,24,28,32,36,40,44,48,52))

d_skin_PsA<-d_step1_PsA %>% select(USUBJID,PSTSCO,AVISITN) %>% na.omit() %>% unique()
dsubjid_skin_ex1_PsA<-d_skin_PsA %>% filter(AVISITN>0) %>% select(USUBJID) %>% unique() %>% pull()
dsubjid_skin_ex2_PsA<-d_skin_PsA %>% filter(AVISITN!=0) %>% select(USUBJID) %>% unique() %>% pull()
d_skin_id_PsA<-d_skin_PsA$USUBJID[!d_skin_PsA$USUBJID %in% c(dsubjid_skin_ex1_PsA,dsubjid_skin_ex2_PsA)]
d_step1_PsA<-d_step1_PsA %>% mutate(PSTSCO=ifelse(USUBJID %in% d_skin_id_PsA, 0, PSTSCO)) 



# subject data table (baseline data)
d_subject_00_PsA <- d1_PsA %>% filter(AVISITN==0) %>% select(USUBJID, TRT01A,ARM, TRT01P, STUDYID, BMI, RACE, SEX,TNFRES) %>% unique()

# subject data table (baseline data)
d_subject_00_RA <- d1_RA %>% filter(AVISITN==0) %>% select(USUBJID, TRT01A,ARM, TRT01P, STUDYID, BMI, RACE, SEX,TNFRES) %>% unique()


## save parameter names
data_vars_all<-d1_PsA   %>% select(PARAMCD,PARAM) %>% unique()

##combine PsA and RA if needed 
d_subject_00<- d_subject_00_PsA %>% full_join(d_subject_00_RA)
d_step1 <-d_step1_PsA  %>% full_join(d_step1_RA)


dir.create("../preproc_data", recursive = T, showWarnings = FALSE)
saveRDS(d_step1, file = "../preproc_data/d_step_1.RDS")
saveRDS(data_vars_all, file = "../preproc_data/data_vars_all.RDS")
saveRDS(d_subject_00, file = "../preproc_data/d_subject_00.RDS")





