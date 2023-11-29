
library(dplyr)
d_time <- as_tibble(d_step1)  # main data frame from preprocessing containing all time series measurements
d_subject <- as_tibble(d_subject_00)  # data frame containing all subject specific information at baseline
all_params <- data_vars_all   # data frame containing all variables (columns) in d00 and further explanation on what they mean

# select which studies to choose
all_PsA_studies <- c("CAIN457F2306", "CAIN457F2312", "CAIN457F2318", "CAIN457F2336", "CAIN457F2342", "CAIN457F2366")  # Note: Do not change order!
all_RA_studies <- c("CAIN457F2201", "CAIN457F2206", "CAIN457F2208", "CAIN457F2302", "CAIN457F2309", "CAIN457F2311")  # Note: Do not change order!

PsA_study <- all_PsA_studies
RA_study <- all_RA_studies
d_time <- d_time %>%  
        mutate(newTRT=ifelse((TRT01P=="adalimumab 40 mg") & STUDYID %in% PsA_study,"Humira",
                ifelse((TRT01P %in% c("AIN457 10mg/kg - 75 mg","AIN457 10mg/kg - 150 mg"))& STUDYID %in% PsA_study,"IV",
                ifelse((TRT01P=="AIN457 300 mg")& STUDYID %in% PsA_study,"300 mg",
                ifelse((TRT01P %in% c("AIN457 150 mg","AIN457 150 mg no load"))& STUDYID %in% PsA_study,"150 mg",
                ifelse((TRT01P %in% c("AIN457 75 mg"))& STUDYID %in% PsA_study,"75 mg",
                ifelse((TRT01P=="placebo")& STUDYID %in% PsA_study,"Placebo",
                ifelse((TRT01P %in% c("AIN457 10mg/kg - 150 mg","AIN457 150 mg","AIN457 150mg","AIN457 150mg s.c. load")) & STUDYID %in% RA_study,"150 mg",
                ifelse((TRT01P %in% c("AIN457 300mg","AIN457 300mg s.c.","AIN457 10mg/kg I.v.")) & STUDYID %in% RA_study,"300 mg",
                ifelse((TRT01P %in% c("AIN457 25mg")) & STUDYID %in% RA_study,"25 mg",
                ifelse((TRT01P %in% c("AIN457 75 mg","AIN457 10mg/kg - 75 mg","AIN457 75mg")) & STUDYID %in% RA_study,"75 mg",
                ifelse((TRT01P %in% c("placebo")) & STUDYID %in% RA_study,"Placebo",
                ifelse((TRT01P %in% c("abatacept")) & STUDYID %in% RA_study,"Abatacept",NA)))))))))))))  %>% 
        mutate(newTRT=as.character(newTRT),STUDYID=as.character(STUDYID),USUBJID=as.character(USUBJID)) %>% 
        select(-c("ARM","TRT01P","TNFRES"))%>% as_tibble()%>% 
        filter(!newTRT=="Abatacept")

d_subject <- d_subject %>% filter(STUDYID  %in% c(PsA_study,RA_study))%>%  
        mutate(newTRT=ifelse((TRT01P=="adalimumab 40 mg") & STUDYID %in% PsA_study,"Humira",
                ifelse((TRT01P %in% c("AIN457 10mg/kg - 75 mg","AIN457 10mg/kg - 150 mg"))& STUDYID %in% PsA_study,"IV",
                ifelse((TRT01P=="AIN457 300 mg")& STUDYID %in% PsA_study,"300 mg",
                ifelse((TRT01P %in% c("AIN457 150 mg","AIN457 150 mg no load"))& STUDYID %in% PsA_study,"150 mg",
                ifelse((TRT01P %in% c("AIN457 75 mg"))& STUDYID %in% PsA_study,"75 mg",
                ifelse((TRT01P=="placebo")& STUDYID %in% PsA_study,"Placebo",
                ifelse((TRT01P %in% c("AIN457 10mg/kg - 150 mg","AIN457 150 mg","AIN457 150mg","AIN457 150mg s.c. load")) & STUDYID %in% RA_study,"150 mg",
                ifelse((TRT01P %in% c("AIN457 300mg","AIN457 300mg s.c.","AIN457 10mg/kg I.v.")) & STUDYID %in% RA_study,"300 mg",
                ifelse((TRT01P %in% c("AIN457 25mg")) & STUDYID %in% RA_study,"25 mg",
                ifelse((TRT01P %in% c("AIN457 75 mg","AIN457 10mg/kg - 75 mg","AIN457 75mg")) & STUDYID %in% RA_study,"75 mg",
                ifelse((TRT01P %in% c("placebo")) & STUDYID %in% RA_study,"Placebo",
                ifelse((TRT01P %in% c("abatacept")) & STUDYID %in% RA_study,"Abatacept",NA))))))))))))) %>% 
        filter(!newTRT=="Abatacept")


non_endpoint_nms <- c("USUBJID", "AVISITN", "STUDYID", "newTRT")

## convert factor to numeric
cols_before <- colnames(d_time)
cols_after <- colnames(d_time)
d_time <- data.frame(d_time[, non_endpoint_nms],
                     lapply(as.data.frame(d_time[, colnames(d_time)[!colnames(d_time) %in% non_endpoint_nms]]),
                            function(x) as.numeric((x)) ))
use_nms<-c("PAP01","CRPSI","ESRSI","DAS28CRP","DAS28ESR","HAQDI","PCS","MCS","PADA01","PADA02","SWJT28","TNJT28")

# ## rename columns 
nms_params_v2 <- all_params %>% 
        filter(PARAMCD %in% use_nms) %>% 
        # filter(!PARAM %in% c("over the last 24 hours","physicians global assessment of disease activity", 
                             # "patients global assessment")) %>% 
        na.omit()
nms_params_v2$PARAM<-gsub("subjects global assessment of disease activity","Patient VAS",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("physicians global assessment of disease activity","Physician VAS",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("erythrocyte sedimentation rate \\(mm\\/h\\)","ESR",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("c\\-reactive protein \\(mg\\/l\\)","CRP",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("disability index score","HAQDI",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("mental component summary","SF36-MCS",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("physical component summary","SF36-PCS",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("rheumatoid arthritis pain today","Pain",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("swollen joint total score for 28 joints","Swollen",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("tender joint total score for 28 joints","Tender",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("psoriatic arthritis pain today","Pain",nms_params_v2$PARAM)
nms_params_v2$PARAM<-gsub("psoriatic arthritis pain today","Pain",nms_params_v2$PARAM)

nms_params_v2 <- na.omit(nms_params_v2[match(use_nms, nms_params_v2$PARAMCD), ] )

d_time$study <- gsub("CAIN457", "", d_time$STUDYID)
d_subject$study <- gsub("CAIN457", "", d_subject$STUDYID)

# Removing space in '300 mg' etc.
d_time$newTRT <- gsub(" ", "", d_time$newTRT)

saveRDS(d_time, file = "../preproc_data/d_time.RDS")
saveRDS(d_subject, file = "../preproc_data/d_subject.RDS")
saveRDS(nms_params_v2, file = "../preproc_data/nms_params_v2")
