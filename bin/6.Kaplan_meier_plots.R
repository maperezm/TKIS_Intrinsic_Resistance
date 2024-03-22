library(survival)
library(survminer)
library(ggthemes)
setwd("~/Dropbox/Trabajo/Repositorios/Intrinsic_CDDP_Resistan_Scripts/Publicacion/bin/")
source("~/Dropbox/Trabajo/Repositorios/Intrinsic_CDDP_Resistan_Scripts/Erlo_vs_Osi_vs_CDDP/Functions/functions.R")
#--------Load survival data----
p <- list()
p[[1]]  <-  kaplan_meier_plot(fit_obj = fit_EGFR,data_surv =  EGFR_CD74,legend_title = "fff")
p[[2]]  <-  kaplan_meier_plot(fit_EPB41L4A.AS1,legend_title = "EPB41L4A.AS1",data_surv =survival_EPB41L4A.AS1 )
p[[3]]  <-  kaplan_meier_plot(fit_SGMS1.AS1,legend_title = "SGMS1.AS1",data_surv =survival_SGMS1.AS1 )
p[[4]]  <-  kaplan_meier_plot(fit_LINC01133,legend_title = "LINC01133",data_surv =survival_LINC01133 )
p[[5]]  <-  kaplan_meier_plot(fit_CERS6.AS1,legend_title = "CERS6.AS1",data_surv =survival_CERS6.AS1 )
p[[6]]  <-  kaplan_meier_plot(fit_Lnc_PCM1.4,legend_title = "Lnc_PCM1.4",data_surv =survival_Lnc_PCM1.4 )
p[[7]]  <-  kaplan_meier_plot(fit_WWP1.AS1,legend_title = "WWP1.AS1",data_surv =survival_WWP1.AS1 )
p[[8]]  <-  kaplan_meier_plot(fit_HSALNG0005948,legend_title = "HSALNG0005948",data_surv =survival_HSALNG0005948 )
p[[9]]  <-  kaplan_meier_plot(fit_LINC00638,legend_title = "LINC00638",data_surv =survival_LINC00638 )
p[[10]]  <-  kaplan_meier_plot(fit_Lnc_ZBTB20_1,legend_title = "Lnc_ZBTB20_1",data_surv =survival_Lnc_ZBTB20_1 )
p[[11]]  <-  kaplan_meier_plot(fit_Lnc_PXDC1_15,legend_title = "Lnc_PXDC1_15",data_surv =survival_Lnc_PXDC1_15 )
p[[12]]  <-  kaplan_meier_plot(fit_Lnc_KDM5A_1,legend_title = "Lnc_KDM5A_1",data_surv =survival_Lnc_KDM5A_1 )
p[[13]]  <-  kaplan_meier_plot(fit_Lnc_PODXL_3,legend_title = "Lnc_PODXL_3",data_surv =survival_Lnc_PODXL_3 )
p[[14]]  <-  kaplan_meier_plot(fit_Lnc_GATA5_9,legend_title = "Lnc_GATA5_9",data_surv =survival_Lnc_GATA5_9 )


pdf("../kapplan_meier.pdf", width = 14, height =9)

kaplan_meier_plot(fit_CERS6.AS1,legend_title = "CERS6-AS1",data_surv =survival_CERS6.AS1 )
kaplan_meier_plot(fit_AGAP2.AS1,legend_title = "AGAP2-AS1",data_surv =survival_AGAP2.AS1 )
dev.off()
s