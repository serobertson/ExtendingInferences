
#source('~/simulated_example_source_code4.R')
DF <- read.csv(file="simulated_data.csv", header=TRUE, sep=",",na.string = "")

DF$Y[is.na(DF$Y)] <- 0
DF$A[is.na(DF$A)] <- 0
DF$epsilon<-NULL
DF$inter<-NULL
DF$d_pop<-NULL
DF$d_out<-NULL
DF$Y_1<-NULL
DF$Y_0<-NULL
weights<-generate_weights(Smod=S~X1+X2+X3,Amod=A~X1+X2+X3, data=DF)

DF2<-weights$dat

OM<-OM_est(data=DF2)
DF2$p1<-OM$p1
DF2$p0<-OM$p0

IOW1<-IOW1_est(data=DF2)
IOW2<-IOW2_est(data=DF2)
  
DR1<-DR1_est(data=DF2)
DR2<-DR2_est(data=DF2)
DR3<-DR3_est(data=DF2)

#Estimates of potential outcome Y^1
results_1<-rbind(OM_1=OM$OM_1, IOW1_1=IOW1$IOW1_1, IOW2_1=IOW2$IOW2_1,
               DR1_1=DR1$DR1_1, DR2_1=DR2$DR2_1,  DR3_1=DR3$DR3_1)
print(results_1) 

#Estimates of potential outcome Y^0
results_0<-rbind(OM_0=OM$OM_0, IOW1_0=IOW1$IOW1_0, IOW2_0=IOW2$IOW2_0,
               DR1_0=DR1$DR1_0, DR2_0=DR2$DR2_0,  DR3_0=DR3$DR3_0)
print(results_0)

#Estimates of ATE
results<-rbind(OM=OM$OM, IOW1=IOW1$IOW1, IOW2=IOW2$IOW2,
                 DR1=DR1$DR1, DR2=DR2$DR2,  DR3=DR3$DR3)
print(results)



#-------sandwich variance------#

library("geex")

#OM
(param_start_OM<-c(coef(OM$OM1mod), coef(OM$OM0mod), 
                   m1=OM$OM_1, m0=OM$OM_0,ate=OM$OM))

OM_mest<-m_estimate(
  estFUN = OM_EE,
  data  = DF,
  root_control = setup_root_control(start = param_start_OM),
  compute_roots = T,
  compute_vcov = T
)

#return potential outcomes means and their difference (ate) and corresponding standard errors
OM_m1<-extractEST(geex_output=OM_mest, est_name="m1",param_start=param_start_OM)
OM_m0<-extractEST(geex_output=OM_mest, est_name="m0",param_start=param_start_OM)
OM_ate<-extractEST(geex_output=OM_mest, est_name="ate",param_start=param_start_OM)


#IOW1 

(param_start_IOW1<-c(coef(weights$Smod) , coef(weights$Amod),
                     m1=IOW1$IOW1_1, m0=IOW1$IOW1_0, ate=IOW1$IOW1) )

IOW1_mest <-m_estimate(
  estFUN = IOW1_EE,
  data  = DF,
  root_control = setup_root_control(start = param_start_IOW1),
  compute_roots = T,
  compute_vcov = T
) 

#save variance + SE
IOW1_m1<-extractEST(geex_output=IOW1_mest, est_name="m1",param_start=param_start_IOW1)
IOW1_m0<-extractEST(geex_output=IOW1_mest, est_name="m0",param_start=param_start_IOW1)
IOW1_ate<-extractEST(geex_output=IOW1_mest, est_name="ate",param_start=param_start_IOW1)


#IOW2 sandwich variance

(param_start_IOW2<-c(coef(weights$Smod),coef(weights$Amod),
                     int1=coef(IOW2$IOW1mod)["(Intercept)"],int0=coef(IOW2$IOW0mod)["(Intercept)"],
                     m1=IOW2$IOW2_1, m0=IOW2$IOW2_0, ate=IOW2$IOW2))

IOW2_mest <-m_estimate(
  estFUN = IOW2_EE,
  data  = DF,
  root_control = setup_root_control(start = param_start_IOW2),
  compute_roots = T,
  compute_vcov = T
) 

#save variance + SE
IOW2_m1<-extractEST(geex_output=IOW2_mest, est_name="m1",param_start=param_start_IOW2)
IOW2_m0<-extractEST(geex_output=IOW2_mest, est_name="m0",param_start=param_start_IOW2)
IOW2_ate<-extractEST(geex_output=IOW2_mest, est_name="ate",param_start=param_start_IOW2)


#DR1

(coef_DR1est<-c(coef(OM$OM1mod), coef(OM$OM0mod), m1=DR1$DR1_1, m0=DR1$DR1_0, ate=DR1$DR1))

param_start_DR1<-c(coef(weights$Smod) , coef(weights$Amod), coef_DR1est)

DR1_mest<-m_estimate(
  estFUN = DR1_EE,
  data  = DF,
  root_control = setup_root_control(start = param_start_DR1),
  compute_roots = T,
  compute_vcov = T
) 

#save variance + SE
DR1_m1<-extractEST(geex_output=DR1_mest, est_name="m1",param_start=param_start_DR1)
DR1_m0<-extractEST(geex_output=DR1_mest, est_name="m0",param_start=param_start_DR1)
DR1_ate<-extractEST(geex_output=DR1_mest, est_name="ate",param_start=param_start_DR1)


#DR2

param_start_DR2<-c(coef(weights$Smod), coef(weights$Amod), 0.5, 0.5, 0.5,
                   coef(OM$OM1mod), coef(OM$OM0mod),
                   m1=DR2$DR2_1, m0=DR2$DR2_0, ate=DR2$DR2)

DR2_mest<-m_estimate(
  estFUN = DR2_EE,
  data  = DF2,
  root_control = setup_root_control(start = param_start_DR2),
  compute_roots = T,
  compute_vcov = T
) 

#save variance + SE
DR2_m1<-extractEST(geex_output=DR2_mest, est_name="m1",param_start=param_start_DR2)
DR2_m0<-extractEST(geex_output=DR2_mest, est_name="m0",param_start=param_start_DR2)
DR2_ate<-extractEST(geex_output=DR2_mest, est_name="ate",param_start=param_start_DR2)


#DR3

param_start_DR3<-c(coef(weights$Smod) , coef(weights$Amod), 
                   coef(DR3$DR1mod), coef(DR3$DR0mod),
                   m1=DR3$DR3_1, m0=DR3$DR3_0, ate=DR3$DR3)

DR3_mest<-m_estimate(
  estFUN = DR3_EE,
  data  = DF,
  root_control = setup_root_control(start = param_start_DR3),
  compute_roots = T,
  compute_vcov = T
) 

#save variance + SE
DR3_m1<-extractEST(geex_output=DR3_mest, est_name="m1",param_start=param_start_DR3)
DR3_m0<-extractEST(geex_output=DR3_mest, est_name="m0",param_start=param_start_DR3)
DR3_ate<-extractEST(geex_output=DR3_mest, est_name="ate",param_start=param_start_DR3)

#Estimates of potential outcome mean Y^1
summary_results1<-data.frame(OM_m1, IOW1_m1, IOW2_m1, DR1_m1,DR2_m1, DR3_m1)
print(summary_results1)

#Estimates of potential outcome mean Y^0
summary_results0<-data.frame(OM_m0, IOW1_m0, IOW2_m0, DR1_m0,DR2_m0, DR3_m0)
print(summary_results0)

#Estimates of ate
summary_results<-data.frame(OM_ate, IOW1_ate, IOW2_ate, DR1_ate,DR2_ate, DR3_ate)
print(summary_results)
