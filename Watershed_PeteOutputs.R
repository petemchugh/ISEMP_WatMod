###################################################
# This File contains two functions:
#
# CreatePeteOut(X1,X2,X3,X4...) which computes summary statistics for each iteration
# WritePeteOut(...) which writes the CSV file containing everything...
# 
# This function is very straightforward, it just does a bunch of post-processing
# It has been moved here for tidiness reasons
# NOTE: at this time this has been designed for a single population and is steelhead centric
###################################################


CreatePeteOut<-function(r,results,header,parameters,outstuff)
  {
      #This is a bunch of output bundling for making sure the model properly reflects the input data
      #It's tailored to the John Day O mykiss case, i.e., 1-3 salt returns and smolt ages 1-3
      #Also, it is currently structured for testing purposes only, i.e., single rep runs (testing/validation)
      #First is an adult check, including returns to bonneville, age comp, sex ratio on spawning grounds, etc.
      #Note this only covers a life history with up to 3 ocean ages, plus repeat spawners
      SimYr<-seq(1:header$Tr)
      REP<-rep(r,length(SimYr))
      Male_Spawners_NR[,,,r] = results$Male_Spawners #for computing sex ratio of spawners
      Female_Spawners_NR[,,,r] = results$Female_Spawners
      Bonneville_NT[,,,,r] = results$NT #for computing age comp of adult spawners
      Candidate_Smolt_ByAge[,,,,r] = results$CandSmoltsByAge
      Resident_Spawners_NR[,,,r] =results$RainbowSpawners
      ResidentSexRatio[,,,r] =results$RainbowFemSpawners/results$RainbowSpawners
      NR[,,,,r]=results$N
      N5R[,,,,r]=results$N5
      Prod[,,,r]=results$p
      Cap[,,,r]=results$c
      ###***End-Pete May 2015 Addition***
      adults<-t(Bonneville_NT[1,1:4,,1,r])
      tot<-adults[,1]+adults[,2]+adults[,3]+adults[,4]
      pct1<-adults[,1]/tot
      pct2<-adults[,2]/tot
      pct3<-(adults[,3]+adults[,4])/tot
      adults<-cbind(adults,pct1,pct2,pct3)
      sexrat<-cbind(Female_Spawners_NR[,,1,r],Male_Spawners_NR[,,1,r],
                    Female_Spawners_NR[,,1,r]/(Male_Spawners_NR[,,1,r]+Female_Spawners_NR[,,1,r]))
      adults<-cbind(REP,SimYr,adults,sexrat,(Male_Spawners_NR[,,1,r]+Female_Spawners_NR[,,1,r])/(tot))
      #Second, some smolt checks/calcs.
      smolty<-t(Candidate_Smolt_ByAge[1,,,1,r])[,1:3]
      totsmSAR<-smolty[,1]+smolty[,2]+smolty[,3]
      smolty<-smolty/t(Prod[1,,,r])[,c(5)] # Must divide by survival from trap to JDA to get it in Sm/Sp units at the trap...
      totsm<-smolty[,1]+smolty[,2]+smolty[,3]
      pct1sm<-smolty[,1]/totsm
      pct2sm<-smolty[,2]/totsm
      pct3sm<-smolty[,3]/totsm
      smolty<-cbind(smolty,totsm,pct1sm,pct2sm,pct3sm)
      #Third, some combined SAR and sm/sp calcs
      SAR<-rep(NA,header$Tr)
      sm.per.sp<-SAR
      for(ot in 2:(header$Tr-3))
      {
        SAR[ot]<-(adults[ot+1,3]+adults[ot+2,4]+adults[ot+3,5])/totsmSAR[ot]  # SAR is JDA to BON
        sm.per.sp[ot]<-sum(smolty[ot+1,1],smolty[ot+2,2],smolty[ot+3,3])/(sexrat[ot,1]+sexrat[ot,2]) 
      }
      stages<-t(NR[1,,,1,r])
      stages<-stages[,1:11]
      stagesN5<-t(N5R[1,,,1,r])
      stagesN5<-stagesN5[,1:3]
      rainbowdat<-cbind(Resident_Spawners_NR[,,1,r],ResidentSexRatio[,,1,r],
                        t(results$Rainbow_N[1,2:3,,1]))
      
      AdRiverSurv<-1-t(parameters$harvest.wild)
      AdOcnSurv<-t(parameters$Sr[1,,])[,8:11]
      SpostSp<-t(parameters$Post_Spawn_Survival_Anadromous_M[1,,,1])[,1:3]
      Surv<-cbind(t(Prod[1,,,r])[,c(2:3)],t(parameters$SR5[1,,])[,1:3],t(Prod[1,,,r])[,c(5)],AdOcnSurv,AdRiverSurv,SpostSp)
      Capac<-cbind(t(Cap[1,,,r])[,1:5],t(parameters$C_ocean[1,,])[,1])
      SmoltProbF<-t(parameters$N5.Psmolt_F[1,,])[,1:3]
      SmoltProbM<-t(parameters$N5.Psmolt_M[1,,])[,1:3]
      Mat<-t(parameters$Mat8Plus_F[1,,])[,1:4]
      Fec<-t(parameters$Female_Fecundity[1,,,1])[,1:4]
      DemParm<-cbind(Surv,SmoltProbF,SmoltProbM,Mat,Fec,Capac)
      
      if(r==1)
      {outstuff<-cbind(adults,smolty,SAR,sm.per.sp,stages,stagesN5,rainbowdat,DemParm)}
      else
      {outstuff<-rbind(outstuff,cbind(adults,smolty,SAR,sm.per.sp,stages,stagesN5,rainbowdat,DemParm))}
      
      return(outstuff)
  }




#This has been moved to this program for tidiness reasons.
WritePeteOut <- function(outstuff)
  {
      namescol<-c("rep","yr","1SaltBON","2SaltBON","3SaltBON","4SaltBON","1S.pct","2S.pct","3S+.pct","FemSp","MaleSp","pct.FemSp","Convert")
      namescolsmolt<-c("1smolt","2smolt","3smolt","totsmolt","1sm.pct","2sm.pct","3sm.pct")
      namestages<-c("Nsp","Negg","Nfry","Nparr","Npsm","Nsm","Nad0","Nad1","Nad2","Nad3","Nad4","sm1","sm2","sm3")
      nameDem<-c("Segg","Sfry","Sps1","Sps2","Sps3","Ssm","Sad1o","Sad2o","Sad3o","Sad4o","SadR","Spost1","Spost2","Spost3")
      nameDem<-c(nameDem,"Psm1f","Psm2f","Psm3f","Psm1m","Psm2m","Psm3m","Pad1","Pad2","Pad3","Pad4")
      nameDem<-c(nameDem,"Fec1","Fec2","Fec3","Fec4","Cegg","Cfy","CParr","Cpresm","Csmolt","Cocean")
      namesnames<-c(namescol,namescolsmolt,"SAR","sm.per.sp",namestages,"RBTsp","RBTpctFem","RBTegg","RBTfry",nameDem)
      write.table(outstuff,col.names=namesnames,
                  row.names=FALSE,sep = ",",file="WritePeteOutFile.csv")
  }
