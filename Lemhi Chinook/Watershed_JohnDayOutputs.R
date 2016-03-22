###################################################
# This File contains two functions:
#
# CreateJohnDayOut(X1,X2,X3,X4...) which computes summary statistics for each iteration
# WriteJohnDayOut(...) which writes the CSV file containing everything...
# 
# This function is very straightforward, it just does a bunch of post-processing
# It has been moved here for tidiness reasons
# NOTE: at this time this has been designed for a single population and is steelhead centric
###################################################

CreateJohnDayOut<-function(r,results,header,parameters,outstuff)
  {
      #This is a bunch of output bundling for making sure the model properly reflects the input data
      #It's tailored to the John Day O mykiss case, i.e., 1-3 salt returns and smolt ages 1-3
      #Also, it is currently structured for testing purposes only, i.e., single rep runs (testing/validation)
      #First is an adult check, including returns to bonneville, age comp, sex ratio on spawning grounds, etc.
      #Note this only covers a life history with up to 3 ocean ages, plus repeat spawners
      
  #Pete Feb2016 multipop, anywhere kk appears
  #gg<-3 kk<-1
  Gname<-c("Natural","H1","N.H2","H2","N-N.H3","N-H3","N.H3-H3","N.H3-H2",
                 "H3-H2", "H3","N.H3-N.H3")
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
  
  for(gg in 1:header$G)
  {
      for(kk in 1:header$K)
      {
        Gtype<-rep(Gname[gg],header$Tr)
        SimYr<-seq(1:header$Tr)
        REP<-rep(r,length(SimYr))
        SubPop<-rep(kk,length(SimYr))
        SubPopName<-rep(header$Site.Names[kk],length(SimYr))
        ###***End-Pete May 2015 Addition***
        adults<-t(Bonneville_NT[kk,1:4,,gg,r]) #Pete Feb2016 Lewis function
        tot<-adults[,1]+adults[,2]+adults[,3]+adults[,4]
        pct1<-adults[,1]/tot
        pct2<-adults[,2]/tot
        pct3<-(adults[,3]+adults[,4])/tot
        adults<-cbind(adults,pct1,pct2,pct3)
        sexrat<-cbind(Female_Spawners_NR[kk,,gg,r],Male_Spawners_NR[kk,,gg,r],
                      Female_Spawners_NR[kk,,gg,r]/(Male_Spawners_NR[kk,,gg,r]+Female_Spawners_NR[kk,,gg,r]))
        adults<-cbind(REP,SimYr,adults,sexrat,(Male_Spawners_NR[kk,,gg,r]+Female_Spawners_NR[kk,,gg,r])/(tot))
        #Second, some smolt checks/calcs.
        smolty<-t(Candidate_Smolt_ByAge[kk,,,gg,r])[,1:3] #Pete Feb2016 Lewis function
        totsmSAR<-smolty[,1]+smolty[,2]+smolty[,3]
        smolty<-smolty/t(Prod[kk,,,r])[,c(5)] # Must divide by survival from trap to JDA to get it in Sm/Sp units at the trap...
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
        adults<-cbind(Gtype,SubPop,SubPopName,adults)
        stages<-t(NR[kk,,,gg,r])
        stages<-stages[,1:11]
        stagesN5<-t(N5R[kk,,,gg,r])
        stagesN5<-stagesN5[,1:3]
        rainbowdat<-cbind(Resident_Spawners_NR[kk,,gg,r],ResidentSexRatio[kk,,gg,r],
                          t(results$Rainbow_N[kk,2:3,,gg]))
        
        AdRiverSurv<-as.matrix(1-parameters$harvest.wild[kk,])
        AdOcnSurv<-t(parameters$Sr[kk,,])[,8:11]
        SpostSp<-t(parameters$Post_Spawn_Survival_Anadromous_M[kk,,,gg])[,1:3]
        Surv<-cbind(t(Prod[kk,,,r])[,c(2:3)],t(parameters$SR5[kk,,])[,1:3],t(Prod[kk,,,r])[,c(5)],AdOcnSurv,AdRiverSurv,SpostSp)
        Capac<-cbind(t(Cap[kk,,,r])[,1:5],t(parameters$C_ocean[kk,,])[,1])
        SmoltProbF<-t(parameters$N5.Psmolt_F[kk,,])[,1:3]
        SmoltProbM<-t(parameters$N5.Psmolt_M[kk,,])[,1:3]
        Mat<-t(parameters$Mat8Plus_F[kk,,])[,1:4]
        Fec<-t(parameters$Female_Fecundity[kk,,,gg])[,1:4]
        DemParm<-cbind(Surv,SmoltProbF,SmoltProbM,Mat,Fec,Capac)
        
        if(r==1&&kk==1&&gg==1)
        {outstuff<-cbind(adults,smolty,SAR,sm.per.sp,stages,stagesN5,rainbowdat,DemParm)} else
        {outstuff<-rbind(outstuff,cbind(adults,smolty,SAR,sm.per.sp,stages,stagesN5,rainbowdat,DemParm))}
        
      }
      
  }
      return(outstuff)
  }




#This has been moved to this program for tidiness reasons.
WriteJohnDayOut <- function(outstuff)
  {
      namescol<-c("Gtype","SubPopN","Name","rep","yr","1SaltBON","2SaltBON","3SaltBON","4SaltBON","1S.pct","2S.pct","3S+.pct","FemSp","MaleSp","pct.FemSp","Convert")
      namescolsmolt<-c("1smolt","2smolt","3smolt","totsmolt","1sm.pct","2sm.pct","3sm.pct")
      namestages<-c("Nsp","Negg","Nfry","Nparr","Npsm","Nsm","Nad0","Nad1","Nad2","Nad3","Nad4","sm1","sm2","sm3")
      nameDem<-c("Segg","Sfry","Sps1","Sps2","Sps3","Ssm","Sad1o","Sad2o","Sad3o","Sad4o","SadR","Spost1","Spost2","Spost3")
      nameDem<-c(nameDem,"Psm1f","Psm2f","Psm3f","Psm1m","Psm2m","Psm3m","Pad1","Pad2","Pad3","Pad4")
      nameDem<-c(nameDem,"Fec1","Fec2","Fec3","Fec4","Cegg","Cfy","CParr","Cpresm","Csmolt","Cocean")
      namesnames<-c(namescol,namescolsmolt,"SAR","sm.per.sp",namestages,"ResSp","ResFem%","ResEgg","ResFry",nameDem)
      Fname<-paste("JohnDayOutputs",Sys.time(),".csv")
      Fname<-gsub("-","",Fname)
      Fname<-gsub(" ","_",Fname)
      Fname<-gsub(":","",Fname)
      write.table(outstuff,col.names=namesnames,
                  row.names=FALSE,sep = ",",file=Fname)
  }
