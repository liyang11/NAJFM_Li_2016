
//FILE CREATED ON 3/3/2011 BY   TRAVIS BRENDEN & KYLE MOLTON & WEIHAI LUI
//MODIFIED BY YANG LI  FOR WHITEFISH MSE SIMULATION PROJECT 08/12/1012
//LAST MODIFIED 3/10/2015

GLOBALS_SECTION
  #include "admodel.h"
  #include <time.h>

  //define constant variable 
  const double MathPI = 3.141592654;        // or using MI_PI
  const double EPS = 1e-30; 
  const double MathE = 2.71828183; 

DATA_SECTION
  !! ad_comm::change_datafile_name("simulation1.dat");
  init_int numSim                           // num sims
  init_int numYear                          // total yrs
  init_int numAge                           // total ages
  init_int numPop                           // num pops/productivity: 3
  init_int numMix                           // num mixing scenaros  : 2
  init_int numTagA                          // num Target mortality : 3
  init_int numRecV                          // num Recruit Variation: 2
  init_int assessYr                         // ass starting yr=21
  init_int assessStepYr                     // ass every 3 year
  init_vector natMort(1,numPop)             // M
  init_matrix targetA(1,numTagA,1,numPop)             // target A (According to control rule)
  init_vector lvbGrowth(1,3)                // LVB growth paras(Linf, K, to)
  init_vector weightPar(1,2)                // weihgt condition paras(alpha, beta)
  init_number propFemale                    // prop. of fem
  init_number eggsKg                        // eggs/kg   
  init_vector logistMat(1,3)                // logistic maturity paras(maturity at inf length,curv param, length at inflect of logistic model) 
  init_matrix recruitPar(1,numPop,1,2)      // #column for(fecundity alpha, density depend. beta, std for recruit.) 
  init_matrix recruitInt(1,numPop,1,numTagA)
  init_matrix recruitVar(1,numRecV,1,2)
  init_vector gammaSel(1,2)                 // selectivity at age
  init_number q                             // catchability for 4 independent pops
  init_vector sdEffort(1,numPop)            // std error for fishing effort
  init_matrix move(1,numMix,1,numPop)      // movement parameters square matrix
  init_vector implSdTAC(1,numPop)           // implement error for TAC(sd)
  init_vector stdObsCatch(1,numPop)         // std error for observed catch for assessment model
  init_vector effSampleSizes(1,2)            // effective sample size 200 or 800
  init_adstring outFile                     // output file name
  init_adstring outPerform1                  // output file name for independent assessment model
  init_adstring outPerform2                  // output file name for pooled assessment model
  init_int useAllYr                         // use all year data for assessment, 1 = all,0 != all(provide yr# later) =0
  init_int nYrAssess                        // use how many yr data for assess =20
  init_int assessFAge                       // ass 1st age =3
  init_int assessLAge                       // ass last age =12
  init_int lastNumYrMort                    // how many last yrs mortality as mean mortality for TAC
  init_int lastNumYrRec                     // number of last yrs on first age as mean recruitment for TAC
  init_int useRnd                           // use fixed or random seed
  init_int seed                             // if using fixed, then use this seed as provided in dat file,ow using system time as seed
  init_int model                                                                // model1,2,3,4 represent four assessment methods
  init_int mixLev
  init_int TarALev
  init_int RecVLev
  init_int forLinux                        // 0 for windows, 1 for linux

  // ***************************************following derived directly from the input data
  vector initAbund(1,numPop)                // initial abundance /population
  number Linf                               // growth para
  number K                                  // growth para
  number t0                                 // growth para
  number wgtAlpha                           // weight condition para
  number wgtBeta                            // weight condition para
  number matInf                             // maturity parameter
  number matK                               // maturity parameter
  number matFlect                           // maturity parameter
  vector targetTotMort(1,numPop)            // vector for Z's
  vector recruAlpha_Pri(1,numPop)               // logistic recruit parameter
  number rec_St_Var                         // stationary var for AR(1) recruitment
  vector recruAlpha(1,numPop)               // logistic recruit parameter
  vector recruBeta(1,numPop)                // logistic recruit parameter
  number recruStd                           // logistic recruit parameter
  number recruAcor                          // logistic recruit parameter
  number effSampleSize
  vector ages(1,numAge)                     // age vector  (12)
  vector years(1,numYear)                   // year vector (100)
  vector len(1,numAge)                      // fish length at age 
  vector weight(1,numAge)                   // weight at length
  vector maturity(1,numAge)                 // maturity at length
  vector select(1,numAge)                   // selectivity at length or age
  vector fishMort(1,numPop)                 // fishing mortality/population
  

  //following holding the sim result data
  5darray subPop(1,numSim,1,numYear,1,numPop,1,numPop,1,numAge)                 // sub population abundance
  5darray subMixPop(1,numSim,1,numYear,1,numPop,1,numPop,1,numAge)              // sub population abundance only after mixing
  5darray subCatch(1,numSim,1,numYear,1,numPop,1,numPop,1,numAge)               // harvest number for sub pop
  4darray areaCatch(1,numSim,1,numYear,1,numPop,1,numAge)                       // catch at age for each area
  4darray pop(1,numSim,1,numYear,1,numPop,1,numAge)                             // population abundance at age (at the end of each yr)
  4darray mixStock(1,numSim,1,numYear,1,numPop,1,numAge)                        // area abundance at age
  4darray subSSB(1,numSim,1,numYear,1,numPop,1,numPop)                          // cross all age, (pop,area)
  4darray subSpawn(1,numSim,1,numYear,1,numPop,1,numPop)                        // sub population abundance
  3darray SSB(1,numSim,1,numYear,1,numPop)                                      // cross all pop
  3darray Spawn(1,numSim,1,numYear,1,numPop)                                    // cross all pop
  3darray Recruit(1,numSim,1,numYear,1,numPop)                                  // cross all pop
  3darray effort(1,numSim,1,numYear,1,numPop)                                   // fishing effort for each area
  4darray popBegYr(1,numSim,1,numPop,1,numYear,1,numAge)                        // pop abundance at the beginning of yr

  4darray mixStock1(1,numSim,1,numPop,1,numYear,1,numAge)                       //switch second and 3rd dimension as mixStock 
  4darray areaCatch1(1,numSim,1,numPop,1,numYear,1,numAge)                      //switch second and 3rd dimension as areaCatch    
  4darray areaYield1(1,numSim,1,numPop,1,numYear,1,numAge)                      //switch second and 3rd dimension as areaYield
  4darray pop1(1,numSim,1,numPop,1,numYear,1,numAge)                            //switch second and 3rd dimension as pop
  
  int assessOK                                                                  // control which yr to do assessment
  int assessLastOK
  vector readInTAC(1,numPop)                                                    // from assessment model
  vector realTAC(1,numPop)                                                      // add the implementation error
  vector newF(1,numPop)                                                         // estimate fishing mort for each area, from Newton-raphson on TAC
  vector warnConverg(1,numPop)                                                  // keep track of convergence warning for assessment run
  vector warnHessian(1,numPop)                                                  // keep track of hessian  warning for assessment run
  int totAssess                                                                 // total assessment running
  3darray annualTAC(1,numSim,1,numYear,1,numPop)
  3darray annualYield(1,numSim,1,numPop,1,numYear)
  matrix indicator(1,numSim,1,numPop)
  matrix hess_ind(1,numSim,1,numPop)
  matrix counter(1,numSim,1,numPop)
  matrix hess_count(1,numSim,1,numPop)
  int totalAssessYrs
 
  
  //added by yang

  number readInTACpooled                                                        // tmp TACpooled value from ass model
  vector OutTACpooled(1,numYear)
  matrix OutTAC(1,numPop,1,numYear)
  //number catchability                                                           
  matrix pooledeffort(1,numSim,1,numYear)
  3darray pooledYield(1,numSim,1,numYear,1,numAge)    
  3darray pooledpopBegYr(1,numSim,1,numYear,1,numAge)
  3darray pooledcatch(1,numSim,1,numYear,1,numAge)
  3darray F(1,numSim,1,numYear,1,numPop)                                        // fishing mortality for each pop each yr
  number pooledwarnConverg                                                      // keep track of convergence warning for assessment run
  number pooledwarnHessian                                                      // keep track of hessian  warning for assessment run
  vector pooledindicator(1,numSim)
  vector pooledhess_ind(1,numSim)
  vector pooledcounter(1,numSim)
  vector pooledhess_count(1,numSim)
  3darray yieldTolarea(1,numSim,1,numYear,1,numPop)
  3darray catchTolarea(1,numSim,1,numYear,1,numPop)
  3darray areaCPE(1,numSim,1,numYear,1,numPop) 
  matrix TotalCPE(1,numSim,1,numYear)   
  3darray CPEpar(1,numSim,1,numYear,1,numPop)
  matrix pooledTAC(1,numSim,1,numYear)
  matrix assess_N(1,numPop,assessFAge,assessLAge)                                //area specific abundance *weight
  vector assess_ssb(1,numPop)
  number pooled_assess_ssb                            //used to cal SSB mre, this is abundance *weight ~SSB
  3darray N_mre(1,numSim,1,numPop,1,numYear)                                //mre abundance
  3darray N_mare(1,numSim,1,numPop,1,numYear)                                //mare abundance
  matrix pooled_ssb_mre(1,numSim,1,numYear)                                    // total abundance *weight for 4 areas mre
  matrix pooled_ssb_mare(1,numSim,1,numYear)                                   // total abundance *weight for 4 areas mare                                                                   
  matrix pooledcatch_mre(1,numSim,1,numYear)
  matrix pooledcatch_mare(1,numSim,1,numYear)  
  3darray ssb_mre(1,numSim,1,numPop,1,numYear)
  3darray ssb_mare(1,numSim,1,numPop,1,numYear)
  3darray catch_mre(1,numSim,1,numPop,1,numYear)                                // readin MRE  from assessment model
  3darray catch_mare(1,numSim,1,numPop,1,numYear)                               // readin MARE from assessment model
  matrix movePar(1,numPop,1,numPop)                                             // movement parameters square matrix
  3darray R_err(1,numSim,1,numYear,1,numPop)

  //add by yang for check result
  //matrix n_est(1,20,3,12) 
  //vector N_rc_est(1,20)  
  //vector new_rc_cch_est(1,20)  
  //vector rc_s_est(3,12)
  //matrix rc_s_est(1,20,3,12)   
  //vector Ftac_est(3,12)
  //vector tacAge_est(1,10) 


 // ********begin loc_cals to filled data in to targetTotMort,  growth paras,   weight condition paras,  maturity paras, logistic recruit paras... 
 LOC_CALCS
    
   for(int i=1;i<=numPop;i++)  
   targetTotMort(i)=-log(1.0-targetA(TarALev,i)); // set the targetTotMort to Z instead of A, ie make it instantaneous
   //following derived directly from the input data
   effSampleSize = effSampleSizes(model);
   Linf       = lvbGrowth(1);         // growth parameter
   K          = lvbGrowth(2);
   t0         = lvbGrowth(3);
   wgtAlpha   = weightPar(1);         // weight condition parameter
   wgtBeta    = weightPar(2);
   matInf     = logistMat(1);         // maturity parameter
   matK       = logistMat(2);
   matFlect   = logistMat(3);
   recruAlpha_Pri = column(recruitPar,1); // read in value was alpha'= alpha * exp(var/2) for E(R)=alpha'*S*exp(-beta*S)
   recruBeta  = column(recruitPar,2);
   initAbund  = column(recruitInt,TarALev); 
   recruAcor  = recruitVar(RecVLev,1);
   recruStd  = recruitVar(RecVLev,2);
   rec_St_Var = square(recruStd)/(1-square(recruAcor));
   recruAlpha = recruAlpha_Pri*mfexp(-rec_St_Var/2);  // read in value was alpha'= alpha * exp(var/2) for E(R)=alpha'*S*exp(-beta*S)

   ages.fill_seqadd(1,1);                                 // fill age vector
   years.fill_seqadd(1,1);                                // fill year vector
   fishMort  = targetTotMort-natMort;
   len       = Linf*(1.-mfexp(-K*(ages-t0)));             // fish length at age, vector
   weight    = wgtAlpha*pow(len,wgtBeta);                 // weight at length
   maturity  = matInf/(1.+mfexp(-matK*(len- matFlect)));  // maturity at length
   maturity(1,3)=0;
   select.initialize();
   for(int i=assessFAge; i<=numAge; i++)  select(i) = pow(i,gammaSel(1))*mfexp(-gammaSel(2)*i); // selectivity at age/len
   select   /= select(10);                               // changed from 8 to 10   
   dvector stay(1,numPop);
   for(int i=1; i<=numPop; i++)   stay(i)= 1- move(mixLev,i);
   dvector stayother(1,numPop);
   for(int i=1; i<=numPop; i++)   stayother(i)= sum(stay)-stay(i);
   for(int i=1; i<=numPop; i++){
    for (int j=1; j<=numPop; j++){
     if (j==i) movePar(i,j)=stay(i);
     if (j!=i) movePar(i,j)=(1-stay(i))* stay(j)/stayother(i);
   }
  } 

  
  
 END_CALCS


  4darray mixBegYr(1,numSim,1,numYear,1,numPop,1,numAge) // mixture abundance for each area at beginning of yr
  !! totalAssessYrs=int(ceil(double(numYear-assessYr+1)/assessStepYr)); // used for defined variables for mre 
  ivector assessYearVec(1,totalAssessYrs)
  !! assessYearVec.fill_seqadd(assessYr,assessStepYr); 
  //added to record grad
  //3darray Grad(1,numSim,1,numPop,1,numYear)        // readin Grad from independent assessment model
  //matrix pooledGrad(1,numSim,1,numYear)            // readin Grad from pooled assessment model
  
  
  
PARAMETER_SECTION
  init_number dummy
  objective_function_value fake


PRELIMINARY_CALCS_SECTION 
  // ******************************************************initialize vectors
  // ******************************************************change random seed based on value from data file
  indicator =0.;
  hess_ind  =0.;
  counter   =0.;
  hess_count=0.;
  if(useRnd == 1) seed=(unsigned)time(0);                 // using system time as seed for total random
  random_number_generator rnd(seed); 
  subPop.initialize();      pop.initialize();       subCatch.initialize();  areaCatch.initialize(); 
  mixStock.initialize();    subSSB.initialize();    SSB.initialize();       warnConverg.initialize();
  warnHessian.initialize(); popBegYr.initialize();  subSpawn.initialize();  Spawn.initialize();
  Recruit.initialize(); F.initialize();  R_err.initialize(); 
  
  pooledpopBegYr.initialize();yieldTolarea.initialize();pooledTAC.initialize();
  catchTolarea.initialize(); areaCPE.initialize();TotalCPE.initialize();CPEpar.initialize();
  N_mre.initialize(); N_mare.initialize();mixBegYr.initialize();
  pooled_ssb_mre.initialize();pooled_ssb_mare.initialize(); catch_mre.initialize(); catch_mare.initialize();subMixPop.initialize();
  pooledcatch_mre.initialize();pooledcatch_mare.initialize();
  OutTAC.initialize();OutTACpooled.initialize(); ssb_mre.initialize();ssb_mare.initialize();
  assess_N.initialize(); assess_ssb.initialize(); 
  pooledindicator   =0.;
  pooledhess_ind    =0.;
  pooledcounter     =0.;
  pooledhess_count  =0.;
  pooledwarnConverg =0.;
  pooledwarnHessian =0.;
  pooled_assess_ssb =0.;


   // *************************************data genetation begins********************************
   // year 1
  totAssess=0;
  for(int i=1;i<=numSim;i++){          // **simulation loop  i
    assessOK=assessYr;                 // set assessment start yr 
    assessLastOK=assessOK-assessStepYr;
    for(int j=1;j<=numYear;j++){       // **year loop j
    //*********************************   first year data simulation 
      if(j == 1){ 
        for(int n=1;n<=numPop;n++){    // **population  loop n
          double  tmpPop=initAbund(n); // initial bundance for each pop at first year, first age
          popBegYr(i,n,j,1)=tmpPop;    // keep track of abundance at beginning of yr, first age hold initial abund,other ages as 0
          for(int k=1;k<=numPop;k++){  // **area  loop k
            for(int m=1;m<=numAge;m++){// **age loop m
              if(m == 1){              // age=1                              
                subMixPop(i,j,n,k,m)=tmpPop*movePar(n,k);                                          // only apply mixing on age 0,1st                
                subPop(i,j,n,k,m)=tmpPop*movePar(n,k)*mfexp(-(natMort(k)+select(m)*fishMort(k)));  // apply mortality on this age0, 1st yr ; after fishing: P'= P * exp(-(M+F))
                subCatch(i,j,n,k,m)=(select(m)*fishMort(k)/(natMort(k)+select(m)*fishMort(k)))*tmpPop*movePar(n,k)*(1.-mfexp(-(natMort(k)+select(m)*fishMort(k))));
                                       // catch= P*(1-exp(-(M+F)))*(F/(M+F)) 
              }
              else {                   // age>=2
                                       // initial other ages(based on previous age),1st yr; assume P'[age2, 0yr]=P'[age1,1yr]
               subPop(i,j,n,k,m)=subPop(i,j,n,k,m-1)*mfexp(-(natMort(k)+select(m)*fishMort(k)));   // eg. P`[age2]= P`[age1] * exp(-(M+F))
               subCatch(i,j,n,k,m)=(select(m)*fishMort(k)/(natMort(k)+select(m)*fishMort(k)))*subPop(i,j,n,k,m-1)*(1.-mfexp(-(natMort(k)+select(m)*fishMort(k)))); 
                                       // eg. catch[age2]= P`[age1]*(1-exp(-(M+F)))*(F/(M+F))
              }
              double tmp=subPop(i,j,n,k,m)*eggsKg*propFemale*maturity(m)*weight(m);  
              subSSB(i,j,n,k)+= tmp;   // accumulate SSB, cross all age 
              subSpawn(i,j,n,k)+=(tmp/eggsKg)/weight(m);
              pop(i,j,n,m)+=subPop(i,j,n,k,m);         // sum up all area to each pop
              if(m>=2) popBegYr(i,n,j,m)=pop(i,j,n,m-1);
              mixBegYr(i,j,k,m)+=subMixPop(i,j,n,k,m); // sum up all population to each area just after mixing
              mixStock(i,j,k,m)+=subPop(i,j,n,k,m);    // sum up all population to each area            
              areaCatch(i,j,k,m)+=subCatch(i,j,n,k,m); // sum up all subgroup catch to each area              
          }                                         // end age m             
            SSB(i,j,n)+=subSSB(i,j,n,k);            // sum up all area to each population SSB
            Spawn(i,j,n)+=subSpawn(i,j,n,k);
                     
        }                                           // end area k 
        R_err(i,j,n)=0;
      }    
      for(int k=1;k<=numPop;k++) {
       effort(i,j,k)=1./q*fishMort(k)*mfexp((sdEffort(k)*randn(rnd))-(0.5*square(sdEffort(k))));      // fishing effort:  E= F*exp(sd-0.5*sigma^2) /q
       F(i,j,k)= fishMort(k);
      }
      if (model==2) {
       for(int k=1;k<=numPop;k++){
        for(int m=1;m<=numAge;m++)  yieldTolarea(i,j,k) += areaCatch(i,j,k,m)*weight(m);     
        yieldTolarea(i,j,k)=yieldTolarea(i,j,k)*mfexp((stdObsCatch(k)*randn(rnd))-(0.5*square(stdObsCatch(k))));     
        areaCPE(i,j,k)=  yieldTolarea(i,j,k)/effort(i,j,k);                                    
        TotalCPE(i,j)+=  areaCPE(i,j,k);  
       }        
       for(int k=1;k<=numPop;k++)  CPEpar(i,j,k)= areaCPE(i,j,k)/ TotalCPE(i,j);
      }
      
      if(model==2){                                 // begin model 2
       for(int m=1;m<=numAge;m++){                  // begin  m
        for(int k=1;k<=numPop;k++){                 // begin  k
         pooledcatch(i,j,m)+= areaCatch(i,j,k,m);   // sum up area catch to pooled catch 
         pooledpopBegYr(i,j,m)+= popBegYr(i,k,j,m); // cal pooledmixBegYr
        }                                           // end k      
        pooledYield(i,j,m)=pooledcatch(i,j,m)*weight(m);
       }                                            // end m
       for(int k=1;k<=numPop;k++)  pooledeffort(i,j)+= effort(i,j,k);
      }                                             // end model 2    
    } 
    // *************************************************end of first year data simulation 
    //year 2-21
    // **************************************************************after first year &  before starting assessment year data simulation
      if(j > 1 && j<=assessYr){ 
        for(int n=1;n<=numPop;n++){    
          dvector tmpPop(1,numAge);
          R_err(i,j,n)= recruStd*randn(rnd)+ recruAcor*R_err(i,j-1,n);
          tmpPop(1)=recruAlpha(n)*SSB(i,j-1,n)*mfexp(-recruBeta(n)*SSB(i,j-1,n))*mfexp(R_err(i,j,n)); // using SSB from previous year
          Recruit(i,j,n)=tmpPop(1);                               // record recruitment
          for(int t=2;t<=numAge;t++) tmpPop(t)=pop(i,j-1,n,t-1);  //initialize abun at age: eg. initial tmpPop[age2,yr2]=Pop[age1, end of yr1] (no lastyr class) 
          tmpPop(numAge)+=pop(i,j-1,n,numAge);                    //initialize last yr class: tmp[age12,yr2] =Pop[age11, end of yr1] + Pop[age12, end of yr1]    
          popBegYr(i,n,j)=tmpPop;                                 //keep track of abundance at all age at beginning of yr
          for(int k=1;k<=numPop;k++){  
            for(int m=1;m<=numAge;m++){ 
              subMixPop(i,j,n,k,m)=tmpPop(m)*movePar(n,k);
              subPop(i,j,n,k,m)=tmpPop(m)*movePar(n,k)*mfexp(-(natMort(k)+select(m)*fishMort(k)));
              subCatch(i,j,n,k,m)=(select(m)*fishMort(k)/(natMort(k)+select(m)*fishMort(k)))*tmpPop(m)*movePar(n,k)*(1.-mfexp(-(natMort(k)+select(m)*fishMort(k))));
              double tmp=subPop(i,j,n,k,m)*eggsKg*propFemale*maturity(m)*weight(m);  
              subSSB(i,j,n,k)+= tmp; 
              subSpawn(i,j,n,k)+=(tmp/eggsKg)/weight(m);
              pop(i,j,n,m)+=subPop(i,j,n,k,m); 
              mixBegYr(i,j,k,m)+=subMixPop(i,j,n,k,m); 
              mixStock(i,j,k,m)+=subPop(i,j,n,k,m);            
              areaCatch(i,j,k,m)+=subCatch(i,j,n,k,m); 
            } 
            SSB(i,j,n)+=subSSB(i,j,n,k);  
            Spawn(i,j,n)+=subSpawn(i,j,n,k);
    
          }            // end area k                
        }              // end pop n 
        for(int k=1;k<=numPop;k++)  {
         effort(i,j,k)=1./q*fishMort(k)*mfexp((sdEffort(k)*randn(rnd))-(0.5*square(sdEffort(k))));      // fishing effort:  E= F*exp(sd-0.5*sigma^2) /q
         F(i,j,k)= fishMort(k);
        }
        if (model==2) {
        for(int k=1;k<=numPop;k++){
         for(int m=1;m<=numAge;m++)  yieldTolarea(i,j,k) += areaCatch(i,j,k,m)*weight(m); 
         yieldTolarea(i,j,k)=yieldTolarea(i,j,k)*mfexp((stdObsCatch(k)*randn(rnd))-(0.5*square(stdObsCatch(k))));        
         areaCPE(i,j,k)=  yieldTolarea(i,j,k)/effort(i,j,k);                                    
         TotalCPE(i,j)+=  areaCPE(i,j,k);  
        }        
        for(int k=1;k<=numPop;k++)  CPEpar(i,j,k)= areaCPE(i,j,k)/ TotalCPE(i,j);
        cout<<"Total Yield"<< yieldTolarea(i,j) <<"  "<<endl;
        cout<<"effort"<< effort(i,j)<<"  "<<endl;
        cout<<"areaCPE"<< areaCPE(i,j)<<"  "<<endl;
        cout<<"CPEpar(for next year use) :"<<CPEpar(i,j)<<"  "<<endl;
       }

        if(model==2){  
         for(int m=1;m<=numAge;m++){ 
          for(int k=1;k<=numPop;k++){  
           pooledcatch(i,j,m)+= areaCatch(i,j,k,m);  
           pooledpopBegYr(i,j,m)+= popBegYr(i,k,j,m);
          }         // end  k
          pooledYield(i,j,m)=pooledcatch(i,j,m)*weight(m);
         }            // end  m
        for(int k=1;k<=numPop;k++)  pooledeffort(i,j)+= effort(i,j,k);
        }             // end model !=1  
      } 
      // end  first year &  before starting assessment year data simulation
      // *********************************************************************
            

      
      
// ****************assessment begin**************
// model=1 :: independent populatation assessment model
// model=2 :: pooled population assessment model (cpe3)
// model=3 :: pooled population assessment model (msy)
// model=4 :: pooled population assessment model (ey)
// Begins in year 22
// ***********************************************

  if(j > assessYr){                        // use new F for year after assessment started
        dmatrix tmpPop(1,numPop,1,numAge); // build the population abundance at start of the yr
        dmatrix tmpPopMove(1,numPop,1,numAge);
        for(int k=1;k<=numPop;k++){ 
          if (model==1) readInTAC(k)= OutTAC(k,j);  //OutTAC has the correct year num defined in assess part of this tpl
          if (model==2) readInTACpooled=OutTACpooled(j);
        }
        for(int n=1;n<=numPop;n++){ 
          R_err(i,j,n)= recruStd*randn(rnd)+ recruAcor*R_err(i,j-1,n);
          tmpPop(n,1)=recruAlpha(n)*SSB(i,j-1,n)*mfexp(-recruBeta(n)*SSB(i,j-1,n))*mfexp(R_err(i,j,n)); // using SSB from previous year
          Recruit(i,j,n)=tmpPop(n,1);
          for(int t=2;t<=numAge;t++) 
          tmpPop(n,t)=pop(i,j-1,n,t-1);    
          tmpPop(n,numAge)+=pop(i,j-1,n,numAge);
          popBegYr(i,n,j)=tmpPop(n);
        }    // end pop n 
        int assessYrnext=assessLastOK+1;
        if(j == assessYrnext){
         int j_correct=j-1;      //we saved the mre/mare as the full assessment year results, but now we are in year full assessment+1
         //cout<<"Here mre CAL"<<endl;
         dvector real_ssb(1,numPop);                                //area specific abundance *weight
         dvector real_N(1,numPop);                               //area specific abundance
         for(int n=1;n<=numPop;n++){ 
          real_ssb(n)=sum(elem_prod(elem_prod(popBegYr(i,n,j),weight),maturity)*propFemale);
          ssb_mre(i,n,j_correct)=(assess_ssb(n)- real_ssb(n))/(real_ssb(n)+EPS);  
          ssb_mare(i,n,j_correct)=abs(assess_ssb(n)- real_ssb(n))/(real_ssb(n)+EPS); 
         }
         double pooled_real_ssb=sum(real_ssb);
         pooled_assess_ssb=sum(assess_ssb); //number of fish * weight~~~SSB assessment result  
         pooled_ssb_mre(i,j_correct)=(pooled_assess_ssb- pooled_real_ssb)/(pooled_real_ssb+EPS); //mre of SSB
         pooled_ssb_mare(i,j_correct)=abs(pooled_assess_ssb- pooled_real_ssb)/(pooled_real_ssb+EPS); //mare of SSB   
        }
	dvector tmp(1,numPop);
        if (model==1){ // begin model 1 : independent population model 1
         for(int k=1;k<=numPop;k++){ 
           tmp(k)=mfexp((implSdTAC(k)*randn(rnd))-(0.5*square(implSdTAC(k))));
           realTAC(k)=readInTAC(k)*tmp(k); // *error
           annualTAC(i,j,k)=readInTAC(k);
         }  // end area k
        }   // end model 1
        
        if (model==2){  // begin model!=1 :     pooled population model 2
          pooledTAC(i,j)= readInTACpooled;
          for(int k=1;k<=numPop;k++){ 
           annualTAC(i,j,k)=readInTACpooled*(CPEpar(i,j-1,k)+CPEpar(i,j-2,k)+CPEpar(i,j-3,k))/3;         
           tmp(k)=mfexp((implSdTAC(k)*randn(rnd))-(0.5*square(implSdTAC(k))));           
           realTAC(k)=annualTAC(i,j,k)*tmp(k); 
          }  // end area k        
        }   // end model !=1
            // update new F each year 
            
        cout<< "######Sim, Year: " <<i<<" "<<j<<endl;
        //cout<< "TAC for each area(without error): " << annualTAC(i,j)<< "  " <<endl;
        //cout<< "TAC for each area(with error): " << realTAC<< "  " <<endl; 
 
        for(int n=1;n<=numPop;n++){      
         for(int k=1;k<=numPop;k++){   
          for(int m=1;m<=numAge;m++){      
           subMixPop(i,j,n,k,m)=tmpPop(n,m)*movePar(n,k); 
           mixBegYr(i,j,k,m)+=subMixPop(i,j,n,k,m);     
           tmpPopMove(k,m) = mixBegYr(i,j,k,m);     
          }       
         }
        }         
        
        NewtonRaphsonForF(tmpPopMove); // get new F based on annualTAC and using the abundance at two year before                  
        cout<<" updated F: "<< newF <<  "  " <<endl;
        
        for(int n=1;n<=numPop;n++){      
          for(int k=1;k<=numPop;k++){   
            for(int m=1;m<=numAge;m++){ 
              subPop(i,j,n,k,m)=tmpPop(n,m)*movePar(n,k)*mfexp(-(natMort(k)+select(m)*newF(k)));
              subCatch(i,j,n,k,m)=(select(m)*newF(k)/(natMort(k)+select(m)*newF(k)))*tmpPop(n,m)*movePar(n,k)*(1.-mfexp(-(natMort(k)+select(m)*newF(k))));
              double tmp=subPop(i,j,n,k,m)*eggsKg*propFemale*maturity(m)*weight(m);  
              subSSB(i,j,n,k)+= tmp; 
              subSpawn(i,j,n,k)+=(tmp/eggsKg)/weight(m);
              pop(i,j,n,m)+=subPop(i,j,n,k,m); 
              mixStock(i,j,k,m)+=subPop(i,j,n,k,m);           
              areaCatch(i,j,k,m)+=subCatch(i,j,n,k,m); 
            } // end age m
            SSB(i,j,n)+=subSSB(i,j,n,k);  
            Spawn(i,j,n)+=subSpawn(i,j,n,k);
          }   // end area k
        }     // end pop n
        for(int k=1;k<=numPop;k++)  {
         effort(i,j,k)=1./q*newF(k)*mfexp((sdEffort(k)*randn(rnd))-(0.5*square(sdEffort(k))));      // fishing effort:  E= F*exp(sd-0.5*sigma^2) /q
          F(i,j,k) = newF(k); 
         }
         
        if (model==2) {
         for(int k=1;k<=numPop;k++){
          for(int m=1;m<=numAge;m++)  yieldTolarea(i,j,k) += areaCatch(i,j,k,m)*weight(m);  
          yieldTolarea(i,j,k)=yieldTolarea(i,j,k)*mfexp((stdObsCatch(k)*randn(rnd))-(0.5*square(stdObsCatch(k))));        
          areaCPE(i,j,k)=  yieldTolarea(i,j,k)/effort(i,j,k);                                    
          TotalCPE(i,j)+=  areaCPE(i,j,k);  
         }        
         for(int k=1;k<=numPop;k++)  CPEpar(i,j,k)= areaCPE(i,j,k)/ TotalCPE(i,j);
        }

        if(model==2){  
         for(int m=1;m<=numAge;m++){ 
          for(int k=1;k<=numPop;k++){  
           pooledcatch(i,j,m)+= areaCatch(i,j,k,m);  
           pooledpopBegYr(i,j,m)+= popBegYr(i,j,k,m);
          }         // end  k
          pooledYield(i,j,m)=pooledcatch(i,j,m)*weight(m);
         }            // end  m
         for(int k=1;k<=numPop;k++)  pooledeffort(i,j)+= effort(i,j,k);
        }             // end model !=1       
     }         // end assessment if branch: j > assessYr
 
 //****************************************************
 //population data simulation ends here
 //****************************************************
  

  //run accessment to get TAC from assessment model
  //the order for following if branch does matters, should not change it 
  //*****assessmentOK ::assessOK+=assessStepYr 
  //**********************assessment on each assessment OK year***********************

      if(j == assessOK){ //start on assessYr, then every 3 yrs after 
                         //check if std file exist or not, if true, then delete it         
       ifstream fexist("sim4assessment.std");
       if(fexist) {
        fexist.close(); 
        if(forLinux==0) int tmp11=system("del sim4assessment.std");
        else int tmp11=system("rm -f sim4assessment.std");
       }
       totAssess++;
       //*****************independent population model begins here***************
       if(model==1){
       
         dvector IntPop(1,numPop);
         dvector RecPop(1,numPop);
         for(int k=1;k<=numPop;k++){ 
          //cout<<endl<<"#==== N  (yr x age) "<<k <<" ===="<<endl; //output catch at age     
          if(useAllYr==1){
            IntPop(k)=mean(popBegYr(i,k,1)(assessFAge+1,assessLAge));
            for(int yr=1;yr<=j-1;yr++) { //use all previous yr data
             dvector TmpPop(1,j-1);
             TmpPop(yr)=popBegYr(i,k,yr,assessFAge);
              //cout<<popBegYr(i,k,yr)<<endl;  //all ages for this yr
             if(yr==j-1)  RecPop(k)=mean(TmpPop);
            }	   
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            IntPop(k)=mean(popBegYr(i,k,startYr)(assessFAge+1,assessLAge));
            for(int yr=startYr;yr<=j-1;yr++){
              dvector TmpPop(startYr,j-1);
              TmpPop(yr)=popBegYr(i,k,yr,assessFAge);
              //cout<<popBegYr(i,k,yr)<<endl;
              if(yr==j-1) RecPop(k)=mean(TmpPop);
            }
           }                                
         }
         IntPop=log(IntPop);
         RecPop=log(RecPop);
         //cout<<"Ini Mean"<<IntPop<<endl;
         //cout<<"Rec Mean"<<RecPop<<endl;
         
         
          for(int k=1;k<=numPop;k++){      
                         //output data file for assessment
          ofstream report("sim4assessment1.dat");
          report<<"# sim,yr,area = "<<i<<", "<<j-1<<", "<<k<<endl; //atucally yr data from
          if(nYrAssess>assessYr-1) nYrAssess=assessYr-1;           //in case user input big number than actual assessment starting yr
          report<<"# num years, first Age, last age "<<endl;
          int tmpYr=(useAllYr==1)? j-1:nYrAssess;
          report<<tmpYr <<" "<<assessFAge<<" "<<assessLAge<<endl;
          report<<endl<<"# natural mort, total target mort "<<endl;
          report<<natMort(k) <<" "<<targetTotMort(k)<<endl;
          report<<endl<<"# weight at age "<<endl;
          report<<weight(assessFAge,assessLAge)<<endl;
          report<<endl<<"# selectivity at age "<<endl;
          report<<select(assessFAge,assessLAge)<<endl; 
          report<<endl<<"# maturity at age "<<endl;
          report<<maturity(assessFAge,assessLAge)<<endl;
          report<<endl<<"# female proportion "<<endl;
          report<<propFemale<<endl;  
          report<<endl<<"#==== catch at age for each area  (yr x age)  ===="<<endl; //output catch at age     
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<areaCatch(i,yr,k)(assessFAge,assessLAge)<<endl;  //all ages for this yr
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<areaCatch(i,yr,k)(assessFAge,assessLAge)<<endl;
          }   
                                        //output fishing effort
          report<<endl<<"#==== fishing effort  (year)  ===="<<endl; 
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<effort(i,yr,k) <<endl;
          }else{                        //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<effort(i,yr,k)<<endl;
          }  
          report<<endl<<"# number of last yrs mortality  as mean mortality for TAC "<<endl<<lastNumYrMort<<endl;
          report<<endl<<"# number of last yrs on first age as mean recruitment for TAC "<<endl<<lastNumYrRec<<endl;
          report<<endl<<"# lognormal std error for observed catch for each area "<<endl<<stdObsCatch(k)<<endl;
          report<<endl<<"# effective sample size "<<endl<<effSampleSize <<endl;
          report<<endl<<"# assess Freq "<<endl<<assessStepYr <<endl;  
          report<<endl<<"# testing "<<endl;
          report<<"-9999 -8888 -7777"<<endl;
          report.close();
          //exit(7);

          if(forLinux==0) int tmp12=system("sim4assessment1 -nox > NUL"); //run assessment program
          else int tmp12=system("./sim4assessment1 -nox > /dev/null");    //run assessment program
          //read in TAC from assessment report file
          ifstream readin("sim4assessment1.rep");  
          for(int q= assessOK+1; q<=assessOK+assessStepYr; q++) readin>>OutTAC(k,q);  //based on full assessment done in year assessOK_start from year 21, we predict TAC for year 22,23,24 or 22-26 based on data of year 1-20, so year 21 is the lag year
          //readin>>readInTAC(k);
          dvariable grad;
          readin>>grad;
          //Grad(i,k,j)= double(grad);
          readin>>assess_N(k);readin>>assess_ssb(k); 
          readin>>catch_mre(i,k,j);readin>>catch_mare(i,k,j); 
          //*******here to check
          //cout<<"sim="<<i<<" yr="<<j<<" area="<<k<<"Fnew"<<newF(k)<<" tac="<<readInTAC(k) <<endl;
          readin.close();
          
          
          //cout<<"sim="<<i<<" yr="<<j<<" area="<<k<<"Fnew"<<newF(k)<<" tac="<<readInTAC(k)<<endl;
          //exit(7);
          
          //if max gradient too big
          if(fabs(grad)>0.001) {
            cout<<"*** Warning: Proper convergence could not be reached, gradient >0.001 !"<<endl;
            warnConverg(k)++;  //convergence warning counter for each area
            indicator(i,k)=1.0;
            counter(i,k)+=1.0;
            if(annualTAC(i,j,k)!=0) readInTAC(k)=annualTAC(i,j,k);
          }
          //if got hessian warning
          ifstream ifs("sim4assessment1.std");
          if(ifs) {
            ifs.close(); 
            if(forLinux==0) int tmp14=system("del sim4assessment1.std");
            else int tmp14=system("rm -f sim4assessment1.std");
          }
          else {
            cout<<"*** Warning: The function maximizer failed !" <<endl;
            warnHessian(k)++;  //convergence warning counter for each area
            hess_ind(i,k)=1.0;
            hess_count(i,k)+=1.0;
            if(annualTAC(i,j,k)!=0) readInTAC(k)=annualTAC(i,j,k);
          }

          if(forLinux==0) int tmp13=system("del sim4assessment1.rep");   
          else int tmp13=system("rm -f sim4assessment1.rep");
                                         
         } //end k area loop        
        } //end if model==1
       
       
       
       //********************pooled population assessment model begins here*****************************
         if( model==2){     
          ofstream report("sim4assessment.dat");         
          report<<"# sim,yr = "<<i<<", "<<j-1<<endl; //atucally yr data from
          if(nYrAssess>assessYr-1) nYrAssess=assessYr-1;  //in case user input big number than actual assessment starting yr
          report<<"# num years, first Age, last age "<<endl;
          int tmpYr=(useAllYr==1)? j-1:nYrAssess;
          report<<tmpYr <<" "<<assessFAge<<" "<<assessLAge<<endl;
          report<<endl<<"# natural mort, total target mort "<<endl;
          report<<natMort(1)<<" "<<targetTotMort(1)<<endl;
          report<<endl<<"# weight at age "<<endl;
          report<<weight(assessFAge,assessLAge)<<endl;
          report<<endl<<"# selectivity at age "<<endl;
          report<<select(assessFAge,assessLAge)<<endl;
          report<<endl<<"# maturity at age "<<endl;
          report<<maturity(assessFAge,assessLAge)<<endl;
          report<<endl<<"# female proportion "<<endl;
          report<<propFemale<<endl; 
          //output catch at age
          report<<endl<<"#==== catch at age for pooled area  (yr x age)  ===="<<endl;
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<pooledcatch(i, yr)(assessFAge,assessLAge)<<endl;
           }else{ //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
             report<<pooledcatch(i, yr)(assessFAge,assessLAge)<<endl;
          }
           
          //output fishing effort
          report<<endl<<"#==== pooled fishing effort  (year)  ===="<<endl;
          if(useAllYr==1){
            for(int yr=1;yr<=j-1;yr++)  //use all previous yr data
              report<<pooledeffort(i,yr) <<endl;
           }else{ //output fixed num yr data
            int startYr=j-nYrAssess;
            for(int yr=startYr;yr<=j-1;yr++)
              report<<pooledeffort(i,yr)<<endl;
          }
          report<<endl<<"# number of last yrs mortality  as mean mortality for TAC "<<endl<<lastNumYrMort<<endl;
          report<<endl<<"# number of last yrs on first age as mean recruitment for TAC "<<endl<<lastNumYrRec<<endl;
          report<<endl<<"# lognormal std error for observed catch for pooled area "<<endl<<stdObsCatch(1)<<endl;
          report<<endl<<"# effective sample size "<<endl<<effSampleSize <<endl;
          report<<endl<<"# assess Freq "<<endl<<assessStepYr <<endl;  
          report<<endl<<"# testing "<<endl;
          report<<"-6666 -5555 -4444"<<endl;          
          report.close();
          //exit(10);
          
          if(forLinux==0) int tmp12=system("sim4assessment -nox > NUL"); //run assessment program
          else int tmp12=system("./sim4assessment -nox > /dev/null"); //run assessment program

          //read in TAC from assessment report file
          ifstream readin("sim4assessment.rep");
          for(int q= assessOK+1; q<=assessOK+assessStepYr; q++) readin>>OutTACpooled(q); //from Year j+1 to j+assessStepYr
          double pooledgrad;
          readin>>pooledgrad;
          readin>>pooled_assess_ssb; 
          readin>>pooledcatch_mre(i, j); readin>>pooledcatch_mare(i, j);
          //readin>>tacAge_est;
          //readin>>Ftac_est;
          //readin>>rc_s_est;
          //for(int q=1; q<=20; q++) {
          //readin>>new_rc_cch_est(q);     
          //readin>>N_rc_est(q);
          //}
          //readin>>n_est;  
          //readin>>catchability;   
          readin.close();
          

         //if max gradient too big
         for(int k=1;k<=numPop;k++){
          if(fabs(pooledgrad)>0.001) {
            cout<<"*** Warning: Proper convergence could not be reached, pooledgradient >0.001 !"<<endl;
            pooledwarnConverg++;  //convergence warning counter for each area
            pooledindicator(i)=1.0;
            pooledcounter(i)+=1.0;
            if(pooledTAC(i,j)!=0) readInTACpooled=pooledTAC(i,j);
          }
         }
          //if got hessian warning
         ifstream ifs("sim4assessment.std");
         if(ifs) {                                                                                                              
            ifs.close(); 
            if(forLinux==0) int tmp14=system("del sim4assessment.std");
            else int tmp14=system("rm -f sim4assessment.std");
          }
         else {
            cout<<"*** Warning: The function maximizer failed !" <<endl;
            pooledwarnHessian++;  //convergence warning counter for each area
            pooledhess_ind(i)=1.0;
            pooledhess_count(i)+=1.0;
            if(pooledTAC(i,j)!=0) readInTACpooled=pooledTAC(i,j);
         }

          if(forLinux==0) int tmp13=system("del sim4assessment.rep");
          else int tmp13=system("rm -f sim4assessment.rep");                                  
         } //end if model!=1
         
         assessLastOK=assessOK;
         assessOK+=assessStepYr; //update to next assess yr
      
     }
     //************************end of each time assessment****************************** 
     if( (j>assessYr) && (j>assessLastOK) && (j<assessOK) ){
      if(model==1){
       for (int k=1;k<=numPop;k++) readInTAC(k)=OutTAC(k,j+1);
      }
      if(model==2){
       readInTACpooled= OutTACpooled(j+1);
      }
     } //end if((j>assessYr)&&(j>assessLastOK) && (j<assessOK)) 
    
   } //end year loop j
  }//end sim loop i

//**********************************end of data generation*************************************
//*********************************************************************************************
  for(int i=1;i<=numSim;i++){
   for(int j=1;j<=numYear;j++){
    for(int k=1;k<=numPop;k++) {
     for(int m=1;m<=numAge;m++){ 
      mixStock1(i,k,j,m)=mixStock(i,j,k,m);
      pop1(i,k,j,m)=pop(i,j,k,m);
      areaCatch1(i,k,j,m)=areaCatch(i,j,k,m);
      areaYield1(i,k,j,m)=areaCatch1(i,k,j,m)*weight(m);
      annualYield(i,k,j)+=areaYield1(i,k,j,m);
     }    
    }
   }
  }


//**********************************result analysis********************************************* 
 if(model==1){
  Saveindependent();
  }
 if (model==2){
  Savepooled();
 } 
  SaveSimFile();
  exit(8);
 
 
 
 
PROCEDURE_SECTION
  exit(9);


// run newton method to estiamte F, based on real TAC
//popN is pop initial abundance at start of yr, with 1st age as recruitment,
// other ages from the N at end of previous yr but excluding the last age(age move up)
FUNCTION int NewtonRaphsonForF(const dmatrix&stockN)
  double criteria=0.0001;
  int maxLoop=100;
  for(int k=1;k<=numPop;k++){ //for each area
    double Fnew=0.1; //starting value, initial as .1
    double holder=0;  //hold Fnew value from last run
    double delta,tmpCatch,tmpCatchH,tmpCatchL;    
    for(int l=1;l<=maxLoop;l++){  //fixed the loop number as 10
      delta=Fnew*0.001; //try make the step small
      tmpCatch=0; tmpCatchH=0;tmpCatchL=0;
      for(int m=1;m<=numAge;m++){ //loop for all ages
       tmpCatch+=(select(m)*Fnew/(natMort(k)+select(m)*Fnew))*stockN(k,m)*(1.-mfexp(-(natMort(k)+select(m)*Fnew)));//use Fnew, fill in catch equation
       tmpCatchH+=(select(m)*(Fnew+delta)/(natMort(k)+select(m)*(Fnew+delta)))*stockN(k,m)*(1.-mfexp(-(natMort(k)+select(m)*(Fnew+delta))));//use Fnew+delta
       tmpCatchL+=(select(m)*(Fnew-delta)/(natMort(k)+select(m)*(Fnew-delta)))*stockN(k,m)*(1.-mfexp(-(natMort(k)+select(m)*(Fnew-delta))));//use Fnew-delta
      }
      //cout<<k<<" "<<l<<"     "<<Fnew<<"     "<<tmpCatch<<"      "<<delta<<endl;
      //following Fratio=f(x0)/f'(x0)
      double Fratio=(tmpCatch-realTAC(k))/((tmpCatchH-tmpCatchL)/(2.*delta));
      Fnew=Fnew-Fratio;  //update Fnew
      if(Fnew<0) Fnew=0.5*(Fnew+Fratio);  //adjust it
      if(Fnew>3) Fnew=3;//At low N, Fnew would sometimes be really high, so limit here
      if(fabs(Fnew-holder)<criteria) break; //if meeting the convergence criteria, then exit
      holder=Fnew;
    }
    newF(k)=Fnew;
  }
  //cout<<endl<<newF<<endl;
  return 0;
 

FUNCTION Saveindependent
  ofstream report(outPerform1);
  /*
  report<<"#============= Number of populations ==================" << endl;
  report << "# numPop" << endl; 
  report<<numPop<<endl; 
  report<<"#=============== Number of iterations =======================" << endl;
  report << "# numSim" << endl; 
  report<<numSim<<endl; 
  */
  

    report<<"#===================== MRE for area-specific ssb (sim )======================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### area_SSB_mre(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_SSB_mre_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<ssb_mre(i,n)(assessYearVec)<<endl<<endl;
    }
  }
  
  report<<"#===================== MARE for area-specific ssb (sim )======================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### area_SSB_mare(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_SSB_mare_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<ssb_mare(i,n)(assessYearVec)<<endl<<endl;
    }
  }  

  
  
  report<<"#===================== MRE for area-specific catch (sim )======================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### area_catch_mre(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_catch_mre_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<catch_mre(i,n)(assessYearVec)<<endl<<endl;
    }
  }
  
  report<<"#===================== MARE for area-specific catch (sim )======================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### area_catch_mare(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_catch_mare_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<catch_mare(i,n)(assessYearVec)<<endl<<endl;
    }
  }

 
  report<<"#============== Convergence Gradient Indicator (sim x area) ===="<<endl;
  report <<"# assess_conv_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<indicator<<endl<<endl;
  report<<"#============== Hession Warning Indicator (sim x area) ===="<<endl;
  report <<"# assess_hess_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<hess_ind<<endl<<endl;
  report<<"#============== Convergence Gradient Counter (sim x area) ===="<<endl;
  report <<"# counter_conv_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<counter<<endl<<endl;
  report<<"#============== Hession Warning Counter (sim x area) ===="<<endl;
  report <<"# counter_hess_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<hess_count<<endl<<endl;
 
  
  report<<"#============== Assessment Convergence Warnings Prop. for each area cross all sim. ===="<<endl;
  report<<"# convergWarn"<<endl; //that is the name when read in R
  report<<setfixed()<<setprecision(6)<<setw(10)<<warnConverg/(totAssess)<<endl<<endl;
  report<<"#============== Assessment Hessian  Warnings Prop. for each area cross all sim. ===="<<endl;
  report<<"# hessianWarn"<<endl; //that is the name when read in R
  report<<setfixed()<<setprecision(6)<<setw(10)<<warnHessian/(totAssess)<<endl<<endl; 
  /*
  report<<"#============== Grad(area x year) ===="<<endl;
  report <<"# Gradient" << endl;
  for(int i=1;i<=numSim;i++)
  {    
   report<<setfixed()<<setprecision(4)<<setw(10)<<Grad(i)<<" ";
   report<<endl;
  }
 report<<"#===================== MRE for area-specific abundance (sim )======================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### area_N_mre(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_N_mre_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<N_mre(i,n)(assessYearVec)<<endl<<endl;
    }
  }

  report<<"#===================== MARE for area-specific abundance (sim )======================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### area_N_mare(Year x population) = "<<n<<endl;
      report<<"### Sim"<<i<<"_N_mare_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(4)<<setw(10)<<N_mare(i,n)(assessYearVec)<<endl<<endl;
    }
  }
  report<<"#======================= MRE for Pooled SSB (sim )  =========================="<<endl;  
  report<<"# Pooled_SSB_mre "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolSSB_mre"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooled_ssb_mre(i)(assessYearVec)<<endl<<endl;
  }
  
  report<<"#======================= MARE for Pooled SSB (sim )  =========================="<<endl;  
  report<<"# Pooled_SSB_mare "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolSSB_mare"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooled_ssb_mare(i)(assessYearVec)<<endl<<endl;
  } 
  
 
  */
  report.close();  
  
  
  
  
  
  
FUNCTION Savepooled
  ofstream report(outPerform2);
  report<<"#============= ***For pooled population populations*** ==================" << endl;
  report<<"#============= *************************************** ==================" << endl;
  /*
  report<<"#=============== Number of iterations =======================" << endl;
  report << "# numSim" << endl; 
  report<<numSim<<endl;
  */
  
  report<<"#======================= MRE for Pooled catch (sim )  =========================="<<endl;  
  report<<"#Pooled_catch_mre "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {   
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolcatch_mre"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooledcatch_mre(i)(assessYearVec)<<endl<<endl;
  }
  
  report<<"#======================= MARE for Pooled catch (sim)  =========================="<<endl;  
  report<<"# Pooled_catch_mare "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolcatch_mare"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooledcatch_mare(i)(assessYearVec)<<endl<<endl;
  }  
  
  report<<"#======================= MRE for Pooled SSB (sim )  =========================="<<endl;  
  report<<"# Pooled_SSB_mre "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolSSB_mre"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooled_ssb_mre(i)(assessYearVec)<<endl<<endl;
  }
  
  report<<"#======================= MARE for Pooled SSB (sim )  =========================="<<endl;  
  report<<"# Pooled_SSB_mare "<<endl; //that is the name when read in R
  for(int i=1;i<=numSim;i++)
  {    
   report<<"## Simulation ( Year)= "<<i<<endl;
   report<<"### Sim"<<i<<"_poolSSB_mare"<<endl; //that is the name when read in R
   report<<setfixed()<<setprecision(4)<<setw(10)<<pooled_ssb_mare(i)(assessYearVec)<<endl<<endl;
  } 
  /*    
  report<<"#================= Area TAC=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation ( Area x Year)= "<<i<<endl;
    report<<"### Sim"<<i<<"_TAC"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(5)<<setw(10)<<annualTAC(i)<<endl<<endl;
  }
   */
  report<<"#============== Convergence Gradient Indicator for pooled population (sim ) ===="<<endl;
  report <<"# Pooled_assess_conv_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<pooledindicator<<endl<<endl;
  report<<"#============== Hession Warning Indicator for pooled population(sim ) ===="<<endl;
  report <<"# Pooled_assess_hess_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<pooledhess_ind<<endl<<endl;
  report<<"#============== Convergence Gradient Counter for pooled population(sim ) ===="<<endl;
  report <<"# Pooled_counter_conv_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<pooledcounter<<endl<<endl;
  report<<"#============== Hession Warning Counter for pooled population(sim ) ===="<<endl;
  report <<"# Pooled_counter_hess_warn" << endl;
  report<<setfixed()<<setprecision(0)<<setw(10)<<pooledhess_count<<endl<<endl;
 
  
  report<<"#============== Assessment Convergence Warnings Prop. for pooled population cross all sim. ===="<<endl;
  report<<"# Pooled_convergWarn"<<endl; //that is the name when read in R
  report<<setfixed()<<setprecision(6)<<setw(10)<<pooledwarnConverg/(totAssess*4)<<endl<<endl;
  report<<"#============== Assessment Hessian  Warnings Prop. for pooled population cross all sim. ===="<<endl;
  report<<"# Pooled_hessianWarn"<<endl; //that is the name when read in R
  report<<setfixed()<<setprecision(6)<<setw(10)<<pooledwarnHessian/(totAssess*4)<<endl<<endl;
  report.close();
  
  
  
// save the simulation  results as file, filename from dat file
FUNCTION SaveSimFile
  ofstream report(outFile);
  report<<"# assessment model"<<endl; //for plot in r
  report<<model<<endl<<endl;
  
  //report<<"# Years"<<endl; //for plot in r
  //report<<years<<endl<<endl;

  report<<endl<<"#====================== SSB(Year x Pop) ========================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    report<<"## Sim"<<i<<"_SSB"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(0)<<setw(15)<<SSB(i)<<endl<<endl;  
  }
 
  report<<endl<<"#====================== R(Year x Pop) ========================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    report<<"## Sim"<<i<<"_Recruits"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(0)<<setw(15)<<Recruit(i)<<endl<<endl; 
  }

  report<<endl<<"#====================== R_error(Year x Pop) ========================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    report<<"## Sim"<<i<<"_SR_ProcessErr"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(4)<<setw(10)<<R_err(i)<<endl<<endl; 
  }
 

  
    report<<"#================= Area Harvest Yield =================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### harvest yield (Year x Age) = "<<n<<endl;
      report<<"### Sim"<<i<<"_Yield_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(2)<<setw(10)<<areaYield1(i,n)<<endl<<endl;
    }
  }
  
    report<<"#==================== Fishing Mortality ====================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation (Year x Area)= "<<i<<endl;
    report<<"### Sim"<<i<<"_F"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(5)<<setw(10)<<F(i)<<endl<<endl;
  }
  
    report<<"#================= Area TAC=================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation ( Area x Year)= "<<i<<endl;
    report<<"### Sim"<<i<<"_TAC"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(5)<<setw(10)<<annualTAC(i)<<endl<<endl;
  }
  /*
    report<<"#==================== Area Fishing Effort ====================="<<endl;
  for(int i=1;i<=numSim;i++)
  {          
    report<<"## Simulation= "<<i<<endl;
    report<<"### Sim"<<i<<"_effort"<<endl; //that is the name when read in R
    report<<setfixed()<<setprecision(5)<<setw(10)<<effort(i)<<endl<<endl;
  }




  
  report<<"#======================= Area Abundance =========================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### mix area abundance (Year x Age) = "<<n<<endl;
      report<<"### Sim"<<i<<"_Abund_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(0)<<setw(10)<<mixStock1(i,n)<<endl<<endl;
    }
  }

  report<<"#===================== Area Catch ======================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### catch N (Year x Age) = "<<n<<endl;
      report<<"### Sim"<<i<<"_Catch_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(1)<<setw(10)<<areaCatch1(i,n)<<endl<<endl;
    }
  }



  report<<"#==================== Population abundance ====================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### population abundance after fishing(Year x Age) = "<<n<<endl;
      report<<"### Sim"<<i<<"_Spawning_Area"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(0)<<setw(10)<<pop1(i,n)<<endl<<endl;
    }
  }

  report<<"#======================= Stocks Population Abundance (beginning of yr) =========================="<<endl;
  for(int i=1;i<=numSim;i++)
  {
    report<<"## Simulation = "<<i<<endl;
    for(int n=1;n<=numPop;n++)
    {
      report<<"### Pop num (Year x Age) = "<<n<<endl;
      report<<"### Sim"<<i<<"_PopBeg_Pop"<<n<<endl; //that is the name when read in R
      report<<setfixed()<<setprecision(0)<<setw(10)<<popBegYr(i,n)<<endl<<endl;
    }
  }
  */
   
  report.close();
