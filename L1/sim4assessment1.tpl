//FILE CREATED ON 08/13/09 BY TRAVIS BRENDEN FOR ADMB 2009 SHORT COURSE
//MODIFIED BY Yang Li FOR WHITEFISH MSE SIMULATION PROJECT
//LAST MODIFIED 3/10/2015
//SIMPLE SCAA MODEL

GLOBALS_SECTION
  //following add by liu, should use system random seed in the final delivery version
  #include <time.h>
  #include <admodel.h>
  //int seed=45456; //fixed random seed
  int seed=(unsigned)time(0); //using system time as random seed
  //the pesudo code I used from
  //http://commons.apache.org/math/apidocs/org/apache/commons/math/stat/descriptive/rank/Percentile.html
  //this match to type=6 in quantile function in R, which is also used by SPSS and minitab,
  // but R by defult use type=7, and R also implemented total 9 types
  //require the index of data start from 1

  double quantile(dvector data, double p){
    int n=data.size();
    dvector sortdata(1,n); 
    sortdata=sort(data);
    //cout<<sortdata<<endl<<endl;

    if(n==1) return sortdata(1);
    double pos=p*(n+1); // estimated percentile position
    double d=pos - floor(pos); //difference
    //cout<<pos<<" "<< d<<endl;

    if(pos<1) return sortdata(1);
    else if(pos>=n) return sortdata(n);
    else{
      double lower=sortdata(floor(pos));
      double upper=sortdata(floor(pos)+1);
      return (lower+ d*(upper-lower));
    }
  }

  double quantile(dvar_vector data, double p){
      return quantile(value( data),p);
    
  }
DATA_SECTION
  int fyear
  !!fyear=1;
  init_int lyear                              //Number of years
  init_int fage
  init_int lage                              //Number of ages
  init_number nat_mort
  init_number tot_targ_Z
  init_vector wt_age(fage,lage)         //weight at age
  init_vector sel_age(fage,lage)        //not being used
  init_vector mat_age(fage,lage)        //maturity at age
  init_number propFemale                 //female proportion
  init_matrix rc_catch(fyear,lyear,fage,lage)      //recreational catch
  init_vector rc_eff(fyear,lyear)		        //observed recreational effort
  //following just added by liu
  init_int lastNumYrMort  //how many last yrs mortality as mean mortality for TAC
  init_int   lastNumYrRec    //the last 5yrs first age as mean recruitment
  init_number stdObsCatch  //lognormal std error for observed catch
  init_int effSampleSize  //effective sample size
  init_vector test(1,3)				//test vector for ensuring data are read in correctly
  
  int i                                         //index for year loop
  int j                                         //index for age loop
  matrix M(fyear,lyear,fage,lage)               //Natural mortality matrix
  vector rc_cch(fyear,lyear)		        //observed recreational catch
  matrix rc_cps(fyear,lyear,fage,lage)          //age composition from the recreational catch

  //following new variables add by liu
  matrix obsResampleCatch(fyear,lyear,fage,lage)  //observed catch based on input actual age composition and effective sample size
  vector new_rc_cch(fyear,lyear)   //observed annual total catch after add in obs.error
  matrix new_rc_catch(fyear,lyear,fage,lage)  //obs catch at age based on observed total catch times new age composition
  matrix new_rc_cps(fyear,lyear,fage,lage)    //age composition from resampling readin catch
  number mre
  number mare
  //!!cout << test << endl;
  //!!exit(88);


INITIALIZATION_SECTION

  log_qrc -13.41005              //Recreational catchability
  log_rc_s_p1 2.570625525          //Recreational selectivity parameter 1, estimated on a log scale
  rc_s_p2 1.26                   //Recreational selectivity parameter 2
  log_initial 12.0
  log_mean_recruit 13.5
  log_sigma_rc -1.0



PARAMETER_SECTION

  //Declaration of parameters to be estimated
  init_bounded_number log_qrc(-15.,-10.,1)                 //Recreational catchability coefficient
  init_bounded_number log_rc_s_p1(-1.,5.,1)                //Recreational selectivity parameter 1
  init_bounded_number rc_s_p2(-1.,4.,1)                    //Recreational selectivity parameter 2
  init_bounded_number log_initial(8.,16.,1)                 //Numbers of age-2 and older fish in the first year
  init_bounded_dev_vector initial_devs(fage+1,lage,-8.,8.,1)
  init_bounded_number log_mean_recruit(8.0,17.5,1)           //Average recruitment over the years - log scale
  init_bounded_dev_vector recruit_devs(fyear,lyear,-10.,10.,1)  //Annual recruitment deviations
  init_bounded_dev_vector log_effort_devs(fyear,lyear,-15.,15.,1)
  init_bounded_number log_sigma_rc(-5.,5.,1)                 //log sd for rec catch

  vector initial(fage+1,lage)                  //Exponent of age-2 and older fish in the first year
  number qrc                                   //Exponent of recreational catchability
  number rc_s_p1                               //Exponent of recreational selectivity parameter
  number rc_s_max                              //Maximum of recrational selectivites for scaling to 1
  vector rc_s_age(fage,lage)                   //Recreational selectivities by age
  matrix rc_s(fyear,lyear,fage,lage)           //Matrix of recreational selectivities
  vector re(fyear,lyear)
  vector abs_re(fyear,lyear)

  matrix F(fyear,lyear,fage,lage)              //Instantaneous fishing mortality matrix
  matrix Z(fyear,lyear,fage,lage)              //Total mortality matrix
  matrix S(fyear,lyear,fage,lage)              //Annual survival matrix
  matrix n(fyear,lyear,fage,lage)	       //Predicted abundance at age
  matrix ssb(fyear,lyear,fage,lage)	       //Predicted abundance weight(~ssb) at age
  number Total_ssb                             //Predicted total abundance weight of the year when TAC calculate
  vector N(fyear,lyear)                        //Predicted total abundance
  matrix n_rc(fyear,lyear,fage,lage)	       //Predicted recreational catch at age
  vector N_rc(fyear,lyear)                     //Predicted total recreational catch
  matrix n_rc_cps(fyear,lyear,fage,lage)       //Predicted recreational catch age composition
  matrix d_M(fyear,lyear,fage,lage)            //Predicted natural deaths at age
  vector effort_devs(fyear,lyear)
  vector rc_diff(fyear,lyear)                  //Differences between observed and predicted recreational catches
  number log_sigma_effort
  number variance_ratio
  matrix abund(lyear,lyear+2,fage,lage);
  number L1  // components of the objective function
  number L2
  number L3
  number rc_cps_spl
  vector tacAge(fage,lage) //TAC at age
  vector Ftac(fage,lage) //fishing mort for tac
  vector Mvec(fage,lage)
  vector Zvec(fage,lage)
  vector Svec(fage,lage)

  objective_function_value negLL;

RUNTIME_SECTION
  maximum_function_evaluations 10000
  
  
PRELIMINARY_CALCS_SECTION 
  random_number_generator rnd(seed);

  //creating the natural mortality matrix; M not estimated as part of the model
  M=nat_mort;
   
  rc_cch=rowsum(rc_catch);                     //Total catch in a year
  for(i=fyear;i<=lyear;i++)
  {
    for(j=fage;j<=lage;j++)
    {
      rc_cps(i,j)=rc_catch(i,j)/rc_cch(i);     //age comp by division
    }   
  }

  //cout<<"Total rec catch"<<endl;   
  //cout<<rc_cch<<endl;
  //cout<<"Age composition"<<endl;
  //cout<<rc_cps<<endl;
  //exit(0);

 
  //following add by liu 
  ivector fish(1,effSampleSize); //resampling fish from readin catch_at_age
  for(i=fyear;i<=lyear;i++){
    fish.fill_multinomial(rnd,rc_cps(i)); //based on readin age composition and effect sample size  
    obsResampleCatch(i)= CountFreqAges(fish,fage,lage);  //total is effSampleSize
    new_rc_cps(i)=obsResampleCatch(i)/double(effSampleSize);  //age composition based on resampling from effective sample size
    //add observation error for read in catch
    new_rc_cch(i)=rc_cch(i)*mfexp((stdObsCatch*randn(rnd))-(0.5*square(stdObsCatch)));  //add  lognormal err   
    new_rc_catch(i)=new_rc_cch(i)*new_rc_cps(i);  //based on new obs. annual catch and new age compos., reallocate fish 
  }
  //cout<<endl<<rc_cps <<endl<<endl<<new_rc_cps <<endl<<endl<<new_rc_cch<<endl;  //exit(9);
  //cout<<endl<<obsResampleCatch<<endl<<endl<<rc_catch <<endl<<endl<<new_rc_catch<<endl; exit(9);

PROCEDURE_SECTION
  get_mortality();
  get_population();
  get_catch();
  get_diffs();
  get_objective();




FUNCTION get_mortality
  qrc=mfexp(log_qrc);
  effort_devs=mfexp(log_effort_devs);

  //Recreational selectivity modeled as a gamma function
  rc_s_p1=mfexp(log_rc_s_p1);

  for (j=fage; j<=lage; j++)
  {
    rc_s_age(j)=pow(j,rc_s_p1)*mfexp(-rc_s_p2*(j));
  }

  for (i=fyear;i<=lyear;i++)
  {
  rc_s(i)=rc_s_age/rc_s_age(10);   // scaling
  }

  //Fishing mortality set proportional to effort
  for(i=fyear;i<=lyear;i++)
   {
    for(j=fage;j<=lage;j++)
     {
     F(i,j)=qrc*rc_eff(i)*rc_s(i,j)*effort_devs(i);
     }
   }
  Z= M + F;
  S = mfexp(-1.0*Z);     //Annual survival




FUNCTION get_population 
  for (j=fage+1;j<=lage;j++)
  {
   initial(j)=mfexp(log_initial+initial_devs(j));
  }

  for (i=fyear;i<=lyear;i++)
  {
  n(i,fage)=mfexp(log_mean_recruit+recruit_devs(i));  // yearling numbers over years
  }

  n(fyear)(fage+1,lage)=initial; //set abundance in the fist year for age 2 and older fish

  for (i=fyear+1; i<=lyear; i++)
   {
     for (j=fage+1; j<=lage; j++)
     {
     // n(i,j) =n(i-1,j-1)*S(i-1,j-1);
     n(i,j) =n(i-1,j-1)*S(i-1,j-1);
     }
     n(i,lage)+=n(i-1,lage)*S(i-1,lage);  // the last age group is a plus age group
   }
   
  for (i=fyear; i<=lyear; i++)
   {
    for (j=fage; j<=lage; j++)   ssb(i,j) =  n(i,j) * wt_age(j)*mat_age(j)*propFemale;
   } 
  
  for (i=fyear; i<=lyear; i++)
   {
    N(i)=sum(n(i));          // total population over years

   }




FUNCTION get_catch
  n_rc=elem_prod(elem_div(F,Z),elem_prod(1.0-S,n));	//Baranov catch equation for recreational harvest
  d_M =elem_prod(elem_div(M,Z),elem_prod(1.0-S,n));     //Deaths due to nat mort(M)

  //Estimates total catch by year
  for (i=fyear; i<=lyear; i++)
   {
     N_rc(i)=sum(n_rc(i));
   }

   //Estimated age composition of recreational catch
  for (i=fyear; i<=lyear; i++)
  {
      n_rc_cps(i)=n_rc(i)/N_rc(i);
  }



FUNCTION get_diffs
//calculate differences in total catch, survvey cpe, and juvenile index
  for(i=fyear;i<=lyear;i++)
  {
      //rc_diff(i)=(log(rc_cch(i)+0.000001)-log(N_rc(i)+0.0000001)); //old way, comment by liu
      rc_diff(i)=(log(new_rc_cch(i)+0.000001)-log(N_rc(i)+0.0000001)); 
      re(i)=(N_rc(i)+0.000001-new_rc_cch(i)+0.000001)/(new_rc_cch(i)+0.000001);
      abs_re(i)=re(i);
      if(re(i)<0) abs_re(i)*=-1;
 }
 mre=quantile(re,0.5);
 mare=quantile(abs_re,0.5);

FUNCTION get_objective
  variance_ratio=0.25;
  log_sigma_effort=log(sqrt((1./variance_ratio)*square(mfexp(log_sigma_rc))));
  L1 = size_count(rc_diff)*log_sigma_rc+1.0/(2.0*square(mfexp(log_sigma_rc)))*norm2(rc_diff);       //Likelihood for recreati
  L2 = -sum(double(effSampleSize)*elem_prod(rc_cps,log(0.0000001+n_rc_cps)));     //age comp for recreational fishery, old way
  //L2 = -sum(100.*elem_prod(new_rc_cps,log(0.0000001+n_rc_cps))); //age comp for recreational fishery, change by liu
  L3 = size_count(log_effort_devs)*log_sigma_effort+1.0/(2.0*square(mfexp(log_sigma_effort)))*norm2(log_effort_devs);
  negLL=L1+L2+L3;


//only being called in report section
FUNCTION get_TAC
    //only used for tac calculation
       
  abund(lyear)=n(lyear); //fill in the last yr abundance from n
    //if not believe lyear recru esti
  abund(lyear,fage) =sum(column(n,fage)(lyear-lastNumYrRec,lyear-1) )/double(lastNumYrRec);    //use average recruitment
  
  //use last-year Z as the Z for year last+1  ,use average recruitment
  for(j=fage+1; j<=lage; j++) abund(lyear+1,j) =abund(lyear,j-1)*mfexp(-Z(lyear,j-1));
  abund(lyear+1,lage)+=abund(lyear,lage)*mfexp(-Z(lyear,lage)); //last age group is a plus age group
  //abund(lyear+1,fage)=sum(column(n,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec);    //use average recruitment
    //if not believe lyear recru esti
  abund(lyear+1,fage) =sum(column(n,fage)(lyear-lastNumYrRec,lyear-1) )/double(lastNumYrRec);    //use average recruitment
  
  
  //use average Z as the Z for year last+1  ,use average recruitment
  dvariable tmpMort; 
  for(j=fage+1; j<=lage; j++)
   {
   tmpMort=mean(column(Z,j-1)(lyear-lastNumYrMort+1,lyear) ); //how many last yrs mortality as mean mortality for TAC
   abund(lyear+2,j) =abund(lyear+1,j-1)*mfexp(-tmpMort);  
   }
  tmpMort=mean(column(Z,lage)(lyear-lastNumYrMort+1,lyear) );
  abund(lyear+2,lage)+=abund(lyear+1,lage)*mfexp(-tmpMort);  // the last age group is a plus age group
  //abund(lyear+2,fage)=sum(column(n,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec);  //the last 5yrs first age as mean recruitment
    //if not believe lyear recru esti
  abund(lyear+2,fage) =sum(column(n,fage)(lyear-lastNumYrRec,lyear-1) )/double(lastNumYrRec);    //use average recruitment


  Ftac=(tot_targ_Z-nat_mort)*rc_s(lyear);  //also the targetF*rc_S(lyear)for later use;  target total mortality and last yr selectivity
  Mvec=M(lyear);
  Zvec=Mvec+Ftac;
  Svec=mfexp(-1.0*Zvec);
  tacAge=elem_prod(elem_div(Ftac,Zvec),elem_prod(1.0-Svec,abund(lyear+2)));
  Total_ssb = sum(elem_prod(abund(lyear+2),elem_prod(wt_age,mat_age))*propFemale);

//return num of fish at same age 
FUNCTION dvector CountFreqAges(const ivector& age,int fage,int lage)
  int numAge=lage-fage+1;
  ivector uniqAge(fage,lage);
  uniqAge.fill_seqadd(fage,1); 
  dvector freq(fage,lage); 
  freq.initialize();    

  //count the number of ages which at the same age
  for(int i=age.indexmin();i<=age.indexmax();i++){
    for(int j=fage;j<=lage;j++){
      if(age(i)==uniqAge(j)){ // found match
        freq(j)+=1;
        break;
      }
    }
  }
  return freq;



REPORT_SECTION
  get_TAC();
  report<<sum(tacAge)<<endl;
  report<<objective_function_value::gmax <<endl;
  report<<abund(lyear+2)<<endl;  //estimated abundance at age for begining of the TAC setting year(the year after the lag year)
  report<<mre<<endl;  //median relative error between observed and predicted catch for the 20 yrs included in this assessment
  report<<mare<<endl;  //median absolute relative error between observed and predicted for the 20yrs included in this assessment
  report<<Total_ssb<<endl; //SSB estimation for the begining of the TAC setting year(the year after the lag year)
  /*
  report << "Relative errir" << endl;
  for (i=fyear;i<=lyear;i++)
  {
  report << re(i) << endl;
  }
  report << "Absolute Relative errir" << endl;
  for (i=fyear;i<=lyear;i++)
  {
  report << abs_re(i) << endl;
  }
  report << "Total Abundance" << endl;
  for (i=fyear;i<=lyear;i++)
  {
  report << N(i) << endl;
  }
  report << "Abundance at Age" << endl;
  report << n << endl;
  report << "Fishing Mortality" << endl;
  report << F << endl;
  report << "Natural Mortality" << endl;
  report << M << endl;
  report << "Total Mortality" << endl;
  report << Z << endl;
  report << "Observed Recreational Catch" << " " << "Predicted Recreational Catch" << endl;
  for (i=fyear;i<=lyear;i++)
      {
          report << new_rc_cch(i) << "                 " << N_rc(i) << endl; //change by liu
      }

  report << "Recreational Selectivity" << endl;
  report << rc_s<< endl;  //changed by liu
  report << "Recreational catch standard deviation" << endl;  
  report << mfexp(log_sigma_rc) << endl;
  report << "Recreational catchability" << endl;
  report << qrc << endl;
  report << "Recreational differences" << endl;
  report << rc_diff << endl;
  report << "F_TAC at Age" << endl;
  report << Ftac << endl;
  report << "TAC at Age" << endl;
  report << tacAge << endl;
  */



