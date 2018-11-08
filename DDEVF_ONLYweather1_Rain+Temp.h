double DDEVF(void *Paramstuff,double *RandNumsPass,size_t dim,int pop,int maxy_t,double (*sim_results)[7])
//double DDEVF2(struct STRUCTURE *Params,int pop,double sim_results[1+Params->MAXT[pop]/7][5])
{
// DDEVF sets up model and calls 0DE_SOLVER, returns Params->sim_results
//printf("DDEVF: population:%d\n",pop);

STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;


int MAXT3=(Params->EXPDATA[pop][Params->MAXT2[pop]][2]+1)*7;

if(pop==6){MAXT3=77;}  ///can't figure out why plot 6 (UMBS 2012) is not working.  says it only has 2 weeks of data, but that's not true at all...

//printf("data set %d has %d days\n",pop,MAXT3);getc(stdin);

//printf("Final exp week: %d\n", MAXT3);		getc(stdin);

int DIM = Params->PARS[9]+3;         //m classes, where PARS[9] is m, plus 1 for S, 1 for C and 1 for cadavers
int m = Params->PARS[9];

double t=h;		double t_next=h;	double t_0=h;	int i;				// time loop and index
double epsilon = pow(10,-6);
double y_ode[DIM];
double rand_nuR[MAXT3];
double rand_nuF[MAXT3];

//double ave_R = Params->PARS[50+pop];
//double ave_R = Params->PARS[2];  //average R(0) for all the sites.  Fit like all other general params
double ave_R = 1.0;  //Set it to 1.0 to eliminate its effect in the differential equations. (weather-only model)
//double ave_R = Params->PARS[13];  //average R(0) for all the sites.  Fit like all other general params
double specific_muF = Params->PARS[6];   //general intercept for MAX TEMP decay function for conidia
double specific_nuF = Params->PARS[3];   //Site-specific infection rate for conidia
//double R_end = Params->PARS[16];   //CK//  Flat time for end of resting spore blooming
double rain_P = Params->PARS[21];  //fit param used to scale accumulating rain.
double rain_P2 = Params->PARS[26];  //fit param used to scale accumulating rain.
double rain_P3 = Params->PARS[29];  //fit param used to scale accumulating rain.
double RH_P = Params->PARS[22];   //CK//  Parameter for RH data
double temp_P = Params->PARS[23];   //CK//  Parameter for temperature data
double fourth_size=Params->PARS[28];	//CK// degree day when the bugs reach 4th instar

//double stop1	= Params->PARS[16];  //CK// param used for starting date
//double stop2	= Params->PARS[19];  //CK// param used for starting date
double DDstart	= Params->PARS[27];  //CK// param used for starting date
double DDstop	= Params->PARS[19];  //CK// param used for starting date

Params->size_C = 1.0;
double C_end=Params->PARS[16];	  //CK// fit param that turns off new conidia production once a specific size has been reached

double temp_now;  //CIK// used simplify decay functions
double total_rainfall;  //used to sum up rainfall
int rain_day;  //used to sum up rainfall
int beta;		//used to dictate how many days to go back when accumulating rain
int theta;    //used to determine lag period before calculating accumulated rainfall
//beta=Params->PARS[40+pop];
beta=Params->PARS[7];
//beta=1;
theta=1;
//theta=Params->PARS[13];
//printf("Beta: %d Theta: %d\n", beta, theta);		getc(stdin);


double DDtemp_now;  //CIK// used simplify decay functions


double nuF2;
double nuR2;
double DD10=0;    //accumulated degree days about 10 degrees C

//double specific_Rend = Params->PARS[40+pop];   //Site-specific decay rate for conidia


//printf("specific_nuF: %e\n", specific_muF);		getc(stdin);


// ----------------------------------- Generate Random Numbers -------------------------------------------- //
for (i=0;i<MAXT3;i++)	{   //JL: The stochasticity to change the transmission rates of conidia and resting spores vary every day.
	rand_nuR[i]=gsl_cdf_gaussian_Pinv(RandNumsPass[i],Params->PARS[11]); //2nd entry is stdev. (or variance)
	rand_nuF[i]=gsl_cdf_gaussian_Pinv(RandNumsPass[i+MAXT3],Params->PARS[12]);
	//rand_nuF[i]=gsl_cdf_gaussian_Pinv(RandNumsPass[i],Params->PARS[12]);
	//printf("i(%d)\t rand pass=%f\t var=%f\t rand parts: nuV=%e\t nuF=%e\n",i,RandNumsPass[i],Params->PARS[11],rand_nuV[i],rand_nuF[i]);
}//getc(stdin);
// ------------------------------------- Initial Conditions ---------------------------------------------- //
// THIS SHOULD BE A NEW FUNCTION (should be in main?? outside an extra loop??)//

//CK//  Changing initial conditions for fungus only model.  Don't need to take out virus infected neonates
//CK//  Assume that all hatched bugs are susceptible to fungus
//CK//  Also going to convert to bugs per meter^2 for units of S(0), rather than egg masses per 1/40th of a hectare

//Params->INITS[0] = Params->PARS[30+pop]*temp_S/(temp_S+temp_V);			// initS
//Params->INITS[1] = Params->PARS[30+pop]*temp_V/(temp_S+temp_V);			// initV
Params->INITS[0] = Params->PARS[30+pop];			// initS
Params->INITS[1] = 0.0;						// initV  SET TO 0 FOR FUNGUS ONLY MODEL!!
//Params->INITS[3] = Params->PARS[40+pop];								// initR
Params->INITS[3] = ave_R;												// initR  //CK// changed to use average R(0), not site-specific

// END OF NEW FUNCTION //
double initS = Params->INITS[0];			// initS
double initV = Params->INITS[1];			// initV
double initR = Params->INITS[3];			// initR

//printf("CHECKING S(0)!! base S(0)=%f\t temp_S=%f\t temp_V = %f\n",Params->PARS[30+pop],temp_S, temp_V);
//printf("Start DDVF pop:%d\t initR =%f\t initialS = %f\t initV = %f\n",pop,initR, initS, initV ); getc(stdin);

// ----------------------------------------- Fixed Parameters ---------------------------------------------- //
int gstepsV		= (int) Params->PARS[8];	int gstepsF	= (int) Params->PARS[9];
double ratio	= Params->PARS[10];
//double neo_v	= Params->PARS[15];		// latent period of neonates (days)
double neo_v	= 7.0;			// latent period of neonates (days) Not used here!

double R_end;   //CK//  Change value for function of latitude
double R_start;   //CK//  Change value for function of latitude

//Params->PARS[0] = 1/(Params->PARS[4]);		// THIS SHOULD BE IN PARHOOD.C AFTER PARAMS ARE SELECTED (CAREFUL FOR MCMC.C VERSION)

Params->PARS[0]=1.0;


// ------------------------------------- initialize model parameters --------------------------------------- //
int FlagWeek;	int FlagV=0;	int FlagR=0;	int FlagR_end=0;		// keep track of end of week

double ConiBefore;    //CK// Thing to store conidia the day before the feral collection
double RestBefore;    //CK// Thing to store resting spores the day before the feral collection

int day = 0; int week = 0;							// keeps track of day and week number

int line_ticker=0;   //CK// Ticker used to associate t in function with numbered days.
int line_ticker2;
int test_day;	//CK// used to find the line in the weather data that corresponds to the starting day of collections
int num_day =  Params->DATA[pop][0][4];  //CK// Starting day number

//printf("starting day number: %d\n", num_day);		getc(stdin);

while(test_day != num_day){
	test_day = Params->WDATA[pop][line_ticker][0];  //The line in the weather data that corresponds to the starting day of collections
	line_ticker++;}

line_ticker=line_ticker-1;
line_ticker2=line_ticker-1;

//printf("corresponding line in WDATA: %d\n", line_ticker);		getc(stdin);

double S,V,F,R;	double IV=0, IF=0;
double E_V[gstepsV+1]; double E_F[gstepsF+1];
int num_weeks=MAXT3/7;

// -----------------------------------//CK// calculating ending blooming times //CK//--------------------------------------- //

DD10=0.0;		R_start = 0.0;
test_day = line_ticker;

while(DD10 <= DDstart){

	DDtemp_now = Params->WDATA[pop][test_day][4]-10.0;  //CK// begin calculation of accumulated Degree Days
	if(DDtemp_now<0.0){DDtemp_now=0.0;}
	DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time
	R_start++;
	test_day++;
//printf("temp now: %f DD10: %f R_start: %f test_day: %d\n", DDtemp_now, DD10, R_start, test_day);		//getc(stdin);
}
//printf("temp now: %f DD10: %f R_start: %f test_day: %d\n", DDtemp_now, DD10, R_start, test_day);		//getc(stdin);


//R_start=0.0;

DD10=0.0;		R_end = 0.0;
test_day = line_ticker;

while(DD10 <= DDstop){

	DDtemp_now = Params->WDATA[pop][test_day][4]-10.0;  //CK// begin calculation of accumulated Degree Days
	if(DDtemp_now<0.0){DDtemp_now=0.0;}
	DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time
	R_end++;
	test_day++;
//printf("temp now: %f DD10: %f R_end: %f test_day: %d\n", DDtemp_now, DD10, R_end, test_day);		//getc(stdin);
}
//printf("temp now: %f DD10: %f R_end: %f test_day: %d\n", DDtemp_now, DD10, R_end, test_day);		//getc(stdin);

DD10=0.0;

/*
//Lat and long for all 6 sites in current order (KF, RC1, UM1, RC2, RC3, CCR)
//42.363523, -85.348499
//44.463390, -84.604086
//45.483875, -84.680951
//44.465764, -84.595857
//44.465764, -84.595857
//45.188782, -84.22861
//42.614460, -85.453187  <- YS 2011

*/
//printf("Checking blooming params: pop=%d\t Rstart= %f\t lat=%f\t r_time= %f\t Rend= %f\n", pop,Params->DAY_F[pop],lats[pop-1], r_time, R_end);		getc(stdin);



//-----------------------------------//CK// Infections on day 0!!!! -------------------------------//
// At week 0, when the bugs hatch, no infections can occur because the bugs are too young.
sim_results[0][0]=1.0;sim_results[0][1]=0.0;sim_results[0][2]=0.0;  //CK// Feral results for week 0.  No infections on day 0, apparently
sim_results[0][3]=0.0;sim_results[0][4]=0.0;  //CK// No conidia on day 0
//sim_results[0][5]=initR; sim_results[0][6]=initR;
//printf("R_start: %f\n", R_start);


if(R_start<1.0){   //JL: Usually R_start>1.0, which means the starting days of collection are in advance of resting spore blooming.


	total_rainfall = 0.0;

	//printf("stochasticity: %f\t initR: %f\n", rand_nuF[0], initR);

	for (rain_day= (line_ticker - beta - theta - 1);rain_day <= line_ticker - theta -1;rain_day++){
		total_rainfall = Params->WDATA[pop][rain_day][1]+total_rainfall;
	//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
	}

	//Params->nuR = (alpha*total_rainfall)*exp(rand_nuF[0]);
	//Params->nuR = exp(rain_P*total_rainfall + rand_nuF[0]);
	//Params->nuR = (total_rainfall/(rain_P + total_rainfall))*exp(rand_nuF[0]);	//CK// Resting Spore Rain Response #3
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2.0)*exp(rand_nuF[0]);	//CK// Resting Spore Rain Response #4

//Deterministic nuR!//
	Params->nuR=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(rand_nuR[0]);
	//Params->nuR = total_rainfall;
	//Params->nuR = exp(rain_P*total_rainfall);
	//Params->nuR = rain_P2*exp(rain_P*total_rainfall*exp(rand_nuR[0]));
	//Params->nuR = (total_rainfall/(rain_P + total_rainfall));	//CK// Resting Spore Rain Response #3
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2.0);	//CK// Resting Spore Rain Response #4
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))*(total_rainfall/(rain_P + total_rainfall)));	//CK// Resting Spore Rain Response #4
//Deterministic nuR!//

	//if(total_rainfall==0.0){Params->nuR=0.0;}

	//Params->nuR = (rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][6])*exp(rand_nuF[0]);
	//Params->nuR = exp(rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[0]);

	//Params->nuR = exp(rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][2] + RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[0]);

	//printf("initR: %f\n", initR);
	//printf("rain: %f\t stochasticity: %f\t initR: %f\n", total_rainfall, rand_nuF[0], initR);

	//sim_results[0][5]= Params->nuR*initR*(size_S*fourth_size);
	sim_results[0][5]= Params->nuR*initR;


	total_rainfall = 0.0;

//printf("Pop:%d\t Day:%d\t Week:%d\t Rain1: %f\t Stoch1: %f\t Spores1: %f\n", pop, day, week, total_rainfall, junk3, junk4);
//printf("Pop:%d\t Day:%d\t Week:%d\t Rain1: %f\t Stoch1: %f\t Spores1: %f\n", pop, day, week, total_rainfall, exp(rand_nuF[0]), initR);

	for (rain_day= (line_ticker2 - beta - theta - 1);rain_day <= line_ticker2 - theta -1;rain_day++){
		total_rainfall = Params->WDATA[pop][rain_day][1]+total_rainfall;
	//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
	}

	//Params->nuR = (alpha*total_rainfall)*exp(rand_nuF[0]);
	//Params->nuR = exp(rain_P*total_rainfall + rand_nuF[0]);
	//Params->nuR = (total_rainfall/(rain_P + total_rainfall))*exp(rand_nuF[0]);	//CK// Resting Spore Rain Response #3
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2)*exp(rand_nuF[0]);	//CK// Resting Spore Rain Response #4

//Deterministic nuR!//
	Params->nuR=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(rand_nuR[0]);
	//Params->nuR = total_rainfall;
	//Params->nuR = exp(rain_P*total_rainfall);
	//Params->nuR = rain_P2*exp(rain_P*total_rainfall*exp(rand_nuR[0]));
	//Params->nuR = (total_rainfall/(rain_P + total_rainfall));	//CK// Resting Spore Rain Response #3
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2.0);	//CK// Resting Spore Rain Response #4
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))*(total_rainfall/(rain_P + total_rainfall)));	//CK// Resting Spore Rain Response #4
//Deterministic nuR!//

	//if(total_rainfall==0.0){Params->nuR=0.0;}

	//Params->nuR = (rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][6])*exp(rand_nuF[0]);
	//Params->nuR = exp(rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[0]);

	//Params->nuR = exp(rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][2] + RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[0]);

	//sim_results[0][6]= Params->nuR*initR*(size_S*fourth_size);
	sim_results[0][6]= Params->nuR*initR;
}
else{sim_results[0][5]=0.0; sim_results[0][6]=0.0;}

//printf("Pop:%d\t Day:%d\t Week:%d\t Rain1: %f\t Stoch1: %f\t Spores1: %f\n", pop, day, week, total_rainfall, exp(rand_nuF[0]), initR); //getc(stdin);

//printf("Pop:%d\t Day:%d\t Week:%d\t Rain: %f\n", pop, day, week, total_rainfall);
//printf("Pop:%d\t Day:%d\t Week:%d\t Conidia1:%f\t Conidia2:%f\t Spores1: %f\t Spores2: %f\n", pop, day, week, sim_results[week][3],sim_results[week][4], sim_results[week][5], sim_results[week][6]); //getc(stdin);

//R=0.0;  //CK//  Make sure R density starts at zero.  Wasn't sure if that was happening....

// ---------------------------- Calculated Parameters  and Population Sizes ------------------------------- //
double Vstart = ratio*initV;						// viral cadavers after infected neonates die

//printf("Starting Virus Conditions: initV=%f\t ratio=%f\t Vstart=%f\n", initV, ratio, Vstart);	getc(stdin);

//double r_germ = Params->DAY_F[pop]-r_time;		// calculated day of resting spore germination
double r_germ = R_start;		//CK// nixing r_time because I made all germ dates start at beginning of the collections
if (r_germ<0)	r_germ=0;
if		(initS>Vstart)	{	S = initS-Vstart;	}	// S(0) (number of healthy neonates)
else					{	S = 0;				}
V=0.0; F=0.0; R=0.0;
//for (i=1;i<=gstepsV;i++)		{	E_V[i]=0;			}
for (i=1;i<=gstepsF;i++)		{	E_F[i]=0;			}   //Define the initial values (all 0) in the exposed classes.
//double timing[6]={neo_v,r_germ,R_end,MAXT3};
double timing[6]={r_germ,R_end,MAXT3};
// ----------------------------------------- initialize results ------------------------------------------- //
double P[4] = {0,0,0,0};				// model results (0,healthy,viral kill,fungal kill)
double Vkill=0, Fkill=0;				// total individuals infected (used in calculating fraction infected)
double VirFI[100], FungFI[100];			// fraction of infected individals (by week number)
VirFI[0] = 0.0; VirFI[1] = 0.0; FungFI[0] = 0.0; FungFI[1] = 0.0;		// fraction infected first two weeks

// -------------------- MAIN LOOP!! (calculate populations as time is increased) -------------------------- //
//printf("pop:%d\t t_0=%f\t neo_v=%f\t r_germ=%f R_end=%f week_end=%d maxt=%f\n",pop,t_0,timing[0],timing[1],timing[2],7*(week+1),timing[3]);
//printf("begin DDEVF main loop\n");		getc(stdin);
while (t_0<MAXT3+h)	{    //CK// change MAXT to MAXT2 to let it go to the end of the experimental datasets?
	if (week>=num_weeks)	break;
	FlagWeek=0;



	//printf("t_0=%f\t day=%d\t week=%d\t num_weeks=%d max_time=%d\n",t_0,day,week,num_weeks,Params->MAXT[pop]);//getc(stdin);
	//printf("neo=%f\t rgerm=%f\t rend=%f\t day+1=%d\n",timing[0],timing[1],timing[2],day+1);	getc(stdin);
	// ------------------------------- Find Stoppage Event ----------------------------------------------- //
	// --------------------- end of day -------------------- //
	if (day+1<timing[0] && day+1<timing[1])	{
		t_next=day+1;
		day++;
		num_day++;
		line_ticker++;

		//printf("end of day: t=%f\n",t_next);

		if (day%7==0)	{
		  //printf("end of week!!: t=%f\n",t_next);
			week++;
			FlagWeek=1;
		}
		//getc(stdin);
	}
	// --------------------- resting spores bloom -------------------- //
	else if (timing[0]<timing[1] && timing[0]<day+1)	{
		FlagR=1;
		t_next=timing[0];
		timing[0]=999.9;
		//printf("Resting spores bloom: t=%f\n",t_next);		//getc(stdin);
	}
	// --------------------- resting spores done -------------------- //
	else if (timing[1]<timing[0] && timing[1]<day+1)	{
		FlagR_end=1;
		t_next=timing[1];
		timing[0]=999.9;
		timing[1]=999.9;
		//printf("Resting spores done: t=%f\n",t_next);		//getc(stdin);
	}
	// --------------------- resting spores bloom and end of day -------------------- //
	else if (abs(day+1-timing[0])<epsilon)	{
		FlagR=1;
		t_next=day+1;
		timing[0]=999.9;
		day++;
		num_day++;
		line_ticker++;

		if (day%7==0)	{
			week++;
			FlagWeek=1;
		}
		//printf("end of day and Resting Spores bloom: t=%f\n",t_next);	//getc(stdin);
	}
	// --------------------- resting spores done and end of day -------------------- //
	else if (abs(day+1-timing[1])<epsilon)	{
		FlagR_end=1;
		t_next=day+1;
		timing[0]=999.9;
		timing[1]=999.9;
		day++;
		num_day++;
		line_ticker++;

		if (day%7==0)	{
			week++;
			FlagWeek=1;
		}
		//printf("end of day and Resting Spores done: t=%f\n",t_next);	//getc(stdin);
	}
	else {
		printf("ERROR: NO EVENT IS NEXT IN TIME???\n");
		getc(stdin);
	}
	// -------------------------- integrate until next stoppage event ---------------------------------- //
	while (t<t_next)	{
		y_ode[0]=S;	y_ode[m+1]=F;	Params->POPS[3]=R;
		//printf("pop:%d t=%f\t S=%4.3e\t V=%4.3e\t F=%4.3e\t R=%4.3e\n",pop,t,S,V,F,R);

		for (i=1;i<=gstepsF;i++)	{
			//y_ode[1+i]=E_V[i];
			y_ode[i]=E_F[i];     //Multiple exposed classes
		}
		y_ode[m+2]=Fkill;        //Infected cadavers
		//printf("DDEVF:\n parm 2=%f parm 3=%f parm 4=%f parm 5=%f parm 6=%f\n",Params->PARS[2],Params->PARS[3],Params->PARS[4],Params->PARS[5],Params->PARS[6]);
		//getc(stdin);

		//Params->nuF = specific_nuF;
		//Params->nuR = initR;
		//Params->nuV = Params->PARS[2];
		Params->nuV = 0.0;
		//Params->muF = specific_muF;
		Params->nuF = 0.0;
		nuF2 = 0.0;
		Params->muF = 0.0;

		DDtemp_now = Params->WDATA[pop][line_ticker - 1][4]-10.0;  //CK// begin calculation of accumulated Degree Days
		if(DDtemp_now<0.0){DDtemp_now=0.0;}
		DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time

		if(DD10>=C_end){Params->size_C=0.0;}	//CK// stops new conidia production once a fit size has been reached.

		//Params->nuF = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6]);
		//Params->nuF = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[(int)t]);
		//Params->nuF = (DD10/fourth_size)*specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[(int)t]);
		//nuF2 = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[(int)t]);
		//Params->nuF = (DD10/fourth_size)*specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] * exp(rand_nuF[(int)t]));
		//nuF2 = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] * exp(rand_nuF[(int)t]));


		temp_now = Params->WDATA[pop][line_ticker - 1][2];  //CK// putting max temp in smaller object

		//Params->muF = specific_muF;	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		//Params->muF = specific_muF*temp_now;	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		//Params->muF = specific_muF*exp(temp_P*temp_now);	//CK// Conidia Decay Response #2.2  BEST SO FAR!!

		total_rainfall = 0.0;
		rain_day= line_ticker - beta - 1;
//printf("Start rain calc! day: %d\n", line_ticker-1);		getc(stdin);
		for (rain_day= (line_ticker - beta - theta - 1);rain_day <= line_ticker - theta -1;rain_day++)	{
			total_rainfall = Params->WDATA[pop][rain_day][1]+total_rainfall;
//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
		}
/*
		if(Params->POPS[3] > 0.0){
			total_rainfall = Params->WDATA[pop][line_ticker - 2][1]+total_rainfall;
		}
		else{total_rainfall = 0.0;}
*/
		//Params->nuR = total_rainfall;
		//Params->nuR = 1.0-exp(rain_P*total_rainfall);
		//Params->nuR = exp(rain_P*total_rainfall);
		//Params->nuR = rain_P2*exp(rain_P*total_rainfall + rand_nuR[(int)t]);
		//Params->nuR = (DD10/fourth_size)*rain_P2*exp(rain_P*total_rainfall*exp(rand_nuR[(int)t]));
		//nuR2 = rain_P2*exp(rain_P*total_rainfall*exp(rand_nuR[(int)t]));

		//nuR2=rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall*exp(rand_nuR[(int)t]))) - rain_P/(rain_P2+1);
		//nuR2=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp((RH_P*Params->WDATA[pop][line_ticker - 1][6] - temp_P*temp_now) * exp(rand_nuR[(int)t]));  //Rain + RH - Temp
		//nuR2=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp((RH_P*Params->WDATA[pop][line_ticker - 1][6]) * exp(rand_nuR[(int)t]));  //Rain + RH
		nuR2=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(( - temp_P*temp_now) * exp(rand_nuR[(int)t]));
		//JL: Calculating the fungal transmission rate, related to rainfall and temperature only. (Best weather-only model)
		if(nuR2> pow(10.0,10.0)){nuR2= pow(10.0,10.0);}


		Params->nuR = (DD10/fourth_size)*nuR2;


		//Params->nuR = (DD10/fourth_size)*rain_P2*exp(rain_P*total_rainfall*exp(rand_nuR[(int)t]));
		//nuR2 = rain_P2*exp(rain_P*total_rainfall*exp(rand_nuR[(int)t]));

		//if(Params->POPS[3] == 0.0 || total_rainfall==0.0){Params->nuR = 0.0; nuR2 = 0.0;}
		if(Params->POPS[3] == 0.0){Params->nuR = 0.0; nuR2 = 0.0;}

		//if(total_rainfall==0.0){Params->nuR=0.0;}

//printf("day: %d\t line: %d\t accumulated rain: %lf\n", num_day, line_ticker-1, total_rainfall);		getc(stdin);


		//Params->nuF = Params->PARS[3]*exp(rand_nuF[(int)t])*total_rainfall*rain_P;
//		//Params->nuF = Params->PARS[3]*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[(int)t]);
		//Params->nuF = Params->PARS[3]*exp(rand_nuF[(int)t]);

		//printf("call fast_ode: \t t=%f\t day=%d\t t_index=%d\t nuV_determ=%f\t nuV=%f\t rand_nuV=%f\n",
		//	t,day,(int)t,Params->PARS[2],Params->nuV,rand_nuV[(int)t]);	//getc(stdin);
		//printf("call fast_ode: \t t=%f\t day=%d\t t_index=%d\t nuF_determ=%f\t nuF=%f\t rand_nuF=%f\n",
		//	t,day,(int)t,Params->PARS[3],Params->nuF,rand_nuF[(int)t]);
		//printf("nuF: %f\n", Params->nuF);
		//printf("resting spores: %e\t nuR: %f\n", Params->POPS[3], Params->nuR );
		//if(Params->nuR > 0.0000000000001){getc(stdin);}

		//CK// Add in calculation of accumulated rain here.
		//Need to link first day of collections with numbered day for weather.  Will allow back calculation here
		//include total rain in Params, which is passed to ODE_Solver.  Make part of PARS.  Need to find end of PARS...


		t=ODE_Solver(t,t_next,Params,y_ode);  //Numerical integration from t to t_next, y_ode holds group densities
		//Update the values for each class//
		S=y_ode[0];	F=y_ode[m+1];  //set group densities equal to ode output
		IF=0;

//printf("pop: %d day: %d conidia: %f\n", pop, day, y_ode[7]*Params->nuF);



//printf("Pop:%d\t Day:%d\t Host:%f\t Virus:%f\t Conidia:%f\t Resting Spores: %f\t nuR: %f\n", pop, day, y_ode[0], y_ode[1], y_ode[7], Params->POPS[3], Params->nuR);
//if(pop == 1 && y_ode[7] >0.0){getc(stdin);}
//if(Params->POPS[3] >0.0){getc(stdin);}
//getc(stdin);

	if ((day+1)%7==0)	{
		ConiBefore=y_ode[m+1]*nuF2;  //CK// Saving the conidia on the beginning day of the week. 0 here.
		RestBefore = Params->POPS[3]*nuR2;  //Saving the conidia on the beginning day of the week
		//ConiBefore=y_ode[m+1]*Params->nuF;  //CK// Saving the conidia 24 hours before feral collection
		//RestBefore = Params->POPS[3]*Params->nuR;
		//FlagConidia=2;
	}
		//getc(stdin);

		for (i=1;i<=gstepsF;i++)	{
			//E_V[i]=y_ode[1+i];
			E_F[i]=y_ode[i];
			//IV += E_V[i];
			IF += E_F[i];
			//printf("E_V(%d)=%e\t",i,E_V[i]);	//printf("IV=%f\n",IV);
		}//printf("\n");
		if(IF>initS)		IF=initS;
		Fkill=y_ode[m+2];			// killed populations
	}
	if (FlagV==1)			{	//printf("VIRAL INFECTED NEONATES DIE: t=%f\n",t_next);		//getc(stdin);
		V=Vstart;
		FlagV=2;
	}
	if (FlagR==1)			{	//printf("Resting spores bloom: t=%f\n",t_next);			//	getc(stdin);
		R=initR;
		FlagR=2;
	}

	else if (FlagR_end==1)	{	//printf("Resting spores done: t=%f\n",t_next);				// getc(stdin);
		R=0;
		FlagR_end=2;
	}

	// ---------------------------- end of the week updates ---------------------------------- //

	if (FlagWeek==1)	{
		//printf("t=%f S=%e IV=%e IF=%e\n",t,S,IV,IF);
		//week++;							// add one to the week if 7 days are complete
		if ((S+IF)==0)	{
			VirFI[week]=0;	FungFI[week]=0;
		}
		else					{
			//VirFI[week]  = IV/(S+IV+IF);		// fraction infected at the end of each week
			FungFI[week] = IF/(S+IF);       // fraction infected by fungus at the end of each week
		}
		P[1] = 1-(FungFI[week]);
		P[2] = 0;  			//CK//	formerly virus infected
		P[3] = FungFI[week];

		//printf("week=%d:\t VirFI=%f\t FungFI=%f\n",week,VirFI[week],FungFI[week]);
		sim_results[week][0]=P[1];	sim_results[week][1]=P[2];	sim_results[week][2]=P[3];
		//JL: Recording the fraction infected or not infected on the ending day of a week. Passing to hood_pops.h to calculate likelihoods from observational data.
		//printf("pop=%d week %d results: P1=%f P2=%f P3=%f\n",pop,week,sim_results[week][0],sim_results[week][1],sim_results[week][2]);	//getc(stdin);

		//STORING THE F AND R AT EACH WEEK
		//JL: Recording the fungus density and transmission rate on the beginning and ending day of a week. Passing to hood_pops.h to calculate likelihoods from experimental data.
		sim_results[week][3]=ConiBefore;  //CK// Saving the conidia on the beginning day of each week. 0 here.

		//sim_results[week][4]=y_ode[m+1]*Params->nuF;
		sim_results[week][4]=y_ode[m+1]*nuF2;   //Saving the conidia on the ending day of each week. 0 here.

		sim_results[week][5]=RestBefore;    //CK// Saving the resting spores (fungus) on the beginning day of each week. Actually only recording nuR2.
		//sim_results[week][6]=Params->POPS[3]*Params->nuR;
		sim_results[week][6]=Params->POPS[3]*nuR2;   //Saving the resting spores (fungus) on the ending day of each week. Actually only recording nuR2.
        //In Weather-only model, no difference between conidia and resting spore. So only having one spore type for the fungus.

//printf("Pop:%d\t Day:%d\t Week:%d\t Conidia1:%f\t Conidia2:%f\t Spores1: %f\t Spores2: %f\n", pop, day, week, sim_results[week][3],sim_results[week][4], sim_results[week][5], sim_results[week][6]); getc(stdin);
//if( y_ode[7] >0.0){getc(stdin);}
//if(Params->POPS[3] >0.0){getc(stdin);}
//getc(stdin);

	}
	t_0=t_next;
}

//printf("Resting spores day -1: %f\n", sim_results[0][5]);

//printf("end DDEVF\n");	getc(stdin);
return 0;
}
