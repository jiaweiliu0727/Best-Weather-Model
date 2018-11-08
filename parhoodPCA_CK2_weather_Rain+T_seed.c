#include "headWonly2.h"
gsl_rng *r;
//////////////////Begin Dot Product////////////////////////////////////////

float DotProduct (int Length, double *Holder, double *PCA)
{
  double answer = 0;

  int i;
  for(i=0;i<Length;i++)
    answer += PCA[i]*Holder[i];   //JL: Used in the reconstruction of parameter values from PCA

  return(answer);


}

//////////////////End dot product/////////////////////////////////////////



int main(int argc, char *argv[])
{
//int test = atoi(argv[1]);                  // test is 99 for checking, 0 otherwise
int test=0; //JL// Using atoi(argv[1]) the code does not run. Valgrind gives the information "Invalid read of size 1 in stdlib.h.GF changes this to 0.
//int test = 99;
//int test = 66; //CK// Second test mode.  Like the full program but less runs and less MISER calls

////////////
///NOTE:  Going to try to change nuF to a site-specific parameter and make R(0) a general parameter
////////// Just going to try making a new parameter for average R(0).  Make it number 23.
////////// Comment out param #3 (general nuF) but leave structure in place if we go back to it
////////// Pipe param #23 into R for all populations
////////// Make it fit nuF, rather than R(0), for each population.  Sounds simple enough....

STRUCTURE Params;

int pro = 1;//atoi(argv[1]);						// pro and argv[1] are the inputs (argv[i] is the i^th input)
printf("Profile Parameter is %d\n",pro);	fflush(stdout);
// ------------------------------------- Adustable accuracy vs. speed ------------------------------------------------ //
int num_runs	 = 20;
double parm_inc, host_inc, initR_inc;	//int inc_gamma_box= 1;

//if (pro==1)	{	parm_inc=200.0;		host_inc=100.0;	initR_inc=100.0;	}
if (pro==1)	{	parm_inc=32.0;		host_inc=10.0;	initR_inc=18.0;	}
//if (pro==1)	{	parm_inc=34.0;		host_inc=18.0;	initR_inc=20.0;	}
else		{	parm_inc=15.0;		host_inc=10.0;	initR_inc=10.0;	}

//if (test==66)	{	printf("for checking CK MODE2!!!\n");        num_runs=5;} //CK// New test mode!

if (test==66)	{	printf("for checking!!!\n");        num_runs=5;	parm_inc=16.0;	host_inc=8.0;	initR_inc=8.0;}

//if (test==66)	{	printf("for checking!!!\n");        num_runs=5;	parm_inc=6.0;	host_inc=4.0;	initR_inc=4.0;	}
printf("runs=%d\t incs: parm=%2.0f\t S_0=%2.0f\t R_0=%2.0f\n",num_runs,parm_inc,host_inc,initR_inc);


double initialR[8+1];   //Try to fit the resting spore density at each site as individual parameters.

//---------------//CK// Initial conditions for S from field observations //CK//-------------------------------------------------------//
double initialS[8]={1.0 ,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};	 //CK// my observed average egg density per circle at each site S(0)

// ------------------------------------------------------------------------------------------------------------------ //
int i=0; int j;int ii; int jj; int k;
int run;	            int changer;	    double index, tot_index;

int num_adj_pars=29;			// number of adjustable parameters

double inner_parm;
double inner_parm2;
int pop;
double log_pop;
//Params.th_id=0;
// -------------------------------------------- MISTER STUFF --------------------------------------------------------- //
inputdata(&Params);				// gets observation, experiment and weather data from inputdata.h
//int calls=2;					//CK// turned down to run without stochasticity
int calls=400;					// number of stochastic simulations for each parameter and IC set
if (test==99)	calls=2;
if (test==66)	calls=2; //CK// second test mode
//if (test==66)	calls=5; //CK// second test mode
size_t dim;
size_t dim2;
// --------------------------------------- Name for Output Files ----------------------------------------------------- //
char strFileName[99];					// from filenames.h
GetString(pro,0,strFileName,98);	printf("file names:\t   %s\n",strFileName);		fflush(stdout);		//getc(stdin);
FILE *fp_results;
// ---------------------------------------- Random Number Stuff ------------------------------------------------------ //
gsl_rng *r_seed;
r_seed=random_setup();
//printf("Random Seed: %f\n", r_seed); getc(stdin);
// -------------------- parameter high/low values and increments and fixed parameter values ------------------------- //
global_fixed_parms(&Params);  // gets Params.PARS[i] for fixed parameters from bounds.h
parm_range_inc(&Params,parm_inc,host_inc,initR_inc,num_adj_pars); // gets Params.parm_set,low,high,R_END from bounds.h. Some of them will be further fixed later in this file.
// ------------------------------------ Declare Likelihood Quanitites ----------------------------------------------- //
double pop_lhood, pop_lhood2, pop_err,post_hood;	// population lhood (and posterior lhood) calculated for each initS and initR
double pop_best_lhood;					// likelihood and error for best initS and initR
double total_lhood;						// sum of pop_best_lhood over all patches
double best_post_hood;	double best_lhood=0;		// best post_hood and lhood
double prior[num_adj_pars];


/////////////////////////////Begin reading in principal component analysis results//////////////////////

	int NumberOfParams=13;			// 13 parameters

	int Realizations=25000;        //JL: Number of realizations in the MCMC step

	double RandomNumber;

	double LogJumpToNew;
	double LogJumpToOld;
	double ProbOfAcceptance;

	double LogOldPosterior;
	double LogNewPosterior;
	double LogNewPrior;

	const gsl_rng_type * T;


	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r, time(NULL)*(int)getpid());   //JL: Include the process id to avoid taking the same random number at the beginning stage (multiple jobs submitted together at the same time).
	srand(time(NULL));

	double delta=0;

	//Some variables for recording the PCA results//
	double SDpca[NumberOfParams];
	double Coefficients[NumberOfParams*NumberOfParams];
	double Center[NumberOfParams];
	double Scale[NumberOfParams];
	double PC[NumberOfParams];
	double Old_PC[NumberOfParams];
	double Holder[NumberOfParams];
	double PCAparams[NumberOfParams];
	double Old_Params[NumberOfParams];

	//double AcceptedVect[NumberOfParams]={0.0};
	//double LoopVect[NumberOfParams]={0.0};

	int a;
	int b;
	int ticker;
	int ticker2;

	int Case;

	int Accepted=0;
	signed int LoopNumber=-1;

	int ParCnt2 = NumberOfParams;

	double SigmaInflation=1.20;  //CK// Dave had 4 as his sigma inflation factor

	double sigma[NumberOfParams];

	run=1;	changer=1;	best_post_hood=-10000000000;

/////////////////////////////Begin reading in principal component analysis results//////////////////////

	FILE *file;

	file=fopen("PCAsdCK2.txt", "r");        //Reading in the standard deviations of PCA.
	for (a=0; a<(NumberOfParams); a++)
	{
		fscanf(file, "%lf\n", &SDpca[a]);
		//printf("%lf\n", SDpca[a]);
	}
	fclose(file);


	file=fopen("PCArotationsCK2.txt", "r");  //Reading in the rotations (coefficients between PC's and parameters)
	for (a=0; a<(NumberOfParams*NumberOfParams); a++)
	{
		fscanf(file, "%lf\n", &Coefficients[a]);
		//printf("%lf\n", Coefficients[a]);
	}
	fclose(file);

	file=fopen("PCAscaleCK2.txt", "r");      //Reading in the scales to reconstruct parameter values from PC's.
	for (a=0; a<NumberOfParams; a++)
	{
		fscanf(file, "%lf\n", &Scale[a]);
		//printf("%lf\n", Scale[a]);
	}
	fclose(file);

	file=fopen("PCAcenterCK2.txt", "r");        //Reading in the centers of parameters. Used in the reconstruction of parameters.

	for (a=0; a<NumberOfParams; a++)
	{
		fscanf(file, "%lf\n", &Center[a]);
		//printf("%lf\n", Center[a]);
	}
	fclose(file);

/////////////////////////////DONE reading in principal component analysis results//////////////////////


//////////////////////////////Start MCMC/////////////////////////////////////////////////////////////////////////////

///Generate first set of PC values////

	for(a=0; a<NumberOfParams; a++){
		sigma[a]=SigmaInflation*SDpca[a];
		//printf("%lf\n", sigma[a]);
	}


	for (a=0; a<NumberOfParams; a++){
		PC[a]=gsl_ran_gaussian (r, sigma[a]);    //CK//  PC contains current set of PC values
		//printf("%lf\n", PC[a]);
	}
	//getc(stdin);

///Begin Loop that goes through and tweaks PCs one at a time////

//while (1==1) {     //CK// INFINITE LOOP START!!!!
while (LoopNumber<=Realizations) {     //CK// BOUND LOOP START!!!!

		LoopNumber=LoopNumber+1;
//if (LoopNumber % NumberOfParams == 1)
//		{
		Case=LoopNumber%NumberOfParams;					//Determines which PC to change

		for (a=0; a<NumberOfParams; a++)
		{
			Old_PC[a]=PC[a];							//Store old PC values:
		}

		PC[Case]=gsl_ran_gaussian (r, sigma[Case]);		//Draw 1 new PC

		for (a=0;a<NumberOfParams; a++)								//Back transform PC's into model parameters
		{
			for (b=0; b<NumberOfParams; b++){
				Holder[b]=Coefficients[a*NumberOfParams+b];        //Store the coefficients between a certain parameter and all the PC's
			}


			PCAparams[a]=exp(DotProduct(ParCnt2, Holder, PC)*Scale[a]+Center[a]);       //Reconstruct the parameter value after tweaking one PC
			Old_Params[a]=exp(DotProduct(ParCnt2, Holder, Old_PC)*Scale[a]+Center[a]);  //Reconstruct the parameter value before tweaking that PC
			//JL: Parameters logged in PCA. So there is an exp() step here in reconstruction


		    //printf("%f\n", PC[a]);
		    //printf("%f\n", exp(PCAparams[a]));
		    //printf("PCAparams: %e\n", PCAparams[a]);  //getc(stdin);

		}

		LogJumpToNew= -log(gsl_ran_gaussian_pdf(PC[Case], sigma[Case]));		//Metropolis sampling step: LATER YOU WILL USE THIS TO CORRECT FOR PROPOSAL
		LogJumpToOld= -log(gsl_ran_gaussian_pdf(Old_PC[Case], sigma[Case]));

		//printf("New: %f Old: %f\n", LogJumpToNew, LogJumpToOld);

		//---------------------LINE UP THE NEW PCAparams WITH THE CURRENT ORDER OF PARAMS --------------------------//

		ticker=0;

		for (k=0;k<=num_adj_pars;k++)	{

			if (k==3||k==4||k==5||k==6||k==8||k==10||k==12||k==15||k==16||k==17||k==22||k==24)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 0.0;
			}
			//else if (k==6)	{	//CK// set conidia parameters to 0. NO CONIDIA
			//	Params.PARS[k] = 0.0;
			//}
			else if (k==0||k==1||k==2)	{	//JL: Something from the old code. Needs to fix some parameters to certain values to keep the model structure similar as before.
				Params.PARS[k] = 1.0;
			}

			else if (k==7)	{	//JL: Something from the old code.
				Params.PARS[k] = 10.0;
			}
			else if (k==9)	{	//CK// Number of exposed groups
				Params.PARS[k] = 50.0;
			}

			else{Params.PARS[k] = PCAparams[ticker];        //JL: Put in the parameter values according to the order.
				ticker++;}

			//printf("NEW Params.PARS: %e\n",Params.PARS[k]);
		}
//getc(stdin);
//JL: In the testing period, setting upper bounds for some parameters to make the MCMC process go smoothly. Finally, without these bounds, the process can go without being stuck.
		//if(Params.PARS[11]>1.0){Params.PARS[11] = 1.0;}
		//if(Params.PARS[12]>1.0){Params.PARS[12] = 1.0;}
		//if(Params.PARS[13]>10.0){Params.PARS[13] = 10.0;}
		//if(Params.PARS[14]>10.0){Params.PARS[14] = 10.0;}
		//if(Params.PARS[21]>0.10){Params.PARS[21] = 0.10;}
		//if(Params.PARS[22]>0.1){Params.PARS[22] = 0.1;}
		//if(Params.PARS[23]>0.251){Params.PARS[23] = 0.251;}
		//if(Params.PARS[26]>0.251){Params.PARS[26] = 0.251;}


//		exit(0);
/*
		for (j=1;j<=DATA_SETS;j++)	{

			ticker2= NumberOfParams - DATA_SETS + j - 1;       //JL: Weather-only model, no need to record the fungus density

			initialR[j] = PCAparams[ticker2];
			//printf("ssNuF: %e\n",initialR[j]);
		}
*/
		//printf("New: %f Old: %f\n", LogJumpToNew, LogJumpToOld);


				total_lhood=0;

				// ----------------------- loop over patch numbers -------------------------------------------- //
				for (Params.pop=1;Params.pop<=DATA_SETS;Params.pop++)	{
					pop=Params.pop;
					pop_best_lhood = -1e9;

					int MAXT3=(Params.EXPDATA[pop][Params.MAXT2[pop]][2]+1)*7;
					if(pop==6){MAXT3=77;}  ///can't figure out why plot 6 (UMBS 2012) is not working.  says it only has 2 weeks of data, but that's not true at all...

					dim = 2*MAXT3;     //CK//  Changed dim to accomodate the longer EXP data sets... hopefully that works

					//dim = 2*Params.MAXT[pop];		// need to let this vary by patch (double because of nuV and nuF)
					gsl_monte_function G = { &Hood_Pops, dim, &Params };	// declares function calling Hood_Pops.h
					double xl[dim];	double xu[dim];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim;jj++)	{
						xl[jj]=0;
						xu[jj]=1;
					}

					Params.PARS[30+pop]=initialS[pop-1];  //JL: Space for recording initial host densities, starting from Params.PARS[31] for population 1.
					Params.PARS[50+pop]=initialR[pop];  //JL: Space for recording initial resting spore densities, starting from Params.PARS[51] for population 1.
					//Params.PARS[50+pop]=initialR[pop-1];  //CK// Needs to be [pop-1].  initialS is where the conditions are read in, so it starts at 0 and needs to be adjusted for

					//printf("S(0): %e\n",Params.PARS[30+pop]);
					//printf("ssNuF: %e\n",Params.PARS[50+pop]);

						// ----------- Use MISER to call function pop_lhood --------------------------------- //
						//JL: Use MISER to take multiple random values for each stochastic term and then average the likelihood values over these calls.
						gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim);
						gsl_monte_miser_integrate (&G,xl,xu,dim,calls,r_seed,s,&pop_lhood,&pop_err);
						gsl_monte_miser_free(s);
							//printf("pop:%d lhood=%f log_error=%f error=%f best_post_hood=%f\n",
								//pop,pop_lhood,log(pop_err),pop_err,pop_best_lhood);
							// ------------------- check to see if these ICs are better  ------------------------ //
						pop_lhood2=log(pop_lhood)-700.0;  //CK// Converting back to log likelihoods for MCMC
							//printf("pop_lhood:%e\t pop_lhood2:%f\n",pop_lhood,pop_lhood2); //getc(stdin);

							//printf("pop_lhood:%e\t pop_lhood2:%f\t pop_best_lhood=%f\n",pop_lhood,pop_lhood2,pop_best_lhood); getc(stdin);
						total_lhood += pop_lhood2;
					//printf("total_lhood=%f\n",total_lhood);	getc(stdin);
				} //CK// end of going through patches

				//prior[5]=log(prior_dist(5,log10(Params.PARS[5])));	//printf("prior[5]=%f\n",prior[5]);
				//prior[6]=log(prior_dist(6,log10(Params.PARS[6])));	//printf("prior[6]=%f\n",prior[6]);
				//post_hood=total_lhood+prior[5]+prior[6];
				LogNewPosterior = 0.0;
				LogNewPosterior=total_lhood;
				//printf("lhood=%f\t param(5)=%f\t prior(5)=%e\t params(6)=%f\t prior(6)=%e\n",LogNewPosterior,log10(Params.PARS[5]),prior[5],log10(Params.PARS[6]),prior[6]);
				//printf("\t parm:%d posthood=%f\t best_post_hood=%f\t prior5=%f\t prior6=%f\n",i,LogNewPosterior,best_post_hood,prior[5],prior[6]);	//getc(stdin);

		//---------------------Done calculating likelihood of NEW PARAM set --------------------------//

		//---------------------LINE UP THE OLD Params WITH THE CURRENT ORDER OF PARAMS --------------------------//

		ticker=0;

		for (k=0;k<=num_adj_pars;k++)	{

			if (k==3||k==4||k==5||k==6||k==8||k==10||k==12||k==15||k==16||k==17||k==22||k==24)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 0.0;
			}
			//else if (k==6)	{	//CK// set conidia parameters to 0. NO VIRUS
			//	Params.PARS[k] = 0.0;
			//}
			else if (k==0||k==1||k==2)	{	//JL: Something from the old code. Needs to fix some parameters to certain values to keep the model structure similar as before.
				Params.PARS[k] = 1.0;
			}

			else if (k==7)	{	//JL: Something from the old code.
				Params.PARS[k] = 10.0;
			}
			else if (k==9)	{	//CK// Number of exposed classes.
				Params.PARS[k] = 50.0;
			}
			else{Params.PARS[k] = Old_Params[ticker];                //CK// PUT IN OLD_PARAMS HERE!
				ticker++;}

			//printf("Params.PARS: %e\n",Params.PARS[k]);
		}
//JL: In the testing period, setting upper bounds for some parameters to make the MCMC process go smoothly. Finally, without these bounds, the process can go without being stuck.
		//if(Params.PARS[11]>1.0){Params.PARS[11] = 1.0;}
		//if(Params.PARS[12]>1.0){Params.PARS[12] = 1.0;}
		//if(Params.PARS[13]>10.0){Params.PARS[13] = 10.0;}
		//if(Params.PARS[14]>10.0){Params.PARS[14] = 10.0;}
		//if(Params.PARS[21]>0.10){Params.PARS[21] = 0.10;}
		//if(Params.PARS[22]>0.1){Params.PARS[22] = 0.1;}
		//if(Params.PARS[23]>0.251){Params.PARS[23] = 0.251;}
/*
		for (j=1;j<=DATA_SETS;j++)	{

			ticker2= NumberOfParams - DATA_SETS + j - 1;

			initialR[j] = Old_Params[ticker2];       //CK// PUT IN OLD_PARAMS HERE!
			//printf("ssNuF: %e\n",initialR[j]);
		}
*/
		//printf("New: %f Old: %f\n", LogJumpToNew, LogJumpToOld);


				total_lhood=0;

				// ----------------------- loop over patch numbers -------------------------------------------- //
				for (Params.pop=1;Params.pop<=DATA_SETS;Params.pop++)	{
					pop=Params.pop;
					pop_best_lhood = -1e9;

					int MAXT3=(Params.EXPDATA[pop][Params.MAXT2[pop]][2]+1)*7;
					if(pop==6){MAXT3=77;}  ///can't figure out why plot 6 (UMBS 2012) is not working.  says it only has 2 weeks of data, but that's not true at all...

					dim = 2*MAXT3;     //CK//  Changed dim to accomodate the longer EXP data sets... hopefully that works

					//dim = 2*Params.MAXT[pop];		// need to let this vary by patch (double because of nuV and nuF)
					gsl_monte_function G = { &Hood_Pops, dim, &Params };	// declares function calling Hood_Pops.h
					double xl[dim];	double xu[dim];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim;jj++)	{
						xl[jj]=0;
						xu[jj]=1;
					}

					Params.PARS[30+pop]=initialS[pop-1];  //JL: Space for recording initial host densities, starting from Params.PARS[31] for population 1.
					Params.PARS[50+pop]=initialR[pop];  //JL: Space for recording initial resting spore densities, starting from Params.PARS[51] for population 1.

					//printf("S(0): %e\n",Params.PARS[30+pop]);
					//printf("ssNuF: %e\n",Params.PARS[50+pop]);

						// ----------- Use MISER to call function pop_lhood --------------------------------- //
						//JL: Use MISER to take multiple random values for each stochastic term and then average the likelihood values over these calls.
						gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim);
						gsl_monte_miser_integrate (&G,xl,xu,dim,calls,r_seed,s,&pop_lhood,&pop_err);
						gsl_monte_miser_free(s);
							//printf("pop:%d lhood=%f log_error=%f error=%f best_post_hood=%f\n",
								//pop,pop_lhood,log(pop_err),pop_err,pop_best_lhood);
							// ------------------- check to see if these ICs are better  ------------------------ //
						pop_lhood2=log(pop_lhood)-700.0;  //CK// Converting back to log likelihoods for MCMC
							//printf("pop_lhood:%e\t pop_lhood2:%f\n",pop_lhood,pop_lhood2); //getc(stdin);

							//printf("pop_lhood:%e\t pop_lhood2:%f\t pop_best_lhood=%f\n",pop_lhood,pop_lhood2,pop_best_lhood); getc(stdin);
						total_lhood += pop_lhood2;
					//printf("total_lhood=%f\n",total_lhood);	getc(stdin);
				} //CK// end of going through patches

				//prior[5]=log(prior_dist(5,log10(Params.PARS[5])));	//printf("prior[5]=%f\n",prior[5]);
				//prior[6]=log(prior_dist(6,log10(Params.PARS[6])));	//printf("prior[6]=%f\n",prior[6]);
				//post_hood=total_lhood+prior[5]+prior[6];
				LogOldPosterior=0.0;
				LogOldPosterior=total_lhood;
				//printf("lhood=%f\t param(5)=%f\t prior(5)=%e\t params(6)=%f\t prior(6)=%e\n",total_lhood,log10(Params.PARS[5]),prior[5],log10(Params.PARS[6]),prior[6]);
				//printf("\t parm:%d posthood=%f\t best_post_hood=%f\t prior5=%f\t prior6=%f\n",i,post_hood,best_post_hood,prior[5],prior[6]);	//getc(stdin);

		//---------------------Done calculating likelihood of OLD PARAM set --------------------------//

		// ----------------------Compare New Param set with Old Param Set ----------------------------------------- //

		//ProbOfAcceptance=exp(-LogNewPosterior-LogJumpToOld + LogOldPosterior+LogJumpToNew);    //Probability of accepting the new PC
		ProbOfAcceptance=exp(LogNewPosterior+LogJumpToOld - LogOldPosterior-LogJumpToNew);    //Probability of accepting the new PC
//		printf("Reasons to stay=%.0f\n", LogNewPosterior+LogJumpToOld);
//		printf("Reasons to leave=%.0f\n",LogOldPosterior+LogJumpToNew);
//		printf("ProbOfAcceptance=%f\n", ProbOfAcceptance);
//The larger the value, the more likely to accept:  -LogNewPosterior	-LogJumpToOld 	+ LogOldPosterior	+LogJumpToNew


		//printf("LogOldPosterior: %f\t LogNewPosterior: %f\t ProbOfAcceptance: %f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance); getc(stdin);

		Params.LoopVect[Case] = Params.LoopVect[Case] + 1;

		if (ProbOfAcceptance>1 || gsl_rng_uniform_pos (r) < ProbOfAcceptance)   //MH-MCMC algorithm
		{
			LogOldPosterior=LogNewPosterior;     //Accept the new parameter values.
			//Accepted=Accepted+1;
			//printf("Accepted\n");
			Params.AcceptedVect[Case] = Params.AcceptedVect[Case] + 1;
		}

		else
		{
			for (a=0; a<NumberOfParams; a++)
			{
				PC[a]=Old_PC[a];             //Keep the old parameter values.
			}

			//printf("Rejected\n");

		}

		//getc(stdin);

		//exit(0);

		//printf("%f\t%f\t%f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance); getc(stdin);


		if (LoopNumber % 10 == 0)   	//CK// output results every 10 loops (probably not best plan but we'll see)
		{
		// ------------------------------------------ output results to file  --------------------------------------- //
			//Params.PARS[0] = Accepted;
			//Params.PARS[1] = LoopNumber;

			printf("%f\t%f\t%f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance);

			fp_results = fopen(strFileName,"a");
			output_file(&Params,fp_results,LogOldPosterior,num_adj_pars,pro); // prints to output file (filenames.h)
			fclose(fp_results);

			fflush(stdout);
		// ------------------------------------------------------------------------------------------------------- //
		}


}   //end of the infinite while loop


index=(Params.PARS[pro]-Params.parm_low[pro])/Params.parm_step[pro];
tot_index=(Params.parm_high[pro]-Params.parm_low[pro])/Params.parm_step[pro];
//printf("INDEX:%3.0f of %3.0f\t VALUE:%f\t LHOOD:%f\n",index+1,tot_index,Params.PARS[pro],best_post_hood);

free_i3tensor(Params.DATA,0,DATA_SETS,0,MAX_WEEKS,0,3);
free_i3tensor(Params.EXPDATA,0,DATA_SETS,0,MAX_WEEKS,0,3);
//free_i3tensor(Params.WDATA,0,DATA_SETS,0,MAX_WEEKS2,0,3);
printf("DONE!!!\n");

return 0;
}
