double pc_rand_bi(int i,double *ratio_left,double *mean_left,double *mean_right,double *stdev_left, double *stdev_right,gsl_rng *rr)
{
double bimodal_rand = gsl_ran_flat(rr,0,1);		
double pc_rand;

if (i==0)	{
	if		(bimodal_rand < ratio_left[i])	{	pc_rand = mean_left[i] +  gsl_ran_gaussian(rr,stdev_left[i]);	}
	else									{	pc_rand = mean_right[i] + gsl_ran_gaussian(rr,stdev_right[i]);	}
}
else		{
	if		(bimodal_rand < ratio_left[i])	{	pc_rand = 0 - gsl_ran_exponential(rr,stdev_left[i]);			}
	else									{	pc_rand = 0 + gsl_ran_exponential(rr,stdev_left[i]);			}
}

return pc_rand;
}
// ---------------------------------------------------------------------------------------------------- //
double pc_prob_prop(int i,double *ratio_left,double *ratio_right,double *mean_left,double *mean_right,double *stdev_left,double *stdev_right,double jump_pc_rand)
{
double prob_prop;
if (i==0)		{
	prob_prop = ratio_left[0] * gsl_ran_gaussian_pdf(jump_pc_rand-mean_left[0], stdev_left[0])+
				ratio_right[0]* gsl_ran_gaussian_pdf(jump_pc_rand-mean_right[0],stdev_right[0]);
}
else			{
	if (jump_pc_rand<0)	{prob_prop = ratio_left[i]  * gsl_ran_exponential_pdf(-jump_pc_rand,stdev_left[i]);	}
	else				{prob_prop = ratio_right[i] * gsl_ran_exponential_pdf(jump_pc_rand,stdev_right[i]);	}
}
return prob_prop;
}
// --------------------------------------------------------------------------------------- //
double prior_dist(int i,double value)
{
double prob=0;
double mean[19];
double stdev[19];

mean[5]=-.42;	stdev[5]=.08;		//mean[5]=.38;	stdev[5]=.675;
mean[6]=-.545;	stdev[6]=.12;		//mean[6]=.285;	stdev[6]=.0725;

prob=gsl_ran_gaussian_pdf(value-mean[i],stdev[i]);

//printf("%d: log10value=%f\t prior=%f\n",i,value,prob);	//getc(stdin);

return prob;
}

// ----------------------------------------------------------------------------------------------------- //
double prior_func(void *Paramstuff)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int i;
double log_prob_prior=0;
const double delta=pow(10,-10);	
double prob_prior[NUM_PARS];

for (i=1;i<=50+DATA_SETS;i++)	{
	if (i==2||i==3||i==4||i==5||i==6||i==10||i==11||i==12||i==14||i==15||i==16||i==17||i==18||i==19)	{	// continuous uniform for prior
		if (i==5||i==6)	{
			//printf("value=%e log_value=%e\n",Params->PARS[i],log10(Params->PARS[i]));
			prob_prior[i]=prior_dist(i,log10(Params->PARS[i]));
			//printf("prior[%d]=%f\n",i,prob_prior[i]);
		}
		else	{		
			prob_prior[i]		= 1;	//gsl_ran_flat_pdf(PARAMS[0].PARS[i],prior_low[i],prior_high[i]);
			//printf("prior[%d]=%f\n",i,prob_prior[i]);
		}
		/*
		else if (i>20 && i<=30+DATA_SETS)	{	// S(0) uniform proposal/prior
			prob_prior[i]		= .1;	//gsl_ran_flat_pdf(PARAMS[0].PARS[i],prior_low[i],prior_high[i]);
		}
		else if (i>50 && i<=50+DATA_SETS)	{	// R(0) gaussian proposal/uniform prior
			prob_prior[i]		 = .1;	//gsl_ran_flat_pdf(PARAMS[0].PARS[i],prior_low[i],prior_high[i]);
		}
		else {
			prob_prior[i] = 1;
		}
		*/
		if (prob_prior[i]<delta)	{
			printf("ZERO PROBABILITY for PRIOR!! i=%d\t prob_prior=%e\n",i,prob_prior[i]);	getc(stdin);
		}
		log_prob_prior += log(prob_prior[i]);
		//printf("parm:%d\t value=%f\t prob_prior=%f\t log_prob=%f\t sum_log_prob=%f\n",
		//	i,Params->PARS[i],prob_prior[i],log(prob_prior[i]),log_prob_prior);
	}//getc(stdin);
}
//printf("prior=%f\n",log_prob_prior);
return log_prob_prior;
}
// ---------------------------------------------------------------------------------------------------- //
void translate_parms_gamma(void *Paramstuff,int size_parm_index,int *parm_index,double *pc_parms)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int i,j;
//printf(" PARMS:\t");
for (i=0;i<=(size_parm_index-1);i++)			{								// intrinsic parameters
	j=parm_index[i];
	Params->PARS[j]=pc_parms[i];
	//printf("translate::%d or %d=%4.3e\n",j,i,pc_parms[i]);
	//printf("Params->PARS[%d]=%4.3e\n",j,pc_parms[i]);
}

for (i=0;i<=(size_parm_index-1);i++)			{								// intrinsic parameters
	j=parm_index[i];
	//printf("Params->PARS[%d]=%4.3e\n",j,Params->PARS[j]);
}//getc(stdin);

}
// ---------------------------------------------------------------------------------------------------- //
void translate_parms(void *Paramstuff,int size_parm_index,int size_S0_index,int size_r0_index,int *parm_index,int *S0_index,int *r0_index,double *pc_parms)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int i,j;
//printf(" PARMS:\t");
for (i=0;i<=(size_parm_index-1);i++)			{								// intrinsic parameters
	j=parm_index[i];
	Params->PARS[j]=pc_parms[i];
	//printf("%d=%4.3e ",j,pc_parms[i]);
}
//printf("\n S(0):\t");
for (i=size_parm_index;i<=(size_parm_index+size_S0_index-1);i++)	{							// S(0)
	j=S0_index[i-size_parm_index];
	Params->PARS[j]=pc_parms[i];
	//printf("%d=%3.2e ",j-20,pc_parms[i]);
}
//printf("\n R(0):\t");
for (i=size_parm_index+size_S0_index;i<=(size_parm_index+size_S0_index+size_r0_index-1);i++){	// R(0)
	j=r0_index[i-(size_parm_index+size_S0_index)];
	Params->PARS[j]=pc_parms[i];
	//printf("%d=%3.2e ",j-50,pc_parms[i]);
	}
//printf("\n");	getc(stdin);
}
// ---------------------------------------------------------------------------------------------------- //
void testing_parms(void *Paramstuff)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int j;
printf("testing:\n");

j=2;	Params->PARS[2]=3.801894e-02;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);	
j=3;	Params->PARS[3]=4.011745e-06;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=4;	Params->PARS[4]=2.800667e+00;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=5;	Params->PARS[5]=3.981072e-01;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=6;	Params->PARS[6]=3.235937e-01;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=10;	Params->PARS[10]=7.943282e-02;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=11;	Params->PARS[11]=0.000000e+00;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=12;	Params->PARS[12]=0.000000e+00;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=14;	Params->PARS[14]=2.618589e+00;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=15;	Params->PARS[15]=2.310000e+00;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=16;	Params->PARS[16]=1.1e+00;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=17;	Params->PARS[17]=5.250000e-02;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=18;	Params->PARS[18]=2.400000e-01;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);
j=19;	Params->PARS[18]=2.400000e-01;	Params->MLE[j] = Params->PARS[j];	printf("Parm[%d]=%4.3e\n",j,Params->PARS[j]);

//getc(stdin);

// Params->PARS[11]=3.333333e-02;

/*
for (i=0;i<=(size_parm_index-1);i++)			{								// intrinsic parameters
	j=parm_index[i];
	printf("testing:Params->PARS[%d]=%4.3e\n",j,Params->PARS[j]);
	Params->MLE[j] = Params->PARS[j];				// necessary if called by parhood
}
*/


}
