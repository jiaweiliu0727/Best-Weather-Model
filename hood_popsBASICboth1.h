double Hood_Pops(double *RandNumsPass,size_t dim,void *Paramstuff)
{
// calls DDEVF, uses results from DDEVF to calculate and return 'hood'
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int pop = Params->pop;
//printf("lhood_pops func! pop=%d\n",pop);
//printf("hood_pops nuV=%f\n",Params->PARS[2]);//getc(stdin);
int i;
int j; //CK//
int day; //CK//
int day2; //CK//

int MAXT3=(Params->EXPDATA[pop][Params->MAXT2[pop]][2]+1)*7;
int MAXT4=(Params->EXPDATA[pop][Params->MAXT2[pop]][2]+1);

if(pop==6){MAXT3=77;	MAXT4=11;}  ///can't figure out why plot 6 (UMBS 2012) is not working.  says it only has 2 weeks of data, but that's not true at all...

//printf("data set %d has %d days and %d weeks.  Test: %d\n",pop,MAXT3, MAXT4, Params->MAXT2[pop]);getc(stdin);


//printf("data set %d has %d feral weeks, %d exp weeks and %d exp points \n",pop,Params->MAXT[pop]/7,MAXT4, MAXT3);getc(stdin);

//for (i=0;i<=(Params->MAXT[pop]/7);i++)	{
//	printf("FERALS: i=%d\t pop:%d\t healthy:%d\t viral:%d\t fungal:%d\t week:%d\t week2:%d\n",i,pop,Params->DATA[pop][i][0],Params->DATA[pop][i][1],Params->DATA[pop][i][2], Params->DATA[pop][i][3], Params->DATA[pop][i][4]);
//}

//for (i=0;i<=(Params->MAXT2[pop]);i++)	{
//	printf("EXPERIMENTALS: i=%d\t wk_number:%d\t healthy:%d\t fungal:%d\t week:%d\t covered:%d\t date2:%d\n",pop,i,Params->EXPDATA[pop][i][0],Params->EXPDATA[pop][i][1],Params->EXPDATA[pop][i][2],Params->EXPDATA[pop][i][3],Params->EXPDATA[pop][i][4]);
//}

// ------------------------------------- likelihood stuff  ---------------------------------------------- //
double delta = 1.000000e-02;  //CK// Delta was previously fit using PARM[13] but it always came out to be 0.01.  Just set it to that
double gamma = Params->PARS[13];
double gamma2 = Params->PARS[14];

//double cover = 1.0;		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double cover_C = Params->PARS[17];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double cover_R = Params->PARS[20];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
//double cover = Params->PARS[20];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double open_C = Params->PARS[24];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double open_R = Params->PARS[25];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals

int    x1,x2,x3,n;			// data pop sizes (n is total caterpillars, x_i are number healthy, viral inf, and fung inf)
int    y1,y2,y3,n2,k2;			//CK// experimental data pop sizes (n2 is total caterpillars, y_i are number healthy, viral inf, and fung inf)
double p1,p2,p3;	        // simulation results (p_i is the probability of being in each of the three classes)
double r1, r2, r3, c1, c2;	            //CK// simulation results (r1 is resting spore density, c1 is condidia density)
//double a1,a2,a3;	        // calculated from size of pops in simulation results
double a,b;	        // calculated from size of pops in simulation results
double a2,b2;	        // calculated from size of pops in simulation results


double week_lhood = 0.0;		    // used to keep track of weekly likelihoods for feral data
double exp_lhood = 0.0;		    //CK// used to keep track of weekly likelihoods for experimental
double lhood = 0.0;  		// total likelihoods over the length of one epizootic (one simulation)
double lhood2 = 0.0;  		//CK// total likelihoods for experiments (one simulation)
double par, temp_par,prob, prob2, prob3;					//CK//

double sim_results[MAXT4][7];

DDEVF(Params,RandNumsPass,dim,pop,MAXT4,sim_results);
//DDEVF2(Params,pop,sim_results,th_id);
//getc(stdin);

//for (day2=0;day2<(MAXT4);day2++)	{
//	printf("%d\t %d\t %f\t %f\t %f\t %e\t %e\t %e\t %e\n",pop, day2, sim_results[day2][0], sim_results[day2][1], sim_results[day2][2], sim_results[day2][3], sim_results[day2][4], sim_results[day2][5], sim_results[day2][6]);
//}
//getc(stdin);

if(pop==6 || pop==7){lhood = 0;}
else{
//Calculating the likelihoods from observational data
//We are confident enough that the observation days are about to be the ending days of each week, so we need the fraction of uninfected from DDEVF at the end of each week.
// ---------------------------------------- Loop over all weeks ----------------------------------------------- //
for(i=2;i<=(Params->MAXT[pop]/7);i++)	{
	p1=sim_results[i][0];	p2=sim_results[i][1];	p3=sim_results[i][2];	// p1+p2+p3=1
	//JL: From DDEVF.h, p1 is the probability of being in S class (fraction uninfected); p3 is that of being in E class.(fraction infected). p2=0.
	//These probabilities are calculated from the model output. Then they will be fitted to the data.
    x1=Params->DATA[pop][i][0];		x2=Params->DATA[pop][i][1];		x3=Params->DATA[pop][i][2];
//    a1=p1*exp(gamma)+delta;	a2=p2*exp(gamma)+delta;	a3=p3*exp(gamma)+delta;
	n = x1+x2+x3;
	x1=n-x3;  //JL: Here we recalculate x1 to be the individuals not infected by the fungus from data. The ones infected by the virus are treated not infected by the fungus.

//p1=.733;

	if(p1 >= 1.0){p1 = 0.99999999;}
	if(p1 <= 0.0){p1 = 0.00000000001;}

	a = p1*exp(gamma) + delta;
	b = (1-p1)*exp(gamma) + delta;

	//printf("coefficients for beta binomial: a=%f b=%f\n",a,b);

//	if (a1<0||a2<0||a3<0)	{   printf("negative results!! error coming!! a1=%f a2=%f a3=%f\n",a1,a2,a3);	}
	if (a<0||b<0)	{   printf("negative results!! error coming!! a=%f b=%f\n",a,b);	}

	if (n==0) {week_lhood = 0.0;}
	else{
//	week_lhood = gsl_sf_lnfact(n)-gsl_sf_lnfact(x1)-gsl_sf_lnfact(x2)-gsl_sf_lnfact(x3)+
//		gsl_sf_lngamma(a1+a2+a3)-gsl_sf_lngamma(a1)-gsl_sf_lngamma(a2)-gsl_sf_lngamma(a3)+
//		gsl_sf_lngamma(a1+x1)+gsl_sf_lngamma(a2+x2)+gsl_sf_lngamma(a3+x3)-gsl_sf_lngamma(a1+x1+a2+x2+a3+x3);

	//week_lhood = (x1*log(p1) + (n-x1)*log(1-p1));

	week_lhood = gsl_sf_lnbeta(x1+a,n-x1+b) - gsl_sf_lnbeta(a,b) - gsl_sf_lnbeta(x1+1,n-x1+1) - log(n+1);  //CK// beta binomial lhood function
	//JL: Here is the logarithm of the pdf of a beta binomial distribution. The distribution describes the probability of observing x1 individuals not infected by the fungus.
	//JL: The probability of being infected by fungus follows a beta distribution with parameter a and b.
	}

	lhood +=week_lhood;

//	printf("%d\t %d\t %f\t %f\t %f\n",pop,i,p1,p2,p3); //getc(stdin);

	//printf("%d\t %d\t %d\t %d\t %d\t %f\t %f\t %f\n",pop,i,x1,x2,x3,p1,p2,p3);
	//printf("pop:%d\n data:\t n=%d x1=%d x2=%d x3=%d\n",pop,n,x1,x2,x3);	printf("results:\t p1=%f p2=%f p3=%f\n",p1,p2,p3);
	//printf("week_hood=%f\n",week_lhood);

	//printf("part1=%f part2=%f part3=%f part4=%f\n",gsl_sf_lnbeta(x1+a,n-x1+b), - gsl_sf_lnbeta(a,b), - gsl_sf_lnbeta(x1+1,n-x1+1), -log(n+1));
	//printf("part1=%f part2=%f\n",x1*log(p1), (n-x1)*log(1-p1));

	//printf("pop:%d week:%d\t gamma=%f\t temp_hood=%f\t gamma_hood=%f\n",pop,i,gamma,week_lhood,lhood);
 }//getc(stdin);
}
	//printf("pop:%d gamma=%f\t gamma_hood=%f\n",pop,gamma,lhood); getc(stdin);

	//CK// Adding in the likelihoods from the experimental data here
	//JL: In the calculation of likelihood from experimental data, we assume that transmission rate of fungus remains the same for each cage. No fungus density effect.
	//JL: The cages are put in the field for only 24 hours on one day near the ending of a week, and then taken back to the lab, and raise each individual separately until death.
	//JL: Thus, we can directly solve the differential equation to get S(t)/S(0), which is exactly the fraction of uninfected.
	//JL: We estimate the weekly "force" of fungal infection by averaging the products of density and transmission rate on the second last and the ending day of that week.


for(j = 0; j<=(Params->MAXT2[pop]); j++){    //loop calculating the likelihood for each cage

	day=Params->EXPDATA[pop][j][2];   //CK// pull out the week number for that cage

	c1=sim_results[day][3];	c2=sim_results[day][4];			//CK//  conidia and resting spores for the second last and ending of that week from DDEVF (from model output)
	r1=sim_results[day][5];	r2=sim_results[day][6];

	y1=Params->EXPDATA[pop][j][0];		y2=Params->EXPDATA[pop][j][1];		//CK//  S and F(infected by fungus) for that cage
	n2=y1+y2;  //summing total bugs in cage
	y1=n2-y2;  //summing total bugs in cage

	if(Params->EXPDATA[pop][j][3] == 1){   //Covered cages

		//par = (cover+cage)*((r1+r2)/2);
		//par = cover*((r1+r2)/2);
		par = cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2);
		//par = ((c1+c2)/2) + ((r1+r2)/2);


		//printf("%d\t %d\n",pop,day);
		//printf("COVER=%lf\t r1=%lf\t r2=%lf\t par=%lf\n",cover, r1, r2, par);  //getc(stdin);
	}

	else{  //Open cages
		par = open_C*((c1+c2)/2) + open_R*((r1+r2)/2);
		//par = ((c1+c2)/2) + ((r1+r2)/2);

		//printf("%d\t %d\n",pop,day);
		//printf("c1=%lf\t c2=%lf\t r1=%lf\t r2=%lf\t par=%lf\n", c1, c2, r1, r2, par); //getc(stdin);
	}

		if(par <= 0.0){par = 0.0000000000001;}

		prob = exp(-par);

		if(prob >= 1.0){prob = 0.99999999;}
		if(prob <= 0.0){prob = 0.00000000001;}

		//printf("prob=%f\n",prob);

		//exp_lhood = (y1*log(prob) + (n2-y1)*log(1-prob));

		a2 = prob*exp(gamma2) + delta;
		b2 = (1-prob)*exp(gamma2) + delta;


		exp_lhood = gsl_sf_lnbeta(y1+a2,n2-y1+b2) - gsl_sf_lnbeta(a2,b2) - gsl_sf_lnbeta(y1+1,n2-y1+1) - log(n2+1);  //CK// beta binomial lhood function
		//JL: Here is the logarithm of the pdf of a beta binomial distribution. The distribution describes the probability of observing y1 individuals not infected by the fungus.
	    //JL: The probability of being infected by fungus follows a beta distribution with parameter a2 and b2.

		//if(y2==0){exp_lhood=0.0;}    //CK// Condition to skip data where no infection occured.

		lhood2 += exp_lhood;

		//printf("c1=%lf\t c2=%lf\t r1=%lf\t r2=%lf\n",c1, c2, r1, r2);
		//printf("exp_hood=%f\n",exp_lhood);  //getc(stdin);


	if(Params->EXPDATA[pop][j][3] == 0){

		//prob2 = ((open_C*((c1+c2)/2))/(open_C*((c1+c2)/2) + open_R*((r1+r2)/2)))*(n2*(1 - prob))/n2;  //fraction of F infection due to conidia
		//prob3 = ((open_R*((r1+r2)/2))/(open_C*((c1+c2)/2) + open_R*((r1+r2)/2)))*(n2*(1 - prob))/n2;  //fraction of F infection due to resting spores
		prob2 = ((((c1+c2)/2))/(((c1+c2)/2) + ((r1+r2)/2)))*(n2*(1 - prob))/n2;  //fraction of F infection due to conidia
		prob3 = ((((r1+r2)/2))/(((c1+c2)/2) + ((r1+r2)/2)))*(n2*(1 - prob))/n2;  //fraction of F infection due to resting spores


		temp_par=((c1+c2)/2 + (r1+r2)/2);

		if(temp_par == 0.0){ prob2 =0.0; prob3 = 0.0;}

		//printf("%d\t %d\t %d\t %d\t %f\n",pop,day,n2,y1,prob); //getc(stdin);
//		printf("%d\t %d\t %f\t %f\t %f\t %f\t\n",pop,day,prob,1.0-prob, prob2, prob3); // getc(stdin);
	}


	if(Params->EXPDATA[pop][j][3] == 1){

		prob2 = 0.0;  //fraction of F infection due to conidia
		//prob3 = (((r1+r2)/2)/((r1+r2)/2))*(n2*(1 - prob))/n2;  //fraction of F infection due to resting spores
		prob3 = 1-prob;  //fraction of F infection due to resting spores\

		//prob2 = ((cover_C*((c1+c2)/2))/(cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2)))*(n2*(1 - prob))/n2;  //fraction of F infection due to conidia
		//prob3 = ((cover_R*((r1+r2)/2))/(cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2)))*(n2*(1 - prob))/n2;  //fraction of F infection due to resting spores
		//prob2 = ((((c1+c2)/2))/(((c1+c2)/2) + ((r1+r2)/2)))*(n2*(1 - prob))/n2;  //fraction of F infection due to conidia
		//prob3 = ((((r1+r2)/2))/(((c1+c2)/2) + ((r1+r2)/2)))*(n2*(1 - prob))/n2;  //fraction of F infection due to resting spores


		temp_par=((r1+r2)/2);
		//temp_par=((c1+c2)/2 + (r1+r2)/2);


		if(temp_par == 0.0){ prob2 =0.0; prob3 = 0.0;}

		//printf("%d\t %d\t %d\t %d\t %f\n",pop,day,n2,y1,prob); //getc(stdin);
//		printf("%d\t %d\t %f\t %f\t %f\t %f\t\n",pop,day,prob,1.0-prob, prob2, prob3); //getc(stdin);
	}

		//printf("%d\t %d\t %d\t %d\t %f\n",pop,day,n2,y1,prob); //getc(stdin);
		//printf("c1=%lf\t c2=%lf\t r1=%lf\t r2=%lf\n",c1, c2, r1, r2);
		//printf("pop:%d\t week:%d\t N:%d\t K:%d\t resting spores:%f\t prob:%f\t totalLike:%f\n",pop,day,n2,y1,(r1+r2)/2,prob,lhood2); getc(stdin);
		//if(c1 || r1 >0.0){getc(stdin);}
		//}
}

	//CK// Adding in the likelihoods from the experimental data here


//lhood2=-log(lhood2);   //CK// not sure if this is right, but it fixs the difference between our likelihood values

//printf("hood1=%f\t hood2=%f\n",lhood, lhood2);

//printf("pop:%d gamma=%f\t hood1=%f\t hood2=%f\n",pop,gamma,lhood, lhood2);

lhood += lhood2 + 700;
//lhood += lhood2 + 550;
//lhood += lhood2 + 1650;


lhood = exp(lhood);

if(isnan(lhood)){printf("ERROR NaN likelihood!!!"); getc(stdin);}
if(isinf(lhood)){printf("ERROR infinite likelihood!!!"); getc(stdin);}

//printf("hood=%e\n",lhood); getc(stdin);

//if(lhood > 0.0){getc(stdin);}
//getc(stdin);

return lhood;
}





double temp_hood(void *Paramstuff,double *jump_pc_parms,int pop)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

//printf("PARM(%d)=%e\n",2,jump_pc_parms[2]);	getc(stdin);

double var[]= {
 1.423383e-03, 6.116464e-03, 1.073083e+02, 8.457159e+03, 1.995897e-04,
 5.296273e-03, 4.286525e+00, 3.373528e-01, 4.043294e-04, 3.337862e-02,
 2.187840e+06, 7.753218e+05, 1.226094e+03, 6.979874e+06, 3.044502e+04,
 1.373654e+04, 1.816982e+06, 6.749283e+06, 2.903288e+07, 3.681558e+08,
 8.906700e+04, 6.779074e+05, 2.485600e+10, 1.416496e+10, 6.709227e-03,
 7.405556e-03, 1.338063e-03, 3.749663e-05, 2.734404e-05, 2.337219e-05,
 3.013087e-05, 1.811350e-05, 7.151468e-06, 6.439309e-05, 1.157808e-04,
 3.554809e-04	};

double hood=0;
int i;
int j;
for(i=0; i<100; i++)
{
	for(j=0; j<36; j++)
	{
		//hood += log(gsl_ran_gaussian_pdf(Params->test_data[i][j]-jump_pc_parms[j], 1));
		hood -= ((Params->test_data[i][j]-jump_pc_parms[j])*(Params->test_data[i][j]-jump_pc_parms[j])/(2*var[j]));
		//printf("hood=%e\n",hood);	getc(stdin);
		//printf("parms(%d)=%e data=%e\n",j,jump_pc_parms[j],Params->test_data[i][j]);
		//printf("marginalChange %d %d =%e\n",i,j, (Params->test_data[i][j]-jump_pc_parms[j])*(Params->test_data[i][j]-jump_pc_parms[j])/(2*var[j]));
	}
}


//printf("hood=%e\n",hood);

return hood;
}

