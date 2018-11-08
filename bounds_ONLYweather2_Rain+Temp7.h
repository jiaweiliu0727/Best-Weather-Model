// ---------------------------------------------------------------------------------------------------------------- //
double bound(int i,int j)				// bounds on parameters for parhood line search
{
double low;
double high;

//if (i==2)		{	low = -6.1;			high = .1;	}		//nuBAR	(log10)
if (i==2)		{	low = -3.0;			high = 0.20;	}		//nuV	(log10)
else if (i==3)	{	low = -4.0;			high = -2.00;	}		//nuF
//else if (i==4)	{	low = 0.001;			high = 3.0;	}		//k
else if (i==4)	{	low = 0.1;			high = 4.0;	}		//k
else if (i==5)	{	low = -1.0;			high = 0.2;}		//muV	(log10)
else if (i==6)	{	low = -2.0;			high = 1.50;}		//muF	(log10)

else if (i==7)	{	low = 1.01;			high = 9.01;}		//CK//  beta = the size of accumulated rainfall window

else if (i==8)	{	low = 5;			high = 50;	}		//gam_stepsV
else if (i==9)	{	low = 10;			high = 200;	}		//gam_stepsF
else if (i==10)	{	low = -4.0;			high = -1.0;}		//ratio	(log10)
else if (i==11)	{	low = 0.000001;			high = 0.600;	}		//sdR	RESTING SPORES
else if (i==12)	{	low = 0.000001;			high = 0.350;	}		//sdC  CONIDIA

else if (i==13)	{	low = -4.0;			high = 1.00;	}		//gamma
//else if (i==13)	{	low = -4.0;			high = 0.5;	}		//CK// Average R(0) for all populations
//else if (i==13)	{	low = .50;			high = 3.5;	}		//CK// cover_C
//else if (i==13)	{	low = .50;			high = 4.5;	}		//CK// Now it is the COVERED CAGE scaling parameter
												//CK// theta = delay for when to start accumulating rainfall

else if (i==14)	{	low = -4.0;			high = 1.00;	}		//gamma

else if (i==15)	{	low = 1.0;			high = 14.1;}		//neonates_v
//else if (i==16)	{	low = 0.51;			high = 5.1;	}		//r_time
//else if (i==16)	{	low = 38.0;			high = 52.0;	}		//CK// STOP1
else if (i==16)	{	low = 350.0;			high = 600.0;	}		//CK// C end
//CK//  I'm changing param 16 from influencing when spores start blooming to when spores stop blooming
//else if (i==17)	{	low = .005;			high = 3.0;	}		//lambdaV
else if (i==17)	{	low = .5;			high = 10.0;	}		//Covered_C
else if (i==18)	{	low = .066667;			high = 0.25;	}		//lambdaF

//else if (i==19)	{	low = -2.0;			high = -0.4;	}		//CK//
else if (i==19)	{	low = 380.0;			high = 550.0;	}		//CK// stop2 R_end.  time when resting spores stop blooming

else if (i==20)	{	low = 4.0;			high = 18.50;	}		//CK// Covered_R
//else if (i==20)	{	low = 0.5;			high = 4.50;	}		//CK// CAGE scaling parameter
//else if (i==20)	{	low = 1E-5;			high = 0.1;	}		//CK// rain scaling parameter

//else if (i==21)	{	low = 0.01;			high = 20.0;	}		//CK// Rain scaling parameter func 3 and 4
else if (i==21)	{	low = -2.50;			high = 0.50;	}		//CK// Rain scaling parameter func 2
//else if (i==21)	{	low = 0.0000001;			high = 0.2;	}		//CK// temp scaling parameter

else if (i==22)	{	low = -3.0;			high = -1.0;	}		//CK// RH scaling parameter

else if (i==23)	{	low = -1.50;			high = -0.200;	}		//CK// temp scaling parameter for muF
//else if (i==23)	{	low = 0.0001;			high = 0.2;	}		//CK// stop2 R_end.  time when resting spores stop blooming

else if (i==24)	{	low = 0.80;			high = 10.00;	}		//CK// open_C

else if (i==25)	{	low = 3.1;			high = 15.50;	}		//CK// open_R

else if (i==26)	{	low = 4.0;			high = 7.0;	}		// rain_P2

else if (i==27)	{	low = 80.0;			high = 150.0;	}		//CK// R_start.  time when resting spores start blooming
//else if (i==27)	{	low = -4.0;			high = -0.7;	}		//CK// stop2 R_end.  time when resting spores stop blooming
//else if (i==26)	{	low = 1E-4;			high = 0.1;	}		//CK// scaler for TEMP affect of muF func 2
//else if (i==26)	{	low = 0.001;			high = 4.0;	}		//CK// scaler for TEMP affect of muF func 3 and 4

//else if (i==28)	{	low = 0.015;			high = 0.07;	}		//CK// size_S. effect of size on susceptibility.
else if (i==28)	{	low = 250.0;			high = 350.0;	}		//CK// size_S. effect of size on susceptibility.

else if (i==29)	{	low = -1.00;			high = 0.30;	}		//CK// rain_P3

else			{	low = 1;			high = 1;	}

if		(j==1)	return low;
else if (j==2)	return high;
else { printf("PROBLEM WITH BOUNDS ON PARAMETERS!!\n");	return 0;	}
}
// ------------------------------------------------------------------------------------ //
double patch_bound(int i,int j)			// bounds on S(0) and R(0) for parhood line search
{
double Slow;	double Shigh;	double Rlow;	double Rhigh;
double eps=-4.0;
double eps2=1.1;
// i=zero isnt used normally
if      (i==0)  {Slow = eps2;                Shigh = 2.1;					Rlow = eps;     Rhigh = .51;	}

else if (i==1)	{Slow = eps2;				Shigh = 10.1;					Rlow = -3.0;		Rhigh = 0.10;	}
else if (i==2)	{Slow = eps2;				Shigh = 10.1;					Rlow = -4.0;		Rhigh = 0.10;	}
else if (i==3)	{Slow = eps2;				Shigh = 10.1;					Rlow = -3.0;		Rhigh = 0.10;	}
else if (i==4)	{Slow = eps2;				Shigh = 10.1;					Rlow = -4.0;		Rhigh = -0.50;	}
else if (i==5)	{Slow = eps2;				Shigh = 10.1;					Rlow = -4.0;		Rhigh = -0.50;	}
else if (i==6)	{Slow = eps2;				Shigh = 10.1;					Rlow = -5.0;		Rhigh = -2.00;	}
else if (i==7)	{Slow = eps2;				Shigh = 10.1;					Rlow = -12.0;		Rhigh = -9.0;	}
else if (i==8)	{Slow = 0.5;				Shigh = 4.0;					Rlow = -12.0;		Rhigh = -9.0;	}

else if (i==9)	{Slow = 1.5;				Shigh = 4.0;					Rlow = eps;		Rhigh = 0.07;	}
else if (i==10)	{Slow = 1.0;				Shigh = 4.0;					Rlow = eps;		Rhigh = 0.03;	}
else if (i==11)	{Slow = 1.0;				Shigh = 4.0;					Rlow = eps;		Rhigh = 0.09;	}
else if (i==12)	{Slow = 1.0;				Shigh = 4.0;					Rlow = eps;		Rhigh = 0.08;	}
else if (i==13)	{Slow = 1.2;				Shigh = 4.0;					Rlow = eps;		Rhigh = 0.08;	}
else if (i==14)	{Slow = 1.5;				Shigh = 4.0;					Rlow = 0.0;		Rhigh = 1.0;	}
else if (i==15)	{Slow = 1.3;				Shigh = 4.0;					Rlow = 0.0;		Rhigh = 1.0;	}
else if (i==16)	{Slow = 1.3;				Shigh = 4.0;					Rlow = 0.0;		Rhigh = 1.0;	}

else			{Slow = 0;				Shigh = 0;				Rlow = 0;		Rhigh = 0;	}

if		(j==1)	return Slow;
else if (j==2)	return Shigh;
else if (j==3)	return Rlow;
else if (j==4)	return Rhigh;
else { printf("PROBLEM WITH BOUNDS ON S(0) or R(0)!!\n");	return 0;	}
}
// --------------------------------------------------------------------------------------- //
double prior_bound(int i,int j)
{
double low;
double high;

if (i==2)		{	low = 0.00001;		high = .2;	}		//nuBAR
else if (i==3)	{	low = 0.00001;			high = 0.01;	}		//nuF
else if (i==4)	{	low = 0.01;			high = .4;	}		//k
else if (i==5)	{	low = .1666;		high = .9;	}		//muV
else if (i==6)	{	low = .01;			high = 11;	}		//muF

else if (i==7)	{	low = 0.01;			high = 20.01;}		//CK//  beta = the size of accumulated rainfall window

else if (i==8)	{	low = 16;			high = 36;	}		//gam_stepsV
else if (i==9)	{	low = 12;			high = 36;	}		//gam_stepsF
else if (i==10)	{	low = .001;			high = .05;	}		//ratio
else if (i==11)	{	low = 0;			high = .7;	}		//sdv
else if (i==12)	{	low = 0;			high = .7;	}		//sdf

else if (i==13)	{	low = 0.01;			high = 4.01;	}		//CK// theta = delay for when to start accumulating rainfall

else if (i==14)	{	low = 1.9;			high = 3.5;	}		//gamma
else if (i==15)	{	low = 1;			high = 14;	}		//neonates_v
//else if (i==16)	{	low = 0;			high = 10;	}		//r_time
else if (i==16)	{	low = 0.1;			high = 20.0;	}		//r_time  //CK// changed!		//CK// to change!
else if (i==17)	{	low = .1;			high = 9.0;	}		//lambdaV
else if (i==18)	{	low = .1;			high = 9.0;	}		//lambdaF

else if (i==19)	{	low = .01;			high = 0.4;	}		//CK// dummy param for testing!!!

else if (i==20)	{	low = 0.1;			high = 5.4;	}		//CK// covered cage scaling parameter
//else if (i==20)	{	low = 0.00001;			high = 0.4;	}		//CK// rain scaling parameter

else if (i==21)	{	low = 0.00001;			high = 0.4;	}		//CK// temp scaling parameter

else if (i==22)	{	low = 0.00001;			high = 0.4;	}		//CK// RH scaling parameter

else if (i==23)	{	low = 0.0000001;			high = 0.5;	}		//CK// Average R(0) for all populations

else if (i==24)	{	low = 0.1;			high = 3.5;	}		//CK// open_R

else if (i==25)	{	low = 0.5;			high = 4.5;	}		//CK// cover_R

//else if (i==26)	{	low = 1E-6;			high = 0.5;	}		//CK// scaler for TEMP affect of muF func 2
else if (i==26)	{	low = 0.01;			high = 20.0;	}		//CK// scaler for TEMP affect of muF func 3 and 4
else if (i==27)	{	low = 1E-6;			high = 0.5;	}		//CK// scaler for TEMP affect of muF func 2

else if (i==28)	{	low = 0.02;			high = 0.05;	}		//CK// open_R

else if (i==29)	{	low = 0.02;			high = 0.05;	}		//CK// cover_R


else			{	low = 1;			high = 1;	}

if		(j==1)	return low;
else if (j==2)	return high;
else { printf("PROBLEM WITH BOUNDS ON PARAMETERS!!\n");	return 0;	}
}
// --------------------------------------------------------------------------------------- //
double prior_patch_bound(int i,int j)
{
double Slow;	double Shigh;	double Rlow;	double Rhigh;

if (i==1)		{Slow = .1;				Shigh = 1;			Rlow = 0;		Rhigh = 60; }
else if (i==2)	{Slow = .1;				Shigh = 1;			Rlow = 0;		Rhigh = 40; }
else if (i==3)	{Slow = .1;				Shigh = 1;			Rlow = 0;		Rhigh = 1;	}
else if (i==4)	{Slow = .1;				Shigh = 1;			Rlow = 0;		Rhigh = .25;	}
else if (i==5)	{Slow = .1;				Shigh = 1;			Rlow = 0;		Rhigh = .1;	}
else if (i==6)	{Slow = .1;				Shigh = 1;			Rlow = 0;		Rhigh = .04;}
else if (i==7)	{Slow = .7*(35148/30);	Shigh = 1.3*(35148/30);	Rlow = 0;		Rhigh = .06;}
else if (i==8)	{Slow = 9;				Shigh = 1500;			Rlow = 0;		Rhigh = .03;	}
else if (i==9)	{Slow = 10;				Shigh = 9966;			Rlow = 0;		Rhigh = .07;	}
else if (i==10)	{Slow = 100;			Shigh = 12900;			Rlow = 0;		Rhigh = .04;}
else if (i==11)	{Slow = 200;			Shigh = 19900;			Rlow = 0;		Rhigh = .15;	}
else if (i==12)	{Slow = 200;			Shigh = 49900;			Rlow = 0;		Rhigh = .1;}
else if (i==13)	{Slow = 50;				Shigh = 1290;			Rlow = 0;		Rhigh = .16;	}
else if (i==14)	{Slow = 10;				Shigh = 1990;			Rlow = 0;		Rhigh = 1;	}
else if (i==15)	{Slow = 100;			Shigh = 49400;			Rlow = 0;		Rhigh = 1;	}
else if (i==16)	{Slow = 100;			Shigh = 69440;			Rlow = 0;		Rhigh = 1;	}

else			{Slow = 0;				Shigh = 0;				Rlow = 0;		Rhigh = 0;	}

if		(j==1)	return Slow;
else if (j==2)	return Shigh;
else if (j==3)	return Rlow;
else if (j==4)	return Rhigh;
else { printf("PROBLEM WITH BOUNDS ON S(0) or R(0)!!\n");	return 0;	}
}
// --------------------------------------------------------------------------------------- //
double prop_dist(int i,int j)
{
double ulow=0, uhigh=0;
double mean=0, sdev=0;
double Slow=0, Shigh=0;
double Rmean=0,Rsdev=0;

//ulow=0;	uhigh=0;	mean=0;	sdev=0;	Slow=0;	Shigh=0; Rmean=0; Rsdev=0;

if (i==2)		{	ulow =.001;	uhigh=.3;	}
else if (i==4)	{	mean=.13;	sdev =.4;	}
else if (i==5)	{	ulow=.1;	uhigh=1;	}
else if (i==6)	{	ulow=.1;	uhigh=11;	}
else if (i==8)	{	ulow=15;	uhigh=40;	}
else if (i==9)	{	ulow=15;	uhigh=40;	}
else if (i==10)	{	mean=.01;	sdev =.8;	}
else if (i==14)	{	mean=2.6;	sdev =.1;	}
else if (i==15)	{	ulow=1;		uhigh=14;	}
else if (i==16)	{	ulow=0;		uhigh=10;	}

else if ((i>20 && i<=30+DATA_SETS) || (i>50 && i<=50+DATA_SETS))	{
	if (i==21||i==51)		{ Slow=100;		Shigh=6000;	Rmean=20;	Rsdev=.4;	}
	else if (i==22||i==52)	{ Slow=200;		Shigh=1500;	Rmean=20;	Rsdev=.4;	}
	else if (i==23||i==53)	{ Slow=20;		Shigh=700;	Rmean=0;	Rsdev=0;	}
	else if (i==24||i==54)	{ Slow=100;		Shigh=5000;	Rmean=.08;	Rsdev=.4;	}
	else if (i==25||i==55)	{ Slow=10;		Shigh=800;	Rmean=.04;	Rsdev=.3;	}
	else if (i==26||i==56)	{ Slow=10;		Shigh=1500;	Rmean=.05;	Rsdev=.3;	}
	else if (i==27||i==57)	{ Slow=10;		Shigh=1400;	Rmean=.02;	Rsdev=.4;	}
	else if (i==28||i==58)	{ Slow=100;		Shigh=9000;	Rmean=.02;	Rsdev=.4;	}
	else if (i==29||i==59)	{ Slow=100;		Shigh=9000;	Rmean=.03;	Rsdev=.35;	}
	else if (i==30||i==60)	{ Slow=100;		Shigh=2000;	Rmean=.015;	Rsdev=.4;	}
	else if (i==31||i==61)	{ Slow=100;		Shigh=2000;	Rmean=.06;	Rsdev=.5;	}
	else if (i==32||i==62)	{ Slow=100;		Shigh=5000;	Rmean=.05;	Rsdev=.4;	}
	else if (i==33||i==63)	{ Slow=100;		Shigh=1200;	Rmean=.05;	Rsdev=.4;	}
	else if (i==34||i==64)	{ Slow=10;		Shigh=2000;	Rmean=0;	Rsdev=0;	}
	else if (i==35||i==65)	{ Slow=100;		Shigh=5000;	Rmean=0;	Rsdev=0;	}
	else if (i==36||i==66)	{ Slow=100;		Shigh=7000;	Rmean=0;	Rsdev=0;	}

	else {		printf("called index that is not a data-set!!\n");}
}
else	{
	printf("CALLED PROPOSAL DISTRIBUTION FOR NON-EXISTING PARAMETER!!\n");
}

if		(j==1)	return ulow;
else if (j==2)	return uhigh;
else if	(j==3)	return mean;
else if	(j==4)	return sdev;
else if (j==5)	return Slow;
else if (j==6)	return Shigh;
else if (j==7)	return Rmean;
else if (j==8)	return Rsdev;
else {printf("PROBLEM WITH MEAN or SDEV in GAUSSIAN PARAMETERES!!\n");	return 0;	}
}

double r_end(int i)
{
double r_day=31+20;
//r_day=32+30+55;
double end_day;

//Lat and long for all 6 sites in current order (KF, RC1, UM1, RC2, RC3, CCR)
//42.363523, -85.348499
//44.463390, -84.604086
//45.483875, -84.680951
//44.465764, -84.595857
//44.465764, -84.595857
//45.188782, -84.22861

if      (i==0)      end_day = 0;
else if (i==1)		end_day = 20;
else if (i==2)		end_day = 14;
else if (i==3)		end_day = 14;
else if (i==4)		end_day = 12;
else if (i==5)		end_day = 7;
else if (i==6)		end_day = 9;
else if (i==7)		end_day = 17;
else if (i==8)		end_day = 11;
else if (i==9)		end_day = 4;
else if (i==10)		end_day = 5;
else if (i==11)		end_day = 4;
else if (i==12)		end_day = 10;
else if (i==13)		end_day = 25;
else if (i==14)		end_day = 12;
else if (i==15)		end_day = 12;
else if (i==16)		end_day = 12;
//R_end[17]=R_final-21;	R_end[18]=R_final-21;	R_end[19]=R_final-21;

else {
	printf("ERROR IN R END DATE!!!!\n");
	return 0;
}
return (r_day - end_day);
}


double fixed_parm(int i)
{
double value;

if (i==2)			value = 0.0001;			//nuV
else if (i==3)		value = 0.0001;			//nuF
else if (i==4)		value = 0.1;			//k
else if (i==5)		value = 0.2;			//muV
else if (i==6)		value = 0.3;			//muF
else if (i==7)		value = 5.01;		//CK//  beta = the size of accumulated rainfall window
else if (i==8)		value = 5;				//gam_stepsV
else if (i==9)		value = 5;				//gam_stepsF
else if (i==10)		value = 0.03;			//ratio
else if (i==11)		value = 0;				//sdv
else if (i==12)		value = 0;				//sdf
else if (i==13)		value=1.1;			//CK// theta = delay for when to start accumulating rainfall
else if (i==14)		value = 1.9;			//gamma
else if (i==15)		value = 8.3;			//neonates_v
//else if (i==16)		value = 2;				//r_time
else if (i==16)		value = 1.0;				//r_time   //CK// Changed!
else if (i==17)		value = 0.2;			//lambdaV
else if (i==18)		value = 0.2;			//lambdaF

else if (i==19)		value = 0.96;			//CK// Dummy param for checking!!

else if (i==20)		value = 0.01;		//CK// rain scaling parameter

else if (i==21)		value = 0.01;		//CK// temp scaling parameter

else if (i==22)		value = 0.01;		//CK// RH scaling parameter

else if (i==23)		value = 0.01;		//CK// average R(0)

else if (i==24)		value = 1.0;		//CK// RH scaling parameter

else if (i==25)		value = 1.0;		//CK// average R(0)

else if (i==26)		value = 30.0;		//CK// average R(0)

else if (i==27)		value = 0.3;		//CK// average R(0)

else if (i==28)		value = 0.035;		//CK// average R(0)

else if (i==29)		value = 0.035;		//CK// average R(0)

else { printf("PROBLEM WITH BOUNDS ON PARAMETERS!!\n");	return 0;	}

return value;
}


double fixed_patch(int i,int j)
{
double S_value=0;	double R_value=0;

if      (i==0)  {S_value = 100;				R_value = .20;	}
else if (i==1)	{S_value = 100;				R_value = .50;	}
else if (i==2)	{S_value = 300;				R_value = .60;	}
else if (i==3)	{S_value = 100;				R_value = 0;	}
else if (i==4)	{S_value = (1293.3/30);		R_value = .3;	}
else if (i==5)	{S_value = 200;				R_value = .1;	}
else if (i==6)	{S_value = 100;				R_value = .05;	}
else if (i==7)	{S_value = (35148/30);		R_value = .08;	}
else if (i==8)	{S_value = 100;				R_value = .05;	}
else if (i==9)	{S_value = 100;				R_value = .07;	}
else if (i==10)	{S_value = 100;				R_value = .08;	}
else if (i==11)	{S_value = 200;				R_value = .2;	}
else if (i==12)	{S_value = 200;				R_value = .1;	}
else if (i==13)	{S_value = 100;				R_value = .25;	}
else if (i==14)	{S_value = 100;				R_value = 0;	}
else if (i==15)	{S_value = 100;				R_value = 0;	}
else if (i==16)	{S_value = 100;				R_value = 0;	}

else { printf("i=%d\t j=%d\t PROBLEM WITH INDEX ON S(0) or R(0)!!\n",i,j);	getc(stdin);    return 0;	}

if		(j==1)	return S_value;
else if (j==3)	return R_value;
else { printf("PROBLEM WITH CHOOSING S(0) or R(0)!!\n");	return 0;	}

}
// ---------------------------------------------------------------------------------------------------------------- //
double start_pc(int i)
{
double pc_value=0;

	 if (i==0)	pc_value= 5.186043e+00;		else if (i==1)	pc_value= 1.601133e+00;
else if (i==2)	pc_value= 6.347635e-01;		else if (i==3)	pc_value= .373140e+00;
else if (i==4)	pc_value=-5.457952e-01;		else if (i==5)	pc_value= 2.441176e-01;
else if (i==6)	pc_value= 4.549583e-01;		else if (i==7)	pc_value= 4.288782e-01;
else if (i==8)	pc_value=-1.874400e-01;		else if (i==9)	pc_value= 2.059179e-01;
else if (i==10)	pc_value=-2.302962e-01;		else if (i==11)	pc_value=-2.235796e-01;
else if (i==12)	pc_value= 1.475449e-01;		else if (i==13)	pc_value=-1.769097e-01;
else if (i==14)	pc_value=-1.212636e-01;		else if (i==15)	pc_value=-1.272351e-01;
else if (i==16)	pc_value= 4.185948e-01;		else if (i==17)	pc_value=-1.060853e-01;
else if (i==18)	pc_value=-3.476973e-01;		else if (i==19)	pc_value= 9.517782e-01;
else if (i==20)	pc_value= 2.936443e-01;		else if (i==21)	pc_value= 1.013198e-01;
else if (i==22)	pc_value= 8.026273e-01;		else if (i==23)	pc_value=-1.248938e-01;
else if (i==24)	pc_value= 6.162489e-01;		else if (i==25)	pc_value= 7.884241e-01;
else if (i==26)	pc_value= 4.769538e-01;		else if (i==27)	pc_value=-8.913772e-01;
else if (i==28)	pc_value= 1.910217e-01;		else if (i==29)	pc_value= 5.366982e-01;
else if (i==30)	pc_value=-3.458366e-01;		else if (i==31)	pc_value=-2.682097e-01;
else if (i==32)	pc_value= 7.909581e-01;		else if (i==33)	pc_value=-3.420851e-01;
else if (i==34)	pc_value= 2.193495e-02;		else if (i==35)	pc_value= 1.685945e-01;

else printf("bad pc!!!\n");

return pc_value;

}

// ------------------------------------------------------------------------------------------- //
void load_fxd_parms(void *Paramstuff,int numb)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int i=1;
double parm_load[numb];
FILE *ftp_parm_load;
ftp_parm_load = fopen("/users/egoldwyn/chicago/gypsy/loadparm.dat","r");

if (ftp_parm_load==0)	{	printf("loadpc file open error\n");	getc(stdin);	}
while (!feof(ftp_parm_load))	{
	fscanf(ftp_parm_load,"%lf \n",&parm_load[i]);
	printf("parm_load(%d)=%e\n",i,parm_load[i]);
	Params->PARS[i]=parm_load[i];
	i++;
}//getc(stdin);*/
return;
}

// ------------------------------------------------------------------------------------------- //
void global_fixed_parms(void *Paramstuff)                //JL: Called in the main .c file, for fixing some parameter values. Some parameters will be further set in the main .c file.
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

Params->PARS[1] = 1;											// dummy parm to differentiate between profile hood and maxhood
//Params->PARS[3] = 1;											// scaled nuV and nuR
//Params->PARS[7] = 3.0;											//CK//  beta = the size of accumulated rainfall window
Params->PARS[9] = 50.0;											//CK//  beta = the size of accumulated rainfall window


Params->PARS[10] = 0.002;		//ck// ratio for virus infection.  Trying to set it and see what happens

// mu_R (wont matter)
Params->PARS[8] = 5;			Params->PARS[9] = 5;			// gamma steps
//Params->PARS[11]= 0;			Params->PARS[12]= 0;			// noise variance
//Params->PARS[13]= 1.01;											// delta (fixed)
//Params->PARS[24]=1293.3/30;		Params->PARS[27]=35148/30;		// constant S(0)

}
// ---------------------------------------------------------------------------------------------------------------- //
void parm_range_inc(void *Paramstuff,double parm_inc,double host_inc,double initR_inc,int num_adj_pars)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int i;

Params->parm_step[1]=1;
Params->parm_step[7]=1;

if (num_adj_pars>0)	{
	for (i=1;i<=num_adj_pars;i++)	{
		Params->parm_low[i]  = bound(i,1);	Params->parm_high[i] = bound(i,2);
		if (i==2||i==3||i==4||i==5||i==6||i==9||i==10||i==11||i==12||i==13||i==14||i==15||i==16||i==17||i==18||i==19||i==20||i==21||i==22||i==23||i==24||i==25||i==26||i==27||i==28||i==29)	{
			Params->parm_step[i]=(Params->parm_high[i]-Params->parm_low[i])/parm_inc;
		}
		//if (i==14)	{
		//	Params->parm_step[i]=(Params->parm_high[i]-Params->parm_low[i])/(parm_inc/2);
		//}
		//printf("parm:%d\t low=%f high=%f step=%f\n",i,Params->parm_low[i],Params->parm_high[i],Params->parm_step[i]);
	}	//getc(stdin);
}

// ----------------------------- host/initR increments and end of spore germination -------------------------------- //
for (i=0;i<=DATA_SETS;i++){
	//Params->parm_low[30+i]  = patch_bound(i,1);		Params->parm_high[30+i]= patch_bound(i,2);	// S(0)
	Params->parm_low[40+i]  = patch_bound(i,1);		Params->parm_high[40+i]= patch_bound(i,2);	//CK// Site-Specific muF
	Params->parm_low[50+i]  = patch_bound(i,3);		Params->parm_high[50+i]= patch_bound(i,4);	// R(0)
	//Params->parm_step[30+i] = (Params->parm_high[30+i]-Params->parm_low[30+i])/host_inc;
	Params->parm_step[40+i] = (Params->parm_high[40+i]-Params->parm_low[40+i])/host_inc;   //CK// muF
	Params->parm_step[50+i] = (Params->parm_high[50+i]-Params->parm_low[50+i])/initR_inc;
	Params->R_END[i]	= r_end(i);	                                                // fungus stops blooming
}
// -------------------------------------- fixed parameter values -------------------------------------------------- //
//Params->parm_step[30+4]=1.8;	Params->parm_step[30+7]= 1.9;
//Params->parm_step[50+3]=1.8;	Params->parm_step[50+14]=1.8;    Params->parm_step[50+15]=1.8;	Params->parm_step[50+16]=1.8;

}

// ---------------------------------------------------------------------------------------------------------------- //
void random_restart_parms(void *Paramstuff,gsl_rng *r_seed,int num_adj_pars,int test,int pro)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int i;
for (i=1;i<=num_adj_pars;i++)	{
	if ((i==2||i==3||i==4||i==5||i==6||i==7||i==9||i==10||i==11||i==12||i==13||i==14||i==15||i==16||i==17||i==18||i==19||i==20||i==21||i==22||i==23||i==24||i==25||i==26||i==27||i==28||i==29) && (i!=pro))	{
		if (i==2||i==3||i==6||i==13||i==14||i==21||i==22||i==23||i==26||i==29)	{	//CK// Logged fungus params ONLY
			Params->PARS[i] = pow(10,gsl_ran_flat(r_seed,Params->parm_low[i],Params->parm_high[i]));
		}
		else if (i==16||i==18||i==19||i==20||i==25||i==27||i==28)	{	//CK// Fungus params ONLY
			Params->PARS[i] = gsl_ran_flat(r_seed,Params->parm_low[i],Params->parm_high[i]);
		}
		else if (i==11||i==12)						{	// start with very low amounts of stochasticity
		//else if (i==12)						{	// start with very low amounts of stochasticity
			Params->PARS[i] = .01;
		}
		else if (i==9)						{	// start with very low amounts of stochasticity
			Params->PARS[i] = 50.0;
		}
		else if (i==7)						{	// start with very low amounts of stochasticity
			Params->PARS[i] = 10.0;
			//Params->PARS[i] = 5.0;
			//Params->PARS[i] = 3.0;

		}
		//else if (i==10)						{	//CK// Set ratio to experimentall measured value
		//	Params->PARS[i] = 0.002;
		//}
		else if (i==4||i==5||i==8||i==10||i==15||i==17||i==24)	{	//CK// set virus parameters to 0. NO VIRUS
			Params->PARS[i] = 0.0;
		}
		else {	printf("in bounds.h random_restart_parms: bad initial parameter!!!\n"); getc(stdin);	}

		if (test==99)	Params->PARS[i] = fixed_parm(i);				// for fixed init parms
		if (pro==1)		{	printf("initial parm(%d)=%e\n",i,Params->PARS[i]);	}
	}
	Params->MLE[i] = Params->PARS[i];
	//printf("init parm(%d)=%f\n",i,Params->PARS[i]);
}

}

