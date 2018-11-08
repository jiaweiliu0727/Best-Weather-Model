// --------------------------------Begin ODE system of White model --------------------------------------------//
int fast_odes(double t, const double y[], double dydt[],void *Paramstuff)
{
//struct STRUCTURE *Params=(struct STRUCTURE *)Paramstuff;
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int m			= Params->PARS[9];  //CK//  the number of exposed classes.  adjustable and fit.  yay!

double C		= y[m+1];
//double C		= Params->PARS[0];		//CK// Something majorly wrong here!!  Think Eli had the conidia turned off or something
//Need to fix because it is a constant right now	(set to Params->PARS[0]).  Going to try y[7] bc that should be C density...
//double nuV		= Params->nuV;		//Params->PARS[2];
double nuV		= 0.0;		//CK//  FOR FUNGUS ONLY MODEL. NO VIRUS INFECTION
double nuF		= Params->nuF;		//Params->PARS[3];
//double nuR		= 1.0;			
double nuR		= Params->nuR;	//param for nuR after it includes rainfall
//double muV		= Params->PARS[5];	//muF=Params->PARS[6];	//muR=PARAMS[0].PARS[7];// breakdown rate (per day)	
double muV		= 0.0;		//CK// FUNGUS ONLY MODEL.  MAKE SURE DECAY IS ZERO SO NEGATIVE VIRUS DOESN'T HAPPEN!!
double lambdaV	= Params->PARS[17];	double lambdaF=Params->PARS[18];		// decay rate of E_V and E_F				

double size_C		= Params->size_C;  //Scaling effect of size on susceptibility over time.

double muF		= Params->muF;  //CK//  Using the site-specific muF that was loaded up in DDEVF14

int i;

double S0 = Params->INITS[0];
double R  = Params->POPS[3];

double Finf = y[0]*nuF*C;
double Rinf = y[0]*nuF*R;

//printf("fast_odes:\t t=%f\t S=%f nuF=%e nuR=%f R=%f C=%e Cinf=%e Rinf=%e\n",t,y[0],nuF,nuR,R,C, Finf, Rinf);	getc(stdin);
//if(R>0.0){getc(stdin);}

//printf("fast_odes:\t t=%f\t S=%f nuF=%f nuR=%f R=%f\n",t, y[0],nuF,nuR,R);	//getc(stdin);

//printf("fast_odes:\t t=%f\t nuV_determ=%e nuV_stoch=%e\n",t,Params->PARS[2],nuV);	

//if(R>0.0){getc(stdin);}

//double hetero = y[1]*nuV*pow((y[0]/S0),(C));	// save time by doing once (heterogeneity term)
// ------------------------------------------ ODEs -------------------------------------------- //
if		(y[0]<.000001)	dydt[0]=0;
else				

dydt[0]  = -y[0]*(nuF*y[m+1] + nuR*R);	
dydt[1]  = nuF*y[m+1]*y[0] + nuR*R*y[0] - m*lambdaF*y[1];

for(i=2; i <= m; i++){
	dydt[i]=m*lambdaF*(y[i-1] -y[i]);
}

dydt[m+1]  = m*lambdaF*y[m]*size_C - muF*y[m+1];  //Conidia class!  Transission from final exposed class (m) to conidia class (m+1)
dydt[m+2] = m*lambdaF*y[m];

//printf("y(0)=%e hetero=%e muV=%e y(1)=%e\n",y[0],hetero,lambdaV,y[2]);
//printf("t=%f\t y[1]=%e nuV=%e y[0]=%e S0=%e C=%e\n",t,y[1],nuV,y[0],S0,C);
	   
/*for (i=0;i<15;i++)	{
	printf("dydt(%d)=%e\n",i,dydt[i]);
}
*/
return GSL_SUCCESS;
}

////////////////////Begin Jacobian of White model///////////////////////////
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *Paramstuff)
{
	printf("here\n");	getc(stdin);
	int i;
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int m			= Params->PARS[9];  //CK//  the number of exposed classes.  adjustable and fit.  yay!

double C		= y[m+1];
//double C		= Params->PARS[0];		//CK// Something majorly wrong here!!  Think Eli had the conidia turned off or something
//Need to fix because it is a constant right now	(set to Params->PARS[0]).  Going to try y[7] bc that should be C density...
//double nuV		= Params->nuV;		//Params->PARS[2];
double nuV		= 0.0;		//CK//  FOR FUNGUS ONLY MODEL. NO VIRUS INFECTION
double nuF		= Params->nuF;		//Params->PARS[3];
//double nuR		= 1.0;			
double nuR		= Params->nuR;	//param for nuR after it includes rainfall
//double muV		= Params->PARS[5];	//muF=Params->PARS[6];	//muR=PARAMS[0].PARS[7];// breakdown rate (per day)	
double muV		= 0.0;		//CK// FUNGUS ONLY MODEL.  MAKE SURE DECAY IS ZERO SO NEGATIVE VIRUS DOESN'T HAPPEN!!
double lambdaV	= Params->PARS[17];	double lambdaF=Params->PARS[18];		// decay rate of E_V and E_F				

double muF		= Params->muF;  //CK//  Using the site-specific muF that was loaded up in DDEVF14
double size_C		= Params->size_C;  //Scaling effect of size on susceptibility over time.

double k		= Params->PARS[4];

int DIM = Params->PARS[9]+3;

double S0 = Params->INITS[0];
double R  = Params->POPS[3];

//int i;
int a;
int b;

for (i=1;i<=18;i++)	{
	printf("parm(%d)=%e\n",i,Params->PARS[i]);
}getc(stdin);

double hetero_pow = pow((y[0]/S0),C);

*(dfdy + 0 * DIM +0) = - nuF*y[m+1] - nuR*R; // \partial (dS/dt) \partial S
for(i=1; i <= m; i++){               //\partial (dS/dt) \partial E_i
	*(dfdy + 0 * DIM +i) = 0;
}
*(dfdy + 0 * DIM + m + 1) = -nuF*y[0]; //\partial (dS/dt) \partial C
*(dfdy + 0 * DIM + m + 2) = 0; //\partial (dS/dt) \partial R


*(dfdy + 1 * DIM +0) = nuF*y[m+1] + nuR*R; // \partial (dE1/dt) \partial S
*(dfdy + 1 * DIM +1) = -m*lambdaF;  // \partial (dE1/dt) \partial E1
for(i=2; i <= m; i++){              // \partial (dE1/dt) \partial Ei, i > 2
	*(dfdy + 1 * DIM +i) = 0;
}
*(dfdy + 1 * DIM + m + 1) = nuF*y[0]; // \partial (dE1/dt) \partial C
*(dfdy + 1 * DIM + m + 2) = 0; // \partial (dE1/dt) \partial R


//This next bit (7 lines) could apparently be merged with the following chunk
*(dfdy + 2 * DIM +0) = 0;  // \partial (dE2/dt) \partial S
*(dfdy + 2 * DIM +1) = m*lambdaF;  // \partial (dE2/dt) \partial E1
*(dfdy + 2 * DIM +2) = -m*lambdaF;  // \partial (dE2/dt) \partial E2
for(i=3; i <= m+1; i++){
	*(dfdy + 2 * DIM +i) = 0;  // \partial (dE2/dt) \partial Ei, i > 3
}
*(dfdy + 2 * DIM + m + 2) = 0; // \partial (dE2/dt) \partial E2


//Now we do the remaining exposed classes
for(a=3; a <= m; a++){
	for(b=0; b <= m+2; b++){
		if(b == (a-1)){*(dfdy + a * DIM +b) =  m*lambdaF;} // \partial (dEi/dt) \partial Ei
		else if (b == a){*(dfdy + a * DIM +b) = -m*lambdaF;}  // \partial (dEi/dt) \partial Ei+1
 		else {*(dfdy + a * DIM +b) = 0;}
	}
}

//Now for conidia 
for(i=0; i <= m-1; i++){
	*(dfdy + m+1 * DIM +i) = 0;  // \partial (dC/dt) \partial S is m = 0, and then \partial (dC/dt) \partial Ei,  for i = 1 through m-1 
}
*(dfdy + m+1 * DIM +m) = m*lambdaF*size_C; // \partial (dC/dt) \partial Em
*(dfdy + m+1 * DIM +m+1) = -muF;  // \partial (dC/dt) \partial C
*(dfdy + m+1 * DIM + m+2) = 0; // \partial (dC/dt) \partial R


//Now for R, the cumulative dead
for(i=0; i <= m-1; i++){
	*(dfdy + m+2 * DIM +i) = 0;  // \partial (dR/dt) \partial S, \partial Ei, for i = 1 through m - 1;
}
*(dfdy + m+2 * DIM + m) = m*lambdaF;   // \partial (dR/dt)  \partial Em
*(dfdy + m+2 * DIM + m+1) = 0; // \partial (dR/dt)  \partial C
*(dfdy + m+2 * DIM + m+2) = 0; // \partial (dR/dt)  \partial R


//Nothing is a continuous function of time (note that transmission rates are discrete functions of time)
for(i=0; i <= DIM-1; i++){
	dfdt[i] = 0;
}

return GSL_SUCCESS;
}

// ------------------------------------------  ODE Solver  ----------------------------------------------- //
double ODE_Solver(double t_ode,double t_end,void *Paramstuff,double *y_ode)
{
int i;
int status_ode;
double h_init=1.0e-5;

STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int DIM = Params->PARS[9]+3;

//printf("ODE_solver:\n parm 2=%f parm 3=%f parm 4=%f parm 5=%f parm 6=%f\n",
//	   Params->PARS[2],Params->PARS[3],Params->PARS[4],Params->PARS[5],Params->PARS[6]);
//getc(stdin);

const gsl_odeiv_step_type *solver_ode	= gsl_odeiv_step_rkf45; // Runge-Kutta Fehlberg (4, 5)
//const gsl_odeiv_step_type *solver_ode = gsl_odeiv_step_rk4; 

// returns pointer to a newly allocated instance of a stepping function of type 'solver_ode' for a system of DIM dimensions //
gsl_odeiv_step *step_ode	= gsl_odeiv_step_alloc(solver_ode, DIM);

gsl_odeiv_control *tol_ode	= gsl_odeiv_control_standard_new(1.0e-10, 1.0e-5, 1.0, 0.2);
gsl_odeiv_evolve *evol_ode	= gsl_odeiv_evolve_alloc(DIM);

gsl_odeiv_system sys_ode;
sys_ode.function  = fast_odes;
sys_ode.jacobian  = jacobian;
sys_ode.dimension = (size_t)(DIM);
sys_ode.params	  = Params;

//double y_err[DIM]; double dydt_in[DIM];	double dydt_out[DIM];

// ----------------------------------- Integrate Over Time ------------------------------------ //
while (t_ode<t_end)	{
	status_ode = gsl_odeiv_evolve_apply(evol_ode, tol_ode, step_ode, &sys_ode, &t_ode, t_end, &h_init, y_ode);
	//status_ode = gsl_odeiv_step_apply(step_ode, t_ode, h_init, y_ode, y_err, dydt_in, dydt_out, &sys_ode);
	
	/*for (i=0;i<DIM+1;i++)	{
		printf("y_ode(%d)=%e\n",i,y_ode[i]);
		//dydt_in[i]=dydt_out[i];
	}*/

	for (i=0;i<=DIM;i++)	{
		if (y_ode[i]>0)		{
			// keep y_ode as is
		}
		else				{
			//printf("y(%d) NEGATIVE or not a number\n",i);
			y_ode[i]=0;
		}
		if (y_ode[i]>Params->INITS[0])	{
			//printf("y(%d) TOO LARGE!!\n",i);
			y_ode[i]=Params->INITS[0];
		}
	}
	//t_ode+=h_init;

}
// -------------------------------------- Clear Memory ----------------------------------------- //
gsl_odeiv_evolve_free(evol_ode);
gsl_odeiv_control_free(tol_ode);
gsl_odeiv_step_free(step_ode);

return (t_end);
}


