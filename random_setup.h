gsl_rng *random_setup(void)
{
const gsl_rng_type *TT;
long seedy;
srand((unsigned) time(NULL));
seedy = time(NULL)*(int)getpid();  //JL: Including the process id to avoid the situation where the jobs submitted at the same time point have identical initial random numbers drawn.
//seedy = -rand();
gsl_rng_env_setup ();
gsl_rng_default_seed = seedy;
TT = gsl_rng_default;

return gsl_rng_alloc(TT);
}




// -------------------------------------------------------------------- //

int ipow(int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = ipow(x, p/2);
  if (p%2 == 0)             return tmp * tmp;
  else                      return x * tmp * tmp;
}
