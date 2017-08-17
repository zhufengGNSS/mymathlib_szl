////////////////////////////////////////////////////////////////////////////////
// File: gauss_chebyshev_82pts.c                                              //
// Routines:                                                                  //
//    double Gauss_Chebyshev_Integration_82pts( double (*f)(double) )         //
//    void   Gauss_Chebyshev_Zeros_82pts( double zeros[] )                    //
//    void   Gauss_Chebyshev_Coefs_82pts( double coef[] )                     //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// The zeros of the Chebyshev polynomial T82(x) = cos(82 * arccos(x)) are     //
// the positive and negative values of the elements in the array x below.     //
// The coefficient of the Gauss-Chebyshev formula is A = PI / 82.             //
////////////////////////////////////////////////////////////////////////////////

static const double x[] = {
    9.99816528431794335197e-01,    9.98349159803241149804e-01,
    9.95416576111944937527e-01,    9.91023081328871215654e-01,
    9.85175123513469957874e-01,    9.77881285350259174256e-01,
    9.69152271552551843146e-01,    9.59000893151812961216e-01,
    9.47442048695704248755e-01,    9.34492702382410969923e-01,
    9.20171859163341760769e-01,    9.04500536850741684836e-01,
    8.87501735271154441605e-01,    8.69200402510005274293e-01,
    8.49623398296845308206e-01,    8.28799454584994521090e-01,
    8.06759133383438156143e-01,    7.83534781902864084870e-01,
    7.59160485081670494373e-01,    7.33672015561618530660e-01,
    7.07106781186547524382e-01,    6.79503770101205675834e-01,
    6.50903493530771234835e-01,    6.21347926325043119320e-01,
    5.90880445354560571787e-01,    5.59545765849064042664e-01,
    5.27389875771729388452e-01,    4.94459968325490245523e-01,
    4.60804372690504858895e-01,    4.26472483094419684230e-01,
    3.91514686319528932792e-01,    3.55982287753223298945e-01,
    3.19927436090259034555e-01,    2.83403046797357173681e-01,
    2.46462724452459162295e-01,    2.09160684072616744766e-01,
    1.71551671545978277379e-01,    1.33690883284648501394e-01,
    9.56338852163422834625e-02,    5.74365312337232426870e-02,
    1.91548812211141054467e-02
};

static const double A = 3.83121055315828443723e-02;

#define NUM_OF_POSITIVE_ZEROS  sizeof(x) / sizeof(double)
#define NUM_OF_ZEROS           NUM_OF_POSITIVE_ZEROS+NUM_OF_POSITIVE_ZEROS

////////////////////////////////////////////////////////////////////////////////
//  double Gauss_Chebyshev_Integration_82pts( double (*f)(double) )           //
//                                                                            //
//  Description:                                                              //
//     Approximate the integral of f(x) / sqrt(1 - x^2) from -1 to 1 using    //
//     the 82 point Gauss-Chebyshev integral approximation formula.           //
//                                                                            //
//  Arguments:                                                                //
//     double *f   Pointer to function of a single variable of type double.   //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) / sqrt(1 - x^2) from -1 to 1.                     //
//                                                                            //
//  Example:                                                                  //
//     {                                                                      //
//        double f(double);                                                   //
//        double integral;                                                    //
//                                                                            //
//        integral = Gauss_Chebyshev_Integration_82pts( f );                  //
//        ...                                                                 //
//     }                                                                      //
//     double f(double x) { define f }                                        //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Gauss_Chebyshev_Integration_82pts( double (*f)(double) ) {

   double integral = 0.0;
   const double *px;

   for (px = &x[NUM_OF_POSITIVE_ZEROS - 1]; px >= x; px--) 
      integral +=   (*f)(*px) + (*f)(- *px);

   return A * integral; 
}


////////////////////////////////////////////////////////////////////////////////
//  void Gauss_Chebyshev_Zeros_82pts( double zeros[] )                        //
//                                                                            //
//  Description:                                                              //
//     Returns the zeros of the Chebyshev polynomial T82 = cos(82 arccos(x)). //
//                                                                            //
//  Arguments:                                                                //
//     double zeros[] Array in which to store the zeros of T82.  This array   //
//                    should be dimensioned 82 in the caller function.        //
//                    The order is from the minimum zero to the maximum.      //
//                                                                            //
//  Return Values:                                                            //
//     none                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double z[82];                                                          //
//     int i;                                                                 //
//                                                                            //
//     Gauss_Chebyshev_Zeros_82pts( z );                                      //
//     printf("The zeros of the Chebyshev polynomial T82 are:");              //
//     for ( i = 0; i < 82; i++) printf("%12.6le\n",z[i]);                    //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Gauss_Chebyshev_Zeros_82pts( double zeros[] ) {
   
   const double *px = &x[NUM_OF_POSITIVE_ZEROS - 1];
   double *pz = &zeros[NUM_OF_ZEROS - 1];

   for (; px >= x; px--)  {
      *(zeros++) = - *px;
      *(pz--) = *px;
   }   
}


////////////////////////////////////////////////////////////////////////////////
//  void Gauss_Chebyshev_Coefs_82pts( double coef[] )                         //
//                                                                            //
//  Description:                                                              //
//     Returns the coefficients for the 82 point Gauss-Chebyshev formula.     //
//                                                                            //
//  Arguments:                                                                //
//     double coef[]  Array in which to store the coefficients of the Gauss-  //
//                    Chebyshev formula.  For Gauss-Chebyshev integration     //
//                    the coefficients are the same for each term, therefore  //
//                    the dimension of coef[] is only 1. I.e. the argument    //
//                    is the address of a double.                             //
//                                                                            //
//  Return Values:                                                            //
//     none                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double a;                                                              //
//     int i;                                                                 //
//                                                                            //
//     Gauss_Chebyshev_Coefs_82pts( &a );                                     //
//     printf("The coefficient for the Gauss-Chebyshev formula is :\          //
//                                                              %12.6lf",a);  //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Gauss_Chebyshev_Coefs_82pts( double *coef) {

  *coef = A;
}
