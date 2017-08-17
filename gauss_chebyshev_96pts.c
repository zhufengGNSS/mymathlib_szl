////////////////////////////////////////////////////////////////////////////////
// File: gauss_chebyshev_96pts.c                                              //
// Routines:                                                                  //
//    double Gauss_Chebyshev_Integration_96pts( double (*f)(double) )         //
//    void   Gauss_Chebyshev_Zeros_96pts( double zeros[] )                    //
//    void   Gauss_Chebyshev_Coefs_96pts( double coef[] )                     //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// The zeros of the Chebyshev polynomial T96(x) = cos(96 * arccos(x)) are     //
// the positive and negative values of the elements in the array x below.     //
// The coefficient of the Gauss-Chebyshev formula is A = PI / 96.             //
////////////////////////////////////////////////////////////////////////////////

static const double x[] = {
    9.99866137909561782863e-01,    9.98795456205172392701e-01,
    9.96655239309180324928e-01,    9.93447779019444395529e-01,
    9.89176509964780973456e-01,    9.83846005927077416097e-01,
    9.77461974943571863372e-01,    9.70031253194543992616e-01,
    9.61561797682961947144e-01,    9.52062677713924257114e-01,
    9.41544065183020778391e-01,    9.30017223684012117049e-01,
    9.17494496447491307935e-01,    9.03989293123443331582e-01,
    8.89516075421856035267e-01,    8.74090341626758851525e-01,
    8.57728610000272069893e-01,    8.40448401094438021037e-01,
    8.22268218989775107855e-01,    8.03207531480644909799e-01,
    7.83286749228650365376e-01,    7.62527203906388096352e-01,
    7.40951125354959091193e-01,    7.18581617779698057173e-01,
    6.95442635009611651129e-01,    6.71558954847018400619e-01,
    6.46956152534857365421e-01,    6.21660573370077408053e-01,
    5.95699304492433343462e-01,    5.69100145878898230586e-01,
    5.41891580574751716148e-01,    5.14102744193221726607e-01,
    4.85763393716340056256e-01,    4.56903875630420676552e-01,
    4.27555093430282094315e-01,    3.97748474527011052041e-01,
    3.67515936594703565403e-01,    3.36889853392220050703e-01,
    3.05903020096553462752e-01,    2.74588618184932341487e-01,
    2.42980179903263889945e-01,    2.11111552358965165921e-01,
    1.79016861276632682038e-01,    1.46730474455361751659e-01,
    1.14286964966846398116e-01,    8.17210741336682237475e-02,
    4.90676743274180142536e-02,    1.63617316264867816422e-02
};

static const double A = 3.27249234748936795667e-02;

#define NUM_OF_POSITIVE_ZEROS  sizeof(x) / sizeof(double)
#define NUM_OF_ZEROS           NUM_OF_POSITIVE_ZEROS+NUM_OF_POSITIVE_ZEROS

////////////////////////////////////////////////////////////////////////////////
//  double Gauss_Chebyshev_Integration_96pts( double (*f)(double) )           //
//                                                                            //
//  Description:                                                              //
//     Approximate the integral of f(x) / sqrt(1 - x^2) from -1 to 1 using    //
//     the 96 point Gauss-Chebyshev integral approximation formula.           //
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
//        integral = Gauss_Chebyshev_Integration_96pts( f );                  //
//        ...                                                                 //
//     }                                                                      //
//     double f(double x) { define f }                                        //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Gauss_Chebyshev_Integration_96pts( double (*f)(double) ) {
   
   double integral = 0.0;
   const double *px;

   for (px = &x[NUM_OF_POSITIVE_ZEROS - 1]; px >= x; px--) 
      integral +=   (*f)(*px) + (*f)(- *px);

   return A * integral; 
}


////////////////////////////////////////////////////////////////////////////////
//  void Gauss_Chebyshev_Zeros_96pts( double zeros[] )                        //
//                                                                            //
//  Description:                                                              //
//     Returns the zeros of the Chebyshev polynomial T96 = cos(96 arccos(x)). //
//                                                                            //
//  Arguments:                                                                //
//     double zeros[] Array in which to store the zeros of T96.  This array   //
//                    should be dimensioned 96 in the caller function.        //
//                    The order is from the minimum zero to the maximum.      //
//                                                                            //
//  Return Values:                                                            //
//     none                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double z[96];                                                          //
//     int i;                                                                 //
//                                                                            //
//     Gauss_Chebyshev_Zeros_96pts( z );                                      //
//     printf("The zeros of the Chebyshev polynomial T96 are:");              //
//     for ( i = 0; i < 96; i++) printf("%12.6le\n",z[i]);                    //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Gauss_Chebyshev_Zeros_96pts( double zeros[] ) {
   
   const double *px = &x[NUM_OF_POSITIVE_ZEROS - 1];
   double *pz = &zeros[NUM_OF_ZEROS - 1];

   for (; px >= x; px--)  {
      *(zeros++) = - *px;
      *(pz--) = *px;
   }   
}


////////////////////////////////////////////////////////////////////////////////
//  void Gauss_Chebyshev_Coefs_96pts( double coef[] )                         //
//                                                                            //
//  Description:                                                              //
//     Returns the coefficients for the 96 point Gauss-Chebyshev formula.     //
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
//     Gauss_Chebyshev_Coefs_96pts( &a );                                     //
//     printf("The coefficient for the Gauss-Chebyshev formula is :\          //
//                                                              %12.6lf",a);  //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Gauss_Chebyshev_Coefs_96pts( double *coef) {

  *coef = A;
}
