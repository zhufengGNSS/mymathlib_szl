////////////////////////////////////////////////////////////////////////////////
// File: gauss_chebyshev_100pts.c                                             //
// Routines:                                                                  //
//    double Gauss_Chebyshev_Integration_100pts( double (*f)(double) )        //
//    void   Gauss_Chebyshev_Zeros_100pts( double zeros[] )                   //
//    void   Gauss_Chebyshev_Coefs_100pts( double coef[] )                    //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// The zeros of the Chebyshev polynomial T100(x) = cos(100 * arccos(x)) are   //
// the positive and negative values of the elements in the array x below.     //
// The coefficient of the Gauss-Chebyshev formula is A = PI / 100.            //
////////////////////////////////////////////////////////////////////////////////

static const double x[] = {
    9.99876632481660598665e-01,    9.98889874961969972630e-01,
    9.96917333733127976206e-01,    9.93960955455179687732e-01,
    9.90023657716557567258e-01,    9.85109326154773918009e-01,
    9.79222810621765780883e-01,    9.72369920397676601822e-01,
    9.64557418457798093668e-01,    9.55793014798330126634e-01,
    9.46085358827545318532e-01,    9.35444030829867325186e-01,
    9.23879532511286756101e-01,    9.11403276635445248216e-01,
    8.98027575760615630915e-01,    8.83765630088693424299e-01,
    8.68631514438191247747e-01,    8.52640164354092221514e-01,
    8.35807361368270258491e-01,    8.18149717425023432148e-01,
    7.99684658487090538673e-01,    7.80430407338329735838e-01,
    7.60405965600030938169e-01,    7.39631094978609697455e-01,
    7.18126297763188830355e-01,    6.95912796592314325517e-01,
    6.73012513509773338712e-01,    6.49448048330183655701e-01,
    6.25242656335705172918e-01,    6.00420225325884049784e-01,
    5.75005252043278565911e-01,    5.49022817998131743543e-01,
    5.22498564715948865002e-01,    4.95458668432407538050e-01,
    4.67929814260573377236e-01,    4.39939169855915140829e-01,
    4.11514358605108774056e-01,    3.82683432365089771723e-01,
    3.53474843779257124726e-01,    3.23917418198149414390e-01,
    2.94040325232303957766e-01,    2.63873049965372896973e-01,
    2.33445363855905411768e-01,    2.02787295356512483442e-01,
    1.71929100279409546607e-01,    1.40901231937582661159e-01,
    1.09734311091045268027e-01,    7.84590957278449450322e-02,
    4.71064507096426609065e-02,    1.57073173118206757541e-02
};

static const double A = 3.14159265358979323838e-02;

#define NUM_OF_POSITIVE_ZEROS  sizeof(x) / sizeof(double)
#define NUM_OF_ZEROS           NUM_OF_POSITIVE_ZEROS+NUM_OF_POSITIVE_ZEROS

////////////////////////////////////////////////////////////////////////////////
//  double Gauss_Chebyshev_Integration_100pts( double (*f)(double) )          //
//                                                                            //
//  Description:                                                              //
//     Approximate the integral of f(x) / sqrt(1 - x^2) from -1 to 1 using    //
//     the 100 point Gauss-Chebyshev integral approximation formula.          //
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
//        integral = Gauss_Chebyshev_Integration_100pts( f );                 //
//        ...                                                                 //
//     }                                                                      //
//     double f(double x) { define f }                                        //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Gauss_Chebyshev_Integration_100pts( double (*f)(double) ) {

   double integral = 0.0;
   const double *px;

   for (px = &x[NUM_OF_POSITIVE_ZEROS - 1]; px >= x; px--) 
      integral +=   (*f)(*px) + (*f)(- *px);

   return A * integral; 
}


////////////////////////////////////////////////////////////////////////////////
//  void Gauss_Chebyshev_Zeros_100pts( double zeros[] )                       //
//                                                                            //
//  Description:                                                              //
//     Returns the zeros of the Chebyshev polynomial                          //
//     T100(x) = cos(100 arccos(x)).                                          //
//                                                                            //
//  Arguments:                                                                //
//     double zeros[] Array in which to store the zeros of T100.  This array  //
//                    should be dimensioned 100 in the caller function.       //
//                    The order is from the minimum zero to the maximum.      //
//                                                                            //
//  Return Values:                                                            //
//     none                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double z[100];                                                         //
//     int i;                                                                 //
//                                                                            //
//     Gauss_Chebyshev_Zeros_100pts( z );                                     //
//     printf("The zeros of the Chebyshev polynomial T100 are:");             //
//     for ( i = 0; i < 100; i++) printf("%12.6le\n",z[i]);                   //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Gauss_Chebyshev_Zeros_100pts( double zeros[] ) {
   
   const double *px = &x[NUM_OF_POSITIVE_ZEROS - 1];
   double *pz = &zeros[NUM_OF_ZEROS - 1];

   for (; px >= x; px--)  {
      *(zeros++) = - *px;
      *(pz--) = *px;
   }   
}


////////////////////////////////////////////////////////////////////////////////
//  void Gauss_Chebyshev_Coefs_100pts( double coef[] )                        //
//                                                                            //
//  Description:                                                              //
//     Returns the coefficients for the 100 point Gauss-Chebyshev formula.    //
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
//     Gauss_Chebyshev_Coefs_100pts( &a );                                    //
//     printf("The coefficient for the Gauss-Chebyshev formula is :\          //
//                                                              %12.6lf",a);  //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Gauss_Chebyshev_Coefs_100pts( double *coef) {

  *coef = A;
}
