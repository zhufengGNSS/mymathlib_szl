////////////////////////////////////////////////////////////////////////////////
// File: hermite_quadrature_1_derivative.c                                    //
// Routines:                                                                  //
//    Hermite_Quadrature_1_Derivative_LR                                      //
//    Hermite_Quadrature_1_Derivative_RL                                      //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Hermite quadrature formula extends the trapezoidal rule by using   //
//     the values of the first derivatives of the integrand at the endpoints  //
//     of the integration interval.  This technique is convenient when the    //
//     derivatives of the integrand are easy to calculate.                    //
//                                                                            //
//     This version of the Hermite quadrature is sometimes called the         //
//     corrected trapezoidal rule.                                            //
//                                                                            //
//     The integral over a fundamental interval [0,h] is                      //
//                  (h/2) [f(0) + f(h) + (h/6) (f'(0) - f'(h) ]               //
//     with a local truncation error of                                       //
//                             h^5 / 720 f''''(xm).                           //
//                                                                            //
//     In order to integrate a function using Hermite quadrature over a       //
//     closed and bounded interval [a,b], divide the interval [a,b] into N    //
//     subintervals of length h = (b-a) / N.  The integral of the function    //
//     over the interval [a,b] is the sum of the integrals of the function    //
//     over the subintervals [a,a+h], [a+h, a+2h],...,[b-h,b].  The integral  //
//     is approximated then by applying Hermite quadrature to each            //
//     subinterval.                                                           //
//                                                                            //
//     Let I(f) be the approximation of the integral of f(x) over [a,b] then  //
//     Hermite quadrature is:                                                 //
//       I(f) = h { f(x[0]) / 2 + f(x[1]) + ... + f(x[n-1]) + f(x[n]) / 2) }  //
//              + h^2/12 (f'(a) - f'(b))                                      //
//     where h is the step size, h = x[i] - x[i-1]; n is the number of        //
//     steps, n = abs(b - a) / h and x[0] = a, x[i] = x[i-1] + h for          //
//     i = 1,...,n-1 so that x[n] = b.                                        //
//                                                                            //
//     The truncation error estimate for the Hermite quadrature is            //
//                           n * h^5  f''''(xm) / 720                         //
//     where f''''(xm) is the value of the fourth derivative evaluated at     //
//     some (unknown) point a <= xm <= b.                                     //
//                                                                            //
//     Note!  If f'(b) = f'(a), Hermite quadrature is equivalent to the       //
//     trapezoidal rule.  In this case, the trapezoidal rule is superior to   //
//     Simpson's rule.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  double  Hermite_Quadrature_1_Derivative_LR( double a, double h, int n,    //
//                              double (*f)(double), double (*df)(double) );  //
//                                                                            //
//  Description:                                                              //
//     This function integrates f(x) from a to a+nh using the corrected       //
//     trapezoidal rule by summing from the left end of the interval to the   //
//     right end.                                                             //
//                                                                            //
//  Arguments:                                                                //
//     double a   The lower limit of the integration interval.                //
//     double h   The length of each subinterval.                             //
//     int    n   The number of steps, the upper limit of integration is      //
//                a + n * h.                                                  //
//     double *f  Pointer to the integrand, a function of a single variable   //
//                of type double.                                             //
//     double *df Pointer to the derivative of the integrand, a function of a //
//                single variable of type double.                             //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
//  Example:                                                                  //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double  Hermite_Quadrature_1_Derivative_LR( double a, double h, int n, 
                                 double (*f)(double), double (*df)(double) ) {

   double upper_limit = a + (n - .5) * h;
   double integral = 0.5 * (*f)(a)
               + 0.0833333333333333333333 * h * ( (*df)(a) - (*df)(a+n*h) );

   for (a += h; a < upper_limit; a += h) integral += (*f)(a);
   integral += 0.5 * (*f)(a);

   return h * integral;
}


////////////////////////////////////////////////////////////////////////////////
//  double  Hermite_Quadrature_1_Derivative_RL( double b, double h, int n,    //
//                              double (*f)(double), double (*df)(double) );  //
//                                                                            //
//  Description:                                                              //
//     This function integrates f(x) from b-nh to b using the corrected       //
//     trapezoidal rule by summing from the left end of the interval to the   //
//     right end.                                                             //
//                                                                            //
//  Arguments:                                                                //
//     double b   The upper limit of the integration interval.                //
//     double h   The length of each subinterval.                             //
//     int    n   The number of steps, the lower limit of integration is      //
//                b - n * h.                                                  //
//     double *f  Pointer to the integrand, a function of a single variable   //
//                of type double.                                             //
//     double *df Pointer to the derivative of the integrand, a function of a //
//                single variable of type double.                             //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from (b - n * h) to b.                            //
//                                                                            //
//  Example:                                                                  //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double  Hermite_Quadrature_1_Derivative_RL( double b, double h, int n, 
                                 double (*f)(double), double (*df)(double) ) {

   double lower_limit = b - (n - .5) * h;
   double integral = 0.5 * (*f)(b) 
               + 0.0833333333333333333333 * h * ( (*df)(b-n*h) - (*df)(b) );

   for (b -= h; b > lower_limit; b -= h) integral += (*f)(b);
   integral += 0.5 * (*f)(b);

   return h * integral;
}
