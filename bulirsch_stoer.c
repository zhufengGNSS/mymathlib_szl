////////////////////////////////////////////////////////////////////////////////
// File: bulirsch_stoer.c                                                     //
// Routines:                                                                  //
//    Gragg_Bulirsch_Stoer                                                    //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Gragg-Bulirsch-Stoer method for approximating the solution of the  //
//     differential equation y'(x) = f(x,y) with initial condition y = c when //
//     x = x0 is a variable step method which uses Graggs modified midpoint   //
//     method for a given step size h to evaluate y(x+h), a rational function //
//     or polynomial approximation to the limit as h -> 0, and an algorithm to//
//     estimate the error based on successive estimates to y(x+h) and adjust  //
//     the step size so that the absolute error is within a preassigned       //
//     tolerance.                                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>

static int number_of_steps[] = { 2,4,6,8,12,16,24,32,48,64,96,128 };

#define ATTEMPTS sizeof(number_of_steps)/sizeof(number_of_steps[0])

static double Graggs_Method( double (*f)(double, double), double y0, double x0,
                                              double x, int number_of_steps );
static int Rational_Extrapolation_to_Zero( double *fzero, double tableau[],
                                               double x[], double f, int n );
static int Polynomial_Extrapolation_to_Zero( double *fzero, double tableau[], 
                                                double x[], double f, int n );

////////////////////////////////////////////////////////////////////////////////
// int Gragg_Bulirsch_Stoer( double (*f)(double, double), double y0,          //
//             double *y1, double x, double h, double *h_new, double epsilon, //
//                                  double yscale, int rational_extrapolate ) //
//                                                                            //
//  Description:                                                              //
//     This function solves the differential equation y'=f(x,y) with the      //
//     initial condition y=y0 at x.  The value returned is y1 which is the    //
//     is the value of y evaluated at x + h, h_new is the predicted step size //
//     so that the accuracy is maintained. The procedure terminates when      //
//     the absolute difference of two successive extrapolated scaled estimates//
//     of the solution are less than epsilon.  The variable                   //
//     rational_extrapolation determines whether rational or polynomial       //
//     extrapolation is used.  If rational_extrapolate is non-zero (true),    //
//     then rational extrapolation is used while if rational_extrapolate is   //
//     zero (false), then polynomial extrapolation is used.                   //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//        Pointer to the function which returns the slope at (x,y) of         //
//        integral curve of the differential equation y' = f(x,y)             //
//        which passes through the point (x0,y0) corresponding to the         //
//        initial condition y(x0) = y0.                                       //
//     double y0                                                              //
//        The initial value of y at x.                                        //
//     double *y1                                                             //
//        The pointer to the value of y at x + h.                             //
//     double x                                                               //
//        Initial value of x.                                                 //
//     double h                                                               //
//        Initial step size, x + h is abscissa for the return value.          //
//     double *h_new                                                          //
//        The pointer to the new step size required to maintain accuracy.     //
//     double epsilon                                                         //
//        The tolerance.                                                      //
//     double yscale                                                          //
//        A non-zero value which normalizes the estimates for y for testing   //
//        for convergence. i.e. if y1 and y2 are two successive estimates     //
//        for y(x), the test for convergence is                               //
//        |y1/yscale - y2/yscale| < epsilon. yscale must be non-zero.         //
//     int    rational_extrapolate                                            //
//        A flag which if non-zero, then rational extrapolation to zero is    //
//        used and if zero, then polynomial extrapolation to zero is used.    //
//                                                                            //
//  Return Values:                                                            //
//     The solution of y' = f(x,y) at x + h starting with y[0] at x is        //
//     returned via *y1.                                                      //
//     The function returns:                                                  //
//         0 if success                                                       //
//        -1 if the process failed to converge                                //
//        -2 if an attempt was made to divide by zero .                       //
//        -3 yscale is zero.                                                  //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Gragg_Bulirsch_Stoer( double (*f)(double, double), double y0, double *y1,
            double x, double h, double *h_new, double epsilon, double yscale,
                                                  int rational_extrapolate  ) {

   double step_size2[ATTEMPTS];
   double tableau[ATTEMPTS+1];
   double dum;
   double est;
   double old_est;

   int (*Extrapolate)(double*,double*,double*,double,int);
   int i;
   int err;

   /* Perform the first estimate of y(x+h), store the step size which */
   /* was used squared and the estimate for subsequent called to the  */
   /* rational function approximation for the value of y(x+h) as the  */
   /* square of the step size used tends to zero.                     */

   if (yscale == 0.0) return -3;

   if (rational_extrapolate) Extrapolate = Rational_Extrapolation_to_Zero;
   else Extrapolate = Polynomial_Extrapolation_to_Zero;
 
   est = Graggs_Method( f, y0, x, x+h, number_of_steps[0] );
   step_size2[0] = (dum = h / (double) number_of_steps[0], dum * dum);
   *y1 = est;

   if (err = Extrapolate(y1, tableau, step_size2, est, 0) < 0) return err-1;

    /* Continue using Gragg's method with smaller step sizes, followed    */
    /* by an estimate of y(x+h) using the rational function approximation */
    /* as the number of steps becomes large or step size used in Gragg's  */
    /* method tends to zero, and then halt if the absolute difference     */
    /* between the current estimate and the previous estimate is less     */
    /* than the user-specified tolerance.  If an attempt to divide by     */
    /* zero was made in the rational function approximation, then halt    */
    /* with a return code of -2. If the procedure fails to converge       */
    /* within ATTEMPTS - 1 iterations, then halt with a return code of -1.*/
    /* If an attempt is made with a step size of 0, then halt with a      */
    /* return code of -3.                                                 */

   for (i = 1; i < ATTEMPTS; i++) {
      old_est = *y1;
      est = Graggs_Method( f, y0, x, x+h, number_of_steps[i] );
      step_size2[i] = (dum = h / (double) number_of_steps[i], dum * dum);

      if (err = Extrapolate(y1, tableau, step_size2, est, i) < 0) return err-1;

      if ( fabs(*y1 / yscale - old_est / yscale) < epsilon ) {
         if (i > 1) *h_new = 8.0 * h / (double) number_of_steps[i-1];
         else *h_new = h;
         return 0;
      }
   }
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
//  double Graggs_Method( double (*f)(double, double), double y0, double x0,  //
//                                          double x, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     Gragg's method, also called the modified midpoint method, for solving  //
//     a differential equation y'(x) = f(x,y) with initial condition y(x0)= c //
//     is a second order method based upon the following theorem due to       //
//     Gragg.                                                                 //
//                                                                            //
//     Thm:  Let f(x,y) have continuous partial derivatives up to order 2N+2, //
//     y(x) be the exact solution to the differential equation y' = f(x,y)    //
//     with y(x0) = y0, let x[i] = x0 + i*h, and let Y(x,h) be defined        //
//     inductively by: Y(x0,h) = y0, Y(x[1],h) = y0 + h*f(x0,y0), and         //
//     Y(x[i],h) = Y(x[i-2],h) + 2h*f(x[i],Y(x[i],h)) for i >= 2.  Then       //
//     Y(x,h) = y(x) + Sum(h^2k [uk(x) + (-1)^(x-x0)/n vk(x)] + h^(2N+2)e(x,h)//
//     where the Sum extends from k = 1 to N, for all x in the domain of      //
//     definition, all h = (x - x0)/n, n = 1,....  The function e(x,h) for    //
//     fixed x is a bounded function of h.                                    //
//                                                                            //
//     The oscillating term (-1)^(x-x0)/n vk(x) can be removed by noticing    //
//     that       Yest(x,h) = 1/2[Y(x,h) + Y(x-h,h) + h*f(x,Y(x))]            //
//     has an expansion using the theorem above                               //
//            Yest(x,h) = y(x) + h^2 [u1(x) + 0.25*y''(x)] + ...              //
//     in which the oscillating term no longer appears in the leading term.   //
//                                                                            //
//     This function solves the differential equation y'=f(x,y) with the      //
//     initial condition y=y0 at x = x0.  The return value is the value of    //
//     the solution y evaluated at x.                                         //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y0) corresponding to the //
//                initial condition y(x0) = y0.                               //
//     double y0  The initial value of y at x0.                               //
//     double x0  Initial value of x.                                         //
//     double x   Final value of x.                                           //
//     int    number_of_steps  The number_of_steps must be a positive even    //
//                integer.  The step size h is (x - x0) / number_of_steps.    //
//                                                                            //
//  Return Values:                                                            //
//     The solution y(x) of the initial value problem  y' = f(x,y) with       //
//     y(x0) = y0.                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static double Graggs_Method( double (*f)(double, double), double y0, double x0,
                                              double x, int number_of_steps ) {

   double h = (x - x0) / (double) number_of_steps;
   double h2 =  h + h;
   double y1 = y0 + h * (*f)(x0,y0);
   double y2;

   while ( --number_of_steps ) {
      x0 += h;
      y2 = y0 + h2 * (*f)(x0,y1);
      y0 = y1;
      y1 = y2;
   } 
   return  0.5 * ( y0 + y1 + h * (*f)(x,y1) );
}


////////////////////////////////////////////////////////////////////////////////
//  static int Rational_Extrapolation_to_Zero( double *fzero,                 //
//                             double tableau[], double *x, double f, int n ) //
//                                                                            //
//  Description:                                                              //
//                                                                            //
//  The algorithm for rational extrapolation to x = 0 given a sequence of     //
//  support points (x[i],f[i]) for i = 0,...,n consists of forming a tableau  //
//  T[row,col], row = 0,...,n and col = -1,0,...,n for which T[row,col]       //
//  is defined inductively as follows:                                        //
//               T[row,0] = f[row], row = 0,...,n,                            //
//               T[row,-1] = 0. for all row = 0,...,n,                        //
//  and T[row,col] = T[row,col-1] +                                           //
//       (T[row,col-1] - T[row-1,col-1] ) / Denominator(row,col); where       //
//  Denominator(row,col) = ( x[row-col] / x[row] ) *                          //
//   (1 - (T[row,col-1]-T[row-1,col-1]) / (T[row,col-1]-T[row-1,col-2])) - 1, //
//  for col = 1,..., n, row=col,...,n, and T[row,col] is indeterminant        //
//  elsewhere.  The result, f(0), is then T[n,n].                             //
//                                                                            //
//  A detailed explanation can be found in J.Stoer & R.Bulirsch Introduction  //
//  to Numerical Analysis, Second Edition (1992).  The equation above can be  //
//  found on page 71 equation (2.2.3.8) with x set to 0.                      //
//                                                                            //
//  Arguments:                                                                //
//     double *fzero     The estimation of f(0).                              //
//     double tableau[]  A working storage array of dimension at least        //
//                       ATTEMPTS + 1.                                        //
//     double x[]        The support abscissa.                                //
//     double f          The support ordinate f at x[n].                      //
//     int    n          The index of the last support point x[] and f[].     //
//                       i.e. there are n+1 support points.                   //
//                                                                            //
//  Return Values:                                                            //
//     If successful, then the value returned is 0.  If an attempt is made    //
//     to divide by zero, then -1 is returned, if x[n] = 0.0, then -2 is      //
//     returned.                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static int Rational_Extrapolation_to_Zero( double *fzero, double tableau[],
                                               double x[], double f, int n ) {
  
   double t, up, across, denominator, dum;
   int col;

   if (n == 0) {  *fzero = f; tableau[0] = f; return 0; }
   if ( x[n] == 0.0 ) { *fzero = f; return -2; }
   
   across = 0.0;                                                        
   up = tableau[0];                                                    
   tableau[0] = f;                                               

   for (col = 1; col <= n; col++) {
      denominator = tableau[col-1] - across;                                  
      if (denominator == 0.0) return -1;
      dum = 1.0 - (tableau[col-1] - up) / denominator;
      denominator = (x[n - col] / x[n]) * dum - 1.0;
      if (denominator == 0.0) return -1;
      t = tableau[col-1] + ( tableau[col-1] - up ) / denominator;
      across = up;
      up = tableau[col];
      tableau[col] = t;
   }
   *fzero = t;
   return 0;
}

 
////////////////////////////////////////////////////////////////////////////////
// static int Polynomial_Extrapolation_to_Zero( double *fzero,                //
//                           double tableau[], double x[], double f, int n )  //
//                                                                            //
//  Description:                                                              //
//                                                                            //
//  The algorithm for polynomial extrapolation to x = 0 given a sequence of   //
//  support points (x[i],f[i]) for i = 0,...,n consists of forming a tableau  //
//  T[row,col], row = 0,...,n and col = -1,0,...,n for which T[row,col] is    //
//  defined inductively as follows:                                           //
//               T[row,-1] = 0, T[row,0] = f[row], row = 0, ..., n,and        //
//    T[row,col] = T[row,col-1] +                                             //
//      ( (T[row,col-1] - T[row-1,col-1] )* (T[row,col-1] - T[row-1,col-2]) ) //
//          /  ( ( x[row-col] / x[row] ) * (T[row-1,col-1] - T[row-1,col-2] ) //
//               - ( T[row,col-1] - T[row-1,col-2] ) )                        //
//                                                                            //
//  A detailed explanation can be found in J.Stoer & R.Bulirsch Introduction  //
//  to Numerical Analysis, Second Edition (1992), Chapters 2 , 3, and 7.      //
//                                                                            //
//  Arguments:                                                                //
//     double *fzero     The estimation of f(0).                              //
//     double tableau[]  A working storage array of dimension at least        //
//                       ATTEMPTS + 1.                                        //
//     double x[]        The support abscissa, x[n].                          //
//     double f          The support ordinate f and x[n].                     //
//     int    n          The index of the support ordinate.                   //
//                                                                            //
//  Return Values:                                                            //
//     If successful, then the value returned is 0.  If an attempt is made    //
//     to divide by zero, then -1 is returned, if x[n] = 0.0, then -2 is      //
//     returned.                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static int Polynomial_Extrapolation_to_Zero( double *fzero, double tableau[],
                                               double x[], double f, int n ) {

   double back_two_columns;    //  T[row,col-2];
   double old_aux;             //  T[row-1,col];
   double new_value;           //  T[row,col];
   double vertical_diff;       //  T[row,col]-T[row-1,col]
   double backslant_diff;      //  T[row,col]-T[row,col-1]
   double forwardslant_diff;   //  T[row,col]-T[row-1,col-1];
   double denominator;        
   int i;

   if (n == 0) { tableau[0] = f; return 0; }
   if ( x[n] == 0.0 ) { *fzero = f; return -2; }

   back_two_columns = 0.0;
   old_aux = tableau[0];
   tableau[0] = f;
   for (i = 0; i < n; i++) {
      vertical_diff = tableau[i] - old_aux;
      backslant_diff = tableau[i] - back_two_columns;
      forwardslant_diff = backslant_diff - vertical_diff;
      denominator = (x[n-i-1]/x[n]) * forwardslant_diff - backslant_diff;
      if (denominator == 0.0) return -1;
      back_two_columns = old_aux;
      old_aux = tableau[i+1];
      tableau[i+1] = tableau[i]
                               + vertical_diff * backslant_diff / denominator;
   }
   *fzero = tableau[n];
   return 0;
}
