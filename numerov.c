////////////////////////////////////////////////////////////////////////////////
// File: numerov.c                                                            //
// Routines:                                                                  //
//    Numerovs_Method                                                         //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     Numerov's method for approximating the solution of the second order    //
//     differential equation y''(x) = f(x,y) with initial conditions          //
//     y(x0) = a, y'(x0) = c is an implicit fourth order method.  Numerov's   //
//     method can be derived by expanding y in a Taylor series about x and    //
//     adding the expansion evaluated at x+h and the expansion evaluated at   //
//     x-h:                                                                   //
//          y(x0+h) = y(x0) + h y'(x0) + (h^2/2) y''(x0) + (h^3/6)y'''(x0)    //
//                       + (h^4/24) y''''(x0) + h^5/120 y'''''(x0) + O(h^6)   //
//                                                                            //
//          y(x0-h) = y(x0) - h y'(x0) + (h^2/2) y''(x0) - (h^3/6)y'''(x0)    //
//                       + (h^4/24) y''''(x0) - h^5/120 y'''''(x0) + O(h^6)   //
//                                                                            //
//          y(x0+h) + y(x0-h) = 2y(x0) + h^2 y''(x0) + (h^4/12) y''''(x0)     //
//                                                                 + O(h^6)   //
//     then substitute f(x,y) for y''(x):                                     //
//          y(x0+h) + y(x0-h) = 2y(x0) + h^2 f(x0,y0) + (h^4/12) Df^2(x0,y0)  //
//                                                                 + O(h^6)   //
//     and finally approximate Df^2(x0,y0) by                                 //
//         Df^2(x0,y0) = (f(x0+h,y(x0+h)) -2f(x0,y0) + f(x0-h),y(x0-h)))/h^2  //
//                                                                            //
//     which yields the following                                             //
//          y(x0+h) + y(x0-h) = 2y(x0) + h^2 f(x0,y0)                         //
//            + (h^2/12) (f(x0+h,y(x0+h)) -2f(x0,y0) + f(x0-h),y(x0-h)))      //
//                                                                            //
//     Note that here f(x,y) does not depend upon y'.                         //
//                                                                            //
//     Letting x[n] = x0 + nh, y[n] = y(x[n]), and f[n] = f(x[n],y[n]), the   //
//     procedure proceeds recursively via the implicit equation for y[n+1]    //
//     as follows:                                                            //
//                                                                            //
//  y[n+1] = 2y[n] - y[n-1] + h^2 ( f[n] + ( f[n+1] - 2f[n] + f[n-1] ) / 12 ) //
//                                                                            //
//     In order to start the recursion, two successive values of y are        //
//     required, one of which is y[0] = a and the other y[1] needs a separate //
//     approximation.  The routine, Numerovs_Start, below uses the Runge-     //
//     Kutta procedure to approximate y[1].  Particular classes of problems   //
//     may have more accurate estimates for y[1].                             //
//                                                                            //
//     If u[n] is defined by u[n] = y[n] - (h^2 / 12) * f[n], then the        //
//     recursion becomes                                                      //
//                                                                            //
//                 u[n+1] = 2u[n] - u[n-1] + h^2 f[n],                        //
//                                                                            //
//     which is numerically evaluated as follows:                             //
//                                                                            //
//                 u[0] = y[0] - h^2 * f[0],                                  //
//                 z[0] = (y[1] - y[0]) / h - h (f[1] - f[0] ) / 12           //
//                                                                            //
//                 u[i] = u[i-1] + h * z[i-1]                                 //
//                 z[i] = z[i-1] + h * f[i]                                   //
//                                                                            //
//     Richardson extrapolation may be used to increase to increase the       //
//     order.                                                                 //
////////////////////////////////////////////////////////////////////////////////

static const double richardson[] = {  1.0 / 15.0, 1.0 / 63.0, 1.0 / 255.0,
1.0 / 1023.0, 1.0 / 4095.0, 1.0 / 16384.0  };

#define MAX_COLUMNS 1+sizeof(richardson)/sizeof(richardson[0])

static double Numerov( double (*f)(double, double),
                       double (*g)(double,double,double), double x,
                       double *z, double *u, double h, int number_of_steps );

static double Numerovs_Start( double (*f)(double, double), double f0, double x,
                                   double y, double c, double *u, double h  );


////////////////////////////////////////////////////////////////////////////////
// void Numerovs_Method( double (*f)(double, double),                         //
//                       double (*g)(double,double,double),                   //
//                       double y[], double x0, double c, double h,           //
//                       int richardson_columns, int number_of_steps );       //
//                                                                            //
//  Description:                                                              //
//     This function approximates the solution of the second order different- //
//     ial equation y'' = f(x,y) with initial condition y(x0) = y[0] and      //
//     y'(x0) = c.  The solution is returned in the array y[], where          //
//     y[n] = y(x0 + n*h), n = 0,...,number_of_steps.                         //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the rate of change of //
//                the slope with respect to x at the point (x,y) of integral  //
//                curve of the second order differential equation y'' = f(x,y)//
//                which passes through the point (x0,y[0]) with initial slope //
//                y'(x0) = c.                                                 //
//     double *g  Pointer to the function g(x,h,u) which returns the value y  //
//                such that u = y - h^2 f(x,y) / 12.                          //
//     double y[] On input, y[0] is the initial value of y.                   //
//                On output, y[n] is the approximate solution of the initial  //
//                value problem y''=f(x,y), y(x0)=y[0], y'(x0)=c where        //
//                y[n] = y(x0+nh).  The array y[] should be dimensioned at    //
//                least number_of_steps + 1 in the calling routine.           //
//     double x0  The initial value of x.                                     //
//     double c   The initial value y'(x0).                                   //
//     double h   The step size.                                              //
//     int richardson_columns  The number of columns to used for Richardson   //
//                extrapolation to the limit as h -> 0.                       //
//     int number_of_steps     The number of steps.                           //
//                                                                            //
//  Return Values:                                                            //
//     The function is of type void and therefore does not return a value.    //
//     The solution of the second order initial problem is returned in the    //
//     the array y[].                                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Numerovs_Method( double (*f)(double, double), 
                      double (*g)(double,double,double),
                      double y[], double x0, double c, double h, 
                      int richardson_columns, int number_of_steps ) {

   double dt[MAX_COLUMNS];         // dt[i] is the last element in column i.
   double z[MAX_COLUMNS];                                        
   double u[MAX_COLUMNS];
   double f0 = f(x0,y[0]);
   double integral;
   double delta;
   double h_old;
   int i,j,k, number_sub_intervals;

         /* Restrict the number of columns to use for Richardson  */
       /* extrapolation to be between 1 and MAX_COLUMNS inclusively */

   if (richardson_columns < 1) richardson_columns = 1;
   if (richardson_columns > MAX_COLUMNS) richardson_columns = MAX_COLUMNS;

       /* Initialize the starting values for u and z for each column */

   h_old = h;
   for (i = 0; i < richardson_columns; i++) {
      z[i] = Numerovs_Start( f, f0, x0, y[0], c, &u[i], h_old );
      h_old *= 0.5;
   }

        /*  For each step, approximate the solution by performing */
               /*  Richardson Extrapolation to the limit. */

   for (i = 0; i < number_of_steps; i++) {
      number_sub_intervals = 1;
      h_old = h;
      for (j = 0; j < richardson_columns; j++) {
         integral = Numerov( f, g, x0, &z[j], &u[j], h_old,
                                                        number_sub_intervals);
         for ( k = 0; k < j; k++) {
            delta = integral - dt[k];
            dt[k] = integral;
            integral += richardson[k] * delta;
         }
         dt[j] = integral;
         h_old *= 0.5;
         number_sub_intervals *= 2;
      }
      y[i+1] = integral;
      x0 += h;
   }
}


////////////////////////////////////////////////////////////////////////////////
//  static double Numerov( double (*f)(double, double),                       //
//                      double (*g)(double,double,double), double x,          //
//                      double *z, double *u, double h, int number_of_steps );//
//                                                                            //
//  Description:                                                              //
//     This function applies the recursion u += h * z, z += h * f(x,y)        //
//     number_of_steps times and returns the estimate for y at x + h *        //
//     number_of_steps.  The values for u and z are updated for the successive//
//     calls.                                                                 //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the rate of change of //
//                the slope with respect to x at the point (x,y) of integral  //
//                curve of the second order differential equation y'' = f(x,y)//
//                which passes through the point (x0,y[0]) with initial slope //
//                y'(x0) = c.                                                 //
//     double *g  Pointer to the function g(x,h,u) which returns the value y  //
//                such that u = y - h^2 f(x,y) / 12.                          //
//     double x   Initial value of x.                                         //
//     double *z  The auxillary variable used to numerically evaluate u       //
//                recursively.                                                //
//     double *u  u = y - h^2 f(x,y) / 12.                                    //
//     double h   Step size.                                                  //
//     int number_of_steps     The number of steps.                           //
//                                                                            //
//  Return Values:                                                            //
//     The solution y''=f(x,y) at x + h * number_of_steps.                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static double Numerov( double (*f)(double, double),
                       double (*g)(double,double,double), double x,
                       double *z, double *u, double h, int number_of_steps ) {

   double y;
   int i;

   for (i = 0; i < number_of_steps; i++) {
      x += h;
      *u += h * *z;
      y = (*g)(x, h, *u); 
      *z += h * (*f)(x,y);
   }
   return y;
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// static double Numerovs_Start( double (*f)(double, double), double f0,      //
//                    double x0, double y0, double c, double *u, double h  ); //
//                                                                            //
//  Description:                                                              //
//     Numerov's procedure for approximating the solution of the second order //
//     differential equation y''(x) = f(x,y) with initial conditions          //
//     y(x0) = y[0] and y'(x0) = c requires two consecutive values, y[0],     //
//     and y[1], to recursively generate y[n] for n >= 2.  This routine       //
//     returns the value z[0] and sets u[0] = y[0] - h^2 f(x0,y0) / 12  which //
//     are used to recursively generate u[n] and z[n], from which y[n] is     //
//     solved using the implicit relation y[n] = u[n] + h^2 f(x[n],y[n]) / 12.//
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the rate of change of //
//                the slope with respect to x at the point (x,y) of integral  //
//                curve of the second order differential equation y'' = f(x,y)//
//                which passes through the point (x0,y[0]) with initial slope //
//                y'(x0) = c.                                                 //
//     double f0  f(x0,y0) i.e. the value of f at the initial point (x0,y0).  //
//     double x0  Initial value of x.                                         //
//     double y0  The initial value of y.                                     //
//     double c   y'(x0).                                                     //
//     double *u  u = y0 - h^2 * f(x0,y0) / 12.                               //
//     double h   The step size.                                              //
//                                                                            //
//  Return Values:                                                            //
//     The routine returns the estimate of z[0] and sets u[0].                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static double Numerovs_Start( double (*f)(double, double), double f0, double x0,
                                   double y0, double c, double *u, double h  ) {

   static const double one_sixth = 1.0 / 6.0;
   static const double one_twelfth = 1.0 / 12.0;
   double h2 = 0.5 * h;
   double k1 = h * f0;
   double k2 = h * (*f)( x0 + h2, y0 + h2 * c );
   double k3 = h * (*f)( x0 + h2, y0 + h2 * c + 0.25 * k1 * h );
   double dy = c + one_sixth * ( k1 + k2 + k3 );
   double f1 = (*f)( x0 + h, y0 + h * dy );

   *u = y0 - one_twelfth * f0 * h * h;
   return dy + h * (f0 - f1) * one_twelfth;
}
