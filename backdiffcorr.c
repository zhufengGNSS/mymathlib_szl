////////////////////////////////////////////////////////////////////////////////
// File: backdiffcorr.c                                                       //
// Routines:                                                                  //
//    Backward_Difference_Correction                                          //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The backward difference correction procedure for approximating the     //
//     solution of the second order differential equation y''(x) = f(x,y)     //
//     with initial conditions y(x0) = a, y'(x0) = c is an explicit third     //
//     order method.                                                          //
//                                                                            //
//     Note that here f(x,y) does not depend upon y'.                         //
//                                                                            //
//     Letting x[n] = x0 + nh, y[n] = y(x[n]), and f[n] = f(x[n],y[n]), the   //
//     procedure proceeds recursively via the explicit equation for y[n+1]    //
//     as follows:                                                            //
//                                                                            //
//  y[n+1] = 2y[n] - y[n-1] + h^2 ( f[n] + ( f[n] - 2f[n-1] + f[n-2] ) / 12 ) //
//                                                                            //
//     In order to start the recursion, three successive values of y are      //
//     required, one of which is y[0] = a and the other two y[1] and y[2]     //
//     need separate approximations.  The routine, Backward_Start, below uses //
//     the Runge-Kutta procedure to approximate y[1] and y[2].  Particular    //
//     classes of problems may have more accurate estimates for both y[1]     //
//     and y[2].                                                              //
//                                                                            //
//     Richardson extrapolation may be used to increase to increase the       //
//     order.                                                                 //
////////////////////////////////////////////////////////////////////////////////

static const double richardson[] = {  1.0 / 7.0, 1.0 / 15.0, 1.0 / 31.0,
1.0 / 63.0, 1.0 / 127.0, 1.0 / 255.0  };

#define MAX_COLUMNS 1+sizeof(richardson)/sizeof(richardson[0])

static void Backward_Start( double (*f)(double, double), double f0, double x0,
    double y0, double c, double *y1, double *y2, double f_hist[], double h  );

static double Backward( double (*f)(double, double), double f_hist[], double x,
                      double *y1, double *y2, double h, int number_of_steps );

////////////////////////////////////////////////////////////////////////////////
// void Backward_Difference_Correction( double (*f)(double, double),          //
//      double y[], double x0, double c, double h, int richardson_columns,    //
//      int number_of_steps );                                                //
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
void Backward_Difference_Correction( double (*f)(double, double), double y[],
 double x0, double c, double h, int richardson_columns, int number_of_steps ) {

   double dt[MAX_COLUMNS];         // dt[i] is the last element in column i.
   double y1[MAX_COLUMNS];                                        
   double y2[MAX_COLUMNS];
   double f0 = f(x0,y[0]);
   double f_hist[MAX_COLUMNS][3];
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
      Backward_Start( f, f0, x0, y[0], c, &y1[i], &y2[i], &f_hist[i][0], 
                                                                       h_old );
      h_old *= 0.5;
   }

        /*  For each step, approximate the solution by performing */
               /*  Richardson Extrapolation to the limit. */

   for (i = 0; i < number_of_steps; i++) {
      number_sub_intervals = 1;
      h_old = h;
      for (j = 0; j < richardson_columns; j++) {
         integral = Backward( f, &f_hist[j][0], x0, &y1[j], &y2[j],
                                                  h_old, number_sub_intervals);
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
// static double Backward( double (*f)(double, double), double f_hist[],      //
//          double x, double *y1, double *y2, double h, int number_of_steps ) //
//                                                                            //
//  Description:                                                              //
//     This function applies the recursion                                    //
//       y[n+2] = 2 y[n] - y[n-1] + h (f[n] + (f[n] - 2 f[n-1] + f[n-2])/12)  //
//     number_of_steps times and returns the estimate for y at x + h * n,     //
//     where n is the number_of_steps.                                        //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the rate of change of //
//                the slope with respect to x at the point (x,y) of integral  //
//                curve of the second order differential equation y'' = f(x,y)//
//                which passes through the point (x0,y[0]) with initial slope //
//                y'(x0) = c.                                                 //
//     double f_hist[] On input, f_hist[0] is f(x,y0), f_hist[1] is           //
//                f(x+h,y1), and f_hist[2] is f(x+2h,y2).  On output,         //
//                f_hist[0] is f(x+nh, y(x+nh)), f_hist[1] is                 //
//                f(x+(n+1)h,y(x+(n+1)h)) and f_hist[2] is                    //
//                f(x+(n+2)h,y(x+n+2)h)), where n is the number_of_steps.     //
//     double x   The value of x at y0, i.e. y1 is the value of y(x) at x+h   //
//                and y2 is the value of x at x + 2h.                         //
//     double *y1 On input y1 is the value of y(x) at x = x+h,                //
//                on output y1 is the value of y(x) at x+(n+1)h.              //
//     double *y2 On input y2 is the value of y(x) at x = x+2h,               //
//                on output y2 is the value of y(x) at x+(n+2)h.              //
//     double h   Step size.                                                  //
//     int number_of_steps     The number of steps.                           //
//                                                                            //
//  Return Values:                                                            //
//     The solution y(x) of y''=f(x,y) at x + h * number_of_steps.   The      //
//     successive values of y(x) at x + h * (number_of_steps + 1) and at      //
//     x + h * (number_of_steps + 2) are set in y1 and y2 respectively and    //
//     the function calls f_hist are set for succeeding calls to Backward.    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static double Backward( double (*f)(double, double), double f_hist[], double x,
                      double *y1, double *y2, double h, int number_of_steps ) {

   static const double one_twelfth = 1.0 / 12.0;
   double y;
   int i;


   for (i = 0; i < number_of_steps; i++) {
      x += h;
      y = *y1;
      *y1 = *y2;
      *y2 = *y1 + *y1 - y + h * h * ( f_hist[2] +
              one_twelfth * (f_hist[2] - f_hist[1] - f_hist[1] + f_hist[0]) );
      f_hist[0] = f_hist[1];
      f_hist[1] = f_hist[2];
      f_hist[2] = f(x+h+h, *y2);
   }
   return y;
}


////////////////////////////////////////////////////////////////////////////////
//  static void Backward_Start( double (*f)(double, double), double f0,       //
//                              double x0, double y0, double c, double *y1,   //
//                              double *y2, double f_hist[], double h  )      //
//                                                                            //
//  Description:                                                              //
//     The backward differential correction procedure for approximating the   //
//     solution of the second order differential equation y''(x) = f(x,y)     //
//     with initial conditions y(x0) = y[0] and y'(x0) = c requires three     //
//     consecutive values, y[0], y[1], and y[2], to recursively generate y[n] //
//     for n >= 3.  This routine uses the Runge-Kutta method to approximate   //
//     y1 = y[1] and y2 = y[2], and to calculate f_hist[0] = f(x0,y0),        //
//     f_hist[1] = f(x0+h,y1) and f_hist[2] = f(x0+2h,y2).                    //
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
//     double *y1 Output, the initial estimate of y(x0+h).                    //
//     double *y2 Output, the initial estimate of y(x0+2h).                   //
//     double f_hist[]         Output, estimates of {f(x0,y0), f(x+h,y1),     //
//                f(x+2h,y2) }.                                               //
//     double h   The step size.                                              //
//                                                                            //
//  Return Values:                                                            //
//     The routine is of type void and does not return a value.  The values   //
//     y1, y2, and f_hist[] are returned via the addresses in the argument    //
//     list.                                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Backward_Start( double (*f)(double, double), double f0, double x0,
    double y0, double c, double *y1, double *y2, double f_hist[], double h  ) {

   static const double one_sixth = 1.0 / 6.0;
   double h2 = 0.5 * h;
   double k1 = h * f0;
   double k2 = h * (*f)( x0 + h2, y0 + h2 * c );
   double k3 = h * (*f)( x0 + h2, y0 + h2 * c + 0.25 * k1 * h );
   double k4 = h * (*f)( x0 + h, y0 + h * c + h2 * k2 );

   f_hist[0] = f0;

   *y1 = y0 + h * (c + one_sixth * ( k1 + k2 + k3 ) );
   c += one_sixth * ( k1 + k2 + k2 + k3 + k3 + k4 );
   f_hist[1] = (*f)(x0+h, *y1);

   x0 += h;
   k1 = h * f_hist[1];
   k2 = h * (*f)( x0 + h2, *y1 + h2 * c );
   k3 = h * (*f)( x0 + h2, *y1 + h2 * c + 0.25 * k1 * h );
   *y2 = *y1 + h * (c + one_sixth * ( k1 + k2 + k3 ) );
   f_hist[2] = (*f)(x0+h, *y2);
}
