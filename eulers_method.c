////////////////////////////////////////////////////////////////////////////////
// File: eulers_method.c                                                      //
// Routines:                                                                  //
//    Eulers_Method                                                           //
//    Eulers_Method_Richardson                                                //
//    Euler_Integral_Curve                                                    //
//    Euler_Richardson_Integral_Curve                                         //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     Euler's method for approximating the solution of the differential      //
//     equation y'(x) = f(x,y) with initial condition y(a) = c numerically    //
//     approximates the derivative y'(x) by ( y(x+h) - y(x) ) / h.            //
//                                                                            //
//     Letting y[n] = y(a + nh), x[n] = a + nh, and f[n] = f(x[n],y[n]), the  //
//     recursion formula for the Euler method is:                             //
//                          y[n+1] = y[n] + h * f[n].                         //
//                                                                            //
//     Expand y(x+h) in a Taylor series about x,                              //
//       y(x+h) = y(x) + h*y'(x) + (h^2 / 2)*y''(x) + (h^3 / 6)*y'''(x) +...  //
//     and substitute f(x,y(x)) for y'(x)                                     //
//       y(x+h)-y(x)-h*f(x,y(x)) = (h^2 / 2) * y''(e), where x <= e <= x+h.   //
//                                                                            //
//     Thus locally Euler's method is a second order and globally a first     //
//     order procedure.                                                       //
//                                                                            //
//     Richardson extrapolation can be used to increase the order and         //
//     accuracy.                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  double Eulers_Method( double (*f)(double, double), double y0, double x0,  //
//                                          double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Euler's method to approximate the solution at x =    //
//     x0 + h * number_of_steps of the initial value problem y'=f(x,y),       //
//     y(x0) = y0.                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            passes through the point (x0,y0).                               //
//     double y0                                                              //
//            The initial value of y at x = x0.                               //
//     double x0                                                              //
//            The initial value of x.                                         //
//     double h                                                               //
//            The step size.                                                  //
//     int    number_of_steps                                                 //
//            The number of steps. Must be a nonnegative integer.             //
//                                                                            //
//  Return Values:                                                            //
//     The solution of the initial value problem y' = f(x,y), y(x0) = y0 at   //
//     x = x0 + number_of_steps * h.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Eulers_Method( double (*f)(double, double), double y0, double x0,
                                             double h, int number_of_steps ) {

   while ( --number_of_steps >= 0 ) { 
      y0 +=  h * (*f)( x0, y0 );
      x0 += h;
   }

   return y0;
}


////////////////////////////////////////////////////////////////////////////////
//  double Eulers_Method_Richardson( double (*f)(double, double), double y0,  //
//         double x0, double h, int number_of_steps, int richardson_columns)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Euler's method together with Richardson extrapolation//
//     to approximate the solution at x = x0 + h * number_of_steps of the     //
//     initial value problem y'=f(x,y), y(x0) = y0.                           //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            passes through the point (x0,y0).                               //
//     double y0                                                              //
//            The initial value of y at x = x0.                               //
//     double x0                                                              //
//            The initial value of x.                                         //
//     double h                                                               //
//            The step size.                                                  //
//     int    number_of_steps                                                 //
//            The number of steps. Must be nonnegative.                       //
//     int    richardson_columns                                              //
//            The maximum number of columns to use in the Richardson          //
//            extrapolation to the limit.                                     //
//                                                                            //
//  Return Values:                                                            //
//     The solution of the initial value problem y' = f(x,y), y(x0) = y0 at   //
//     x = x0 + number_of_steps * h.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static const double richardson[] = {  
  1.0, 1.0 / 3.0, 1.0 / 7.0, 1.0 / 15.0, 1.0 / 31.0, 1.0 / 63.0, 1.0 / 127.0
};

#define MAX_COLUMNS 1+sizeof(richardson)/sizeof(richardson[0])

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Eulers_Method_Richardson( double (*f)(double, double), double y0,
           double x0, double h, int number_of_steps, int richardson_columns ) {

   double dt[MAX_COLUMNS];         // dt[i] is the last element in column i.
   double integral, delta, h_used;
   int j,k, number_sub_intervals;

   richardson_columns = max(1, min(MAX_COLUMNS, richardson_columns));
   while ( --number_of_steps >= 0 ) {
      h_used = h;
      number_sub_intervals = 1;
      for (j = 0; j < richardson_columns; j++) {
         integral = Eulers_Method( f, y0, x0, h_used, number_sub_intervals);
         for ( k = 0; k < j; k++) {
            delta = integral - dt[k];
            dt[k] = integral;
            integral += richardson[k] * delta;
         }
         dt[j] = integral;
         h_used *= 0.5;
         number_sub_intervals += number_sub_intervals;
      }
      y0 = integral;
      x0 += h;
   }
   return y0;
}


////////////////////////////////////////////////////////////////////////////////
//  void Euler_Integral_Curve( double (*f)(double, double), double y[],       //
//                    double x0, double h, int number_of_steps_per_interval,  //
//                    int number_of_intervals );                              //
//                                                                            //
//  Description:                                                              //
//     This routine uses Euler's method to approximate the solution of the    //
//     differential equation y'=f(x,y) with the initial condition y = y[0] at //
//     x = x0.  The values are returned in y[n] which is the value of y       //
//     evaluated at x = x0 + n * m * h, where m is the number of steps per    //
//     interval and n is the interval number, 0 <= n <= number_of_intervals.  //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            which passes through the point (x0,y[0]).                       //
//     double y[]                                                             //
//            On input y[0] is the initial value of y at x = x0. On output    //
//            y[n+1] = y[n] + m*h*f(a + n*m*h, y[n] ), where m is the number  //
//            of steps per interval and n is the interval number.             //
//     double x0                                                              //
//            Initial value of x.                                             //
//     double h                                                               //
//            Step size                                                       //
//     int    number_of_steps_per_interval                                    //
//            The number of steps of length h used to calculate y[i+1]        //
//            starting from y[i].                                             //
//     int    number_of_intervals                                             //
//            The number of intervals,  y[] should be dimensioned at least    //
//            number_of_intervals + 1.                                        //
//                                                                            //
//  Return Values:                                                            //
//     This routine is of type void and hence does not return a value.        //
//     The solution of y' = f(x,y) from x = x0 to x = x0 + n * m * h,         //
//     where n is the number of intervals and m is the number of steps per    //
//     interval, is stored in the input array y[].                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Euler_Integral_Curve( double (*f)(double, double), double y[], double x0,
      double h, int number_of_steps_per_interval, int number_of_intervals ) {

   int i;

   while ( --number_of_intervals >= 0 ) {
      *(y+1) = *y;
      y++;
      for (i = 0; i < number_of_steps_per_interval; x0 += h, i++) 
         *y += h * (*f)( x0, *y );
   }
}


////////////////////////////////////////////////////////////////////////////////
//  void Euler_Richardson_Integral_Curve( double (*f)(double, double),        //
//        double y[], double x0, double h, int number_of_steps_per_interval,  //
//                         int number_of_intervals, int richardson_columns )  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Euler's method together with Richardson extrapolation//
//     to the limit (as h -> 0) to approximate the solution of the differen-  //
//     tial equation y'=f(x,y) with the initial condition y = y[0] at x = x0. //
//     The values are returned in y[], in which y[n] is the value of y        //
//     evaluated at x = x0 + n * m * h, where m is the number of steps per    //
//     interval and n is the interval number, 0 <= n <= number_of_intervals.  //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            which passes through the point (x0,y[0]).                       //
//     double y[]                                                             //
//            On input y[0] is the initial value of y at x = x0. On output    //
//            y[n+1] = y[n] + m*h*f(a + n*m*h, y[n] ), where m is the number  //
//            of steps per interval and n is the interval number.             //
//     double x0                                                              //
//            Initial value of x.                                             //
//     double h                                                               //
//            Step size                                                       //
//     int    number_of_steps_per_interval                                    //
//            The number of steps of length h used to calculate y[i+1]        //
//            starting from y[i].                                             //
//     int    number_of_intervals                                             //
//            The number of intervals,  y[] should be dimensioned at least    //
//            number_of_intervals + 1.                                        //
//     int    richardson_columns                                              //
//            The maximum number of columns to use in the Richardson          //
//            extrapolation to the limit.                                     //
//                                                                            //
//  Return Values:                                                            //
//     This routine is of type void and hence does not return a value.        //
//     The solution of y' = f(x,y) from x = x0 to x = x0 + n * m * h,         //
//     where n is the number of intervals and m is the number of steps per    //
//     interval, is stored in the input array y[].                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Euler_Richardson_Integral_Curve( double (*f)(double, double), double y[],
                         double x0, double h, int number_of_steps_per_interval,
                            int number_of_intervals, int richardson_columns ) {

   double mh = (double) number_of_steps_per_interval * h;

   while ( --number_of_intervals >= 0 ) {
      *(++y) = Eulers_Method_Richardson( f, *y, x0, h,
                            number_of_steps_per_interval, richardson_columns );
      x0 += mh;
   }
}
