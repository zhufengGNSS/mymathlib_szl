////////////////////////////////////////////////////////////////////////////////
// File: explicit_central_difference.c                                        //
// Routines:                                                                  //
//    Explicit_Central_Difference_Method                                      //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The explicit central difference method for approximating the solution  //
//     of the second order differential equation y''(x) = f(x,y) with initial //
//     conditions y(x0) = a, y'(x0) = c is a second order method based upon   //
//     the following central difference substitutions for y' and y'':         //
//                    y'(x) = (y(x+h) - y(x-h) ) / 2h, and                    //
//                    y''(x) = (y(x+h) - 2y(x) + y(x-h)) / h^2.               //
//     Note that here f(x,y) does not depend upon y'.                         //
//                                                                            //
//     Letting x[n] = x0 + nh, y[n] = y(x[n]), and f[n] = f(x[n],y[n]), the   //
//     algorithm proceeds recursively as follows:                             //
//                                                                            //
//            y[n+1] = 2y[n] - y[n-1] + h^2 f(x[n],y[n])                      //
//                                                                            //
//     The starting procedure is y[0] = a, y[1] = hc + h^2 f(x[0],y[0]) / 2.  //
//                                                                            //
//     Richardson extrapolation may be used to increase to increase the       //
//     order.                                                                 //
////////////////////////////////////////////////////////////////////////////////

static const double richardson[] = {  
  1.0 / 3.0, 1.0 / 7.0, 1.0 / 15.0, 1.0 / 31.0, 1.0 / 63.0
};

#define MAX_COLUMNS 1+sizeof(richardson)/sizeof(richardson[0])

static double Central_Difference_Method( double (*f)(double, double), double y,
                        double x, double *z, double h, int number_of_steps );

////////////////////////////////////////////////////////////////////////////////
//  void Explicit_Central_Difference_Method( double (*f)(double, double),     //
//              double y[], double x0, double c, double h, int max_columns,   //
//                                                    int number_of_steps )   //
//                                                                            //
//  Description:                                                              //
//     This function solves the second order differential equation y''=f(x,y) //
//     with the initial conditions y(x0) = y0 and y'(x0) = c.  The solution   //
//     is returned in the array y[].                                          //
//     Note that the function f does not depend explicitly on y'.             //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the rate of change of //
//                the slope with respect to x at the point (x,y) of integral  //
//                curve of the second order differential equation y'' = f(x,y)//
//                which passes through the point (x0,y[0]) with initial slope //
//                y'(x0) = c.                                                 //
//     double y[] On input y[0] is the initial value of y at x = x0, on output//
//                for i >= 1,  y[i] is given by the recursion described above //
//     double x0  Initial value of x.                                         //
//     double c   The initial value of y'.                                    //
//     double h   The step size.                                              //
//     int    max_columns      The number of columns used for Richardson      //
//                extrapolation. 1 <= max_columns <= MAX_COLUMNS.             //
//     int    number_of_steps  The number of steps.  y[] should be dimensioned//
//                at least number_of_steps + 1.                               //
//                                                                            //
//  Return Values:                                                            //
//     The type of the function is void and therefore does not return a value.//
//     The approximation to the solution to the initial value problem         //
//     y''(x) = f(x,y,y') where y(x0) = y[0] and y'(x0) = c is stored in the  //
//     input array y[].                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Explicit_Central_Difference_Method( double (*f)(double, double), 
                    double y[], double x0, double c, double h, int max_columns,
                                                        int number_of_steps ) {

   double dt[MAX_COLUMNS];         // dt[i] is the last element in column i.
   double z[MAX_COLUMNS];          // approximation of y'.
   double x = x0;                  // first argument to f.
   double integral;
   double delta;
   double h_old;
   int i,j,k, number_sub_intervals;

   if (max_columns < 1) max_columns = 1;
   if (max_columns > MAX_COLUMNS) max_columns = MAX_COLUMNS;

   integral = 0.5 * h * (*f)(x0,y[0]);
   for (i = 0; i < max_columns; i++) {
      z[i] = c + integral;
      integral *= 0.5;
   }

   for (i = 0; i < number_of_steps; i++) {
      number_sub_intervals = 1;
      h_old = h;
      for (j = 0; j < max_columns; j++) {
         integral = Central_Difference_Method( f, y[i], x, &z[j], h_old,
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
      x += h;
   }
   return;
}


////////////////////////////////////////////////////////////////////////////////
//  static double Central_Difference_Method( double (*f)(double, double),     //
//              double y, double x, double *z, double h, int number_of_steps )//
//                                                                            //
//  Description:                                                              //
// Given the initial condition y(x)=y, y'(x)=z approximate y(x+nh) using the  //
// explicit central difference method.                                        //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the rate of change of //
//                the slope with respect to x at the point (x,y) of integral  //
//                curve of the second order differential equation y'' = f(x,y)////                which passes through the point (x0,y[0]) with initial slope //
//                y'(x0) = c.                                                 //
//     double y   The initial value of y at x.                                //
//     double x   Initial value of x.                                         //
//     double *z  On input, the estimate of y' evaluated at x+h/2.  On output //
//                the estimate of y' evaluated at x+h*number_of_steps + h/2.  //
//     double h   The step size.                                              //
//     int    number_of_steps  The number of steps.  The value returned is    //
//                y(x + h*number_of_step).                                    //
//                                                                            //
//  Return Values:                                                            //
//     The approximation at x+h*number_of_steps of the solution of the initial//
//      value problem y''(x) = f(x,y,y') where y(x0) = y and y'(x0+h/2) = *z. //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static double Central_Difference_Method( double (*f)(double, double), double y,
                        double x, double *z, double h, int number_of_steps ) {

   int i;

   for (i = 0; i < number_of_steps; x += h, i++) {
      y += h * *z;
      *z += h * (*f)(x+h, y); 
   }
   return y;
}
