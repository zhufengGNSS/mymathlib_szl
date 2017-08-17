////////////////////////////////////////////////////////////////////////////////
// File: trapezoidal_method.c                                                 //
// Routines:                                                                  //
//    Trapezoidal_Method                                                      //
//    Trapezoidal_Integral_Curve                                              //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The trapezoidal method is an implicit method for approximating the     //
//     solution of the differential equation y'(x) = f(x,y) with initial      //
//     condition y(a) = c.  The trapezoidal method can be derived by          //
//     expanding y(x+h) in a Taylor series about x,                           //
//       y(x+h) = y(x) + h*y'(x) + (h^2 / 2)*y''(x) + (h^3 / 6)*y'''(x) +...  //
//     and substituting (y'(x+h)-y'(x))/h for y''(x) and y'(x) = f(x,y(x)),   //
//          y(x+h) = y(x) + (h/2) * (f(x+h,y(x+h) - f(x,h)) + O(h^3).         //
//                                                                            //
//     Let y[n] = y(a + nh), x[n] = a + nh, then the recursion formula for    //
//     the trapezoidal method is:                                             //
//          y[n+1] = y[n] + (h/2) * ( f(x[n],y[n]) + f(x[n+1],y[n+1]) )       //
//                                                                            //
//     I.e. y[n+1] is the solution y to:                                      //
//              y - (h/2)*(f(x[n+1],y) = y[n] + (h/2)*f(x[n],y[n]).           //
//                                                                            //
//     Thus locally the trapezoidal method is a third order method  and       //
//     globally a second order method.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  double Trapezoidal_Method( double (*f)(double, double),                   //
//                 double (*g)(double,double,double),  double y0, double x0,  //
//                                          double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the trapezoidal method to approximate the solution   //
//     at x = x0 + h * number_of_steps of the initial value problem y'=f(x,y),//
//     y(x0) = y0.                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            passes through the point (x0,y0).                               //
//     double *g                                                              //
//            Pointer to the function g(x,h,u) which returns the value y      //
//            such that y - h *  f(x,y) = u.                                  //
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
double Trapezoidal_Method( double (*f)(double, double),
                     double (*g)(double,double,double), double y0, double x0,
                                             double h, int number_of_steps ) {

   double u;
   double h2 = 0.5 * h;

   while ( --number_of_steps >= 0 ) { 
      u = y0 + h2 * f(x0,y0);
      x0 += h;
      y0 = g(x0, h2, u);
   }

   return y0;
}


////////////////////////////////////////////////////////////////////////////////
//  void Trapezoidal_Integral_Curve( double (*f)(double, double),             //
//                double (*g)(double,double,double),  double y[], double x0,  //
//                double h, int number_of_steps_per_interval,                 //
//                                                int number_of_intervals );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the trapezoidal method to approximate the solution   //
//     of the differential equation y'=f(x,y) with the initial condition      //
//     y = y[0] at x = x0.  The values are returned in y[n] which is the      //
//     value of y evaluated at x = x0 + n * m * h, where m is the number of   //
//     steps per interval and n is the interval number,                       //
//     0 <= n <= number_of_intervals.                                         //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            which passes through the point (x0,y[0]).                       //
//     double *g                                                              //
//            Pointer to the function g(x,h,u) which returns the value y      //
//            such that y - h/2 f(x,y) = u.                                   //
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
void Trapezoidal_Integral_Curve( double (*f)(double, double), 
          double (*g)(double,double,double), double y[], double x0, double h,
                int number_of_steps_per_interval, int number_of_intervals ) {

   double u;
   double h2 = 0.5 * h;
   int i;

   while ( --number_of_intervals >= 0 ) {
      *(y+1) = *y;
      y++;
      for (i = 0; i < number_of_steps_per_interval; i++) {
         u = *y + h2 * f(x0,*y);
         x0 += h;
         *y = g(x0, h2, u);
      }
   }
}
