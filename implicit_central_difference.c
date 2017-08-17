////////////////////////////////////////////////////////////////////////////////
// File: implicit_central_difference.c                                        //
// Routines:                                                                  //
//    Implicit_Central_Difference_Method                                      //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The implicit central difference method for approximating the solution  //
//     of the second order differential equation y''(x) = f(x,y,y') with      //
//     initial conditions y(x0) = a, y'(x0) = c is a second order method      //
//     based upon the following central difference substitutions for y' and   //
//     y'':                                                                   //
//                    y'(x) = (y(x+h) - y(x-h) ) / 2h, and                    //
//                    y''(x) = (y(x+h) - 2y(x) + y(x-h)) / h^2.               //
//                                                                            //
//     Letting x[n] = x0 + nh, y[n] = y(x[n]), and                            //
//                 f[n] = f( x[n], y[n], (y[n+1]-y[n-1])/2h ), then the       //
//     algorithm proceeds recursively as follows:                             //
//                                                                            //
//                      y[n+1] = 2y[n] - y[n-1] + h^2 f[n]                    //
//                                                                            //
//     for which the value y[n+1] is given implicitly.                        //
//     The starting procedure is y[0] = a, y[1] = hc + h^2 f(x0,a,c) / 2.     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  void Implicit_Central_Difference_Method( double f0,                       //
//          double (*g)(double,double,double,double), double y[], double x0,  //
//                                 double c, double h, int number_of_steps ); //
//                                                                            //
//  Description:                                                              //
//     This function solves the differential equation y'=f(x,y,y') with the   //
//     initial condition y(x0) = y0 and y'(x0) = c.  The solution is returned //
//     in the array y[].                                                      //
//                                                                            //
//  Arguments:                                                                //
//     double f0  The value of f(x,y,y') evaluated at the initial values,     //
//     i.e. f0 = f(x0,y[0],c).                                                //
//     double *g  Pointer to the user-supplied function g(x,y[n],y[n-1],h)    //
//                which returns the solution y[n+1] given implicitly by the   //
//                recursion y[n+1] = 2*y[n] - y[n-1]                          //
//                                  + h^2*f(x[n],y[n],(y[n+1)-y[n-1])/(2*h)). //
//     double y[] On input y[0] is the initial value of y at x = x0, on output//
//                for i >= 1,  y[i] is given by the recursion described above //
//     double x0  Initial value of x.                                         //
//     double c   The initial value of y'.                                    //
//     double h   The step size.                                              //
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
void Implicit_Central_Difference_Method( double f0,
          double (*g)(double,double,double,double), double y[], double x0,
                                    double c, double h, int number_of_steps ) {

   int i;

   y[1] = y[0] + h * ( c + 0.5 * h * f0 );

   for (i = 1; i < number_of_steps; x0 += h, i++) {
      y[i+1] = (*g)(x0+h, y[i], y[i-1], h);
   }
   return;
}
