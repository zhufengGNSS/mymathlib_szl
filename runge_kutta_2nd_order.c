////////////////////////////////////////////////////////////////////////////////
// File: runge_kutta_2nd_order.c                                              //
// Routines:                                                                  //
//    Runge_Kutta_2nd_Order                                                   //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     A differential equation of order greater than one can be transformed   //
//     to a multidimensional first order differential equation.  In           //
//     particular, given a second order differential equation y'' = f(x,y,y') //
//     define z = y', then the pair (y,z)' = (z, f(x,y,z) ) is a first order  //
//     differential equation.  The initial conditions y(x0) = y[0],           //
//     y'(x0) = c becomes the initial condition (y,z)(x0) = (y[0],c).         //
//     when x = x0 numerically evaluates f(x,y) four times per step. Then for //
//     step i+1,                                                              //
//           y[i+1] = y[i] + h * ( yp[i] + 1/6 (k1 + k2 + k3 ) and            //
//           yp[i+1] = yp[i] + 1/6 * ( k1 + 2 k2 + 2 k3 + k4 ) where          //
//         k1 = h * f(x[i], y[i], yp[i]),                                     //
//         k2 = h * f(x[i]+h/2, y[i]+(h/2)*yp[i], yp[i]+k1/2),                //
//         k3 = h * f(x[i]+h/2, y[i]+(h/2)*yp[i]+(h/4)*k1, yp[i]+k2/2),       //
//     x[i] = x0 + i * h.                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static const double one_sixth = 1.0 / 6.0;

////////////////////////////////////////////////////////////////////////////////
//  void Runge_Kutta_2nd_Order( double (*f)(double, double, double),          //
//          double x0, double y[], double c, double h, int number_of_steps ); //
//                                                                            //
//  Description:                                                              //
//     This function solves the 2nd order differential equation y''=f(x,y,y') //
//     with the initial conditions y=y[0] and y'= c at x = x0.  The values    //
//     returned are y[n] which is the value of y evaluated at x = x0 + n*h.   //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the rate of change of //
//                the slope with respect to x at the point (x,y) with slope   //
//                y'(x) of integral curve of the second order differential    //
//                equation y'' = f(x,y,y') which passes through the point     //
//                (x0,y[0]) with initial slope y'(x0) = c.                    //
//     double y[] On input y[0] is the initial value of y at x = x0, on output//
//                for i >= 1,  y[i] is given by the algorithm described above.//
//     double x0  Initial value of x.                                         //
//     double h   Step size                                                   //
//     int    number_of_steps  The number of steps, y[] should be dimensioned //
//                             number_of_steps + 1.                           //
//                                                                            //
//  Return Values:                                                            //
//     This routine is of type void and hence does not return a value.  The   //
//     solution of y'' = f(x,y,y') from x = x0 to x = x0 + number_of_steps*h  //
//     is stored in the input array y[].                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Runge_Kutta_2nd_Order( double (*f)(double, double, double), double x0,
                        double y[], double c, double h, int number_of_steps ) {

   double k1, k2, k3, k4;
   double h2 = 0.5 * h;
   double ych2;
   int i;

   for (i = 0; i < number_of_steps; x0 += h, i++) {
      ych2 = y[i] + c * h2;
      k1 = h * (*f)(x0, y[i], c);
      k2 = h * (*f)(x0+h2, ych2, c + 0.5 * k1);
      k3 = h * (*f)(x0+h2, ych2 + 0.25 * k1 * h, c + 0.5 * k2);
      k4 = h * (*f)(x0+h, y[i] + c * h + h2 * k2 , c + k3);
      y[i+1] = y[i] + ( c + one_sixth * (k1 + k2 + k3) ) * h;
      c += one_sixth * ( k1 + k2 + k2 + k3 + k3 + k4 );
   }
}
