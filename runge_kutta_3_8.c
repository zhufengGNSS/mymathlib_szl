////////////////////////////////////////////////////////////////////////////////
// File: runge_kutta_3_8.c                                                    //
// Routines:                                                                  //
//    Runge_Kutta_3_8                                                         //
//    Runge_Kutta_3_8_Richardson                                              //
//    Runge_Kutta_3_8_Integral_Curve                                          //
//    Runge_Kutta_3_8_Richardson_Integral_Curve                               //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     A Runge-Kutta method based upon the 3/8ths quadrature rule used for    //
//     approximating the solution of the differential equation y'(x) = f(x,y) //
//     with initial condition y = c when x = x0 numerically evaluates f(x,y)  //
//     four times per step.                                                   //
//     For step i+1,                                                          //
//          y[i+1] = y[i] + 1/8 * ( k1 + 3 * k2 + 3 * k3 + k4 ) where         //
//     k1 = h * f(x[i],y[i]),                                                 //
//     k2 = h * f(x[i]+h/3,y[i]+k1/3),                                        //
//     k3 = h * f(x[i]+2h/3, y[i] + (k2 - k1 / 3) ),                          //
//     k4 = h * f(x[i]+h, y[i]+k1-k2+k3),                                     //
//     x[i] = x0 + i * h.                                                     //
//                                                                            //
//     This version of the Runge-Kutta method is fourth order.                //
//     Richardson extrapolation can be used to increase the order and         //
//     accuracy.                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static const double one_eighth = 1.0 / 8.0;
static const double one_third = 1.0 / 3.0;
static const double two_thirds = 2.0 / 3.0;

////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_3_8( double (*f)(double, double), double y0, double x0,//
//                                          double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the 4th order Runge-Kutta method described above to  //
//     approximate the solution at x = x0 + h * number_of_steps of the initial//
//     value problem y'=f(x,y), y(x0) = y0.                                   //
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
double Runge_Kutta_3_8( double (*f)(double, double), double y0, double x0,
                                             double h, int number_of_steps ) {

   double k1, k2, k3, k4;
   double h13 = one_third * h;
   double h23 = two_thirds * h;
   double h18 = one_eighth * h;

   while ( --number_of_steps >= 0 ) {
      k1 =  (*f)(x0,y0);
      k2 =  (*f)(x0+h13, y0 + h13 * k1);
      k3 =  (*f)(x0+h23, y0 + h * (k2 - one_third * k1) );
      x0 += h;
      k4 = (*f)(x0, y0 + h * ( k1 - k2 + k3 ) );
      y0 += h18 * ( k1 + 3.0 * k2 + 3.0 * k3 + k4 );
   }

   return y0;
}


////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_3_8_Richardson( double (*f)(double, double), double y0,//
//         double x0, double h, int number_of_steps, int richardson_columns)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the 4th order Runge-Kutta method described above     //
//     together with Richardson extrapolation to approximate the solution at  //
//     x = x0 + h * number_of_steps of the initial value problem y'=f(x,y),   //
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
  1.0 / 15.0, 1.0 / 31.0, 1.0 / 63.0, 1.0 / 127.0, 1.0 / 255.0, 1.0 / 511.0
};

#define MAX_COLUMNS 1+sizeof(richardson)/sizeof(richardson[0])

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Runge_Kutta_3_8_Richardson( double (*f)(double, double), double y0,
           double x0, double h, int number_of_steps, int richardson_columns ) {

   double dt[MAX_COLUMNS];         // dt[i] is the last element in column i.
   double integral, delta, h_used;
   int j,k, number_sub_intervals;

   richardson_columns = max(1, min(MAX_COLUMNS, richardson_columns));
   while ( --number_of_steps >= 0 ) {
      h_used = h;
      number_sub_intervals = 1;
      for (j = 0; j < richardson_columns; j++) {
         integral = Runge_Kutta_3_8( f, y0, x0, h_used, number_sub_intervals);
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
//  void Runge_Kutta_3_8_Integral_Curve( double (*f)(double, double),         //
//        double y[], double x0, double h, int number_of_steps_per_interval,  //
//                                                int number_of_intervals );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the 4th order Runge-Kutta method described above to  //
//     approximate the solution of the differential equation y'=f(x,y) with   //
//     the initial condition y = y[0] at x = x0.  The values are returned in  //
//     y[n] which is the value of y evaluated at x = x0 + n * m * h, where m  //
//     is the number of steps per interval and n is the interval number,      //
//     0 <= n <= number_of_intervals.                                         //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            which passes through the point (x0,y[0]).                       //
//     double y[]                                                             //
//            On input y[0] is the initial value of y at x = x0. On output    //
//            for n >= 1,  y[n] is the appproximation of the solution y(x) of //
//            the initial value problem y' = f(x,y), y(x0) = y[0],  at        //
//            x = x0 + nmh, where m is the number of steps per interval and n //
//            is the interval number.                                         //
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
void Runge_Kutta_3_8_Integral_Curve( double (*f)(double, double), double y[],
                double x0, double h, int number_of_steps_per_interval,
                                                    int number_of_intervals ) {

   double k1, k2, k3, k4;
   double h13 = one_third * h;
   double h23 = two_thirds * h;
   double h18 = one_eighth * h;
   int i;

   while ( --number_of_intervals >= 0 ) {
      *(y+1) = *y;
      y++;
      for (i = 0; i < number_of_steps_per_interval; i++) {
         k1 = (*f)(x0,*y);
         k2 = (*f)(x0+h13, *y + h13 * k1);
         k3 = (*f)(x0+h23, *y + h * (k2 - one_third * k1) );
         x0 += h;
         k4 = (*f)(x0, *y + h * (k1 - k2 + k3));
         *y += h18 * ( k1 + 3.0 * k2 + 3.0 * k3 + k4 );
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
//  void Runge_Kutta_3_8_Richardson_Integral_Curve(                           //
//             double (*f)(double, double),  double y[], double x0, double h, //
//                int number_of_steps_per_interval,  int number_of_intervals, //
//                                                  int richardson_columns )  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the 4th order Runge-Kutta method described above     //
//     together with Richardson extrapolation to the limit (as h -> 0) to     //
//     approximate the solution of the differential equation y'=f(x,y) with   //
//     the initial condition y = y[0] at x = x0.                              //
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
//            for n >= 1,  y[n] is the appproximation of the solution y(x) of //
//            the initial value problem y' = f(x,y), y(x0) = y[0],  at        //
//            x = x0 + nmh, where m is the number of steps per interval and n //
//            is the interval number.                                         //
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
void Runge_Kutta_3_8_Richardson_Integral_Curve( double (*f)(double, double),
             double y[], double x0, double h, int number_of_steps_per_interval,
                            int number_of_intervals, int richardson_columns ) {

   double mh = (double) number_of_steps_per_interval * h;

   while ( --number_of_intervals >= 0 ) {
      *(++y) = Runge_Kutta_3_8_Richardson( f, *y, x0, h,
                            number_of_steps_per_interval, richardson_columns );
      x0 += mh;
   }
   return;
}
