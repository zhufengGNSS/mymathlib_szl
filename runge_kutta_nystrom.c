////////////////////////////////////////////////////////////////////////////////
// File: runge_kutta_nystrom.c                                                //
// Routines:                                                                  //
//    Runge_Kutta_Nystrom                                                     //
//    Runge_Kutta_Nystrom_Richardson                                          //
//    Runge_Kutta_Nystrom_Integral_Curve                                      //
//    Runge_Kutta_Nystrom_Richardson_Integral_Curve                           //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The 6 stage 5th order Runge-Kutta method for solving a differential    //
//     equation y'(x) = f(x,y) with initial condition y = c when x = x0       //
//     evaluates f(x,y) six times per step. For step i+1,                     //
//     y[i+1] = y[i] + ( 23 * k1 + 125 * k3 - 81 * k5 + 125 * k6 ) / 192      //
//     where,                                                                 //
//     k1 = h * f( x[i], y[i] ),                                              //
//     k2 = h * f( x[i]+h/3, y[i]+k1/3 ),                                     //
//     k3 = h * f( x[i]+2h/5, y[i] + ( 4 k1 + 6 k2 ) / 25 ),                  //
//     k4 = h * f( x[i]+h, y[i] + ( k1 - 12 k2 + 15 k3 ) / 4),                //
//     k5 = h * f( x[i]+2h/3, y[i] + (6 k1 + 90 k2 - 50 k3 + 8 k4) / 81),     //
//     k6 = h * f( x[i]+4h/5, y[i] + (6 k1 + 36 k2 + 10 k3 + 8 k4) / 75),     //
//     and x[i+1] = x[i] + h.                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>

static const double one_twentififth = 1.0 / 25.0;
static const double one_fourth = 1.0 / 4.0;
static const double one_eightyfirst = 1.0 / 81.0;
static const double one_seventyfifth = 1.0 / 75.0;
static const double one_third = 1.0 / 3.0;
static const double two_fifths = 2.0 / 5.0;
static const double two_thirds = 2.0 / 3.0;
static const double four_fifths = 4.0 / 5.0;
static const double one_oneninetytwo = 1.0 / 192.0;


////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Nystrom( double (*f)(double, double), double y0,       //
//                               double x0, double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Nystrom method described above to    //
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
double Runge_Kutta_Nystrom( double (*f)(double, double), double y0, double x0,
                                             double h, int number_of_steps ) {

   double k1, k2, k3, k4, k5, k6;
   double h3 = one_third * h;
   double h4 = one_fourth * h;
   double h25 = one_twentififth * h;
   double h81 = one_eightyfirst * h;
   double h2_5 = two_fifths * h;
   double h2_3 = two_thirds * h;
   double h4_5 = four_fifths * h;
   double h75 = one_seventyfifth * h;
   double h192 = one_oneninetytwo * h;

   while ( --number_of_steps >= 0 ) {
      k1 = (*f)(x0,y0);
      k2 = (*f)(x0 + h3, y0 + h3 * k1);
      k3 = (*f)(x0 + h2_5, y0
                             +  h25 * ( 4.0 * k1 + 6.0 * k2 ) );
      k4 = (*f)(x0 + h, y0 + h4 * ( k1 - 12.0 * k2 + 15.0 * k3 ) );
      k5 = (*f)(x0 + h2_3, y0 
                     + h81 * ( 6.0 * k1 + 90.0 * k2 - 50.0 * k3 + 8.0 * k4 ) );
      k6 = (*f)(x0 + h4_5, y0 + h75 * ( 6.0 * k1 + 36.0 * k2 + 10.0 * k3
                                                                + 8.0 * k4 ) );
      y0 += h192 * ( 23.0 * k1 + 125.0 * k3 - 81.0 * k5 + 125.0 * k6 );
      x0 += h;
   }
   return y0;
}


////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Nystrom_Richardson( double (*f)(double, double),       //
//                      double y0, double x0, double h, int number_of_steps,  //
//                                                   int richarson_columns);  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Nystrom method described above       //
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
  1.0 / 31.0, 1.0 / 63.0, 1.0 / 127.0, 1.0 / 255.0, 1.0 / 511.0, 1.0 / 1023.0
};

#define MAX_COLUMNS 1+sizeof(richardson)/sizeof(richardson[0])

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Runge_Kutta_Nystrom_Richardson( double (*f)(double, double), double y0,
           double x0, double h, int number_of_steps, int richardson_columns ) {

   double dt[MAX_COLUMNS];         // dt[i] is the last element in column i.
   double integral, delta, h_used;
   int i,j,k, number_sub_intervals;

   richardson_columns = max(1, min(MAX_COLUMNS, richardson_columns));
   while ( --number_of_steps >= 0 ) {
      h_used = h;
      number_sub_intervals = 1;
      for (j = 0; j < richardson_columns; j++) {
         integral = Runge_Kutta_Nystrom( f, y0, x0, h_used,
                                                          number_sub_intervals);
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
//  void Runge_Kutta_Nystrom_Integral_Curve( double (*f)(double, double),     //
//        double y[], double x0, double h, int number_of_steps_per_interval,  //
//                                                 int number_of_intervals ); //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Nystrom method described above to    //
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
void Runge_Kutta_Nystrom_Integral_Curve( double (*f)(double, double),
            double y[], double x0, double h, int number_of_steps_per_interval,
                                                    int number_of_intervals ) {

   double k1, k2, k3, k4, k5, k6;
   double h3 = one_third * h;
   double h4 = one_fourth * h;
   double h25 = one_twentififth * h;
   double h81 = one_eightyfirst * h;
   double h2_5 = two_fifths * h;
   double h2_3 = two_thirds * h;
   double h4_5 = four_fifths * h;
   double h75 = one_seventyfifth * h;
   double h192 = one_oneninetytwo * h;
   int i;

   while ( --number_of_intervals >= 0 ) {
      *(y+1) = *y;
      y++;
      for (i = 0; i < number_of_steps_per_interval; i++) {
         k1 = (*f)(x0,*y);
         k2 = (*f)(x0 + h3, *y + h3 * k1);
         k3 = (*f)(x0 + h2_5, *y
                             +  h25 * ( 4.0 * k1 + 6.0 * k2 ) );
         k4 = (*f)(x0 + h, *y + h4 * ( k1 -12.0 * k2 + 15.0 * k3 ) );
         k5 = (*f)(x0 + h2_3, *y 
                     + h81 * ( 6.0 * k1 + 90.0 * k2 - 50.0 * k3 + 8.0 * k4 ) );
         k6 = (*f)(x0 + h4_5, *y + h75 * ( 6.0 * k1 + 36.0 * k2 + 10.0 * k3
                                                                + 8.0 * k4 ) );
         *y += h192 * ( 23.0 * k1 + 125.0 * k3 - 81.0 * k5 + 125.0 * k6 );
         x0 += h;
      }   
   }
}



////////////////////////////////////////////////////////////////////////////////
//  void Runge_Kutta_Nystrom_Richardson_Integral_Curve(                       //
//         double (*f)(double, double),  double y[], double x0, double h,     //
//         int number_of_steps_per_interval, int number_of_intervals,         //
//                                                  int richardson_columns )  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Nystrom method described above       //
//     together with Richardson extrapolation to the limit (as h -> 0) to     //
//     approximate the solution of the differential equation y'=f(x,y) with   //
//     the initial condition y = y[0] at x = x0.                              //
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
void Runge_Kutta_Nystrom_Richardson_Integral_Curve( double (*f)(double, double),
             double y[], double x0, double h, int number_of_steps_per_interval,
                            int number_of_intervals, int richardson_columns ) {

   double mh = (double) number_of_steps_per_interval * h;

   while ( --number_of_intervals >= 0 ) {
      *(++y) = Runge_Kutta_Nystrom_Richardson( f, *y, x0, h,
                            number_of_steps_per_interval, richardson_columns );
      x0 += mh;
   }
   return;
}
