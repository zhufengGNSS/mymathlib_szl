////////////////////////////////////////////////////////////////////////////////
// File: runge_kutta_verner.c                                                 //
// Routines:                                                                  //
//    Runge_Kutta_Verner                                                      //
//    Runge_Kutta_Verner_Richardson                                           //
//    Runge_Kutta_Verner_Integral_Curve                                       //
//    Runge_Kutta_Verner_Richardson_Integral_Curve                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The 11 stage 8th order Runge-Kutta method for solving a differential   //
//     equation y'(x) = f(x,y) with initial condition y = c when x = x0       //
//     evaluates f(x,y) eleven times per step. For step i+1,                  //
//     y[i+1] = y[i] + ( 9*k1 + 49*k8 + 64*k9 + 49*k10 + 9*k11 ) / 180        //
//     where, after setting s21 = sqrt(21),                                   //
//     k1 = h * f( x[i], y[i] ),                                              //
//     k2 = h * f( x[i]+h/2, y[i]+k1/2 ),                                     //
//     k3 = h * f( x[i]+h/2, y[i]+k1/4+k2/4 ),                                //
//     k4 = h * f( x[i]+(7+s21)/14*h, y[i]+k1/7+(-7-3*s21)/98*k2+             //
//                                                        (21+5*s21)/49*k3 )  //
//     k5 = h * f( x[i]+(7+s21)/14*h, y[i]+(11+s21)/84*k1+ (18+4*s21)/63*k3+  //
//                                                         (21-s21)/252*k4 )  //
//     k6 = h * f( x[i]+h/2, y[i]+(5+s21)/48*k1+ (9+s21)/36*k3+               //
//                                   (-231+14*s21)/360*k4+(63-7*s21)/80*k5 )  //
//     k7 = h * f( x[i]+(7-s21)/14*h, y[i]+(10-s21)/42*k1+                    //
//                                  (-432+92*s21)/315*k3+(633-145s21)/90*k4+  //
//                                   (-504+115*s21)/70*k5+(63-13s21)/35*k6 )  //
//     k8 = h * f( x[i]+(7-s21)/14*h, y[i]+k1/14+(14-3*s21)/126*k5+           //
//                                                   (13-3*s21)/63*k6+k7/9 )  //
//     k9 = h * f( x[i]+h/2, y[i]+(k1/32+(91-21*s21)/576*k5+                  //
//                     (11/72)*k6+(-385-75*s21)/1152*k7+(63+13*s21)/128*k8 )  //
//     k10 = h * f( x[i]+(7+s21)/14*h, y[i]+(k1/14+k5/9+                      //
//                              (-733-147*s21)/2205*k6+(515+111*s21)/504*k7+  //
//                                  (-51-11*s21)/56*k8+(132+28*s21)/245*k9 )  //
//     k11 = h * f( x[i]+h, y[i]+(-42+7*s21)/18*k5+(-18+28s21)/45*k6          //
//                                   +(-273-53*s21)/72*k7+(301+53s21)/72*k8+  //
//                                      (28-28*s21)/45*k9+(49-7*s21)/18*k10 ) //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>

#define sqrt21 4.58257569495584000680

static const double c1 = 1.0 / 2.0;    
static const double c2 = (7.0 + sqrt21 ) / 14.0;
static const double c3 = (7.0 - sqrt21 ) / 14.0;

static const double a21 =  1.0 / 2.0;
static const double a31 =  1.0 / 4.0;
static const double a32 =  1.0 / 4.0;
static const double a41 =  1.0 / 7.0;
static const double a42 = -(7.0 + 3.0 * sqrt21) / 98.0;
static const double a43 =  (21.0 + 5.0 * sqrt21) / 49.0;
static const double a51 =  (11.0 + sqrt21) / 84.0;
static const double a53 =  (18.0 + 4.0 * sqrt21) / 63.0;
static const double a54 =  (21.0 - sqrt21) / 252.0;
static const double a61 =  (5.0 + sqrt21) / 48.0;
static const double a63 =  (9.0 + sqrt21) / 36.0;
static const double a64 =  (-231.0 + 14.0 * sqrt21) / 360.0;
static const double a65 =  (63.0 - 7.0 * sqrt21) / 80.0;
static const double a71 =  (10.0 - sqrt21) / 42.0;
static const double a73 =  (-432.0 + 92.0 * sqrt21) / 315.0;
static const double a74 =  (633.0 - 145.0 * sqrt21) / 90.0;
static const double a75 =  (-504.0 + 115.0 * sqrt21) / 70.0;
static const double a76 =  (63.0 - 13.0 * sqrt21) / 35.0;
static const double a81 =  1.0 / 14.0;
static const double a85 =  (14.0 - 3.0 * sqrt21) / 126.0;
static const double a86 =  (13.0 - 3.0 * sqrt21) / 63.0;
static const double a87 =  1.0 / 9.0;
static const double a91 =  1.0 / 32.0;
static const double a95 =  (91.0 - 21.0 * sqrt21) / 576.0;
static const double a96 =  11.0 / 72.0;
static const double a97 = -(385.0 + 75.0 * sqrt21) / 1152.0;
static const double a98 =  (63.0 + 13.0 * sqrt21) / 128.0;
static const double a10_1 =  1.0 / 14.0;
static const double a10_5 =  1.0 / 9.0;
static const double a10_6 = -(733.0 + 147.0 * sqrt21) / 2205.0;
static const double a10_7 =  (515.0 + 111.0 * sqrt21) / 504.0;
static const double a10_8 = -(51.0 + 11.0 * sqrt21) / 56.0;
static const double a10_9 =  (132.0 + 28.0 * sqrt21) / 245.0;
static const double a11_5 = (-42.0 + 7.0 * sqrt21) / 18.0;
static const double a11_6 = (-18.0 + 28.0 * sqrt21) / 45.0;
static const double a11_7 = -(273.0 + 53.0 * sqrt21) / 72.0;
static const double a11_8 =  (301.0 + 53.0 * sqrt21) / 72.0;
static const double a11_9 =  (28.0 - 28.0 * sqrt21) / 45.0;
static const double a11_10 = (49.0 - 7.0 * sqrt21) / 18.0;

static const double  b1  = 9.0 / 180.0;
static const double  b8  = 49.0 / 180.0;
static const double  b9  = 64.0 / 180.0;

////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Verner( double (*f)(double, double), double y0,        //
//                               double x0, double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta_Verner method described above to     //
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
double Runge_Kutta_Verner( double (*f)(double, double), double y0, double x0,
                                             double h, int number_of_steps ) {

   double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11;
   double c1h = c1 * h, c2h = c2 * h, c3h = c3 * h;

   while ( --number_of_steps >= 0 ) {
      k1 = h * (*f)(x0,y0);
      k2 = h * (*f)(x0 + c1h, y0 + a21 * k1);
      k3 = h * (*f)(x0 + c1h, y0 + ( a31 * k1 + a32 * k2 ) );
      k4 = h * (*f)(x0 + c2h, y0 + ( a41 * k1 + a42 * k2 + a43 * k3 ) );
      k5 = h * (*f)(x0 + c2h, y0 + ( a51 * k1 + a53 * k3 + a54 * k4 ) );
      k6 = h * (*f)(x0 + c1h, y0 + ( a61 * k1 + a63 * k3 + a64 * k4
                                                                + a65 * k5 ) );
      k7 = h * (*f)(x0 + c3h, y0 + ( a71 * k1 + a73 * k3 + a74 * k4
                                                     + a75 * k5 + a76 * k6 ) );
      k8 = h * (*f)(x0 + c3h, y0 + ( a81 * k1 + a85 * k5 + a86 * k6 
                                                                + a87 * k7 ) );
      k9 = h * (*f)(x0 + c1h, y0 + ( a91 * k1 + a95 * k5 + a96 * k6
                                                     + a97 * k7 + a98 * k8 ) );
      k10 = h * (*f)(x0 + c2h, y0 + ( a10_1 * k1 + a10_5 * k5 + a10_6 * k6 
                                    + a10_7 * k7 + a10_8 * k8 + a10_9 * k9 ) );
      x0 += h;
      k11 = h * (*f)(x0, y0 + ( a11_5 * k5 + a11_6 * k6 + a11_7 * k7
                                  + a11_8 * k8 + a11_9 * k9 + a11_10 * k10 ) );
      y0 += (b1 * k1 + b8 * k8 + b9 * k9 + b8 * k10 + b1 * k11);
   }
   return y0;
}


////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Verner_Richardson( double (*f)(double, double),        //
//                      double y0, double x0, double h, int number_of_steps,  //
//                                                   int richarson_columns);  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Verner method described above        //
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
  1.0 / 255.0, 1.0 / 511.0, 1.0 / 1023.0, 1.0 / 2047.0, 1.0 / 4095.0
};

#define MAX_COLUMNS 1+sizeof(richardson)/sizeof(richardson[0])

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Runge_Kutta_Verner_Richardson( double (*f)(double, double), double y0,
           double x0, double h, int number_of_steps, int richardson_columns ) {

   double dt[MAX_COLUMNS];         // dt[i] is the last element in column i.
   double integral, delta, h_used;
   int i,j,k, number_sub_intervals;

   richardson_columns = max(1, min(MAX_COLUMNS, richardson_columns));
   while ( --number_of_steps >= 0 ) {
      h_used = h;
      number_sub_intervals = 1;
      for (j = 0; j < richardson_columns; j++) {
         integral = Runge_Kutta_Verner( f, y0, x0, h_used,
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
//  void Runge_Kutta_Verner_Integral_Curve( double (*f)(double, double),      //
//        double y[], double x0, double h, int number_of_steps_per_interval,  //
//                                                 int number_of_intervals ); //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Verner method described above to     //
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
void Runge_Kutta_Verner_Integral_Curve( double (*f)(double, double),
            double y[], double x0, double h, int number_of_steps_per_interval,
                                                    int number_of_intervals ) {

   double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11;
   double c1h = c1 * h, c2h = c2 * h, c3h = c3 * h;
   int i;

   while ( --number_of_intervals >= 0 ) {
      *(y+1) = *y;
      y++;
      for (i = 0; i < number_of_steps_per_interval; i++) {
         k1 = h * (*f)(x0,*y);
         k2 = h * (*f)(x0 + c1h, *y + a21 * k1);
         k3 = h * (*f)(x0 + c1h, *y + ( a31 * k1 + a32 * k2 ) );
         k4 = h * (*f)(x0 + c2h, *y + ( a41 * k1 + a42 * k2 + a43 * k3 ) );
         k5 = h * (*f)(x0 + c2h, *y + ( a51 * k1 + a53 * k3 + a54 * k4 ) );
         k6 = h * (*f)(x0 + c1h, *y + ( a61 * k1 + a63 * k3 + a64 * k4
                                                                + a65 * k5 ) );
         k7 = h * (*f)(x0 + c3h, *y + ( a71 * k1 + a73 * k3 + a74 * k4
                                                     + a75 * k5 + a76 * k6 ) );
         k8 = h * (*f)(x0 + c3h, *y + ( a81 * k1 + a85 * k5 + a86 * k6 
                                                                + a87 * k7 ) );
         k9 = h * (*f)(x0 + c1h, *y + ( a91 * k1 + a95 * k5 + a96 * k6
                                                     + a97 * k7 + a98 * k8 ) );
         k10 = h * (*f)(x0 + c2h, *y + ( a10_1 * k1 + a10_5 * k5 + a10_6 * k6 
                                    + a10_7 * k7 + a10_8 * k8 + a10_9 * k9 ) );
         x0 += h;
         k11 = h * (*f)(x0, *y + ( a11_5 * k5 + a11_6 * k6 + a11_7 * k7
                                  + a11_8 * k8 + a11_9 * k9 + a11_10 * k10 ) );
         *y += (b1 * k1 + b8 * k8 + b9 * k9 + b8 * k10 + b1 * k11);
      }   
   }
}


////////////////////////////////////////////////////////////////////////////////
//  void Runge_Kutta_Verner_Richardson_Integral_Curve(                        //
//         double (*f)(double, double),  double y[], double x0, double h,     //
//         int number_of_steps_per_interval, int number_of_intervals,         //
//                                                  int richardson_columns )  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Verner method described above        //
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
void Runge_Kutta_Verner_Richardson_Integral_Curve( double (*f)(double, double),
             double y[], double x0, double h, int number_of_steps_per_interval,
                            int number_of_intervals, int richardson_columns ) {

   double mh = (double) number_of_steps_per_interval * h;

   while ( --number_of_intervals >= 0 ) {
      *(++y) = Runge_Kutta_Verner_Richardson( f, *y, x0, h,
                            number_of_steps_per_interval, richardson_columns );
      x0 += mh;
   }
   return;
}
