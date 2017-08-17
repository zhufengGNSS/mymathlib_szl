////////////////////////////////////////////////////////////////////////////////
// File: simpson_simpson.c                                                    //
// Routines:                                                                  //
//    Simpson_Simpson_Adapative                                               //
////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>                            // required for malloc()
#include <math.h>                              // required for fabs()

struct Subinterval { 
   double upper_limit;
   double lower_limit;
   double function[5];
   struct Subinterval *interval;
};

static double s1, s2;
static struct Subinterval interval;
static struct Subinterval *pinterval;
static double (*fnc)(double);

static void Simpsons_Rule_Update();

////////////////////////////////////////////////////////////////////////////////
//  double Simpson_Simpson_Adaptive( double a, double b, double tolerance,    //
//                             double (*f)(double), double min_h, int *err ); //
//                                                                            //
//  Description:                                                              //
//                                                                            //
//    Starting at the left-end point, a, find the min power of 2, m, so that  //
//    the difference between using Simpson's rule and the composite Simpson's //
//    rule on the interval [a, a+(b-a)/2^m] is less than twice the tolerance  //
//    * (length of the subinterval)/(b-a).  Then repeat the process for       //
//    integrating over the interval [(a+(b-a)/2^m,b] until the right end      //
//    point, b, is finally reached.                                           //
//    The integral is then the sum of the integrals of each subinterval.  If  //
//    at any time, the length of the subinterval for which the estimates      //
//    based on Simpson's rule and the composite Simpson's rule is less than   //
//    min_h the process is terminated.                                        //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double a          The lower limit of the integration interval.         //
//     double b          The upper limit of integration.                      //
//     double tolerance  The acceptable error estimate of the integral.       //
//                       The integral of a subinterval is accepted when the   //
//                       magnitude of the estimated error falls below the     //
//                       pro-rated tolerance of the subinterval.              //
//     double *f         Pointer to the integrand, a function of a single     //
//                       variable of type double.                             //
//     int    min_h      The minimum subinterval length.  If no subinterval   //
//                       of length > min_h is found for which the estimated   //
//                       error falls below the pro-rated tolerance, the       //
//                       process terminates after setting *err to -1.         //
//     int    *err       0 if the process terminates successfully; -1 if no   //
//                       subinterval of length > min_h was found for which    //
//                       the estimated error was less that the pro-rated      //
//                       error, -2 if memory could not be allocated to        //
//                       proceed with a new subinterval.                      //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to a +  h.                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Simpson_Simpson_Adaptive(double a, double b, double tolerance, 
                              double (*f)(double), double min_h, int *err) {

   double integral = 0.0;
   double epsilon_density = 2.0 * tolerance / ( b - a );
   double epsilon;

   struct Subinterval *qinterval;

      // Create the initial level, with lower_limit = a, upper_limit = b,  //   
      // and f(x) evaluated at a, b, and (a + b) / 2.                      //

   interval.interval = NULL;
   interval.upper_limit = b;
   interval.lower_limit = a;
   interval.function[0] = (*f)(interval.lower_limit);
   interval.function[2] = (*f)( 0.5 * ( interval.upper_limit 
                                        + interval.lower_limit ) );
   interval.function[4] = (*f)(interval.upper_limit);
   pinterval = &interval;

            // Calculate the tolerance for the current interval.  //
            // calculate the single subinterval Simpson rule,     //
            // and the two subintervals composite Simpson rule.   //

   fnc = f;
   *err = 0;
   epsilon = epsilon_density * (b - a);
   Simpsons_Rule_Update();

   while ( pinterval->upper_limit - pinterval->lower_limit > min_h ) {
      if ( fabs( s1 - s2 ) < epsilon ) {

            // If the two estimates are close, then increment the    //
            // integral and if we are not at the right end, set the  // 
            // left end of the new interval to the right end of the  //
            // old interval and the right end of the new interval    //
            // remains the same (as the previous right end for this  //
            // interval.                                             //

         integral += s2;
         if (pinterval->interval == NULL) { return integral; }
         qinterval = pinterval->interval;
         qinterval->lower_limit = pinterval->upper_limit;
         qinterval->function[0] = qinterval->function[2];
         qinterval->function[2] = qinterval->function[3];
         free(pinterval);
         pinterval = qinterval;
      }
      else {
            // If the two estimates are not close, then create a new //
            // interval with same left end point and right end point // 
            // at the midpoint of the current interval.              //
         
         qinterval = (struct Subinterval*) malloc( sizeof(struct Subinterval) );
         if ( qinterval == NULL ) { *err = -2; break; }
         qinterval->interval = pinterval;
         qinterval->lower_limit = pinterval->lower_limit;
         qinterval->upper_limit = 0.5 * (pinterval->upper_limit
                                           + pinterval->lower_limit);
         qinterval->function[0] = pinterval->function[0];
         qinterval->function[2] = pinterval->function[1];
         qinterval->function[4] = pinterval->function[2];
         pinterval = qinterval;
      }
      Simpsons_Rule_Update();
      epsilon = epsilon_density * (pinterval->upper_limit
                                                    - pinterval->lower_limit);
   }

            // The process failed, free all allocated memory.  //

   while (pinterval != NULL) {
      qinterval = pinterval->interval;
      free(pinterval); 
      pinterval = qinterval;
   }
   if ( *err == 0 ) { *err = -1; return 0.0; }

   return 0.0;
};



static void Simpsons_Rule_Update( ) {

   double h = pinterval->upper_limit - pinterval->lower_limit;
   double h4 = 0.25 * h;

   pinterval->function[1] = (*fnc)(pinterval->lower_limit + h4);
   pinterval->function[3] = (*fnc)(pinterval->upper_limit - h4);

   s1 = pinterval->function[0] + 4.0 * pinterval->function[2]
         + pinterval->function[4];
   s1 *= 0.166666666666666666666667 * h;
   s2 = pinterval->function[0] + 4.0 * pinterval->function[1] 
         + 2.0 * pinterval->function[2] + 4.0 * pinterval->function[3]
         + pinterval->function[4];
   s2 *= 0.0833333333333333333333333 * h;
}
