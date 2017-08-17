////////////////////////////////////////////////////////////////////////////////
// File: testbulirsch.c                                                       //
// Purpose:                                                                   //
//    Test the Gragg_Bulirsch_Stoer method for solving a differential         //
//    equation in the file bulirsch_stoer.c                                   //
//                                                                            //
// Solve the initial value problem, y' = xy, y(0) = 1.0.                      //
// From x = 0.0 to 1.0, epsilon = 1.0e-10.                                    //
////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>

int Gragg_Bulirsch_Stoer( double (*f)(double, double), double y0, double *y1,
double x, double h, double *h_new, double epsilon, double yscale, 
int rational_extrapolate  );

// y' = f(x,y) = xy, y(0) = 1
double f(double x, double y) { return x*y; }

// The actual solution
double If(double x) { return exp(0.5*x*x); }


double tolerance = 1.e-10;          
double a = 1.0;               // y(x0).
double x0 = 0.0;              // x0 the initial condition.
double h = 0.1;               // step size
int number_of_steps = 10;     // and number of steps, i.e. 
                              // solve for x = x0 to x0 + h * number_of_steps.


FILE *out;

int print_line;

void Print_Header() {
   fprintf(out,"Prog: testbulirsch.c\n");
   fprintf(out,"Gragg-Bulirsch-Stoer Method for the solution of a diffeq\n");
   print_line = 2;
}

void Print_Polynomial_Extrapolation_Integral_Test() {
   double exact;
   double error;
   double h_next;
   double h_used;
   double x;
   double x_start;
   double x_stop;
   double y0 = a;
   double y1;
   int i;
   int err;
 
   fprintf(out,"\n\n\nGragg_Bulirsch_Stoer Method\n\n");
   fprintf(out,"Problem: Solve y' = xy,\n");
   fprintf(out,"Initial condition x = %4.1lf, y = %4.1lf\n",x0,y0);
   fprintf(out,"Number of steps %d and step size %4.2lf\n",number_of_steps,h);
   fprintf(out,"Polynomial Extrapolation\n\n");
   fprintf(out,"  n     x[n]           Estimate                 Exact");
   fprintf(out,"             Error\n");
   print_line += 10;

               /* Load initial conditions, (y0 is already set) */

   x_start = x0;
   x = x_start;

          /* Call Gragg_Bulirsh_Stoer method 'number_of_steps' times */
          /* with a step size of h. */

   for (i = 0; i < number_of_steps; i++) {

               /* Load initial conditions, (y0 is already set) */

      h_used = h;
      x_stop = x_start + h;

               /* Continue calling Gragg_Bulirsch_Stoer until */
            /* we have finished solving for x = xstart to x = stop */

      do {

         /* If the Gragg_Bulirsch_Stoer Method fails, divide the */
         /* step size used by 4 and try again.                   */

         /* NOTE.  In general you would want to put a counter    */
         /* here and terminate the routine if it fails too many  */
         /* times.                                               */
         
         while (Gragg_Bulirsch_Stoer( f, y0, &y1, x, h_used, &h_next, 1.0,
                                             tolerance, 0)) h_used /= 4.0;
         x += h_used;
         y0 = y1;
         if ( x + h_next > x_stop ) h_used = x + h_next - x_stop;
         else h_used = h_next;
      } while ( x < x_stop - 1.e-10 );
                  
      exact = If(x);
      error = exact - y1;
      fprintf(out,"%3d   %8.2le  %20.15le", i+1, x, y1);
      fprintf(out,"   %20.15le  %+9.4le\n", exact,error);
      print_line++;
      x_start = x;
   }
}

void Print_Rational_Extrapolation_Integral_Test() {
   double exact;
   double error;
   double h_next;
   double h_used;
   double x;
   double x_start;
   double x_stop;
   double y0 = a;
   double y1;
   int i;
   int err;
 
   fprintf(out,"\n\n\nGragg_Bulirsch_Stoer Method\n\n");
   fprintf(out,"Problem: Solve y' = xy,\n");
   fprintf(out,"Initial condition x = %4.1lf, y = %4.1lf\n",x0,y0);
   fprintf(out,"Number of steps %d and step size %4.2lf\n",number_of_steps,h);
   fprintf(out,"Rational Extrapolation\n\n");
   fprintf(out,"  n     x[n]           Estimate                 Exact");
   fprintf(out,"             Error\n");
   print_line += 10;

               /* Load initial conditions, (y0 is already set) */

   x_start = x0;
   x = x_start;

          /* Call Gragg_Bulirsh_Stoer method 'number_of_steps' times */
          /* with a step size of h. */

   for (i = 0; i < number_of_steps; i++) {

               /* Load initial conditions, (y0 is already set) */

      h_used = h;
      x_stop = x_start + h;

               /* Continue calling Gragg_Bulirsch_Stoer until */
            /* we have finished solving for x = xstart to x = stop */

      do {

         /* If the Gragg_Bulirsch_Stoer Method fails, divide the */
         /* step size used by 4 and try again.                   */

         /* NOTE.  In general you would want to put a counter    */
         /* here and terminate the routine if it fails too many  */
         /* times.                                               */
         
         while (Gragg_Bulirsch_Stoer( f, y0, &y1, x, h_used, &h_next, 1.0,
                                               tolerance,1)) h_used /= 4.0;
         x += h_used;
         y0 = y1;
         if ( x + h_next > x_stop ) h_used = x + h_next - x_stop;
         else h_used = h_next;
      } while ( x < x_stop - 1.e-10 );
                  
      exact = If(x);
      error = exact - y1;
      fprintf(out,"%3d   %8.2le  %20.15le", i+1, x, y1);
      fprintf(out,"   %20.15le  %+9.4le\n", exact,error);
      print_line++;
      x_start = x;
   }
}
            
int main() 
{
   out = fopen("Bulirsch_Stoer.txt","w");
 
   Print_Header(); 
   Print_Polynomial_Extrapolation_Integral_Test();
   fprintf(out,"\n\n\n");
   Print_Rational_Extrapolation_Integral_Test();
   fclose(out);

   return 0;
}
