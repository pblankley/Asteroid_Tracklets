/*

This code is a simple implementation of the p-iteration approach of Herrick and Liu, as described in Danby's Fundamentals of 
Celestial Mechanics.

The purpose is to solve for the orbit that passes from vector r0 at time t0 to vector r1 at time t1.

 */

#include "heliocentric.h"

State p0, p1;

double GMsun;

double nb_time;

double machine_epsilon;

double dt;
double xi, yi, zi, r0;
double xf, yf, zf, r1;
double dc, dc0;
double c, s, cp;
double f, g;
double y=1.0;
double param, dp;

#include "machine-epsilon.c"

int universal_kepler(double gm, double dt, State *s0, State *s);

int main(int argc, char **argv)
{
  int nmax;
  double test_p(double param);
  double zbrent(double (*func)(double), double xf, double x2, double tol);
  double p_iteration(double xi, double yi, double zi, double xf, double yf, double zf, double dt, double param0, double param1);

  if(argc != 11) {
    printf("args: xi yi zi xf yf zf dt p dp nmax\n");
    exit(-1);
  }

  machine_epsilon = determine_machine_epsilon();

  GMsun = 2.959122082322128e-04;
  GMsun = 1.0;

  sscanf(argv[1], "%lf", &xi);
  sscanf(argv[2], "%lf", &yi);
  sscanf(argv[3], "%lf", &zi);
  sscanf(argv[4], "%lf", &xf);
  sscanf(argv[5], "%lf", &yf);
  sscanf(argv[6], "%lf", &zf);
  sscanf(argv[7], "%lf", &dt);
  sscanf(argv[8], "%lf", &param);
  sscanf(argv[9], "%lf", &dp);
  sscanf(argv[10], "%d", &nmax);


  r0 = sqrt(xi*xi + yi*yi + zi*zi);
  r1 = sqrt(xf*xf + yf*yf + zf*zf);

  p0.x = xi;
  p0.y = yi;
  p0.z = zi;

  p1.x = xf;
  p1.y = yf;
  p1.z = zf;

  c = (xi*xf + yi*yf + zi*zf)/(r0*r1);
  s = y*sqrt(1.0-c*c);

  dc = 0.0; dc0 = 0.0;

  int i = 0;
  while(dc*dc0 >= 0.0 && i<nmax){

    f = (r1/param)*(c-1.0) + 1.0;
    g = (r0*r1)/sqrt(GMsun*param) * s;

    p0.xd = (xf - f*xi)/g;
    p0.yd = (yf - f*yi)/g;
    p0.zd = (zf - f*zi)/g;
  
    universal_kepler(GMsun, dt, &p0, &p1);

    double r1p = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);

    cp = (p0.x*p1.x + p0.y*p1.y + p0.z*p1.z)/(r0*r1p);

    printf("%d %lf %lf %lf %lf\n", i, c, cp, param, cp-c);
    printf("%lf %lf %lf %lf %lf %lf\n", p0.x, p0.y, p0.z, p0.xd, p0.yd, p0.zd);
    printf("%lf %lf %lf %lf %lf %lf\n", p1.x, p1.y, p1.z, p1.xd, p1.yd, p1.zd);
    printf("%lf %lf %lf\n", xf, yf, zf);
    printf("\n");

    dc0 = dc;
    dc = cp-c;
    param += dp;
    i++;

  }

  if(i==nmax){
    printf("No bracket!\n");
  }else{
    printf("%lf %lf %lf %lf\n", param-2.*dp, dc0, param-dp, dc);
  }

  double result = zbrent(test_p, param-2.*dp, param-dp, 1e-12);

  printf("%lf\n", result);

  printf("%lf\n", p_iteration(xi, yi, zi, xf, yf, zf, dt, param-2.*dp, param-dp));
  
}

double p_iteration(double xi, double yi, double zi, double xf, double yf, double zf, double dt, double param0, double param1){

  double test_p(double param);
  double zbrent(double (*func)(double), double xf, double x2, double tol);

  r0 = sqrt(xi*xi + yi*yi + zi*zi);
  r1 = sqrt(xf*xf + yf*yf + zf*zf);

  c = (xi*xf + yi*yf + zi*zf)/(r0*r1);
  s = y*sqrt(1.0-c*c);

  p0.x = xi;
  p0.y = yi;
  p0.z = zi;

  return zbrent(test_p, param0, param1, 1e-12);

}

double test_p(double param){

  f = (r1/param)*(c-1.0) + 1.0;
  g = (r0*r1)/sqrt(GMsun*param) * s;

  p0.xd = (xf - f*p0.x)/g;
  p0.yd = (yf - f*p0.y)/g;
  p0.zd = (zf - f*p0.z)/g;
  
  universal_kepler(GMsun, dt, &p0, &p1);
  double r1p = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);

  cp = (p0.x*p1.x + p0.y*p1.y + p0.z*p1.z)/(r0*r1p);

  return(cp-c);
}


