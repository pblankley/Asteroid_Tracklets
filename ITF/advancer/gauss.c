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
  double ra[3];
  double dec[3];
  double rh[3][3];
  double rm[3][3];

  if(argc != 7) {
    printf("args: ra0 dec0 ra1 dec1 ra2 dec2\n");
    exit(-1);
  }

  machine_epsilon = determine_machine_epsilon();

  GMsun = 2.959122082322128e-04;
  GMsun = 1.0;

  sscanf(argv[1], "%lf", &ra[0]);
  sscanf(argv[2], "%lf", &dec[0]);
  sscanf(argv[3], "%lf", &ra[1]);
  sscanf(argv[4], "%lf", &dec[1]);
  sscanf(argv[5], "%lf", &ra[2]);
  sscanf(argv[6], "%lf", &dec[2]);

  for(int i=0; i<3; i++){
    rh[i][0] = cos(dec[i])*cos(ra[i]);
    rh[i][1] = cos(dec[i])*sin(ra[i]);
    rh[i][2] = sin(dec[i]);
  }

  double w1 = 0.0;
  for(int j=0; j<3; j++){
    rm[0][j] = rh[0][j];
    w1 = w1 + rh[0][j]*rh[2][j];
  }

  double w2 = sqrt(1.0 - w1*w1);
  for(int j=0; j<3; j++){
    rm[1][j] = (rh[2][j] - w1*rh[0][j])/w2;   
  }

  rm[2][0] = rm[0][1]*rm[1][2] - rm[0][2]*rm[1][1];
  rm[2][1] = rm[0][2]*rm[1][0] - rm[0][0]*rm[1][2];
  rm[2][2] = rm[0][0]*rm[1][1] - rm[0][1]*rm[1][0];

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


