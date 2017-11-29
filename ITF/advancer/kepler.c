/* There must be a better version of this around somewhere in the Stark directory */

#include <string.h>
#include <stdio.h>
#include <math.h>

typedef struct {
  double x, y, z, xd, yd, zd;
} State;

extern double machine_epsilon;

#define MAX_KEPLER_COUNT 10

int kepler_step(double gm, double dt, State *s0, State *s)
{
  double r0, v0s, u, a, n, ecosE0, esinE0;
  double dM, x, sx, cx, f, fp, fpp, fppp, dx;
  double fdot, g, gdot;
  double sx2, cx2, x2;
  int count, flag;
  double xx, yy, xx1, yy1, omx, h;
  double k0x, k0y, k1x, k1y, k2x, k2y, k3y;
  double ecosE;

  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v0s = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  u = s0->x*s0->xd + s0->y*s0->yd + s0->z*s0->zd;
  a = 1.0/(2.0/r0 - v0s/gm);
  if(a<0.0) {
    flag = 1;
  } else {
    flag = 0;
    n = sqrt(gm/(a*a*a));
    ecosE0 = 1.0 - r0/a;
    esinE0 = u/(n*a*a);
  
    dM = n*dt;
  
    xx = ecosE0;
    yy = esinE0;
    
    h = dM;

    /* RK4 step for initial guess */
    omx = h/(1.0 - xx);
    k0x = - yy * omx;
    k0y = xx * omx;
    xx1 = xx + k0x/2.0;
    yy1 = yy + k0y/2.0;
    omx = h/(1.0 - xx1);
    k1x = - yy1 * omx;
    k1y = xx1 * omx;
    xx1 = xx + k1x/2.0;
    yy1 = yy + k1y/2.0;
    omx = h/(1.0 - xx1);
    k2x = - yy1 * omx;
    k2y = xx1 * omx;
    xx1 = xx + k2x;
    yy1 = yy + k2y;
    omx = h/(1.0 - xx1);
    k3y = xx1 * omx;
    yy += (k0y + 2.0*(k1y + k2y) + k3y)/6.0;
    
    x = dM - esinE0 + yy;
    
    count = 0;
    do {
      x2 = x/2.0;
      sx2 = sin(x2); cx2 = cos(x2);
      sx = 2.0*sx2*cx2; cx = cx2*cx2 - sx2*sx2;
      f = x + 2.0*sx2*(sx2*esinE0 - cx2*ecosE0) - dM;
      ecosE = cx*ecosE0 - sx*esinE0;
      fp = 1.0 - ecosE;
      fpp = (sx*ecosE0 + cx*esinE0)/2.0;
      fppp = ecosE/6.0;
      dx = -f/fp;
      dx = -f/(fp + dx*fpp);
      dx = -f/(fp + dx*(fpp + dx*fppp));
      x += dx;
      if(count++ > MAX_KEPLER_COUNT) return(3); /* this needs to be fixed */
    } while(fabs(f) > machine_epsilon);
    
    /* compute f and g */
    x2 = x/2.0;
    sx2 = sin(x2); cx2 = cos(x2);
    f = 1.0 - (a/r0)*2.0*sx2*sx2;
    sx = 2.0*sx2*cx2; cx = cx2*cx2 - sx2*sx2;
    g = (2.0*sx2*(esinE0*sx2 + cx2*r0/a))/n;
    fp = 1.0 - cx*ecosE0 + sx*esinE0;
    fdot = -(a/(r0*fp))*n*sx;
    gdot = (1.0 + g*fdot)/f;
    
    /* compute new position and velocity */
    s->x = f*s0->x + g*s0->xd;
    s->y = f*s0->y + g*s0->yd;
    s->z = f*s0->z + g*s0->zd;
    s->xd = fdot*s0->x + gdot*s0->xd;
    s->yd = fdot*s0->y + gdot*s0->yd;
    s->zd = fdot*s0->z + gdot*s0->zd;
    
  }
  return(flag);

}
