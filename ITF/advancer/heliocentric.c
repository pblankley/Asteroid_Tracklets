/*

This code is based on the heliocentric tangent map code that Kat Deck and I developed, which is in turn based on a code that Jack Wisdom
and I wrote.

This version is designed to read a long series of barycentric positions for the sun, planets, and the moon, obtained from the DE405 ephemeris.
It then uses an n-body map in canonical heliocentric coordinates to integrate the orbits of the bunch of asteroids, some with masses and 
some as test particles.

The goal is to rapidly but accurately integrate orbits for the MPC.

I am going first strip out a bunch of stuff, just to make the code clearer.  Then I will add it back in as needed.  I will strip out:

1) The symplectic correctors.
2) The relativity stuff.
3) The tangent equations.
4) The sphere of influence stuff.

There are three types of bodies here:

1) The bodies represented by the DE405, i.e. the sun and planets.
2) The "big three" asteroids that have significant masses.
3) The other asteroids that will be treated as test particles.

Need to read the "big three" from the MPCORB.DAT file.
Need to read the small asteroids also from the MPCORB.DAT file.

7 January 2015:

This is really coming along.  I can integrate ~30 years of Ceres and Juno, as tests, and match the JPL positions to ~50 or ~200km, respectively.
And it is very fast.  

At this point I think I need to first make a few improvements:

1) Add back in the solar GR terms.  Note: this is going to be more work that I anticipated, as I confirmed that the precession approximation
I normally use is not adequate when the mean longitudes are important.

2) Start building the capacity to handle close encounters with the Earth-Moon system.

3) Put in DE430 and its largest asteroids.

 */

#include "heliocentric.h"
#include "chealpix.h"

//#define ECL	(84381.448*(1./3600)*PI/180.) /*Obliquity of ecliptic at J2000*/
#define ECL	(84381.4118*(1./3600)*PI/180.) /*Obliquity of ecliptic at J2000*/

typedef struct {
  double GMbig;
  char desig[10];
  int idx;
} BigAsteroid;

typedef struct {
  char desig[10];
  int idx;
} SmallAsteroid;

// These should be dynamically allocated.
#define MAX_N_PLANETS 20
#define MAX_N_BIG 20
#define MAX_N_PARTICLES 10
#define MAX_N_ASTEROIDS 720000

double GMsun;
double over_GMsun;

double GM[MAX_N_PLANETS];
double GMbig[MAX_N_BIG];
State p[MAX_N_PLANETS];    // States of large asteroids
//State tp[MAX_N_PARTICLES]; // States of small asteroids
State h[MAX_N_PLANETS];    // positions of planets

double influence[MAX_N_PLANETS];

State *tp;

State tmp_state, tmp_state_eq;

BigAsteroid big_asteroids[MAX_N_BIG];
SmallAsteroid *small_asteroids;

double machine_epsilon;
int n_planets, n_big, n_particles;
int n_steps;

Vector *positions;
Vector *velocities;
double *jd_tt;

double speed_of_light = 173.144483; /* AU/day,  beware of fixed units! */
double relativity_constant;

int main(int argc, char **argv)
{
  double dt;
  int nplot, ntotal, planet;

  char ic_filename[100];
  FILE *ic_file;
  char positions_filename[100];
  FILE *positions_file;

  char elements_filename[100];
  FILE *elements_file;

  char big_asteroids_filename[100];
  FILE *big_asteroids_file;

  int i, j, k;
  int n_bodies;

  int step;

  int i_big, i_small;

  double determine_machine_epsilon();

  void read_positions(FILE *f, int n_planets, int n_steps, Vector *positions, Vector *velocities, double *jd_tt);
  void read_masses(FILE *f, int n_bodies, double *GMsun, double *GM, double *influence);
  void read_big_asteroids(FILE *f, int n_big, double *GMbig, BigAsteroid *big_aseteroids);
  void kepler_steps(double dt);
  void direct_kicks(State p[], State tp[], Vector h[], double dt, double jd);
  void coriolis_kicks(State p[], State tp[], Vector tmp);

  double julian_date (short int year, short int month, short int day,
		      double hour);

  int read_elements(FILE *big_asteroids_file, FILE *elements_file, char *desig, double *H, double *G, double *epoch_tdt,
		    double *meananom_epoch, double *argperi, double *longnode, double *incl, double *e, double *meanmotion, double *a);

  char desig[100], Epoch_str[100];
  double H, G, phase_angle;
  double meanmotion, n;
  double epoch_tdt, meananom_epoch;
  double a, e, incl, longnode, argperi, meananom;
  int big;

  void cartesian(double gm, 
		 double a, double e, double i, double longnode, double argperi, double meananom, 
		 State *state);

  void keplerian(double gm, State state, 
		 double *a, double *e, double *i, double *longnode, double *argperi, double *meananom);

  void xyz_ec_to_eq(double x_ec, double y_ec, double z_ec,
		    double *x_eq, double *y_eq, double *z_eq);

  void xyz_eq_to_ec(double x_eq, double y_eq, double z_eq,
		    double *x_ec, double *y_ec, double *z_ec);

  void copy_state(State *s1, State *s2);

  double principal_value(double theta);

  Vector tmp, sun_pos, sun_vel, dsun;

  int kepler_step(double gm, double dt, State *s0, State *s);  
  int kepler_step_array(double gm, double *dt, int ntimes, State *s0, State *s);

  double chebeval(double a, double b, double *c, int m, double x);

  long   nside, nside_outer;
  long  ipix, npix, dpix, ip2, ip1, ip3;

  double theta, phi;
  double vec[3];

  double raDeg0, decDeg0;
  double x0, y0, z0;
  double raDeg, decDeg;

  double ang, radiusDeg;
  
  double *hpx, *hpy, *hpz;

  tmp.x = 0.0; tmp.y = 0.0; tmp.z = 0.0;
  
  if(argc != 7) {
    printf("args: ic_filename positions_filename elements_filename big_asteroids_filename dt ntotal\n");
    exit(-1);
  }

  machine_epsilon = determine_machine_epsilon();

  nside = 32;
  npix = nside2npix(nside);

  hpx = malloc(npix*sizeof(double));
  hpy = malloc(npix*sizeof(double));
  hpz = malloc(npix*sizeof(double));

  for (ipix = 0; ipix < npix; ipix++) {
    pix2ang_nest(nside, ipix, &theta, &phi);

    raDeg = phi*180./PI;
    decDeg = 90.-theta*180./PI;

    hpx[ipix] = cos(raDeg*PI/180.)*cos(decDeg*PI/180.);
    hpy[ipix] = sin(raDeg*PI/180.)*cos(decDeg*PI/180.);
    hpz[ipix] = sin(decDeg*PI/180.);

  }
  

  sscanf(argv[1], "%s", ic_filename);
  sscanf(argv[2], "%s", positions_filename);
  sscanf(argv[3], "%s", elements_filename);
  sscanf(argv[4], "%s", big_asteroids_filename);
  sscanf(argv[5], "%lf", &dt);
  sscanf(argv[6], "%d", &ntotal);

  ic_file = fopen(ic_filename, "r");
  if(ic_file == NULL){
    printf("ic_file does not exist.\n");
    exit(-1);
  }
  positions_file = fopen(positions_filename, "r");
  if(positions_file == NULL){
    printf("positions_file does not exist.\n");
    exit(-1);
  }
  elements_file = fopen(elements_filename, "r");
  if(elements_file == NULL){
    printf("elements_file does not exist.\n");
    exit(-1);
  }
  big_asteroids_file = fopen(big_asteroids_filename, "r");
  if(big_asteroids_file == NULL){
    printf("big_asteroids_file does not exist.\n");
    exit(-1);
  }


  fscanf(positions_file, "%d %d", &n_steps, &n_bodies);
  read_masses(ic_file, n_bodies, &GMsun, GM, influence);

  double GMtotal = GMsun;
  for(int i=0; i<n_planets; i++){
    GMtotal += GM[i];
  }

  positions =  malloc((n_steps + 1) * n_bodies * sizeof(Vector));
  velocities = malloc((n_steps + 1) * n_bodies * sizeof(Vector));
  jd_tt      = malloc((n_steps + 1) * sizeof(double));
  
  read_positions(positions_file, n_bodies, n_steps, positions, velocities, jd_tt);

  sun_pos = positions[n_bodies - 1];
  sun_vel = velocities[n_bodies - 1];

  n_planets = n_bodies - 1;

  fscanf(big_asteroids_file, "%d", &n_big);
  read_big_asteroids(big_asteroids_file, n_big, GMbig, big_asteroids);

  tp              = malloc(MAX_N_ASTEROIDS * sizeof(State));
  small_asteroids = malloc(MAX_N_ASTEROIDS * sizeof(SmallAsteroid));
  
  i_big=0;
  i_small=0;
  big = read_elements(big_asteroids_file, elements_file, desig, &H, &G, &epoch_tdt, &meananom_epoch, &argperi, &longnode, &incl, &e, &meanmotion, &a);
  while(big != 0){

    cartesian(GMsun, a, e, incl, longnode, argperi, meananom_epoch, &tmp_state);
    //cartesian(GMtotal, a, e, incl, longnode, argperi, meananom_epoch, &tmp_state);    

    xyz_ec_to_eq(tmp_state.x, tmp_state.y, tmp_state.z, &tmp_state_eq.x, &tmp_state_eq.y, &tmp_state_eq.z);
    xyz_ec_to_eq(tmp_state.xd, tmp_state.yd, tmp_state.zd, &tmp_state_eq.xd, &tmp_state_eq.yd, &tmp_state_eq.zd);

    tmp_state_eq.xd += sun_vel.x;
    tmp_state_eq.yd += sun_vel.y;
    tmp_state_eq.zd += sun_vel.z;

    // Add sun's velocity to make the velocities barycentric for mapping

    if(big==2){
      copy_state(&tmp_state_eq, &p[i_big]);
      i_big++;
    }else{
      copy_state(&tmp_state_eq, &tp[i_small]);
      strncpy(small_asteroids[i_small].desig, desig, 8);
      i_small++;
    }

    big = read_elements(big_asteroids_file, elements_file, desig, &H, &G, &epoch_tdt, &meananom_epoch, &argperi, &longnode, &incl, &e, &meanmotion, &a);
    
  }
  n_particles= i_small;  

  relativity_constant = 6.0*(GMsun/speed_of_light)*(GMsun/speed_of_light);
  relativity_constant = 0.0;


  // Integrate from the initial epoch to the final epoch,
  // using a canonical heliocentric integrator
  
  sun_vel = velocities[(1)*n_bodies - 1];

  step=0;
  sun_pos = positions[(step+1)*n_bodies - 1];

  printf("%d %lf\n", step, jd_tt[step]);
  sun_vel = velocities[(step+2)*n_bodies - 1];

  for(step = 0; step<ntotal; step +=2) {

    // 1/2 step kick
    direct_kicks(p, tp, (positions+step*n_bodies), dt, jd_tt[step]);

    // 1/2 step drift
    tmp = positions[(step+2)*n_bodies - 1];
    dsun.x = +(sun_pos.x-tmp.x); dsun.y = +(sun_pos.y-tmp.y); dsun.z = +(sun_pos.z-tmp.z);
    coriolis_kicks(p, tp, dsun);
    sun_pos = tmp;

    // full kepler step
    kepler_steps(dt*2.0);

    // 1/2 step drift
    tmp = positions[(step+3)*n_bodies - 1];
    dsun.x = +(sun_pos.x-tmp.x); dsun.y = +(sun_pos.y-tmp.y); dsun.z = +(sun_pos.z-tmp.z);
    coriolis_kicks(p, tp, dsun);
    sun_pos = tmp;

    // 1/2 step kick
    direct_kicks(p, tp, (positions+(step+2)*n_bodies), dt, jd_tt[step+2]);

    sun_vel = velocities[(step+2)*n_bodies - 1];

  }

  // Convert from UTC to TT.
  double leaps = 36.0; // Becomes 37 at the end of 2016.  
  double jd_utc = jd_tt[step]; // This is the desired time of the observation in UTC.
  double jd_tt_final = jd_utc + (32.184 + leaps)/(24.0*60.*60.); // This is the corresponding time of the ephemeris

  // The ephemeris is in JPL's T_eph, which is essentially TDB.  We may need to account for the difference between
  // TT and TDB for highest precision.

  // Now we need to get the position of the geocenter at jd_utc/jd_tt.
  // We have it in the positions array, but for a different time.
  // Let's assume that the step we ended on is the closest in time to what we want,
  // and we can do a kepler_step to the time we want.
  // This could come directly from a call to a JPL ephemeris.

  State geocenter, geocenter_final;
  geocenter.x = positions[step*n_bodies+2].x - positions[(step+1)*n_bodies-1].x;
  geocenter.y = positions[step*n_bodies+2].y - positions[(step+1)*n_bodies-1].y;
  geocenter.z = positions[step*n_bodies+2].z - positions[(step+1)*n_bodies-1].z;

  geocenter.xd = velocities[step*n_bodies+2].x - velocities[(step+1)*n_bodies-1].x;
  geocenter.yd = velocities[step*n_bodies+2].y - velocities[(step+1)*n_bodies-1].y;
  geocenter.zd = velocities[step*n_bodies+2].z - velocities[(step+1)*n_bodies-1].z;

  // Set up some constants for the chebyshev polynomials
  int cheby_order = 5;
  double cheby_dts[cheby_order];
  for(int k=0; k<cheby_order; k++){
    double y = cos(PI*(k+0.5)/cheby_order);
    cheby_dts[k]=y*dt;
  }
  
  double final_dt = (jd_tt_final - jd_utc); 

  // Let's try to do a set of times
  double geo_dts[2*cheby_order+1];
  for(int i=0; i<(2*cheby_order+1); i++){
    geo_dts[i] = -dt + i*2.0*dt/(2*cheby_order) + final_dt;
    printf("%d %lf\n", i, geo_dts[i]);
  }

  kepler_step(GMsun+GM[2], geo_dts[cheby_order], &geocenter, &geocenter_final);
  
  printf("%d %lf\n", step, jd_tt[step]);  
  for(int i=0; i<n_particles; i++){
    tmp_state_eq.x  = tp[i].x;
    tmp_state_eq.y  = tp[i].y;
    tmp_state_eq.z  = tp[i].z;
    tmp_state_eq.xd = tp[i].xd-sun_vel.x;
    tmp_state_eq.yd = tp[i].yd-sun_vel.y;
    tmp_state_eq.zd = tp[i].zd-sun_vel.z;

    State cheby_states[cheby_order];
    double cheby_coeffs_x[cheby_order];
    double cheby_coeffs_y[cheby_order];
    double cheby_coeffs_z[cheby_order];        

    kepler_step_array(GMsun, cheby_dts, cheby_order, &tmp_state_eq, cheby_states);
    for(int j=0; j<cheby_order; j++){
      double sumx = 0.0;
      double sumy = 0.0;
      double sumz = 0.0;
      for(int k=0; k<cheby_order; k++){
	double cosfactor = cos(PI*j*(k+0.5)/cheby_order);
	sumx += cheby_states[k].x * cosfactor;
	sumy += cheby_states[k].y * cosfactor;
	sumz += cheby_states[k].z * cosfactor;
      }
      cheby_coeffs_x[j] = sumx*2.0/cheby_order;
      cheby_coeffs_y[j] = sumy*2.0/cheby_order;
      cheby_coeffs_z[j] = sumz*2.0/cheby_order;      
    }

    // Now that we have the chebyshev coefficients and can quickly evaluate the heliocentric position of
    // the asteroid, let's see if we can do the geocentric light time correction.

    double dtp = 0.0;
    double distance = 0.0;
    double distance_prev;
    double au_km = 149597870.700;
    double speed_of_light = 2.99792458e5 * 86400./au_km;
    double dx, dy, dz;
    double x_mp, y_mp, z_mp;
    do{
      distance_prev = distance;
      x_mp = chebeval(-dt, dt, cheby_coeffs_x, cheby_order, final_dt+dtp);
      y_mp = chebeval(-dt, dt, cheby_coeffs_y, cheby_order, final_dt+dtp);
      z_mp = chebeval(-dt, dt, cheby_coeffs_z, cheby_order, final_dt+dtp);

      dx = x_mp - geocenter_final.x;
      dy = y_mp - geocenter_final.y;
      dz = z_mp - geocenter_final.z;
      distance = sqrt(dx*dx + dy*dy + dz*dz);

      double ltt = distance/speed_of_light; // C is speed of light
      dtp = -ltt;
      
    }while(fabs((distance-distance_prev)/distance)>1e-8);

    vec[0] = dx;
    vec[1] = dy;
    vec[2] = dz;

    vec2pix_nest(nside, vec, &ip2);

    double ltt = distance/speed_of_light; // C is speed of light
    printf("%s ", small_asteroids[i].desig);
    printf("%13.6lf %13.6lf ", principal_value(atan2(dy, dx))*180./PI, asin(dz/distance)*180./PI);
    printf("%8ld ", ip2);
    printf("%le %le %le ", dx, dy, dz);
    printf("\n");

    // Calculate the nominal healpix location of the minor planet, as viewed from the geocenter.  DONE
    // Calculate the set of healpix locations of the minor planet, as viewed from the geocenter, over a time range (a night)
    // Calculate the set of healpix locations of the minor planet, as viewed from the geocenter, over a time range (a night),
    // and over a range of longitudinal uncertainty

  }

}

void direct_kicks(State p[], State tp[], Vector h[], double dt, double jd)
{
  int i, j;
  double dx, dy, dz, rij2, fij, fj, tx, ty, tz;
  double rp2, fr, r2;
  Vector acc[MAX_N_BIG], tp_acc, indirect;
  Vector sun;

  sun = h[n_planets];


  for(i=0; i<n_big; i++) {
    /* relativity */
    rp2 = p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z;
    fr = relativity_constant / (rp2*rp2); 
    acc[i].x = - fr*p[i].x; acc[i].y = - fr*p[i].y; acc[i].z = - fr*p[i].z;
    //acc[i].x = 0.0; acc[i].y = 0.0; acc[i].z = 0.0;
  }

  /* Calculate the accelation of each large asteroid due to the other
     large asteroids. */
  for(i=0; i<n_big; i++) {
    for(j=i+1; j<n_big; j++) {
      dx = p[i].x - p[j].x; dy = p[i].y - p[j].y; dz = p[i].z - p[j].z; 
      rij2 = dx*dx + dy*dy + dz*dz;
      fij = -1.0/(rij2*sqrt(rij2));
      tx = dx*fij; ty = dy*fij; tz = dz*fij;
      acc[i].x += GMbig[j]*tx; acc[i].y += GMbig[j]*ty; acc[i].z += GMbig[j]*tz;
      acc[j].x -= GMbig[i]*tx; acc[j].y -= GMbig[i]*ty; acc[j].z -= GMbig[i]*tz;
    }
  }

  /* Calculate the accelation of each large asteroid due to the planets */
  for(i=0; i<n_big; i++) {
    for(j=0; j<n_planets; j++) {
      //printf("%le %lf %lf %lf\n", GM[j]/GMsun, h[j].x, h[j].y, h[j].z);
      dx = p[i].x - (h[j].x-sun.x); dy = p[i].y - (h[j].y-sun.y); dz = p[i].z - (h[j].z-sun.z); 
      rij2 = dx*dx + dy*dy + dz*dz;
      fij = -1.0/(rij2*sqrt(rij2));
      tx = dx*fij; ty = dy*fij; tz = dz*fij;
      acc[i].x += GM[j]*tx; acc[i].y += GM[j]*ty; acc[i].z += GM[j]*tz;
    }
  }

  /* kick the large asteroids */
  for(i=0; i<n_big; i++) {
    p[i].xd += (dt*acc[i].x);  
    p[i].yd += (dt*acc[i].y);  
    p[i].zd += (dt*acc[i].z);  
  }

  /* Deal with the small asteroids */
  for(i=0; i<n_particles; i++){

    /* relativity */
    rp2 = tp[i].x*tp[i].x + tp[i].y*tp[i].y + tp[i].z*tp[i].z;
    fr = relativity_constant / (rp2*rp2); 
    tp_acc.x = - fr*tp[i].x; tp_acc.y = - fr*tp[i].y; tp_acc.z = - fr*tp[i].z;
    //tp_acc.x = 0.0; tp_acc.y = 0.0; tp_acc.z = 0.0;

    /* Calculate the accelation due to the large asteroids */
    for(j=0; j<n_big; j++){
      dx = tp[i].x - p[j].x;
      dy = tp[i].y - p[j].y;
      dz = tp[i].z - p[j].z;
      r2 = dx*dx + dy*dy + dz*dz;
      fj = GMbig[j]/(r2*sqrt(r2));
      tp_acc.x -= fj*dx;
      tp_acc.y -= fj*dy;
      tp_acc.z -= fj*dz;
    }

    /* Calculate the accelation due to the planets */
    for(j=0; j<n_planets; j++){
      dx = tp[i].x - (h[j].x-sun.x);
      dy = tp[i].y - (h[j].y-sun.y);
      dz = tp[i].z - (h[j].z-sun.z);
      r2 = dx*dx + dy*dy + dz*dz;
      if(strncmp(small_asteroids[i].desig, "D4340 ", 5) !=0 || j!=8){ // Avoid calculating perturbation of Pluto on itself 'D4340'
	if(r2 < influence[j]) {
	  //printf("Close encounter between asteroid %s and planet %d %le at JD %lf\n", small_asteroids[i].desig, j, sqrt(r2/influence[j]), jd);
	  fflush(stdout);
	}
	fj = GM[j]/(r2*sqrt(r2));
	tp_acc.x -= fj*dx;
	tp_acc.y -= fj*dy;
	tp_acc.z -= fj*dz;
      }
    }

    tp[i].xd += (dt*tp_acc.x);  
    tp[i].yd += (dt*tp_acc.y);  
    tp[i].zd += (dt*tp_acc.z);  
  }
}

/* This needs to be altered so that the drift comes directly from 
   the change in position of the sun with respect to the barycenter,
   which comes from the ephemeris.
*/
void coriolis_kicks(State p[], State tp[], Vector tmp)
{
  double dt_over_GMsun;
  int i;

  for(i=0; i<n_big; i++){
    p[i].x += tmp.x;
    p[i].y += tmp.y;
    p[i].z += tmp.z;
  }

  for(i=0; i<n_particles; i++){
    tp[i].x += tmp.x;
    tp[i].y += tmp.y;
    tp[i].z += tmp.z;
  }
}

void kepler_steps(double dt)
{
  State tmp;
  int i; 
  void copy_state(State *s1, State *s2);
  int kepler_step(double gm, double dt, State *s0, State *s);

  for(i=0; i<n_big; i++) {
    kepler_step(GMsun, dt, &p[i], &tmp);
    copy_state(&tmp, &p[i]);
  }
  for(i=0; i<n_particles; i++) {
    kepler_step(GMsun, dt, &tp[i], &tmp);
    copy_state(&tmp, &tp[i]);
  }
}

void copy_state(State *s1, State *s2)
{
  s2->x = s1->x;
  s2->y = s1->y;
  s2->z = s1->z;
  s2->xd = s1->xd;
  s2->yd = s1->yd;
  s2->zd = s1->zd;
}

void read_masses(FILE *f, int n_bodies, double *GMsun, double *GM, double *influence)
{
  int i, j;
  int n_planets = n_bodies-1;
  double Msun_over_M;
  char tmpstr[100];
  double GMsun_;
  double a, tmp;

  fscanf(f, "%s %lf", tmpstr, &GMsun_);
  *GMsun = GMsun_;
  
  for(i=0; i<n_planets; i++){
    fscanf(f, "%d %s %lf %lf", &j, tmpstr, &Msun_over_M, &a);
    fflush(stdout);
    GM[i] = GMsun_/Msun_over_M;
    tmp = a*pow(GM[i]/(3.0*GMsun_), 1.0/3.0);
    influence[i] = tmp*tmp;  /* the square of sphere of influence */
    
    //printf("%d %lf %le %le\n", i, a, GM[i], sqrt(influence[i]));
    
  }

  return;
}  

void read_big_asteroids(FILE *f, int n_big, double *GMbig, BigAsteroid *big_asteroids)
{
  int i, j;
  double M_over_Msun;
  char tmpstr[100];
  char desig[100];

  for(i=0; i<n_big; i++){
    if(fscanf(f, "%s %s %lf", desig, tmpstr, &M_over_Msun) == 0){
      printf("Empty!\n");
      fflush(stdout);
    }
    GMbig[i] = GMsun*M_over_Msun;

    big_asteroids[i].GMbig = GMbig[i];
    big_asteroids[i].idx = i;
    strcpy(big_asteroids[i].desig, desig);

    //printf("%d %s %s %le\n", i, big_asteroids[i].desig, tmpstr, GMbig[i]);

    
  }

}

void read_positions(FILE *f, int n_bodies, int n_steps, Vector *positions, Vector *velocities, double *jd_tt)
{
  int i, j;
  int i_tmp;
  double jdt;

  i = 0;
  while(fscanf(f, "%lf", &jdt) != EOF) {

    jd_tt[i] = jdt;
    for(int k=0; k<n_bodies; k++){
      fscanf(f, "%d %lf %lf %lf", &i_tmp, &positions[i*n_bodies+k].x, &positions[i*n_bodies+k].y, &positions[i*n_bodies+k].z);
      fscanf(f, "%lf %lf %lf",  &velocities[i*n_bodies+k].x, &velocities[i*n_bodies+k].y, &velocities[i*n_bodies+k].z);
    }

    i++;
  }

}


void read_planetary_system(char *filename)
{
  FILE *f;
  int i, j;
  double solar_mass, GMplanet, a, GM_total, over_GM;
  State s, sun_i;

  f = fopen(filename, "r");
 
  fscanf(f, "%lf", &GMsun);

  over_GMsun = 1.0/GMsun;

  i = 0; j = 0;
  while(fscanf(f, "%lf", &GMplanet) != EOF) {
    GM[i] = GMplanet;
    if(GMplanet != 0.0){
      fscanf(f, "%lf %lf %lf", &p[i].x, &p[i].y, &p[i].z);
      fscanf(f, "%lf %lf %lf", &p[i].xd, &p[i].yd, &p[i].zd);
      i++;
    }else{
      fscanf(f, "%lf %lf %lf", &tp[j].x, &tp[j].y, &tp[j].z);
      fscanf(f, "%lf %lf %lf", &tp[j].xd, &tp[j].yd, &tp[j].zd);
      j++;
    }
  }
  n_planets = i;
  n_particles = j;

  fclose(f);
}


void print_planetary_system()
{
  int i;

  printf("%.16lf\n", GMsun);
  for(i=0; i<n_planets; i++) {
    printf("%.16lf\n", GM[i]);
    printf("%.16lf %.16lf %.16lf\n", p[i].x, p[i].y, p[i].z);
    printf("%.16lf %.16lf %.16lf\n", p[i].xd, p[i].yd, p[i].zd);
  }
  for(i=0; i<n_particles; i++) {
    printf("%.16lf\n", 0.0);
    printf("%.16lf %.16lf %.16lf\n", tp[i].x, tp[i].y, tp[i].z);
    printf("%.16lf %.16lf %.16lf\n", tp[i].xd, tp[i].yd, tp[i].zd);
  }
}


void inertial_heliocentric(State *inertial, State sun, State *helio)
{
  int i;
  double mass;

  for(i=0; i<n_planets; i++){
    (helio+i)->x = (inertial+i)->x - sun.x;
    (helio+i)->y = (inertial+i)->y - sun.y;
    (helio+i)->z = (inertial+i)->z - sun.z;
    (helio+i)->xd = (inertial+i)->xd - sun.xd;
    (helio+i)->yd = (inertial+i)->yd - sun.yd;
    (helio+i)->zd = (inertial+i)->zd - sun.zd;
  
  }
}


/* Read a line from specified input stream, skipping blank or
 * commented (#) lines.  These are echoed to fpout if it's non-null.
 * Returns NULL at EOF.
 */
char *
fgets_nocomment(char *inbuff, int length, 
		FILE *fpin, FILE *fpout)
{
  int needmore=0;
  char test[10];
  while (1) {
    if (fgets(inbuff,length,fpin)==NULL) return(NULL);

    /* Skip blank lines or comments */
    if (sscanf(inbuff,"%8s",test)>0
	&& test[0]!='#') return(inbuff);

    /* Echo comments to output if desired */
    if (fpout!=NULL) 
      fputs(inbuff,fpout);
  }
}

int read_elements(FILE *big_asteroids_file, FILE *elements_file, char *desig, double *H, double *G, double *epoch_tdt,
		  double *meananom_epoch, double *argperi, double *longnode, double *incl, double *e, double *meanmotion, double *a)
{

  static int first=1;
  char *fgets_nocomment(char *inbuff, int length, FILE *fpin, FILE *fpout);
  int interpret_epoch(char *s, int *year, int *month, int *day);
  char line[1000], the_rest[1000];
  char Epoch_str[100];
  int year, month, day;
  double a_tmp;

  double gconst = GMsun;

  double julian_date (short int year, short int month, short int day,
		      double hour);

  int big = 0;

  if(first == 1){
    do{
      fgets_nocomment(line,1000,elements_file,NULL);
    }while(strncmp(line, "-----", 5) != 0);
    first = 0;
  }

  if(fgets_nocomment(desig,8,elements_file,NULL) !=0){

    for(int k=0; k<n_big; k++){
      if(strncmp(desig, big_asteroids[k].desig, 5)==0){
	big = 1;
      }
    }

    fgets(line,13,elements_file);

    if(sscanf(line, "%lf %lf", H, G) == EOF){
      *H = 99.0;
      *G = 0.15;
    }
    
    fscanf(elements_file, "%s %lf %lf %lf %lf %lf %lf %lf", 
	   Epoch_str, meananom_epoch, argperi, longnode, incl, e, meanmotion, a);
    fgets_nocomment(the_rest,1000,elements_file,NULL);

    interpret_epoch(Epoch_str, &year, &month, &day);
    *epoch_tdt = julian_date((short int)year, (short int)month, (short int)day, 0.0);  /* epoch_tdt is a TDB Julian date */

    *incl *= PI/180.;  *longnode *= PI/180.;  *argperi *= PI/180.;  *meananom_epoch *= PI/180.;  *meanmotion *= PI/180.;

    a_tmp = *a;
    /* Compute mean motion from definition */
    *meanmotion = sqrt(gconst/(a_tmp*a_tmp*a_tmp));

    if(big==1){
      return(2);
    }else{
      return(1);
    }
  }else{
    return(0);
  }
}

int interpret_epoch(char *s, int *year, int *month, int *day)
{

  char subyear[3];

  switch (s[0]) {
  case 'I':
    *year = 1800;
    break;
  case 'J':
    *year = 1900;
    break;
  case 'K':
    *year = 2000;
    break;
  default:
    return 0;
  }

  subyear[0] = s[1];
  subyear[1] = s[2];
  subyear[2] = '\0';

  /*printf("s: %s sy: %s atoi: %d\n", s, subyear, atoi(subyear));*/

  *year += atoi(subyear);

  switch (s[3]) {
  case '1':
    *month = 1;
    break;
  case '2':
    *month = 2;
    break;
  case '3':
    *month = 3;
    break;
  case '4':
    *month = 4;
    break;
  case '5':
    *month = 5;
    break;
  case '6':
    *month = 6;
    break;
  case '7':
    *month = 7;
    break;
  case '8':
    *month = 8;
    break;
  case '9':
    *month = 9;
    break;
  case 'A':
    *month = 10;
    break;
  case 'B':
    *month = 11;
    break;
  case 'C':
    *month = 12;
    break;
  default:
    return 0;
  }

  switch (s[4]) {
  case '1':
    *day = 1;
    break;
  case '2':
    *day = 2;
    break;
  case '3':
    *day = 3;
    break;
  case '4':
    *day = 4;
    break;
  case '5':
    *day = 5;
    break;
  case '6':
    *day = 6;
    break;
  case '7':
    *day = 7;
    break;
  case '8':
    *day = 8;
    break;
  case '9':
    *day = 9;
    break;
  case 'A':
    *day = 10;
    break;
  case 'B':
    *day = 11;
    break;
  case 'C':
    *day = 12;
    break;
  case 'D':
    *day = 13;
    break;
  case 'E':
    *day = 14;
    break;
  case 'F':
    *day = 15;
    break;
  case 'G':
    *day = 16;
    break;
  case 'H':
    *day = 17;
    break;
  case 'I':
    *day = 18;
    break;
  case 'J':
    *day = 19;
    break;
  case 'K':
    *day = 20;
    break;
  case 'L':
    *day = 21;
    break;
  case 'M':
    *day = 22;
    break;
  case 'N':
    *day = 23;
    break;
  case 'O':
    *day = 24;
    break;
  case 'P':
    *day = 25;
    break;
  case 'Q':
    *day = 26;
    break;
  case 'R':
    *day = 27;
    break;
  case 'S':
    *day = 28;
    break;
  case 'T':
    *day = 29;
    break;
  case 'U':
    *day = 30;
    break;
  case 'V':
    *day = 31;
    break;

  default:
    return 0;
  }

  return 1;
}

/********julian_date */

double julian_date (short int year, short int month, short int day,
                    double hour)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function will compute the Julian date for a given calendar
      date (year, month, day, hour).

   REFERENCES: 
      Fliegel & Van Flandern, Comm. of the ACM, Vol. 11, No. 10, October
      1968, p. 657.

   INPUT
   ARGUMENTS:
      year (short int)
         Year.
      month (short int)
         Month number.
      day (short int)
         Day-of-month.
      hour (double)
         Hour-of-day.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
         Julian date.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/06-98/JAB (USNO/AA)

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'juldat'.
      2. This function makes no checks for a valid input calendar
      date.
------------------------------------------------------------------------
*/
{
   long int jd12h;

   double tjd;

   jd12h = (long) day - 32075L + 1461L * ((long) year + 4800L
      + ((long) month - 14L) / 12L) / 4L
      + 367L * ((long) month - 2L - ((long) month - 14L) / 12L * 12L)
      / 12L - 3L * (((long) year + 4900L + ((long) month - 14L) / 12L)
      / 100L) / 4L;
   tjd = (double) jd12h - 0.5 + hour / 24.0;

   return (tjd);
}

/* And transform x,y,z from ecliptic to eq */
void xyz_ec_to_eq(double x_ec, double y_ec, double z_ec,
		  double *x_eq, double *y_eq, double *z_eq)
{
  double se,ce;

  se = sin(-ECL);
  ce = cos(ECL);

  *x_eq = x_ec;
  *y_eq = ce*y_ec + se*z_ec;
  *z_eq = -se*y_ec + ce*z_ec;

  return;
}

/* And transform x,y,z from eq to ecliptic */
void
xyz_eq_to_ec(double x_eq, double y_eq, double z_eq,
	     double *x_ec, double *y_ec, double *z_ec)
{
  double	se,ce;

  se = sin(ECL);
  ce = cos(ECL);

  *x_ec = x_eq;
  *y_ec = ce*y_eq + se*z_eq;
  *z_ec = -se*y_eq + ce*z_eq;

  return;
}

double principal_value(double theta)
{
  theta -= 2*PI*floor(theta/(2*PI));
  return(theta);
}

double chebeval(double a, double b, double *c, int m, double x){

  double d=0.0, dd=0.0, sv, y, y2;
  int j;

  // Put in check for range
  //if((x-a)*(x-b) > 0.0)

  y = (2.0*x-a-b)/(b-a);
  y2 = 2.0*y;
  for(j=m-1; j>=1; j--){
    sv = d;
    d = y2*d-dd+c[j];
    dd = sv;
  }

  return y*d - dd + 0.5*c[0];

}
