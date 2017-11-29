/*

2017 November 22
M. Holman

This code reads a set of time, position vector, velocity vector tuples
and uses a kepler stepper in universal variables to advance, or retard,
the position and velocity vectors to a common reference time.

 */

#include "heliocentric.h"

#define ECL	(84381.4118*(1./3600)*PI/180.) /*Obliquity of ecliptic at J2000*/
double speed_of_light = 173.144483; /* AU/day,  beware of fixed units! */

double GMsun;

State state, ref_state;
double machine_epsilon;

int main(int argc, char **argv)
{
  double t, tref, dt;

  char filename[100];
  FILE *file;

  double determine_machine_epsilon();

  int read_tuple(FILE *file, char *desig, double *t, State *state);

  char desig[100];

  void cartesian(double gm, 
		 double a, double e, double i, double longnode, double argperi, double meananom, 
		 State *state);

  void keplerian(double gm, State state, 
		 double *a, double *e, double *i, double *longnode, double *argperi, double *meananom);

  void xyz_ec_to_eq(double x_ec, double y_ec, double z_ec,
		    double *x_eq, double *y_eq, double *z_eq);

  void xyz_eq_to_ec(double x_eq, double y_eq, double z_eq,
		    double *x_ec, double *y_ec, double *z_ec);

  double principal_value(double theta);

  int kepler_step(double gm, double dt, State *s0, State *s);  
  int universal_kepler(double gm, double dt, State *s0, State *s);

  if(argc != 4) {
    printf("args: GMsun filename tref\n");
    exit(-1);
  }

  machine_epsilon = determine_machine_epsilon();

  sscanf(argv[1], "%lf", &GMsun);  
  sscanf(argv[2], "%s", filename);
  sscanf(argv[3], "%lf", &tref);

  file = fopen(filename, "r");
  if(file == NULL){
    printf("%s does not exist.\n", filename);
    exit(-1);
  }
  
  int result = read_tuple(file, desig, &t, &state);

  while(result != 0){
    universal_kepler(GMsun, t-tref, &state, &ref_state);
    printf("%s %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n", desig, tref, ref_state.x, ref_state.y, ref_state.z, ref_state.xd, ref_state.yd, ref_state.zd);

    result = read_tuple(file, desig, &t, &state);      
  }

  return 1;
  
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

int read_tuple(FILE *file, char *desig, double *t, State *state){

  char *fgets_nocomment(char *inbuff, int length, FILE *fpin, FILE *fpout);
  char line[1000];

  if(fgets_nocomment(line, 1000, file, NULL) !=0){

    sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf", desig, t, &(state->x), &(state->y), &(state->z), &(state->xd), &(state->yd), &(state->zd));

    return(1);

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
