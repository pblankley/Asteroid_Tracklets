
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double x, y, z, xd, yd, zd;
} State;

typedef struct {
  double x, y, z, xd, yd, zd;
} Tangent;

typedef struct {
  double x, y, z;
} Vector;

typedef struct {
  double time;
  int body;
  State state;
} BinaryData;

#define FAILURE -1
#define SUCCESS 0
#define PI 3.14159265358979323846
