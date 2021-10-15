#include "mts.h"
#include "vecfunc.h"
double vec_dot (struct Point *P1, struct Point *P2) { return (((*P1).x * (*P2).x) + ((*P1).y *
(*P2).y) + ((*P1).z * (*P2).z)); }
double vec_norm (struct Point *v) { return (sqrt(vec_dot(v, v))); }
void vec_mid (struct Point *P1, struct Point *P2, struct Point *scatPt) {
(*scatPt).x = ((*P1).x+(*P2).x)/2;
(*scatPt).y = ((*P1).y+(*P2).y)/2,
(*scatPt).z = ((*P1).z+(*P2).z)/2;
return;
}
void vec_sub (struct Point *P1, struct Point *P2, struct Point *vec) {
(*vec).x = ( (*P1).x - (*P2).x );
(*vec).y = ( (*P1).y - (*P2).y );
(*vec).z = ( (*P1).z - (*P2).z );
return;
}
void vec_add (struct Point *P1, struct Point *P2, struct Point *vec) {
(*vec).x = ( (*P1).x + (*P2).x );
(*vec).y = ( (*P1).y + (*P2).y );
(*vec).z = ( (*P1).z + (*P2).z );
return;
}
void vec_mult (double sc, struct Point *P1, struct Point *vec) {
(*vec).x = sc * (*P1).x;
(*vec).y = sc * (*P1).y;
(*vec).z = sc * (*P1).z;
return;
}
void vec_div (double sc, struct Point *P1, struct Point *vec) {
(*vec).x = (*P1).x / sc;
(*vec).y = (*P1).y / sc;
(*vec).z = (*P1).z / sc;
return;
}
void vec_fit (struct Point** points, struct Line* muonTrack, int detectors, int fitX) {
double vecXY[MAX_DIMENSION][MAX_DETECTORS];
double vecZ[MAX_DETECTORS];
double normal[MAX_DETECTORS][MAX_DIMENSION];
double lhs[MAX_DIMENSION][MAX_DIMENSION];
double rhs[MAX_DIMENSION];
double augmented[MAX_DIMENSION+1][MAX_DIMENSION];
double lSum,rSum,newXY1,newZ1,newXY2,newZ2,factor,c0,c1;
int i,j,k;
struct Point* dummy;
if (fitX==FIT_NONE) {
dummy = (*muonTrack).P1;
(*dummy).x = (*points[0]).x;
(*dummy).y = (*points[0]).y;
(*dummy).z = (*points[0]).z;
dummy = (*muonTrack).P2;
(*dummy).x = (*points[detectors-1]).x;
(*dummy).y = (*points[detectors-1]).y;
(*dummy).z = (*points[detectors-1]).z;
} else {
if ((fitX == FIT_X) && ((*points[0]).x == (*points[detectors-1]).x)) {
dummy = (*muonTrack).P1;
(*dummy).x = (*points[0]).x;
dummy = (*muonTrack).P2;
(*dummy).x = (*points[detectors-1]).x;
return;
}
else if ((fitX ==FIT_Y) && ((*points[0]).y == (*points[detectors-1]).y)) {
dummy = (*muonTrack).P1;
(*dummy).y = (*points[0]).y;
dummy = (*muonTrack).P2;
(*dummy).y = (*points[detectors-1]).y;
return;
}
for (i = 0; i<detectors; i++) {
vecXY[0][i] = 1;
if (fitX == FIT_X) vecXY[MAX_DIMENSION-1][i] = (*points[i]).x;
else vecXY[MAX_DIMENSION-1][i] = (*points[i]).y;
vecZ[i] = (*points[i]).z;
}
for (i = 0; i < MAX_DIMENSION; i++) {
for (j = 0; j < detectors; j++) {
normal[j][i] = vecXY[i][j];
}
}
for (i = 0; i < MAX_DIMENSION; i++) {
for (j = 0; j < MAX_DIMENSION; j++) {
lSum = 0;
rSum = 0;
for (k = 0; k <detectors; k++) {
lSum += normal[k][i] * vecXY[j][k];
rSum += normal[k][i] * vecZ[k];
}
lhs[i][j] = lSum;
rhs[i] = rSum;
}
}
//create augmented matrix to solve 2x2 system
for (i = 0; i < MAX_DIMENSION; i++) {
for (j = 0; j < MAX_DIMENSION; j++) {
augmented[i][j] = lhs[i][j];
}
augmented[j][i] = rhs[i];
}
//create augmented matrix to solve 2x2 system
augmented[0][0] = lhs[0][0];
augmented[0][1] = lhs[0][1];
augmented[1][0] = lhs[1][0];
augmented[1][1] = lhs[1][1];
augmented[2][0] = rhs[0];
augmented[2][1] = rhs[1];
factor = augmented[0][1]/augmented[0][0];
for (i = 0; i < 3; i++) {
augmented[i][1] = augmented[i][1] - (factor*augmented[i][0]);
}
c1 = augmented[2][1] / augmented[1][1];
c0 = (augmented[2][0] - augmented[1][0]*c1) / augmented[0][0];
newZ1 = vecZ[0];
newXY1 = (newZ1-c0)/c1;
newZ2 = vecZ[detectors-1];
newXY2 = (newZ2-c0)/c1;
dummy = (*muonTrack).P1;
if (fitX==FIT_X) (*dummy).x = newXY1;
else (*dummy).y = newXY1;
(*dummy).z = newZ1;
dummy = (*muonTrack).P2;
if (fitX==FIT_X) (*dummy).x = newXY2;
else (*dummy).y = newXY2;
(*dummy).z = newZ2;
}
return;
}
double vec_rad_to_deg (double angle) { return (RAD * angle); }
//vec_angle uses a.b=|a||b|acos(theta) to compute the angle between two vectors
double vec_angle (struct Line* L1, struct Line* L2, int component) {
double scatAng, distance, dotUV, normU, normV;
struct Point *u,*v;
if ((u = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((v = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
vec_sub((*L1).P2, (*L1).P1, u);
vec_sub((*L2).P2, (*L2).P1, v);
if (component!=ALL_COMPONENTS) {
if (component==X_COMPONENT) (*u).x = (*v).x = 0;
if (component==Y_COMPONENT) (*u).y = (*v).y = 0;
}
normU = vec_norm(u);
normV = vec_norm(v);
dotUV = vec_dot(u, v);
if (((dotUV) / ((normU) * (normV))) >= 1) scatAng = 0;
else scatAng = acos(dotUV / (normU * normV));
if (dotUV<0) {
if (scatAng > 0) scatAng+=M_PI;
else scatAng-=M_PI;
}
free(v);
free(u);
return scatAng;
}
//vecCopy copies the contents of L1 to L2
void vec_copy(struct Line* L1, struct Line* L2) {
struct Point *P1, *P2;
P1 = (*L1).P1;
P2 = (*L2).P1;
(*P2).x = (*P1).x;
(*P2).y = (*P1).y;
(*P2).z = (*P1).z;
P1 = (*L1).P2;
P2 = (*L2).P2;
(*P2).x = (*P1).x;
(*P2).y = (*P1).y;
(*P2).z = (*P1).z;
return;
}
void vec_copy_point(struct Point* P1, struct Point* P2) {
(*P2).x = (*P1).x;
(*P2).y = (*P1).y;
(*P2).z = (*P1).z;
return;
}