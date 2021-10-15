#include "mts.h"
#include "preprocessing.h"
#include "mtserr.h"
#include "mtsio.h"
#include "vecfunc.h"
int preprocessing(double* lambda, double* M, struct muon* head, double** params, FILE** fps) {
fprintf(ERR_OUT, "%s POCA RECONSTRUCTION %s\n\n", BANNER, BANNER);
int i, j, voxel, xVox, yVox, zVox, errCode = CONTINUE;
double *voxel_std, *voxel_avg, *voxel_n, *doca, *lambdaLocal, *lambdaMedLocal;
//pointers to incoming and outgoing muon tracks
struct Line *muonInc, *muonOut;
//create pointers to detector points on incoming/outgoing muon tracks
struct Point *incPoints[MAX_DETECTORS], *outPoints[MAX_DETECTORS];
struct Point *scatPt, *tempP;
struct muon *curMuon;
struct voxel *trackHead, *dummy = (struct voxel*) malloc(sizeof(struct voxel));
//set aside memory for the different structs
if ((muonInc = (struct Line*) malloc(sizeof(struct Line)))==NULL) return memError();
if (((*muonInc).P1 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if (((*muonInc).P2 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((muonOut = (struct Line*) malloc(sizeof(struct Line)))==NULL) return memError();
if (((*muonOut).P1 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if (((*muonOut).P2 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((scatPt = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
for (i = 0; i < MAX_DETECTORS; i++) {
if ((incPoints[i] = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((outPoints[i] = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
}
if (header(lambda, M, params, fps) != CONTINUE) return 0;
if ((doca = (double*) malloc(sizeof(double) * *params[PARAM_EVENTS]/5))==NULL) return
memError();
if (*params[PARAM_STD]) {
if ((voxel_std = (double*) malloc(sizeof(double) *
*params[PARAM_ALL_VOXELS]))==NULL) return memError();
if ((voxel_avg = (double*) malloc(sizeof(double) *
*params[PARAM_ALL_VOXELS]))==NULL) return memError();
if ((voxel_n = (double*) malloc(sizeof(double) *
*params[PARAM_ALL_VOXELS]))==NULL) return memError();
memset(voxel_std, '\0', sizeof(double) * *params[PARAM_ALL_VOXELS]);
memset(voxel_avg, '\0', sizeof(double) * *params[PARAM_ALL_VOXELS]);
memset(voxel_n, '\0', sizeof(double) * *params[PARAM_ALL_VOXELS]);
}
curMuon = head;
while(errCode!=END) {
//get muon tracks from either input file or manual input
errCode = get_geant_input(incPoints, outPoints, params, fps);
if (*params[PARAM_MOM_CUT]) {
if ((*params[PARAM_MOMENTUM] > *params[PARAM_MOM_HIGH_CUT]) &&
(*params[PARAM_MOM_HIGH_CUT])) continue;
if ((*params[PARAM_MOMENTUM] < *params[PARAM_MOM_LOW_CUT]) &&
(*params[PARAM_MOM_LOW_CUT])) continue;
}
if (errCode==END) {
if (*params[PARAM_STD]==1) {
(*params[PARAM_STD])++;
*params[PARAM_CUR_EVENT]=0;
rewind(fps[FP_IN]);
for (i=0;i<*params[PARAM_ALL_VOXELS];i++) voxel_avg[i] = voxel_avg[i] / voxel_n[i];
errCode=CONTINUE;
continue;
} else break;
}
(*params[PARAM_CUR_EVENT])++;
if (errCode==ERR_CODE_INVALID_INPUT) continue;
vec_fit(incPoints, muonInc, (int) *params[PARAM_INC_DECT], FIT_X);
vec_fit(incPoints, muonInc, (int) *params[PARAM_INC_DECT], FIT_Y);
vec_fit(outPoints, muonOut, (int) *params[PARAM_OUT_DECT], FIT_X);
vec_fit(outPoints, muonOut, (int) *params[PARAM_OUT_DECT], FIT_Y);
//find point of closest approach between incoming/outgoing muon tracks and return doca
*params[PARAM_DOCA] = pocaLtoL(muonInc, muonOut, scatPt, params);
//check if tracks were parallel and if the scatter point is inside the detector volume and set flag to
0 if not
if (*params[PARAM_DOCA]==PARALLEL)*params[PARAM_CONTINUE] = 0;
if (!(in_volume(scatPt, params))) *params[PARAM_CONTINUE] = 0;
*params[PARAM_SCAT_ANG] = vec_angle(muonInc, muonOut, ALL_COMPONENTS);
if (*params[PARAM_STD] && *params[PARAM_CONTINUE]) {
xVox = (int) floor((((*scatPt).x + *params[PARAM_X_MAX]) /
*params[PARAM_X_VOXEL_SIZE]));
yVox = (int) floor((((*scatPt).y + *params[PARAM_Y_MAX]) /
*params[PARAM_Y_VOXEL_SIZE]));
zVox = (int) floor((((*scatPt).z + *params[PARAM_Z_MAX]) /
*params[PARAM_Z_VOXEL_SIZE]));
//voxel number determined in z direction first, then y, then x
voxel = (xVox * *params[PARAM_Y_VOXEL_TOTAL] *
*params[PARAM_Z_VOXEL_TOTAL]) + (yVox * *params[PARAM_Z_VOXEL_TOTAL]) +
zVox;
if (*params[PARAM_STD]==1) {
voxel_avg[voxel] = voxel_avg[voxel] + *params[PARAM_SCAT_ANG];
voxel_n[voxel]++;
}
if (*params[PARAM_STD]==2) voxel_std[voxel] = voxel_std[voxel] +
pow((*params[PARAM_SCAT_ANG] - voxel_avg[voxel]), 2);
}
//if em analysis is to be done, print the appropriate data to file or store it in the appropriate data
structures
if (((*curMuon).nextMuon=(struct muon*) malloc(sizeof(struct muon)))==NULL) return
memError();
curMuon = (*curMuon).nextMuon;
if ((errCode = em_data(muonInc, muonOut, curMuon, params, fps))!=CONTINUE) return
errCode;
if ((*curMuon).event==7699) fprintf(stderr, "angle %f dx %f\n", (*curMuon).dX,
(*curMuon).dtX);
if (*params[PARAM_OUT]) write_optional((*muonInc).P2, scatPt, (*muonOut).P1, params,
fps);
if (((*curMuon).muonTrack=(struct voxel*) malloc(sizeof(struct voxel)))==NULL) return
memError();
//the following block of code takes care of track analysis; if the scatter point was valid then the
track is calculate along the POCA
//path, if not then the path between the entering and exiting tracks is used
trackHead = dummy;
(*trackHead).nextVoxel = (*curMuon).muonTrack;
if (*params[PARAM_CONTINUE]) {
if ((trackHead = track((*muonInc).P2, scatPt, (*muonOut).P1, curMuon, trackHead, M,
params, fps))==NULL) return memError();
if (((*trackHead).nextVoxel=(struct voxel*) malloc(sizeof(struct voxel)))==NULL) return
memError();
if ((trackHead = track(scatPt, (*muonOut).P1, NULL, curMuon, trackHead, M, params,
fps))==NULL) return memError();
} else {
if ((trackHead = track((*muonInc).P2, (*muonOut).P1, NULL, curMuon, trackHead, M,
params, fps))==NULL) return memError();
}
if (*params[PARAM_EM_ONLINE]) {
if (fmod(*params[PARAM_CUR_EVENT], *params[PARAM_EM_ONLINE])==0)
em(lambdaLocal, lambdaMedLocal, head, params, fps);
}
*params[PARAM_PREV_VOXEL] = -1;
*params[PARAM_CONTINUE]=1;
//break;
}
(*curMuon).nextMuon = NULL;
if (*params[PARAM_EM_ONLINE]) {
if (fmod(*params[PARAM_CUR_EVENT], *params[PARAM_EM_ONLINE])!=0)
em(lambdaLocal, lambdaMedLocal, head, params, fps);
free(lambdaLocal);
free(lambdaMedLocal);
}
fprintf(stderr, "\n\nPARALLEL PARAM: %.10f\n\n", SMALL_NUM);
fprintf(stderr, "PARALLEL TRACKS: %f\n", *params[PARAM_PARALLEL]);
fprintf(stderr, "TOTAL TRACKS: %f\n", *params[PARAM_CUR_EVENT]);
double max = 0;
if (*params[PARAM_STD]) {
for (i=0;i<*params[PARAM_ALL_VOXELS];i++) {
voxel_std[i] = pow((voxel_std[i] / voxel_n[i]), 0.5);
if (voxel_std[i]>max) max=voxel_std[i];
//if (*params[PARAM_EM]) lambda[i] = voxel_std[i];
}
fprintf(stderr, "max = %f\n\n", max);
if (fps[FP_STD_OUT]!=NULL) write_lambda(voxel_std, voxel_n, NULL, params,
fps[FP_STD_OUT]);
free(voxel_avg);
free(voxel_n);
free(voxel_std);
}
free(dummy);
free(doca);
free((*muonInc).P1);
free((*muonInc).P2);
free((*muonOut).P1);
free((*muonOut).P2);
free(muonInc);
free(muonOut);
for (i = 0; i <MAX_DETECTORS; i++) {
free(incPoints[i]);
free(outPoints[i]);
}
return errCode;
}
//header processes the first 2 lines from the input file which containts the number of events in
//the run and the lenght of the detector volume in x, y and z and then determines other paramaters
//based on available information
int header(double* lambda, double* M, double** params, FILE** fps) {
int i;
char numEvents[20]; //string which contains number of *params[PARAM_EVENTS] in
simulation
//gets the number of total *params[PARAM_EVENTS] from the second *params[PARAM_LINE]
in the input file
if ((fgets(numEvents, 100, fps[FP_IN])==NULL)) return emptyError(INPUT_FILE);
(*params[PARAM_LINE])++;
if (numEvents[0]!=EVENTS_START) return formatError(INPUT_FILE, numEvents,
(*params[PARAM_LINE])++);
if ((fgets(numEvents, 100, fps[FP_IN])==NULL)) return emptyError(INPUT_FILE);
if ((*params[PARAM_EVENTS] = atof(strtok(numEvents, DELIMS))) == 0) return
formatError(INPUT_FILE, numEvents, (*params[PARAM_LINE])++);
*params[PARAM_X_LENGTH] = atof(strtok(NULL, DELIMS));
if (*params[PARAM_X_LENGTH]==0) return formatError(INPUT_FILE, numEvents,
(*params[PARAM_LINE])++);
*params[PARAM_Y_LENGTH] = atof(strtok(NULL, DELIMS));
if (*params[PARAM_Y_LENGTH]==0) return formatError(INPUT_FILE, numEvents,
(*params[PARAM_LINE])++);
*params[PARAM_Z_LENGTH] = atof(strtok(NULL, DELIMS));
if (*params[PARAM_Z_LENGTH]==0) return formatError(INPUT_FILE, numEvents,
(*params[PARAM_LINE])++);
*params[PARAM_X_MAX] = *params[PARAM_X_LENGTH]/2;
*params[PARAM_Y_MAX] = *params[PARAM_Y_LENGTH]/2;
*params[PARAM_Z_MAX] = *params[PARAM_Z_LENGTH]/2;
*params[PARAM_X_MIN] = *params[PARAM_X_LENGTH]/2 * -1;
*params[PARAM_Y_MIN] = *params[PARAM_Y_LENGTH]/2 * -1;
*params[PARAM_Z_MIN] = *params[PARAM_Z_LENGTH]/2 * -1;
*params[PARAM_X_VOXEL_TOTAL] = *params[PARAM_X_LENGTH] /
*params[PARAM_X_VOXEL_SIZE];
*params[PARAM_Y_VOXEL_TOTAL] = *params[PARAM_Y_LENGTH] /
*params[PARAM_Y_VOXEL_SIZE];
*params[PARAM_Z_VOXEL_TOTAL] = *params[PARAM_Z_LENGTH] /
*params[PARAM_Z_VOXEL_SIZE];
*params[PARAM_ALL_VOXELS] = *params[PARAM_X_VOXEL_TOTAL] *
*params[PARAM_Y_VOXEL_TOTAL] * *params[PARAM_Z_VOXEL_TOTAL];
return CONTINUE;
}
// pocaLtoL():
// Input: two lines L1 and L2:
// Return: the shortest distance between L1 and L2
double pocaLtoL(struct Line *L1, struct Line *L2, struct Point* scatPt, double** params) {
double a, b, c, d, e, D, sc, tc, doca;
struct Point *u, *v, *w, *dPVector, *scXu, *tcXv, *scXuMtcXv, *cp1, *cp2;
if ((u = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((v = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((w = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((dPVector = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((scXu = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((tcXv = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((scXuMtcXv = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((cp1 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((cp2 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
vec_sub((*L1).P2, (*L1).P1, u);
vec_sub((*L2).P2, (*L2).P1, v);
vec_sub((*L1).P1, (*L2).P1, w);
a = vec_dot(u, u);
b = vec_dot(u, v);
c = vec_dot(v, v);
d = vec_dot(u, w);
e = vec_dot(v, w);
D = (a)*(c) - (b)*(b);
// compute the line parameters of the two closest points
if (D < SMALL_NUM) { // the lines are almost parallel
(*params[PARAM_PARALLEL])++;
(*scatPt).x = (*((*L1).P2)).x;
(*scatPt).y = (*((*L1).P2)).y;
(*scatPt).z = (*((*L1).P2)).z;
return PARALLEL;
} else {
sc = ((b)*(e) - (c)*(d)) / D;
tc = ((a)*(e) - (b)*(d)) / D;
}
vec_mult(sc, u, scXu);
vec_mult(tc, v, tcXv);
//next three lines are for distance of closest approach
vec_sub(scXu, tcXv, scXuMtcXv);
vec_add(w, scXuMtcXv, dPVector);
doca = vec_norm(dPVector);
//next three lines are for point of closest approach
vec_add((*L1).P1, scXu, cp1);
vec_add((*L2).P1, tcXv, cp2);
vec_mid(cp1, cp2, scatPt);
free(cp2);
free(cp1);
free(scXuMtcXv);
free(scXu);
free(tcXv);
free(dPVector);
free(w);
free(v);
free(u);
return doca;
}
double in_volume(struct Point* p, double** params) {
if (((*p).x > *params[PARAM_X_MAX]) || ((*p).x < *params[PARAM_X_MIN])) return 0;
if (((*p).y > *params[PARAM_Y_MAX]) || ((*p).y < *params[PARAM_Y_MIN])) return 0;
if (((*p).z > *params[PARAM_Z_MAX]) || ((*p).z < *params[PARAM_Z_MIN])) return 0;
return 1;
}
int em_data(struct Line* muonIn, struct Line* muonOut, struct muon* mu, double** params,
FILE** fps) {
int i, j;
double thetaX, thetaX0, thetaX1, thetaY, thetaY0, thetaY1, difX, difY, deltaX, deltaY, pr2, Lxy;
struct Line *vertIn, *vertOut;
struct Point *dvIn, *p1, *distance;
if ((vertIn = (struct Line*) malloc(sizeof(struct Line)))==NULL) return memError();
if (((*vertIn).P1 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if (((*vertIn).P2 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((vertOut = (struct Line*) malloc(sizeof(struct Line)))==NULL) return memError();
if (((*vertOut).P1 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if (((*vertOut).P2 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((dvIn = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
if ((p1 = (struct Point*) malloc(sizeof(struct Point)))==NULL) return memError();
//dvIn is the position vector of the incoming muon track
vec_sub((*muonIn).P2, (*muonIn).P1, dvIn);
//setting up vertical vectors for the top and bottom outgoing tracks
vec_copy(muonIn, vertIn);
(*(*vertIn).P1).x = (*(*muonIn).P2).x;
(*(*vertIn).P1).y = (*(*muonIn).P2).y;
vec_copy(muonOut, vertOut);
(*(*vertOut).P2).x = (*(*muonOut).P1).x;
(*(*vertOut).P2).y = (*(*muonOut).P1).y;
thetaX0 = vec_angle(muonIn, vertIn, X_COMPONENT);
thetaX1 = vec_angle(muonOut, vertOut, X_COMPONENT);
if (*params[PARAM_EM_3D]) thetaX = *params[PARAM_SCAT_ANG];
else thetaX = thetaX1 - thetaX0;
thetaY0 = vec_angle(muonIn, vertIn, Y_COMPONENT);
thetaY1 = vec_angle(muonOut, vertOut, Y_COMPONENT);
if (*params[PARAM_EM_3D]) thetaY = 0;
else thetaY = thetaY1 - thetaY0;
Lxy = sqrt(1 + pow(tan(thetaX0), 2) + pow(tan(thetaY0), 2));
*params[PARAM_L] = (*params[PARAM_X_VOXEL_SIZE] * Lxy) / 10;
(*p1).x = (*(*muonIn).P2).x;
(*p1).y = (*(*muonIn).P2).y;
(*p1).z = (*(*muonIn).P2).z;
travel(p1, dvIn, Z_COMPONENT, (*(*muonOut).P1).z);
difX = ((*(*muonOut).P1).x - (*p1).x);
difY = ((*(*muonOut).P1).y - (*p1).y);
deltaX = difX * cos(thetaX0) * Lxy * (cos(thetaX + thetaX0) / cos(thetaX));
deltaY = difY * cos(thetaY0) * Lxy * (cos(thetaY + thetaY0) / cos(thetaY));
(*mu).event = (int) *params[PARAM_CUR_EVENT];
if (*params[PARAM_DOCA]==PARALLEL || isnan(thetaX) || isnan(thetaY)) {
(*mu).dtX = 0;
(*mu).dtY = 0;
(*mu).dX = 0;
(*mu).dY = 0;
} else {
(*mu).dtX = fabs(thetaX * *params[PARAM_MILLIRADIANS]);
(*mu).dtY = fabs(thetaY * *params[PARAM_MILLIRADIANS]);
(*mu).dX = fabs(deltaX / *params[PARAM_UNITS_LENGTH]);
(*mu).dY = fabs(deltaY / *params[PARAM_UNITS_LENGTH]);
}
(*mu).pr2 = pow(*params[PARAM_NOM_MOMENTUM] /
*params[PARAM_MOMENTUM], 2);
*params[PARAM_DTX] = (*mu).dtX;
*params[PARAM_DTY] = (*mu).dtY;
*params[PARAM_DX] = (*mu).dX;
*params[PARAM_DY] = (*mu).dY;
(*mu).a = 0;
(*mu).b = -1;
free(p1);
free(dvIn);
free((*vertIn).P1);
free((*vertIn).P2);
free((*vertOut).P1);
free((*vertOut).P2);
free(vertIn);
free(vertOut);
return CONTINUE;
}
//travel moves a point a long a vector
// p is the Point to move
// V is the position vector
// xyz is what component is being solved for
// end is the value of that component
void travel(struct Point* p, struct Point* v, int xyz, double end_point) {
double t=0, initial_point=0, direction_vec=0;
if (xyz == X_COMPONENT) {
initial_point = (*p).x;
direction_vec = (*v).x;
} else if (xyz == Y_COMPONENT) {
initial_point = (*p).y;
direction_vec = (*v).y;
} else if (xyz == Z_COMPONENT) {
initial_point = (*p).z;
direction_vec = (*v).z;
}
//the t parameter will give the units needed to move along the vector until the desired end_point is
reached
t = (end_point - initial_point) / direction_vec;
(*p).x = (*p).x + (*v).x * t;
(*p).y = (*p).y + (*v).y * t;
(*p).z = (*p).z + (*v).z * t;
return;
}
struct voxel* track(struct Point* start, struct Point* end, struct Point* volume_end, struct muon*
mu, struct voxel* path, double* M, double** params, FILE** fps) {
int xVox, yVox, zVox, curVoxel;
double vx, vy, vz;
double tNew, tOld=0, tLast, L, T;
double xOut = (*start).x, yOut = (*start).y, zOut = (*start).z, depth = 0;
struct Point *new, *old, *v;
struct voxel *track;
track = path;
if ((new = (struct Point*) malloc(sizeof(struct Point)))==NULL) return NULL;
if ((old = (struct Point*) malloc(sizeof(struct Point)))==NULL) return NULL;
if ((v = (struct Point*) malloc(sizeof(struct Point)))==NULL) return NULL;
vec_copy_point(start, new);
vx = (*end).x - (*start).x;
vy = (*end).y - (*start).y;
vz = (*end).z - (*start).z;
while (zOut > (*end).z) {
xVox = (int) floor(((xOut - *params[PARAM_X_MIN]) /
*params[PARAM_X_VOXEL_SIZE]));
yVox = (int) floor(((yOut - *params[PARAM_Y_MIN]) /
*params[PARAM_Y_VOXEL_SIZE]));
zVox = (int) floor(((zOut - *params[PARAM_Z_MIN]) /
*params[PARAM_Z_VOXEL_SIZE]));
//voxel number determined in z direction first, then y, then x
curVoxel = xVox * ((int) *params[PARAM_Y_VOXEL_TOTAL]) * ((int)
*params[PARAM_Z_VOXEL_TOTAL]) + yVox * ((int) *params[PARAM_Z_VOXEL_TOTAL])
+ zVox;
tLast = tOld;
tOld = 100;
//calculations for minimum border of voxel in X direction
xOut = *params[PARAM_X_MIN] + xVox * *params[PARAM_X_VOXEL_SIZE];
tNew = (xOut-(*start).x)/vx;
cliv
if ((tNew < tOld) && (tNew > tLast) && (xOut >= *params[PARAM_X_MIN])) tOld=tNew;
//calculations for maximum border of voxel in X direction
xOut = *params[PARAM_X_MIN] + (xVox+1) * *params[PARAM_X_VOXEL_SIZE];
tNew = (xOut-(*start).x)/vx;
if ((tNew < tOld) && (tNew > tLast) && (xOut <= *params[PARAM_X_MAX])) tOld=tNew;
//calculations for minimum border of voxel in Y direction
yOut = *params[PARAM_Y_MIN] + yVox * *params[PARAM_Y_VOXEL_SIZE];
tNew = (yOut-(*start).y)/vy;
if ((tNew < tOld) && (tNew > tLast) && (yOut >= *params[PARAM_Y_MIN])) tOld=tNew;
//calculations for maximum border of voxel in Y direction
yOut = *params[PARAM_Y_MIN] + (yVox+1)**params[PARAM_Y_VOXEL_SIZE];
tNew = (yOut-(*start).y)/vy;
if ((tNew < tOld) && (tNew > tLast) && (yOut <= *params[PARAM_Y_MAX])) tOld=tNew;
//calculations for minimum border of voxel in Z direction
zOut = *params[PARAM_Z_MIN] + zVox**params[PARAM_Z_VOXEL_SIZE];
tNew = (zOut-(*start).z)/vz;
if ((tNew <tOld) && (tNew > tLast) && (zOut <= *params[PARAM_Z_MAX])) tOld=tNew;
tNew = tOld + TRACK_PUSH;
if (*params[PARAM_PRECISE_L] || *params[PARAM_PRECISE_T]) {
vec_copy_point(new, old);
(*new).x = (*start).x + vx*tOld;
(*new).y = (*start).y + vy*tOld;
(*new).z = (*start).z + vz*tOld;
}
if ((zVox < *params[PARAM_Z_VOXEL_TOTAL]) && (zOut < (*start).z) && (zVox >= 0)
&& (curVoxel!= *params[PARAM_PREV_VOXEL])) {
track=(*track).nextVoxel;
(*track).ID = curVoxel;
if (volume_end!=NULL) (*mu).b = (*mu).b + 1;
else (*mu).a = (*mu).a + 1;
//L calculation
if (*params[PARAM_PRECISE_L]) {
vec_sub(new, old, v);
L = vec_norm(v);
} else L = *params[PARAM_L];
//T calculation
if (*params[PARAM_PRECISE_T]) {
if (volume_end==NULL) {
vec_sub(end, new, v);
T = vec_norm(v);
} else {
vec_sub(end, new, v);
T = vec_norm(v);
vec_sub(volume_end, new, v);
T = T + vec_norm(v);
}
} else T = *params[PARAM_L] * zVox;
L = L / *params[PARAM_UNITS_LENGTH];
T = T / *params[PARAM_UNITS_LENGTH];
(*track).wt = L;
(*track).wtX = (pow(L, 2)/2) + L*T;
(*track).wX = (pow(L, 3)/3) + (pow(L, 2)*T) + (L*pow(T, 2));
if (((*track).nextVoxel = (struct voxel*) malloc(sizeof(struct voxel)))==NULL) return NULL;
}
xOut = (*start).x + vx*tNew;
yOut = (*start).y + vy*tNew;
zOut = (*start).z + vz*tNew;
*params[PARAM_PREV_VOXEL] = curVoxel;
}
free((*track).nextVoxel);
(*track).nextVoxel= NULL;
free(v);
free(old);
free(new);
return track;
}
