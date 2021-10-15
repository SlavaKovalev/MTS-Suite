//define constants for markers in poca input file
#define TOP_START 'a'
#define TOP_END 'b'
#define BOTTOM_START 'c'
#define BOTTOM_END 'd'
#define EVENTS_START 'e'
int get_geant_input(struct Point**, struct Point**, double**, FILE**);
void write_lambda(double*, double*, struct muon*, double**, FILE*);
void write_optional(struct Point*, struct Point*, struct Point*, double**, FILE**);
//mtsio.c
#include "mts.h"
#include "mtsio.h"
#include "mtserr.h"
#include "vecfunc.h"
//takes input from a file and parses it into the incoming/outgoing points
//and the scattering angle
int get_geant_input(struct Point** incPoints, struct Point** outPoints, double** params, FILE**
fps) {
//coordinates from file
char *token, input[MAX_LINE];
int i = 0;
struct Point* dummy;
while (1) {
//get next params[PARAM_LINE] from file
if ((fgets(input, MAX_LINE, fps[FP_IN])==NULL)) return END;
(*params[PARAM_LINE])++;
//if the next points aren't from top detector disregard this event
if (input[0]!=TOP_START) {
if (input[0]==EVENTS_START) {
if ((fgets(input, MAX_LINE, fps[FP_IN])==NULL)) return END;
}
continue;
}
//take care of momentum
if ((fgets(input, MAX_LINE, fps[FP_IN])==NULL)) return END;
(*params[PARAM_LINE])++;
if (!isdigit(input[0]) && input[0]!='-') break;
*params[PARAM_MOMENTUM] = atof(strtok(input, DELIMS));
i = 0;
//read in incoming muon tracks
while (input[0]!=TOP_END) {
if ((fgets(input, MAX_LINE, fps[FP_IN])==NULL)) return END;
(*params[PARAM_LINE])++;
if (!isdigit(input[0]) && input[0]!='-') break;
dummy = incPoints[i];
if ((token=strtok( input, DELIMS ))!=NULL) (*dummy).x = atof(token);
else return formatError("inputFile", input, *params[PARAM_LINE]);
if ((token=strtok( NULL, DELIMS ))!=NULL) (*dummy).y = atof(token);
else return formatError("inputFile", input, *params[PARAM_LINE]);
if ((token=strtok( NULL, DELIMS ))!=NULL) (*dummy).z = atof(token);
else return formatError("inputFile", input, *params[PARAM_LINE]);
i++;
}
*params[PARAM_INC_DECT] = i;
//if less than 2 detectors were hit disregard this event
if (i<2) continue;
if (input[0]!=TOP_END) return formatError("inputFile", input, *params[PARAM_LINE]);
if ((fgets(input, MAX_LINE, fps[FP_IN])==NULL)) return END;
(*params[PARAM_LINE])++;
//if the next points aren't from bottom detector disregard this event
if (input[0]!=BOTTOM_START) continue;
i=0;
//read in outgoing muon tracks
while (input[0]!=BOTTOM_END) {
if ((fgets(input, MAX_LINE, fps[FP_IN])==NULL)) return END;
(*params[PARAM_LINE])++;
if (!isdigit(input[0]) && input[0]!='-') break;
dummy = outPoints[i];
if ((token=strtok( input, DELIMS ))!=NULL) (*dummy).x = atof(token);
else return formatError("inputFile", input, *params[PARAM_LINE]);
if ((token=strtok( NULL, DELIMS ))!=NULL) (*dummy).y = atof(token);
else return formatError("inputFile", input, *params[PARAM_LINE]);
if ((token=strtok( NULL, DELIMS ))!=NULL) (*dummy).z = atof(token);
else return formatError("inputFile", input, *params[PARAM_LINE]);
i++;
}
*params[PARAM_OUT_DECT] = i;
//if less than 2 detectors were hit disregard this event
if (i<2) continue;
if (input[0]!=BOTTOM_END) return formatError("inputFile", input, *params[PARAM_LINE]);
break;
}
return CONTINUE;
}
void write_lambda(double* lambda, double* M, struct muon* mu, double** params, FILE* out) {
int i, x, y, z, startVoxel=0;
struct voxel* tempVoxel;
for (i=0; i < *params[PARAM_ALL_VOXELS]; i++) {
if (M[i]==0) lambda[i]=0;
if (*params[PARAM_MILLIRADIANS]==1) lambda[i] = lambda[i]*1000000;
if (*params[PARAM_UNITS_LENGTH]!=10) lambda[i] = lambda[i] * (10 /
*params[PARAM_UNITS_LENGTH]);
x = (int) floor(i/(*params[PARAM_Y_VOXEL_TOTAL] *
*params[PARAM_Z_VOXEL_TOTAL]));
y = (int) floor((i - (x * *params[PARAM_Y_VOXEL_TOTAL] *
*params[PARAM_Z_VOXEL_TOTAL]))/ *params[PARAM_Z_VOXEL_TOTAL]);
z = (int) floor(i - (x * *params[PARAM_Y_VOXEL_TOTAL] *
*params[PARAM_Z_VOXEL_TOTAL]) - (y * *params[PARAM_Z_VOXEL_TOTAL]));
x = x * *params[PARAM_X_VOXEL_SIZE] + *params[PARAM_X_MIN] +
(*params[PARAM_X_VOXEL_SIZE]/2);
y = y * *params[PARAM_Y_VOXEL_SIZE] + *params[PARAM_Y_MIN] +
(*params[PARAM_Y_VOXEL_SIZE]/2);
z = z * *params[PARAM_Z_VOXEL_SIZE] + *params[PARAM_Z_MIN] +
(*params[PARAM_Z_VOXEL_SIZE]/2);
fprintf(out, "%d %d %d %f %d %f\n", x, y, z, lambda[i], i, M[i]);

}
return;
}
void write_optional(struct Point* in, struct Point* scat, struct Point* out, double** params, FILE**
fps) {
if (*params[PARAM_OP_OUT]) {
fprintf(fps[FP_OP_OUT], "Event: %f\n", *params[PARAM_CUR_EVENT]);
fprintf(fps[FP_OP_OUT], "IN: %f %f %f\n", (*in).x, (*in).y, (*in).z);
fprintf(fps[FP_OP_OUT], "POCA & DOCA: %f %f %f %f\n", (*scat).x, (*scat).y, (*scat).z,
*params[PARAM_DOCA]);
fprintf(fps[FP_OP_OUT], "OUT: %f %f %f\n", (*out).x, (*out).y, (*out).z);
fprintf(fps[FP_OP_OUT], "ANGLE(radians/degrees): %f / %f\n",
*params[PARAM_SCAT_ANG], vec_rad_to_deg(*params[PARAM_SCAT_ANG]));
fprintf(fps[FP_OP_OUT], "EM (dtX, dtY, dX, dY, mom, L): %f %f %f %f %f %f\n\n",
*params[PARAM_DTX], *params[PARAM_DTY], *params[PARAM_DX],
*params[PARAM_DY], *params[PARAM_MOMENTUM], *params[PARAM_L]);
}
if (*params[PARAM_POCA])
fprintf(fps[FP_POCA_OUT], "%f %f %f %f\n", (*scat).x, (*scat).y, (*scat).z,
vec_rad_to_deg(*params[PARAM_SCAT_ANG]));
if (*params[PARAM_DIST] && (*params[PARAM_DOCA]!=PARALLEL)) {
fprintf(fps[FP_DIST_OUT], "%f %f", *params[PARAM_SCAT_ANG],
vec_rad_to_deg(*params[PARAM_SCAT_ANG]));
fprintf(fps[FP_DIST_OUT], " %f %f %f %f\n", *params[PARAM_DTX],
vec_rad_to_deg(*params[PARAM_DTX]), *params[PARAM_DTY],
vec_rad_to_deg(*params[PARAM_DTY]));
}
return;
}