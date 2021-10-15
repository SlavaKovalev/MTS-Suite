#include "mts.h"
#include "mtserr.h"
#include "preprocessing.h"
#include "em.h"
#include <getopt.h>
void set_default_parameters(double**);
int get_opts(int, double**, char**, FILE**);
int getopt_file(FILE*, struct option*);
FILE* fopen_ext(char*, char*, char*);
char* create_file_ext(char*, int, int, double**);
int main (int argc, char **argv) {
int i, errCode;
double *params[MAX_PARAMS], *lambda, *lambdaMed, *M;
FILE *filepointers[MAX_FILEPOINTERS];
struct muon* head;
for (i=0;i<MAX_PARAMS;i++) if ((params[i] = (double*) malloc(sizeof(double)))==NULL)
return memError();
for (i=0;i<MAX_FILEPOINTERS;i++) if ((filepointers[i] = (FILE*)
malloc(sizeof(FILE)))==NULL) return memError();
if ((head=(struct muon*) malloc(sizeof(struct muon)))==NULL) return memError();
set_default_parameters(params);
if (get_opts(argc, params, argv, filepointers)!=CONTINUE) return
ERR_CODE_COMMAND_LINE;
free(optarg);
preprocessing(lambda, M, head, params, filepointers);
if (*params[PARAM_EM] && !*params[PARAM_EM_ONLINE]) em(lambda, lambdaMed,
head, params, filepointers);
for (i=0;i<MAX_PARAMS;i++) free(params[i]);
for (i=0;i<MAX_FILEPOINTERS;i++) free(filepointers[i]);
fprintf(stderr, "\n");
return errCode;
}
void set_default_parameters (double** params) {
int i;
//all paramaters are intially zero unless otherwise set in this modle
for (i=0;i<MAX_PARAMS;i++) *params[i] = 0;
*params[PARAM_X_LENGTH] = 4000;
*params[PARAM_Y_LENGTH] = 4000;
*params[PARAM_Z_LENGTH] = 3000;
*params[PARAM_X_VOXEL_SIZE] = 100;
*params[PARAM_Y_VOXEL_SIZE] = 100;
*params[PARAM_Z_VOXEL_SIZE] = 100;
*params[PARAM_X_VOXEL_TOTAL] = *params[PARAM_X_LENGTH] /
*params[PARAM_X_VOXEL_SIZE];
*params[PARAM_Y_VOXEL_TOTAL] = *params[PARAM_Y_LENGTH] /
*params[PARAM_Y_VOXEL_SIZE];
*params[PARAM_Z_VOXEL_TOTAL] = *params[PARAM_Z_LENGTH] /
*params[PARAM_Z_VOXEL_SIZE];
*params[PARAM_ALL_VOXELS] = *params[PARAM_X_VOXEL_TOTAL] *
*params[PARAM_Y_VOXEL_TOTAL] * *params[PARAM_Z_VOXEL_TOTAL];
*params[PARAM_X_MIN] = *params[PARAM_X_LENGTH] / -2;
*params[PARAM_Y_MIN] = *params[PARAM_Y_LENGTH] / -2;
*params[PARAM_Z_MIN] = *params[PARAM_Z_LENGTH] / -2;
*params[PARAM_X_MAX] = *params[PARAM_X_LENGTH] / 2;
*params[PARAM_Y_MAX] = *params[PARAM_Y_LENGTH] / 2;
*params[PARAM_Z_MAX] = *params[PARAM_Z_LENGTH] / 2;
*params[PARAM_CONTINUE] = 1;
*params[PARAM_ITERATIONS] = 100;
*params[PARAM_INIT_LAMBDA] = 0.1;
*params[PARAM_NOM_MOMENTUM] = 3;
*params[PARAM_MIN_MOMENTUM] = 1;
*params[PARAM_MILLIRADIANS] = 1;
*params[PARAM_EM] = 0;
*params[PARAM_EM_BIN_SIZE] = 10000;
*params[PARAM_EM_BINS] = 100;
*params[PARAM_UNITS_LENGTH] = 10;
return;
}
int get_opts (int argc, double** params, char** argv, FILE** fps) {
char fnOut[MAX_FILENAME], extAvg[MAX_FILENAME], extMed[MAX_FILENAME];
int option, argvIndex=0, i=0, places=0, zeros;
FILE* config = NULL;
static struct option long_options[] = {
{LONG_BINS, required_argument, NULL, SHORT_BINS},
{LONG_BIN_SIZE, required_argument, NULL, SHORT_BIN_SIZE},
{LONG_C_PRINT, required_argument, NULL, SHORT_C_PRINT},
{LONG_CUTOFF, required_argument, NULL, SHORT_CUTOFF},
{LONG_DETAILS, no_argument, NULL, SHORT_DETAILS},
{LONG_DIST, required_argument, NULL, SHORT_DIST},
{LONG_EM, optional_argument, NULL, SHORT_EM},
{LONG_EM_WEIGHTED, no_argument, NULL, SHORT_EM_WEIGHTED},
{LONG_HELP, no_argument, NULL, SHORT_HELP},
{LONG_INPUT, required_argument, NULL, SHORT_INPUT},
{LONG_ITERATIONS, required_argument, NULL, SHORT_ITERATIONS},
{LONG_LAMBDA, required_argument, NULL, SHORT_LAMBDA},
{LONG_MILLIRADIANS, no_argument, NULL, SHORT_MILLIRADIANS},
{LONG_MOM_HIGH_CUT, required_argument, NULL, SHORT_MOM_HIGH_CUT},
{LONG_MOM_LOW_CUT, required_argument, NULL, SHORT_MOM_LOW_CUT},
{LONG_NOMINAL, required_argument, NULL, SHORT_NOMINAL},
{LONG_ONLINE, required_argument, NULL, SHORT_ONLINE},
{LONG_OUTPUT, required_argument, NULL, SHORT_OUTPUT},
{LONG_POCA, required_argument, NULL, SHORT_POCA},
{LONG_PRECISE_L, no_argument, NULL, SHORT_PRECISE_L},
{LONG_PRECISE_T, no_argument, NULL, SHORT_PRECISE_T},
{LONG_STD, optional_argument, NULL, SHORT_STD},
{LONG_UNITS_LENGTH, required_argument, NULL, SHORT_UNITS_LENGTH},
{LONG_X, required_argument, NULL, SHORT_X},
{LONG_Y, required_argument, NULL, SHORT_Y},
{LONG_Z, required_argument, NULL, SHORT_Z},
{0, 0, 0, 0}
};
if (argc<=1) return commandError(SHORT_HELP, SHORT_HELP);
if (argv[1][0]!='-') if ((config=fopen(argv[1], "r"))==NULL) return
commandError(SHORT_HELP, SHORT_HELP);
while (1) {
if (config==NULL) option = getopt_long (argc, argv, SHORT_OPTIONS, long_options,
&argvIndex);
else option = getopt_file (config, long_options);
if (option == -1) break;
switch (option) {
case SHORT_BINS:
*params[PARAM_EM_BINS] = atof(optarg);
break;
case SHORT_BIN_SIZE:
*params[PARAM_EM_BIN_SIZE] = atof(optarg);
break;
case SHORT_C_PRINT:
*params[PARAM_C_PRINT] = 1;
break;
case SHORT_CUTOFF:
*params[PARAM_CUTOFF_ANGLE] = atof(optarg);
break;
case SHORT_DETAILS:
*params[PARAM_DETAILS] = 1;
break;
case SHORT_DIST:
if (optarg!=NULL) fps[FP_DIST_OUT] = fopen(optarg, "w");
else fps[FP_DIST_OUT] = fopen("default.dist", "w");
*params[PARAM_OUT] = *params[PARAM_DIST] = 1;
break;
case SHORT_EM:
*params[PARAM_EM] = 2;
if (optarg==NULL) *params[PARAM_EM_AVERAGE] = 1;
else if (strcmp(optarg, LONG_EM_AVERAGE)==0) *params[PARAM_EM_AVERAGE] =
1;
else if (strcmp(optarg, LONG_EM_MEDIAN)==0) *params[PARAM_EM_MEDIAN] = 1;
else if (strcmp(optarg, LONG_EM_3D)==0) *params[PARAM_EM_3D] = 1;
break;
case SHORT_EM_WEIGHTED:
*params[PARAM_EM_WEIGHTED] = 1;
break;
case SHORT_HELP:
return (commandError(option, option));
break;
case SHORT_INPUT:
if ((fps[FP_IN] = fopen(optarg, "r"))==NULL) return fileError(optarg);
strcpy(fnOut, optarg);
break;
case SHORT_ITERATIONS:
*params[PARAM_ITERATIONS] = atof(optarg);
break;
case SHORT_LAMBDA:
*params[PARAM_INIT_LAMBDA] = atof(optarg);
for (i=0; optarg[i]!='.' && optarg[i]!='\0'; i++);
for (i=i+1, places=0, zeros=0; optarg[i]=='0'; places++, zeros++, i++);
for (;optarg[i]!='\0'; places++, i++);
break;
case SHORT_MILLIRADIANS:
*params[PARAM_MILLIRADIANS] = 1000;
break;
case SHORT_MOM_HIGH_CUT:
*params[PARAM_MOM_HIGH_CUT] = atof(optarg);
*params[PARAM_MOM_CUT] = 1;
break;
case SHORT_MOM_LOW_CUT:
*params[PARAM_MOM_LOW_CUT] = atof(optarg);
*params[PARAM_MOM_CUT] = 1;
break;
case SHORT_NOMINAL:
*params[PARAM_NOM_MOMENTUM] = atof(optarg);
break;
case SHORT_ONLINE:
*params[PARAM_EM_ONLINE] = atof(optarg);
break;
case SHORT_OUTPUT:
if (optarg!=NULL) fps[FP_OP_OUT] = fopen(optarg, "w");
else fps[FP_OP_OUT] = fopen("default.opout", "w");
*params[PARAM_OUT] = *params[PARAM_OP_OUT] = 1;
break;
case SHORT_POCA:
if (optarg!=NULL) fps[FP_POCA_OUT] = fopen(optarg, "w");
else fps[FP_POCA_OUT] = fopen("default.poca", "w");
*params[PARAM_OUT] = *params[PARAM_POCA] = 1;
break;
case SHORT_PRECISE_L:
*params[PARAM_PRECISE_L] = 1;
break;
case SHORT_PRECISE_T:
*params[PARAM_PRECISE_T] = 1;
break;
case SHORT_STD:
if (optarg!=NULL) fps[FP_STD_OUT] = fopen(optarg, "w");
else fps[FP_STD_OUT] = NULL;
*params[PARAM_STD] = 1;
break;
case SHORT_UNITS_LENGTH:
*params[PARAM_UNITS_LENGTH] = atof(optarg);
break;
case SHORT_X:
*params[PARAM_X_VOXEL_SIZE] = atof(optarg);
break;
case SHORT_Y:
*params[PARAM_Y_VOXEL_SIZE] = atof(optarg);
break;
case SHORT_Z:
*params[PARAM_Z_VOXEL_SIZE] = atof(optarg);
break;
case '?':
return (commandError(optopt, optopt));
break;
default:
printf ("\n\nIf you are seeing this then quantum mechanics is for real: %c", option);
}
free(optarg);
}
if (config!=NULL) fclose(config);
if (*params[PARAM_EM_3D]) *params[PARAM_EM] = 1;
if (*params[PARAM_EM_AVERAGE]) {
sprintf(extAvg, "avg");
if ((fps[FP_OUT_AVG] = fopen_ext(fnOut, create_file_ext(extAvg, places, zeros, params),
"w"))==NULL) return fileError(fnOut);
}
if (*params[PARAM_EM_MEDIAN]) {
sprintf(extMed, "med%dbins%dsize", (int) *params[PARAM_EM_BINS], (int)
*params[PARAM_EM_BIN_SIZE]);
if ((fps[FP_OUT_MED] = fopen_ext(fnOut, create_file_ext(extMed, places, zeros, params),
"w"))==NULL) return fileError(fnOut);
}
if (*params[PARAM_MILLIRADIANS]==1) *params[PARAM_INIT_LAMBDA] =
*params[PARAM_INIT_LAMBDA] / 1000000;
if (*params[PARAM_UNITS_LENGTH]!=10) *params[PARAM_INIT_LAMBDA] =
*params[PARAM_INIT_LAMBDA] / (10 / *params[PARAM_UNITS_LENGTH]);
return CONTINUE;
}
int getopt_file(FILE* config, struct option* options) {
int i;
char line[MAX_LINE], *token;
optarg = (char*) malloc(MAX_LINE);
while (1) {
if (fgets(line, MAX_LINE-1, config)==NULL) return -1;
if ((token = strtok(line, DELIMS))==NULL) continue;
if (token[0]=='#') continue;
for (i=0; options[i].val!=0; i++) {
if (strcmp(token, options[i].name)==0) {
if ((token = strtok(NULL, DELIMS))!=NULL) strcpy(optarg, token);
else {
free(optarg);
optarg=NULL;
}
return options[i].val;
}
}
return '?';
}
}
FILE* fopen_ext(char* fn, char* ext, char* param) {
char newFn[200];
int i, j;
strcpy(newFn, fn);
j = 0;
i = strcspn(fn, ".");
while (ext[j]!='\0') newFn[++i] = ext[j++];
newFn[++i]='\0';
return fopen(newFn, param);
}
char* create_file_ext(char* ext, int places, int zeros, double** params) {
int i, after_zero = (int) (ceil((((double) (*params[PARAM_INIT_LAMBDA] - ((int)
*params[PARAM_INIT_LAMBDA]))) * pow(10, ((double) places)))));
if (*params[PARAM_EM_3D]) sprintf(ext, "%s-3D", ext);
if (*params[PARAM_X_VOXEL_SIZE]) sprintf(ext, "%s-vox%dcm", ext, (int)
(*params[PARAM_X_VOXEL_SIZE]/10));
if (*params[PARAM_ITERATIONS]) sprintf(ext, "%s-itr%d", ext, (int)
*params[PARAM_ITERATIONS]);
if (*params[PARAM_INIT_LAMBDA]) {
sprintf(ext, "%s-lam%i.", ext, (int) *params[PARAM_INIT_LAMBDA]);
for (i=0; i<zeros; i++) sprintf(ext, "%s0", ext);
sprintf(ext, "%s%i", ext, after_zero);
}
if (*params[PARAM_CUTOFF_ANGLE]) sprintf(ext, "%s-angcut%d", ext, (int)
*params[PARAM_CUTOFF_ANGLE]);
if (*params[PARAM_EM_WEIGHTED]) sprintf(ext, "%s-weight%d", ext, (int)
*params[PARAM_EM_WEIGHTED]);
if (*params[PARAM_EM_ONLINE]) sprintf(ext, "%s-inc%d", ext, (int)
*params[PARAM_EM_ONLINE]);
if (*params[PARAM_PRECISE_L]) sprintf(ext, "%s-pL", ext);
if (*params[PARAM_PRECISE_T]) sprintf(ext, "%s-pT", ext);
if (*params[PARAM_MOM_CUT] = 1) sprintf(ext, "%s-momcutHi%dLo%d", ext, (int)
*params[PARAM_MOM_HIGH_CUT], (int) *params[PARAM_MOM_LOW_CUT]);
if (*params[PARAM_UNITS_LENGTH]) sprintf(ext, "%s-units%fmm", ext,
*params[PARAM_UNITS_LENGTH]);
if (*params[PARAM_MILLIRADIANS]==1000) sprintf(ext, "%s-mrad", ext);
return ext;
}