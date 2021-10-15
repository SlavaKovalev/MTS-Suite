#include "mts.h"
#include "em.h"
int em (double* lambdaMed, double* lambda, struct muon* head, double** params, FILE** fps) {
fprintf(stderr, "\n\n%s EM %s\n\n", BANNER, BANNER);
int i, j, k, bin, allVoxels, iterations, done=0;
double lambdaTemp, *M, *C, **Cbin, **Cn, bCount;
time_t startLocal, endLocal;
struct muon *tempMuon;
allVoxels = *params[PARAM_ALL_VOXELS];
if ((M=(double*) malloc(allVoxels*sizeof(double)))==NULL) return memError();
if (*params[PARAM_EM_AVERAGE]) {
if (!*params[PARAM_EM_ONLINE] || ((*params[PARAM_CUR_EVENT] /
*params[PARAM_EM_ONLINE])==1)) {
if ((lambda=(double*) malloc(allVoxels*sizeof(double)))==NULL) return memError();
}
if ((C=(double*) malloc(allVoxels*sizeof(double)))==NULL) return memError();
}
if (*params[PARAM_EM_MEDIAN]) {
if (!*params[PARAM_EM_ONLINE] || ((*params[PARAM_CUR_EVENT] /
*params[PARAM_EM_ONLINE])==1)) {
if ((lambdaMed=(double*) malloc(allVoxels*sizeof(double)))==NULL) return memError();
}
if ((Cbin=(double**) malloc(allVoxels*sizeof(double)))==NULL) return memError();
if ((Cn=(double**) malloc(allVoxels*sizeof(double)))==NULL) return memError();
for (i=0; i<*params[PARAM_ALL_VOXELS]; i++) {
if ((Cbin[i]=(double*) malloc(*params[PARAM_EM_BINS]*sizeof(double)))==NULL)
return memError();
if ((Cn[i]=(double*) malloc(*params[PARAM_EM_BINS]*sizeof(double)))==NULL) return
memError();
}
}
for (i=0;i<*params[PARAM_ALL_VOXELS];i++) {
M[i]=0;
if (!*params[PARAM_EM_ONLINE] || ((*params[PARAM_CUR_EVENT] /
*params[PARAM_EM_ONLINE])==1)) {
if (*params[PARAM_EM_AVERAGE]) lambda[i]=*params[PARAM_INIT_LAMBDA];
}
if (*params[PARAM_EM_MEDIAN]) {
if (!*params[PARAM_EM_ONLINE] || ((*params[PARAM_CUR_EVENT] /
*params[PARAM_EM_ONLINE])==1)) {
lambdaMed[i]=*params[PARAM_INIT_LAMBDA];
}
for (j=0; j<*params[PARAM_EM_BINS]; j++) {
Cbin[i][j]=0;
Cn[i][j]=0;
}
}
}
char fnVoxel[100];
fprintf(stderr, "\nEntering EM Main Loop...\n");
time(&startLocal);
iterations = *params[PARAM_ITERATIONS];
for (i=0; i<iterations; i++) {
fprintf(stderr, "\n\tIteration %d of %d ", (i+1), iterations);
if (*params[PARAM_EM_AVERAGE]) memset(C, '\0', sizeof(double) *
*params[PARAM_ALL_VOXELS]);
tempMuon = head;
while ((tempMuon=(*tempMuon).nextMuon)!=NULL) {
compute_v(tempMuon, M, lambda, lambdaMed, i, params);
compute_c(tempMuon, C, lambdaMed, i, Cbin, Cn, params);
}
for (j=0; j<allVoxels; j++) {
if (M[j]!=0) {
if (*params[PARAM_EM_AVERAGE]) {
lambdaTemp = lambda[j];
lambda[j] = lambda[j] + pow(lambda[j], 2) * (1/M[j]) * C[j];
if (lambda[j]<=0) lambda[j]=*params[PARAM_INIT_LAMBDA];
}
if (*params[PARAM_EM_MEDIAN]) {
bin = 0;
for (k=0, bCount=0; k<*params[PARAM_EM_BINS]; k++) {
bCount = bCount + Cn[j][k];
if (bCount > (M[j]/2)) {
bin = k;
break;
}
}
lambdaMed[j] = 0.5 * (Cbin[j][bin]/Cn[j][bin]);
if (lambdaMed[j]<=0) lambdaMed[j]=*params[PARAM_INIT_LAMBDA];
for (k=0; k<*params[PARAM_EM_BINS]; k++) {
Cbin[j][k]=0;
Cn[j][k]=0;
}
}
}
}
}
time(&endLocal);
fprintf(stderr, "\n\nEM Program ran for %f seconds\n\n", difftime(endLocal, startLocal));
if (!*params[PARAM_EM_ONLINE] || (fmod(*params[PARAM_CUR_EVENT],
*params[PARAM_EM_ONLINE])!=0)) {
if (*params[PARAM_EM_AVERAGE]) write_lambda(lambda, M, NULL, params,
fps[FP_OUT_AVG]);
if (*params[PARAM_EM_MEDIAN]) write_lambda(lambdaMed, M, NULL, params,
fps[FP_OUT_MED]);
}
free(M);
if (*params[PARAM_EM_AVERAGE]) free(C);
if (!*params[PARAM_EM_ONLINE] || (fmod(*params[PARAM_CUR_EVENT],
*params[PARAM_EM_ONLINE])!=0)) {
if (*params[PARAM_EM_AVERAGE]) free(lambda);
if (*params[PARAM_EM_MEDIAN]) free(lambdaMed);
}
if (*params[PARAM_EM_MEDIAN]) {
for (i=0;i<*params[PARAM_ALL_VOXELS];i++) {
free(Cbin[i]);
free(Cn[i]);
}
}
free(Cbin);
free(Cn);
return;
}
void compute_v (struct muon* mu, double* M, double* lambda, double* lambdaMed, int iteration,
double** params) {
double det, sigma0, sigma1, sigma2, I[4], c=0, weight=1, noInv=0, ID0, ID1, ID2, ID3;
struct voxel *track;
if (*params[PARAM_EM_AVERAGE]) {
(*mu).sigma[0]=0;
(*mu).sigma[1]=0;
(*mu).sigma[2]=0;
}
if (*params[PARAM_EM_MEDIAN]) {
(*mu).sigmaMed[0]=0;
(*mu).sigmaMed[1]=0;
(*mu).sigmaMed[2]=0;
}
track = (*mu).muonTrack;
do {
c = c + 1;
if (*params[PARAM_EM_WEIGHTED]) weight = calc_voxel_weight(mu, (*track).ID, c,
params);
if (iteration==0) M[(*track).ID] = M[(*track).ID] + 1;
if (*params[PARAM_EM_AVERAGE]) {
(*mu).sigma[0] = (*mu).sigma[0] + weight * (*track).wt * lambda[(*track).ID];
(*mu).sigma[1] = (*mu).sigma[1] + weight * (*track).wtX * lambda[(*track).ID];
(*mu).sigma[2] = (*mu).sigma[2] + weight * (*track).wX * lambda[(*track).ID];
}
if (*params[PARAM_EM_MEDIAN]) {
(*mu).sigmaMed[0] = (*mu).sigmaMed[0] + weight * (*track).wt * lambdaMed[(*track).ID];
(*mu).sigmaMed[1] = (*mu).sigmaMed[1] + weight * (*track).wtX *
lambdaMed[(*track).ID];
(*mu).sigmaMed[2] = (*mu).sigmaMed[2] + weight * (*track).wX *
lambdaMed[(*track).ID];
}
} while ((track=(*track).nextVoxel)!=NULL);
if (*params[PARAM_EM_AVERAGE]) {
(*mu).sigma[0] = (*mu).sigma[0] * (*mu).pr2;
(*mu).sigma[1] = (*mu).sigma[1] * (*mu).pr2;
(*mu).sigma[2] = (*mu).sigma[2] * (*mu).pr2;
if ((det = ((*mu).sigma[0] * (*mu).sigma[2]) - ((*mu).sigma[1] * (*mu).sigma[1]))==0) {
fprintf(stderr, "\nError Singular Matrix: Event %d\n\n", (*mu).event);
det = 0;
}
sigma0 = (*mu).sigma[0];
sigma1 = (*mu).sigma[1];
sigma2 = (*mu).sigma[2];
(*mu).sigma[0] = (1/det) * sigma2;
(*mu).sigma[1] = (1/det) * sigma1 * -1;
(*mu).sigma[2] = (1/det) * sigma0;
}
if (*params[PARAM_EM_MEDIAN]) {
(*mu).sigmaMed[0] = (*mu).sigmaMed[0] * (*mu).pr2;
(*mu).sigmaMed[1] = (*mu).sigmaMed[1] * (*mu).pr2;
(*mu).sigmaMed[2] = (*mu).sigmaMed[2] * (*mu).pr2;
if ((det = ((*mu).sigmaMed[0] * (*mu).sigmaMed[2]) - ((*mu).sigmaMed[1] *
(*mu).sigmaMed[1]))==0) {
fprintf(stderr, "\nError Singular Matrix: Event %d\n\n", (*mu).event);
det = 0;
}
sigma0 = (*mu).sigmaMed[0];
sigma1 = (*mu).sigmaMed[1];
sigma2 = (*mu).sigmaMed[2];
(*mu).sigmaMed[0] = (1/det) * sigma2;
(*mu).sigmaMed[1] = (1/det) * sigma1 * -1;
(*mu).sigmaMed[2] = (1/det) * sigma0;
}
return;
}
double calc_voxel_weight(struct muon* mu, unsigned int ID, double c, double** params) {
double a, b, d;
if ((*mu).b==0) return 1;
a = (double) (*mu).a; //voxels after scatter voxel
b = (double) (*mu).b; //voxels before scatter voxel
d = b + 1; //voxel of scattering
if (*params[PARAM_EM_WEIGHTED]==1) {
if (c!=d) return 0;
else return 1;
}
if (*params[PARAM_EM_WEIGHTED]==2) {
if (a==0 || b==0) return 1;
if (c<d) return ((b - (d - c)) / b);
else if (c==d) return 1;
else return ((a - (c - d)) / a);
}
return 1;
}
void compute_c (struct muon* mu, double* C, double* lambdaMed, int iteration, double** Cbin,
double** Cn, double** params) {
int ID, ibin, i;
double dtX, dtY, dX, dY, pr2, fbin, bin_size = 0, neg_bin_size = 2, big_bin_size = 10000;
double wt, wtX, wX, v11, v12, v22, v11m, v12m, v22m, a, b, c, ax, bx, cx, ay, by, cy, oldC, newC;
struct voxel *track, *trackTemp, *trackMed;
dtX = (*mu).dtX;
dtY = (*mu).dtY;
dX = (*mu).dX;
dY = (*mu).dY;
pr2 = (*mu).pr2;
v11 = (*mu).sigma[0];
v12 = (*mu).sigma[1];
v22 = (*mu).sigma[2];
v11m = (*mu).sigmaMed[0];
v12m = (*mu).sigmaMed[1];
v22m = (*mu).sigmaMed[2];
track = (*mu).muonTrack;
do {
ID = (*track).ID;
wt = (*track).wt;
wtX = (*track).wtX;
wX = (*track).wX;
if (*params[PARAM_EM_AVERAGE]) {
ax = dtX * v11 + dX * v12;
bx = dtX * v12 + dX * v22;
cx = ax * bx;
ay = dtY * v11 + dY * v12;
by = dtY * v12 + dY * v22;
cy = ay * by;
a = (pow(ax, 2) + pow(ay, 2)) / *params[PARAM_EM];
b = (pow(bx, 2) + pow(by, 2)) / *params[PARAM_EM];
c = (cx + cy) / *params[PARAM_EM];
oldC = C[ID];
newC = ((pr2 * (a - v11) * wt) + (pr2 * 2 * (c - v12) * wtX) + (pr2 * (b - v22) * wX));
C[ID] = C[ID] + newC;
}
if (*params[PARAM_EM_MEDIAN]) {
ax = dtX * v11m + dX * v12m;
bx = dtX * v12m + dX * v22m;
cx = ax * bx;
ay = dtY * v11m + dY * v12m;
by = dtY * v12m + dY * v22m;
cy = ay * by;
a = (pow(ax, 2) + pow(ay, 2)) / *params[PARAM_EM];
b = (pow(bx, 2) + pow(by, 2)) / *params[PARAM_EM];
c = (cx + cy) / *params[PARAM_EM];
newC = ((pr2 * (a - v11m) * wt) + (pr2 * 2 * (c - v12m) * wtX) + (pr2 * (b - v22m) * wX));
bin_size = *params[PARAM_EM_BIN_SIZE];
fbin = (newC / bin_size) + (*params[PARAM_EM_BINS]/2) + 1;
if (fbin < 0) fbin = 0;
else if (fbin >= *params[PARAM_EM_BINS]) fbin = *params[PARAM_EM_BINS] - 1;
ibin = ((int) floor(fbin));
Cbin[ID][ibin] = Cbin[ID][ibin] + (2*lambdaMed[ID] + pow(lambdaMed[ID], 2) * newC);
Cn[ID][ibin] = Cn[ID][ibin] + 1;
}
} while ((track=(*track).nextVoxel)!=NULL);
return;
}