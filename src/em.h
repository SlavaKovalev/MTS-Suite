//define standard extensions
#define EXT_EM "em"
#define EXT_LT "lt"
#define EXT_OUT "out"
#define EXT_MED "med"
#define EXT_AVG "avg"
#define EXT_MUON_DATA "mu"
#define EXT_VOXEL_DATA "vox"
#define EXT_V "vvals"
#define EXT_C "cvals"
#define EXT_CONVERGE "conv"
#define EXT_LAMBDAS "lam"
#define EXT_SAMP_VOX "smp"
//define standard constants
#define MAX_ITER 2 //for debugging purposes, not computational
#define RAD_CONVERT 1000 //standard (=1) is radians
#define LENGTH_CONVERT 10 //standard (=1) is millimeters
#define X 0
#define Y 1
#define Z 2
//#define PRINT_ITERATION i>0
//#define PRINT_ITERATION ((i+1)%5)==0
#define PRINT_ITERATION ((i+1)%10)==0
int em(double*, double*, struct muon*, double**, FILE**);
void compute_c(struct muon*, double*, double*, int, double**, double**, double**);
void compute_v(struct muon*, double*, double*, double*, int, double**);
double calc_voxel_weight(struct muon*, unsigned int, double, double**);