//define constants for xyz componets
#define ALL_COMPONENTS -1
#define X_COMPONENT 0
#define Y_COMPONENT 1
#define Z_COMPONENT 2
//define constants for fit function
#define FIT_X 1
#define FIT_Y 0
#define FIT_NONE -1
#define RAD 180/(4.0*atan(1.0))
#define MAX_DETECTORS 5
#define MAX_DIMENSION 2
double vec_dot (struct Point*, struct Point*);
double vec_norm (struct Point*);
void vec_mid (struct Point*, struct Point*, struct Point*);
void vec_sub (struct Point*, struct Point*, struct Point*);
void vec_add (struct Point*, struct Point*, struct Point*);
void vec_mult (double, struct Point*, struct Point*);
void vec_div (double, struct Point*, struct Point*);
void vec_fit (struct Point**, struct Line*, int, int);
double vec_rad_to_deg (double);
double vec_angle (struct Line*, struct Line*, int);
void vec_copy(struct Line*, struct Line*);
