#define TRACK_PUSH 0.00000001
#define SMALL_NUM 0.1
#define INPUT_FILE "geant_input_file"
int preprocessing(double*, double*, struct muon*, double**, FILE**);
int header(double*, double*, double**, FILE**);
int em_data(struct Line*, struct Line*, struct muon*, double**, FILE**);
double pocaLtoL(struct Line*, struct Line*, struct Point*, double**);
void travel(struct Point*, struct Point*, int, double);
double in_volume(struct Point*, double**);
struct voxel* track(struct Point*, struct Point*, struct Point*, struct muon*, struct voxel*, double*,
double**, FILE**);