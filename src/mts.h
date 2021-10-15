#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//define various program constants
#define DELIMS " \n\t"
#define ERR_PREC 0.000001
//define various maximum constants
#define MAX_COMMAND 100
#define MAX_FILENAME 200
#define MAX_LINE 1000
//define constants used for validation
#define CONTINUE 1
#define END 0
#define PARALLEL -1
//define file pointer constants
#define FP_IN 0
#define FP_OUT 1
#define FP_OP_OUT 2
#define FP_STD_OUT 3
#define FP_DIST_OUT 4
#define FP_POCA_OUT 5
#define FP_OUT_AVG 6
#define FP_OUT_MED 7
#define MAX_FILEPOINTERS 8
//define constants for error return values
#define ERR_CODE_COMMAND_LINE 1
#define ERR_CODE_EMPTY_FILE 2
#define ERR_CODE_INVALID_COMMAND 3
#define ERR_CODE_INVALID_FILE 4
#define ERR_CODE_INVALID_INPUT 5
#define ERR_CODE_MEM 6
#define ERR_CODE_NO_FILE 7
#define ERR_CODE_OPTIONAL_REQUIRED 8
#define ERR_CODE_UNIX 9
#define ERR_CODE_UNUSED_FILE 10
//define constant for where to print error messages
#define ERR_OUT stderr
//define array indexes for various command line option flags
#define PARAM_DETAILS 0
#define PARAM_DEPENDENT 1
#define PARAM_DIST 2
#define PARAM_EM_AVERAGE 3
#define PARAM_EM_MEDIAN 4
#define PARAM_EM_3D 5
#define PARAM_EM_WEIGHTED 6
#define PARAM_END_OF_PATH 7
#define PARAM_HIT_TARGET 8
#define PARAM_INDEPENDENT 9
#define PARAM_ITERATIONS 10
#define PARAM_INIT_LAMBDA 11
#define PARAM_PRECISE_L 12
#define PARAM_PRECISE_T 13
#define PARAM_STD 14
#define PARAM_X_LENGTH 15
#define PARAM_Y_LENGTH 16
#define PARAM_Z_LENGTH 17
#define PARAM_X_VOXEL_SIZE 18
#define PARAM_Y_VOXEL_SIZE 19
#define PARAM_Z_VOXEL_SIZE 20
#define PARAM_X_VOXEL_TOTAL 21
#define PARAM_Y_VOXEL_TOTAL 22
#define PARAM_Z_VOXEL_TOTAL 23
#define PARAM_ALL_VOXELS 24
#define PARAM_X_MIN 25
#define PARAM_Y_MIN 26
#define PARAM_Z_MIN 27
#define PARAM_X_MAX 28
#define PARAM_Y_MAX 29
#define PARAM_Z_MAX 30
#define PARAM_EVENTS 31
#define PARAM_INC_DECT 32
#define PARAM_OUT_DECT 33
#define PARAM_MOMENTUM 34
#define PARAM_LINE 35
#define PARAM_HITS 36
#define PARAM_PARALLEL 37
#define PARAM_SCAT_ANG 38
#define PARAM_L 39
#define PARAM_EM 40
#define PARAM_IN_VOLUME 41
#define PARAM_CUR_EVENT 42
#define PARAM_CONTINUE 43
#define PARAM_NOM_MOMENTUM 44
#define PARAM_MILLIRADIANS 45
#define PARAM_MIN_MOMENTUM 46
#define PARAM_CUTOFF_ANGLE 47
#define PARAM_OP_OUT 48
#define PARAM_POCA 49
#define PARAM_DOCA 50
#define PARAM_DTX 51
#define PARAM_DTY 52
#define PARAM_DX 53
#define PARAM_DY 54
#define PARAM_OUT 55
#define PARAM_C_PRINT 56
#define PARAM_EM_BINS 57
#define PARAM_EM_BIN_SIZE 58
#define PARAM_EM_ONLINE 59
#define PARAM_UNITS_LENGTH 60
#define PARAM_MOM_CUT 61
#define PARAM_MOM_HIGH_CUT 62
#define PARAM_MOM_LOW_CUT 63
#define PARAM_PREV_VOXEL 64
#define MAX_PARAMS 65
//define constants for command line options
#define LONG_BINS "bins"
#define LONG_BIN_SIZE "bin_size"
#define LONG_C_PRINT "c_print"
#define LONG_CUTOFF "cutoff"
#define LONG_DETAILS "details"
#define LONG_DIST "dist"
#define LONG_EM "em"
#define LONG_EM_AVERAGE "average"
#define LONG_EM_MEDIAN "median"
#define LONG_EM_3D "3D"
#define LONG_EM_WEIGHTED "weight"
#define LONG_HELP "help"
#define LONG_INPUT "input"
#define LONG_ITERATIONS "iterations"
#define LONG_LAMBDA "lambda"
#define LONG_MILLIRADIANS "milliradians"
#define LONG_MOM_HIGH_CUT "mom_cut_high"
#define LONG_MOM_LOW_CUT "mom_cut_low"
#define LONG_NOMINAL "nominal"
#define LONG_ONLINE "online"
#define LONG_OUTPUT "output"
#define LONG_POCA "poca"
#define LONG_PRECISE_L "precise_l"
#define LONG_PRECISE_T "precise_t"
#define LONG_STD "std"
#define LONG_UNITS_LENGTH "units"
#define LONG_X "x_voxel_size"
#define LONG_Y "y_voxel_size"
#define LONG_Z "z_voxel_size"
#define SHORT_BINS 'B'
#define SHORT_BIN_SIZE 'b'
#define SHORT_C_PRINT 'C'
#define SHORT_CUTOFF 'c'
#define SHORT_DETAILS 'd'
#define SHORT_DIST 'D'
#define SHORT_EM 'e'
#define SHORT_EM_WEIGHTED 'w'
#define SHORT_HELP 'h'
#define SHORT_INPUT 'i'
#define SHORT_ITERATIONS 'I'
#define SHORT_LAMBDA 'l'
#define SHORT_MILLIRADIANS 'm'
#define SHORT_MOM_HIGH_CUT 'K'
#define SHORT_MOM_LOW_CUT 'k'
#define SHORT_NOMINAL 'n'
#define SHORT_ONLINE 'O'
#define SHORT_OUTPUT 'o'
#define SHORT_POCA 'p'
#define SHORT_PRECISE_L 'L'
#define SHORT_PRECISE_T 'T'
#define SHORT_STD 's'
#define SHORT_UNITS_LENGTH 'u'
#define SHORT_X 'x'
#define SHORT_Y 'y'
#define SHORT_Z 'z'
#define SHORT_OPTIONS "B:b:c:dD::e::whi:I:k:K:l:Lmn:O:o::p::s::Tu:x:y:z:"
#define SAMPLE_VOXELS {24603, 23384, 23385, 23386, 23413, 23414, 23415, 24584, 24585,
24613, 24614, 24615, 24626, -1}
//#define SAMPLE_VOXELS {14374, 17994, 20638, 22735, -1}
#define NUMOF_SAMP_VOX 13
//cosmetic constants
#define BANNER "======================"
struct Point {
double x, y, z;
};