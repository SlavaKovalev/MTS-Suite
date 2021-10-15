#include "mts.h"
#include "mtsio.h"
//constants defining the error messages to display for the command line
#define ERR_MSG_TITLE "Muon Tomography Suite Command Line Options"
#define ERR_MSG_OPTION "Illegal Command Line Option"
#define ERR_MSG_COMBO "The following commands may not be selected together"
#define ERR_MSG_COVERAGE1 "run coverage analysis after reconstruction"
#define ERR_MSG_COVERAGE2 "(input file must be provided unless poca/em is being run or
stdin option is chosen)"
#define ERR_MSG_DETAILS1 "provides detailed information on inner processes"
#define ERR_MSG_DETAILS2 0
#define ERR_MSG_EM1 "run maximum likelihood algorithm"
#define ERR_MSG_EM2 "(input file must be provided unless stdin option is selected)"
#define ERR_MSG_HELP1 "list command line options"
#define ERR_MSG_HELP2 0
#define ERR_MSG_NO_FIT1 "run reconstrunction without line fitting the data"
#define ERR_MSG_NO_FIT2 0
#define ERR_MSG_NORM1 "run reconstrunction with normalized vectors"
#define ERR_MSG_NORM2 0
#define ERR_MSG_POCA1 "run poca algorithm"
#define ERR_MSG_POCA2 "(input file must be provided unless stdin option is selected)"
#define ERR_MSG_ROOT1 "run root analysis after reconstruction (cannot be selected with
stdout option)"
#define ERR_MSG_ROOT2 0
#define ERR_MSG_STDIN1 "accept input from stdin instead of file"
#define ERR_MSG_STDIN2 "(all input must come from stdin; provided input files are ignored)"
#define ERR_MSG_STDOUT1 "print output to stdout instead of file"
#define ERR_MSG_STDOUT2 0
#define ERR_MSG_VALIDATE1 "run validation analysis after reconstruction"
#define ERR_MSG_VALIDATE2 "(input file must be provided unless poca/em is being run or stdin
option is chosen)"
#define ERR_MSG_X1 "extended details, prints algorithmic computations step by step"
#define ERR_MSG_X2 0
//constants defining the error messages for internal MTS runtime errors
#define ERR_MSG_EMPTY "is empty!"
#define ERR_MSG_FILE "is not a valid file!"
#define ERR_MSG_FORMAT1 "is incorrectly formatted!"
#define ERR_MSG_FORMAT2 "Invalid Line"
#define ERR_MSG_MEM "insufficient memory available"
#define ERR_MSG_OPTIONAL "An input file must be provided for the following options"
#define ERR_MSG_UNIX "unix system command failed:"
#define ERR_MSG_UNUSED "unused: Input from stdin"
char commandError(char option, char option2) {
if (option!=option2) printf("\n%s: -%c and -%c\n", ERR_MSG_COMBO, option, option2);
else if (option!=SHORT_HELP) printf("\n%s: |%c|\n", ERR_MSG_OPTION, option);
printf("\n%s:\n\n", ERR_MSG_TITLE);
printf("\t-%c or --%s: \n\n\t\t%s \n\n", SHORT_DETAILS, LONG_DETAILS,
ERR_MSG_DETAILS1, ERR_MSG_DETAILS2);
printf("\t-%c or --%s: \n\n\t\t%s \n\t\t%s \n\n", SHORT_EM, LONG_EM,
ERR_MSG_EM1, ERR_MSG_EM2);
printf("\t-%c or --%s: \n\n\t\t%s \n\n", SHORT_HELP, LONG_HELP,
ERR_MSG_HELP1, ERR_MSG_HELP2);
return option;
}
int emptyError(char* file) {
fprintf(ERR_OUT, "\n\n\t%s %s\n\n", file, ERR_MSG_EMPTY);
return ERR_CODE_EMPTY_FILE;
}
int fileError(char* file) {
fprintf(ERR_OUT, "\n\n\t%s %s!\n\n", file, ERR_MSG_FILE);
return ERR_CODE_INVALID_FILE;
}
int formatError(char* file, char* input, int line) {
fprintf(ERR_OUT, "\n\n\t%s %s\n\t\t%s %d: %s\n\n", file, ERR_MSG_FORMAT1,
ERR_MSG_FORMAT2, line, input);
return ERR_CODE_INVALID_INPUT;
}
int memError() {
fprintf(ERR_OUT, "\n\n\t%s\n\n", ERR_MSG_MEM);
return ERR_CODE_MEM;
}
int optionalError(char* longOp, char shortOp) {
fprintf(ERR_OUT, "\n\n\t%s: -%c and --%s\n\n", ERR_MSG_OPTIONAL, shortOp, longOp);
return ERR_CODE_OPTIONAL_REQUIRED;
}
int unixError(char* command) {
fprintf(ERR_OUT, "\n\n\t%s %s\n\n", ERR_MSG_UNIX, command);
return ERR_CODE_UNIX;
}
int unusedError(char* file) {
fprintf(ERR_OUT, "\n\t%s %s\n\n", file, ERR_MSG_UNUSED);
return ERR_CODE_UNUSED_FILE;
}