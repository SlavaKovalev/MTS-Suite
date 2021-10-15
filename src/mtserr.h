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
char commandError(char, char);
int emptyError(char*);
int fileError(char*);
int formatError(char*, char*, int);
int memError();
int optionalError(char*, char);
int unixError(char*);
int unusedError(char*);