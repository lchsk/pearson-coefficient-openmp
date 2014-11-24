//------------------------------------------
// definitions.h
// Lists constants used across other files
//------------------------------------------

// how many digits should be printed as an output
#define FLOAT_PRECISION 10

// length of the output
#define RESULT_LENGTH 13

// rank of the master processor
#define ROOT_PROC 0

// lenght of input values
//#define INPUT_LENGTH 1000000

// tags (needed when using send and recv)
#define TAG_SET_A 0
#define TAG_SET_B 1

// whether to use Bcast (1) or send/recv (0)
#define USE_BROADCAST 1

// if set to 1, calculating mean and stddev (in parallel mode)
// will involve setting indices in different places of the array
// (which should be faster),
// otherwise (0) will use scatter/reduce
#define USE_PARALLEL_INDICES 1

#define THREADS 10