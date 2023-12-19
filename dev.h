#ifndef __DEV__
#define __DEV__

#define DEV 0

#if 1 == DEV

// ----------------------------------------------------------------------------->
// DEV MODE: ON
// ----------------------------------------------------------------------------->

/* -> DEVELOPMENT/DEBUGGING MACROS
 * --------------------------------------------*/
 // enabled if `make DEV=1` is used


 /*
  * Macro:[DEVASSERT]
  * Works the same way as ASSERT in shared.h, but only run if dev mode is on
  */
#define DEVASSERT_EXPAND(x) x
#define DEVASSERT(expr)                                                        \
do {                                                                    \
	if (!(DEVASSERT_EXPAND(expr))){ \
		fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n*******\n", \
				__FILE__, __FUNCTION__, __LINE__, #expr);               \
		exit(1);                                                        \
	}                                                                   \
} while (0);



 /*
  * Macro:[DEVRUN]
  * @param expr - an expression to run
  * runs a given expression if DEV==1
  */
#define DEVRUN_EXPAND(x) x
#define DEVRUN(expr)                                                            \
    do {                                                                       \
DEVRUN_EXPAND(expr);       \
    } while (0);

 /*
  * Macro:[DEVRUNV]
  * @param expr - an expression to run
  * same as macro:DEVRUN but verbose (prints what it is running)
  */
#define DEVRUNV(expr)                                                            \
    do {                                                                       \
        fprintf(stderr, "\n\n\n*****************************************\n");  \
        fprintf(stderr, "*** [DEV]<%s:%d>\n\n-> Running:\t%s\n\n-> Output:\t", \
                __FILE__, __LINE__, #expr);                                     \
DEVRUN_EXPAND(expr);       \
        fprintf(stderr, "\n*****************************************\n\n\n");  \
    } while (0);

  /*
   * Macro:[DEVRUNVX]
   * @param expr - an expression to run
   * same as DEVRUNV, but exits after running
   */
#define DEVRUNVX(expr)                                                            \
    do {                                                                       \
        fprintf(stderr, "\n\n\n*****************************************\n");  \
        fprintf(stderr, "*** [DEV]<%s:%d>\n\n-> Running:\t%s\n\n-> Output:\t", \
                __FILE__, __LINE__, #expr);                                     \
DEVRUN_EXPAND(expr);       \
        fprintf(stderr, "\n*****************************************\n\n\n");  \
        exit(1);                                                               \
    } while (0);

   /*
    *
    * Macro:[DEVPRINT]
    * @param msg - a message to print
    * @details uses variable arguments to print a message to stderr with additional
    * info
    * @usage DEVPRINT("Hello %s", "World");
    * @note VAARGS may not work on all compilers, so only defined in dev mode
    */
#define DEVPRINT(msg, ...)                                                    \
    do {                                                                      \
        fprintf(stderr, "\n\n\n*****************************************\n"); \
        fprintf(stderr, "*** [DEV]<%s:%d>\n\n-> Message:\t", __FILE__,        \
                __LINE__);                                                    \
        fprintf(stderr, msg, ##__VA_ARGS__);                                  \
        fprintf(stderr, "\n*****************************************\n\n\n"); \
    } while (0);

    /*
     *
     * Macro:[DEVPRINTX]
     * @param msg - a message to print
     * @details same as macro:DEVPRINT, but exits after printing
     */
#define DEVPRINTX(msg, ...)                                                   \
    do {                                                                      \
        fprintf(stderr, "\n\n\n*****************************************\n"); \
        fprintf(stderr, "*** [DEV]<%s:%d>\n\n-> Message:\t", __FILE__,        \
                __LINE__);                                                    \
        fprintf(stderr, msg, ##__VA_ARGS__);                                  \
        fprintf(stderr, "\n*****************************************\n\n\n"); \
        exit(1);                                                              \
    } while (0);

#else

// ----------------------------------------------------------------------------->
// DEV MODE: OFF
// ----------------------------------------------------------------------------->

#define DEVASSERT(expr)

#define DEVRUN(expr) 

#define DEVRUNX(expr)

#define DEVPRINT(msg, ...)

#define DEVPRINTX(msg, ...)


#endif

/* --------------------------------------------------------------------------*/

#endif  // __DEV__
