#ifndef __DEV__
#define __DEV__

#if 1 == DEV

// ----------------------------------------------------------------------------->
// DEV MODE: ON
// ----------------------------------------------------------------------------->

/* -> DEVELOPMENT/DEBUGGING MACROS
 * --------------------------------------------*/
// enabled if `make dev` is used

/*
 * Macro:[DEVRUN]
 * @param run - a function to run
 * prints the output of a function to stderr with additional info
 */
#define DEVRUN(run)                                                            \
    do {                                                                       \
        fprintf(stderr, "\n\n\n*****************************************\n");  \
        fprintf(stderr, "*** [DEV]<%s:%d>\n\n-> Running:\t%s\n\n-> Output:\t", \
                __FILE__, __LINE__, #run);                                     \
        run;                                                                   \
        fprintf(stderr, "\n*****************************************\n\n\n");  \
    } while (0);

/*
 * Macro:[DEVRUNX]
 * @param run - a function to run
 * just like DEVRUN, but exits after running
 */
#define DEVRUNX(run)                                                           \
    do {                                                                       \
        fprintf(stderr, "\n\n\n*****************************************\n");  \
        fprintf(stderr, "*** [DEV]<%s:%d>\n\n-> Running:\t%s\n\n-> Output:\t", \
                __FILE__, __LINE__, #run);                                     \
        run;                                                                   \
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
 * @details same as DEVPRINT, but exits after printing
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

#define DEVRUN(run) \
    do {            \
    } while (0)

#define DEVRUNX(run) \
    do {             \
    } while (0)

#define DEVPRINT(msg, ...) \
    do {                   \
    } while (0)

#define DEVPRINTX(msg, ...) \
    do {                    \
    } while (0)

#endif

/* --------------------------------------------------------------------------*/

#endif  // __DEV__
