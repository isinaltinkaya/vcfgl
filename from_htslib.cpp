#include "from_htslib.h"

/**********************
 ***     Memory     ***
 **********************/

/* For use with hts_expand macros *only* */
// HTSLIB_EXPORT
size_t hts_realloc_or_die(size_t n, size_t m, size_t m_sz, size_t size,
                          int clear, void **ptr, const char *func) {
    /* If new_m and size are both below this limit, multiplying them
       together can't overflow */
    const size_t safe = (size_t) 1 << (sizeof(size_t) * 4);
    void *new_ptr;
    size_t bytes, new_m;

    new_m = n;
    kroundup_size_t(new_m);

    bytes = size * new_m;

    /* Check for overflow.  Both ensure that new_m will fit in m (we make the
       pessimistic assumption that m is signed), and that bytes has not
       wrapped around. */
    if (new_m > (((size_t) 1 << (m_sz * 8 - 1)) - 1)
        || ((size > safe || new_m > safe)
            && bytes / new_m != size)) {
        errno = ENOMEM;
        goto die;
    }

    new_ptr = realloc(*ptr, bytes);
    if (new_ptr == NULL) goto die;

    if (clear) {
        if (new_m > m) {
            memset((char *) new_ptr + m * size, 0, (new_m - m) * size);
        }
    }

    *ptr = new_ptr;

    return new_m;

 die:
    hts_log_error("%s", strerror(errno));
    exit(1);
}
