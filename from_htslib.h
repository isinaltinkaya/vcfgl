#include <cstddef>
#include <cerrno>
#include <cstring>
#include <stddef.h>
#include <stdint.h>
#include <inttypes.h>

#include <htslib/hts_defs.h>
#include <htslib/hts_log.h>
#include <htslib/kstring.h>
#include <htslib/kroundup.h>

// source: htslib/hts.c
size_t hts_realloc_or_die(size_t n, size_t m, size_t m_sz, size_t size,int clear, void **ptr, const char *func);

