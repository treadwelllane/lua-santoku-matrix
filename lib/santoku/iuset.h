#ifndef TK_IUSET_H
#define TK_IUSET_H

#include <santoku/klib.h>

KHASH_INIT(tk_iuset, int64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal)
typedef khash_t(tk_iuset) tk_iuset_t;

#define tk_iuset_put(...) kh_put(tk_iuset, __VA_ARGS__)
#define tk_iuset_get(...) kh_get(tk_iuset, __VA_ARGS__)
#define tk_iuset_del(...) kh_del(tk_iuset, __VA_ARGS__)
#define tk_iuset_exist(...) kh_exist(__VA_ARGS__)
#define tk_iuset_key(...) kh_key(__VA_ARGS__)
#define tk_iuset_value(...) kh_value(__VA_ARGS__)
#define tk_iuset_begin(...) kh_begin(__VA_ARGS__)
#define tk_iuset_end(...) kh_end(__VA_ARGS__)
#define tk_iuset_size(...) kh_size(__VA_ARGS__)
#define tk_iuset_resize(...) kh_resize(tk_iuset, __VA_ARGS__)
#define tk_iuset_clear(...) kh_clear(tk_iuset, __VA_ARGS__)
#define tk_iuset_destroy(...) kh_destroy(tk_iuset, __VA_ARGS__)
#define tk_iuset_create() kh_init(tk_iuset)
#define tk_iuset_contains(h, v) (tk_iuset_get(h, v) != tk_iuset_end(h))
// Note: this changes the default behavior so that a value var isn't required
#define tk_iuset_foreach(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

static inline void tk_iuset_union (tk_iuset_t *a, tk_iuset_t *b)
{
  int kha;
  int64_t x;
  tk_iuset_foreach(b, x, ({
    tk_iuset_put(a, x, &kha);
  }));
}

static inline double tk_iuset_jaccard (tk_iuset_t *a, tk_iuset_t *b)
{
  uint64_t intersection = 0;
  uint64_t union_count = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (tk_iuset_contains(b, x))
      intersection ++;
  }));
  union_count = tk_iuset_size(a) + tk_iuset_size(b) - intersection;
  if (union_count == 0)
    return 0.0;
  return (double) intersection / (double) union_count;
}

static inline void tk_iuset_intersect (tk_iuset_t *a, tk_iuset_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (!tk_iuset_contains(b, x))
      kh_del(tk_iuset, a, i);
  }))
}

static inline void tk_iuset_difference (tk_iuset_t *a, tk_iuset_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (tk_iuset_contains(b, x))
      kh_del(tk_iuset, a, i);
  }))
}

// Operations with iumap (using map keys as the secondary set)
#include <santoku/iumap.h>

static inline void tk_iuset_union_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  int kha;
  for (khint_t i = tk_iumap_begin(b); i < tk_iumap_end(b); i++) {
    if (tk_iumap_exist(b, i)) {
      int64_t key = tk_iumap_key(b, i);
      tk_iuset_put(a, key, &kha);
    }
  }
}

static inline void tk_iuset_intersect_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    khint_t khi = tk_iumap_get(b, x);
    if (khi == tk_iumap_end(b))
      kh_del(tk_iuset, a, i);
  }))
}

static inline void tk_iuset_difference_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    khint_t khi = tk_iumap_get(b, x);
    if (khi != tk_iumap_end(b))
      kh_del(tk_iuset, a, i);
  }))
}

// TODO
// #define tk_iuset_dup(a)
// #define tk_iuset_union(a, b)
// #define tk_iuset_intersect(a, b)
// #define tk_iuset_difference(a, b)
// #define tk_iuset_jaccard(a, b)

#endif
