#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#ifndef tk_vec_lt
#define tk_vec_lt(a, b) ((a) < (b))
#endif

#ifndef tk_vec_gt
#define tk_vec_gt(a, b) ((a) > (b))
#endif

#ifndef tk_vec_mt
#define tk_vec_mt tk_pp_xstr(tk_vec_pfx(t))
#endif

#define tk_vec_pfx(name) tk_pp_strcat(tk_vec_name, name)
#define tk_vec_err(L, name, n, ...) \
  tk_lua_verror((L), ((n) + 1), tk_pp_xstr(tk_vec_pfx(name)), __VA_ARGS__)

typedef kvec_t(tk_vec_base) tk_vec_pfx(t);
#define tk_vec_ksort(...) KSORT_INIT(__VA_ARGS__)
tk_vec_ksort(tk_vec_pfx(asc), tk_vec_base, tk_vec_lt)
tk_vec_ksort(tk_vec_pfx(desc), tk_vec_base, tk_vec_gt)

#define tk_vec_shuffle(...) ks_shuffle(__VA_ARGS__)
#define tk_vec_introsort(...) ks_introsort(__VA_ARGS__)
#define tk_vec_ksmall(...) ks_ksmall(__VA_ARGS__)

static inline tk_vec_pfx(t) tk_vec_name (tk_vec_base *a, uint64_t n)
{
  return (tk_vec_pfx(t)) { .a = a, .n = n, .m = n };
}

static inline tk_vec_pfx(t) *tk_vec_pfx(peek) (lua_State *L, int i, const char *name)
{
  return (tk_vec_pfx(t) *) tk_lua_checkuserdata(L, i, tk_vec_mt, name);
}

static inline tk_vec_pfx(t) *tk_vec_pfx(peekopt) (lua_State *L, int i)
{
  return (tk_vec_pfx(t) *) tk_lua_testuserdata(L, i, tk_vec_mt);
}

static inline void tk_vec_pfx(destroy) (tk_vec_pfx(t) *r)
{
#ifdef tk_vec_destroy_item
  for (uint64_t i = 0; i < r->n; i ++)
    tk_vec_destroy_item(r->a[i]);
#endif
  kv_destroy(*r);
  memset(r, 0, sizeof(tk_vec_pfx(t)));
}

static inline int tk_vec_pfx(destroy_lua) (lua_State *L)
{
  tk_vec_pfx(destroy)(tk_vec_pfx(peek)(L, 1, "vector"));
  return 0;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(create) (lua_State *, size_t, tk_vec_base *, tk_vec_pfx(t) *);

static inline void tk_vec_pfx(resize) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  size_t m,
  bool setn
) {
  tk_vec_pfx(t) v = *m0;
  if (m == 0) {
    free(v.a);
    v.a = NULL;
    v.m = v.n = 0;
    *m0 = v;
  } else {
    kv_resize(tk_vec_base, v, m);
    if (m0->a == NULL)
      tk_vec_err(L, resize, 1, "resize failed");
    v.m = m;
    if (setn)
      v.n = m;
    *m0 = v;
  }
}

static inline void tk_vec_pfx(setn) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  size_t n
) {
  if (n > m0->m)
    tk_vec_pfx(resize)(L, m0, n, true);
  else
    m0->n = n;
}

static inline void tk_vec_pfx(ensure) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  size_t m
) {
  if (m0->m >= m)
    return;
  tk_vec_pfx(resize)(L, m0, m, false);
}

static inline void tk_vec_pfx(copy) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  int64_t start,
  int64_t end,
  int64_t dest
) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_vec_pfx(ensure)(L, m0, m);
  if (m0 == m1)
    memmove(m0->a + dest, m1->a + start, sizeof(tk_vec_base) * (uint64_t) (end - start));
  else
    memcpy(m0->a + dest, m1->a + start, sizeof(tk_vec_base) * (uint64_t) (end - start));
  if (m0->n < m)
    m0->n = m;
}

static inline void tk_vec_pfx(transpose) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  uint64_t cols
) {
  tk_vec_pfx(ensure)(L, 0, m1->n);
  uint64_t rows = m1->n / cols;
  for (uint64_t r = 0; r < rows; r ++)
    for (uint64_t c = 0; c < cols; c ++)
      m0->a[c * rows + r] = m1->a[r * cols + c];
}

static inline uint64_t tk_vec_pfx(capacity) (tk_vec_pfx(t) *m0)
{
  return m0->m;
}

static inline uint64_t tk_vec_pfx(size) (tk_vec_pfx(t) *m0)
{
  return m0->n;
}

static inline void tk_vec_pfx(shrink) (
  lua_State *L,
  tk_vec_pfx(t) *v
) {
  tk_vec_pfx(resize)(L, v, v->n, true);
}

static inline void tk_vec_pfx(clear) (
  tk_vec_pfx(t) *v
) {
  v->n = 0;
}

static inline void tk_vec_pfx(zero) (
  tk_vec_pfx(t) *v
) {
  memset(v->a, 0, v->m * sizeof(tk_vec_base));
}

static inline tk_vec_base tk_vec_pfx(get) (
  lua_State *L,
  tk_vec_pfx(t) *v,
  uint64_t i
) {
  if (i >= v->n)
    tk_vec_err(L, get, 1, "index out of bounds");
  return v->a[i];
}

static inline void tk_vec_pfx(push) (
  tk_vec_pfx(t) *v,
  tk_vec_base x
) {
  tk_vec_pfx(t) v0 = *v;
  kv_push(tk_vec_base, v0, x);
  *v = v0;
}

static inline void tk_vec_pfx(set) (
  lua_State *L,
  tk_vec_pfx(t) *v,
  uint64_t i,
  tk_vec_base x
) {
  tk_vec_pfx(ensure)(L, v, i + 1);
  v->a[i] = x;
  if (i + 1 > v->n)
    v->n = i + 1;
}

#ifdef tk_vec_pushbase

static inline void tk_vec_pfx(table) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t start,
  uint64_t end
) {
  lua_newtable(L); // vals
  for (uint64_t i = start; i < end; i ++) {
    lua_pushinteger(L, (int64_t) i + 1); // vals idx
    tk_vec_pushbase(L, m0->a[i]); // vals idx v
    lua_settable(L, -3); // vals
  }
}

static inline void tk_vec_pfx(ctable) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols,
  uint64_t start,
  uint64_t end
) {
  bool done = false;
  lua_newtable(L); // cols
  for (uint64_t c = start; !done && c < end; c ++) {
  // for (uint64_t c = 0; !done && c < (m0->n + cols) / cols - 1; r ++) {
    lua_pushinteger(L, (int64_t) c + 1); // cols c
    lua_newtable(L); // cols c rows
    for (uint64_t r = 0; r < m0->n / cols; r ++) {
      if (r * cols + c >= m0->n)
        continue;
      lua_pushinteger(L, (int64_t) r + 1); // cols c rows r
      tk_vec_pushbase(L, m0->a[r * cols + c]); // ols c rows r v
      lua_settable(L, -3); // cols c rows
    }
    lua_settable(L, -3); // cols
  }
}

static inline void tk_vec_pfx(rtable) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols,
  uint64_t start,
  uint64_t end
) {
  bool done = false;
  lua_newtable(L); // rows
  for (uint64_t r = start; !done && r < end; r ++) {
    lua_pushinteger(L, (int64_t) r + 1); // rows r
    lua_newtable(L); // rows r cols
    for (uint64_t c = 0; c < cols; c ++) {
      if (r * cols + c >= m0->n) {
        done = true;
        break;
      }
      lua_pushinteger(L, (int64_t) c + 1); // rows r cols c
      tk_vec_pushbase(L, m0->a[r * cols + c]); // rows r cols c v
      lua_settable(L, -3); // rows r cols
    }
    lua_settable(L, -3); // rows
  }
}

#endif

static inline void tk_vec_pfx(shuffle) (tk_vec_pfx(t) *v) {
  tk_vec_shuffle(tk_vec_pfx(asc), v->n, v->a);
}

static inline void tk_vec_pfx(asc) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  tk_vec_introsort(tk_vec_pfx(asc), e - s, v->a + s);
}

static inline void tk_vec_pfx(desc) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  tk_vec_introsort(tk_vec_pfx(desc), e - s, v->a + s);
}

static inline void tk_vec_pfx(kasc) (tk_vec_pfx(t) *v, size_t k, uint64_t s, uint64_t e) {
  tk_vec_ksmall(tk_vec_pfx(asc), e - s, v->a + s, k);
}

static inline void tk_vec_pfx(kdesc) (tk_vec_pfx(t) *v, size_t k, uint64_t s, uint64_t e) {
  tk_vec_ksmall(tk_vec_pfx(desc), e - s, v->a + s, k);
}

#ifndef tk_vec_limited

static inline double tk_vec_pfx(magnitude) (
  tk_vec_pfx(t) *m0
) {
  tk_vec_base sum = 0.0;
  for (size_t i = 0; i < m0->n; i ++) {
    tk_vec_base val = m0->a[i];
    sum += val * val;
  }
  return sqrt(sum);
}

static inline tk_vec_base tk_vec_pfx(sum) (
  tk_vec_pfx(t) *m0
) {
  tk_vec_base sum = 0.0;
  for (size_t i = 0; i < m0->n; i ++) {
    tk_vec_base val = m0->a[i];
    sum += val;
  }
  return sum;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(csums) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    tk_vec_base sum = 0;
    for (uint64_t r = 0; r < m0->n / cols; r ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[c] = sum;
  }
  return out;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(rsums) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
  }
  return out;
}

static inline tk_vec_base tk_vec_pfx(max) (
  tk_vec_pfx(t) *m0,
  uint64_t *colp
) {
  if (!m0->n)
    return 0;
  uint64_t maxcol = 0;
  tk_vec_base maxval = m0->a[0];
  for (size_t i = 1; i < m0->n; i ++) {
    if (m0->a[i] > maxval) {
      maxcol = i;
      maxval = m0->a[i];
    }
  }
  *colp = maxcol;
  return maxval;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(cmaxs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    tk_vec_base maxv = m0->a[0 * cols + c];
    for (size_t r = 1; r < m0->n / cols; r ++) {
      if (m0->a[r * cols + c] > maxv) {
        maxv = m0->a[r * cols + c];
      }
    }
    out->a[c] = maxv;
  }
  return out;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(rmaxs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
    tk_vec_base maxv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      if (m0->a[r * cols + c] > maxv) {
        maxv = m0->a[r * cols + c];
      }
    }
    out->a[r] = maxv;
  }
  return out;
}

static inline tk_vec_base tk_vec_pfx(min) (
  tk_vec_pfx(t) *m0,
  uint64_t *colp
) {
  if (!m0->n)
    return 0;
  uint64_t mincol = 0;
  tk_vec_base minval = m0->a[0];
  for (size_t i = 1; i < m0->n; i ++) {
    if (m0->a[i] < minval) {
      mincol = i;
      minval = m0->a[i];
    }
  }
  *colp = mincol;
  return minval;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(cmins) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    tk_vec_base minv = m0->a[0 * cols + c];
    for (size_t r = 1; r < m0->n / cols; r ++) {
      if (m0->a[r * cols + c] < minv) {
        minv = m0->a[r * cols + c];
      }
    }
    out->a[c] = minv;
  }
  return out;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(rmins) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
    tk_vec_base minv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      if (m0->a[r * cols + c] < minv) {
        minv = m0->a[r * cols + c];
      }
    }
    out->a[r] = minv;
  }
  return out;
}

static inline void tk_vec_pfx(scale) (
  tk_vec_pfx(t) *m0,
  tk_vec_base scale,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  for (size_t i = start; i < end; i ++)
    m0->a[i] *= scale;
}

static inline void tk_vec_pfx(add) (
  tk_vec_pfx(t) *m0,
  tk_vec_base add,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  for (size_t i = start; i < end; i ++)
    m0->a[i] += add;
}

#ifdef tk_vec_abs
static inline void tk_vec_pfx(abs) (
  tk_vec_pfx(t) *m0,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  for (size_t i = start; i < end; i ++)
    m0->a[i] = tk_vec_abs(m0->a[i]);
}
#endif

static inline double tk_vec_pfx(dot) (
  lua_State *L,
  tk_vec_pfx(t) *a,
  tk_vec_pfx(t) *b
) {
  if (a->n != b->n)
    tk_vec_err(L, dot, 1, "dimension mismatch");
  size_t n = a->n;
  double sum = 0.0;
  for (size_t i = 0; i < n; i ++)
    sum += a->a[i] * b->a[i];
  return sum;
}

static inline void tk_vec_pfx(fill) (tk_vec_pfx(t) *v, tk_vec_base x)
{
  for (uint64_t i = 0; i < v->n; i ++)
    v->a[i] = x;
}

static inline void tk_vec_pfx(multiply) (
  lua_State *L,
  tk_vec_pfx(t) *a,
  tk_vec_pfx(t) *b,
  tk_vec_pfx(t) *c,
  uint64_t k,
  bool transpose_a,
  bool transpose_b
) {
  size_t m = transpose_a ? k : a->n / k;
  size_t n = transpose_b ? k : b->n / k;
  tk_vec_pfx(ensure)(L, a, transpose_a ? k * m : m * k);
  tk_vec_pfx(ensure)(L, b, transpose_b ? k * n : k * n);
  tk_vec_pfx(ensure)(L, c, m * n);
  for (size_t i = 0; i < m; i ++) {
    for (size_t j = 0; j < n; j ++) {
      double sum = 0.0;
      for (size_t l = 0; l < k; l ++) {
        size_t a_idx = transpose_a ? l * m + i : i * k + l;
        size_t b_idx = transpose_b ? j * k + l : l * n + j;
        sum += a->a[a_idx] * b->a[b_idx];
      }
      c->a[i * n + j] = sum;
    }
  }
}

static inline void tk_vec_pfx(exp) (
  tk_vec_pfx(t) *m0,
  double exp,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start > end || start >= m0->n)
    return;
  if (end > m0->n - 1)
    end = m0->n - 1;
  for (size_t i = start; i <= end; i ++)
    m0->a[i] = pow(m0->a[i], exp);
}

static inline const char *tk_vec_pfx(raw) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  const char *fmt
) {
  if (fmt == NULL) {
    lua_pushlstring(L, (char *) m0->a, sizeof(tk_vec_base) * m0->n);
    return lua_tostring(L, -1);
  }
  void *out = NULL;
  size_t out_size = 0;
  if (!strcmp(fmt, "u32")) {
    out_size = sizeof(uint32_t);
    out = malloc(out_size * m0->n);
    if (!out) goto err_mem;
    for (size_t i = 0; i < m0->n; i ++) {
      tk_vec_base r = m0->a[i];
      if (r < 0)
        goto err_negative;
      if (r > (tk_vec_base) UINT32_MAX)
        goto err_toobig;
      ((uint32_t *) out)[i] = (uint32_t) r;
    }
  } else if (!strcmp(fmt, "u64")) {
    out_size = sizeof(uint64_t);
    out = malloc(out_size * m0->n);
    if (!out) goto err_mem;
    for (size_t i = 0; i < m0->n; i ++) {
      tk_vec_base r = m0->a[i];
      if (r < 0)
        goto err_negative;
      if (r > (tk_vec_base) UINT64_MAX)
        goto err_toobig;
      ((uint64_t *) out)[i] = (uint64_t) r;
    }
  } else if (!strcmp(fmt, "f32")) {
    out_size = sizeof(float);
    out = malloc(out_size * m0->n);
    if (!out)
      goto err_mem;
    for (size_t i = 0; i < m0->n; i ++)
      ((float *) out)[i] = (float) m0->a[i];
  } else if (!strcmp(fmt, "f64")) {
    out_size = sizeof(double);
    out = malloc(out_size * m0->n);
    if (!out)
      goto err_mem;
    for (size_t i = 0; i < m0->n; i ++)
      ((double *) out)[i] = (double) m0->a[i];
  } else {
    goto err_format;
  }
  lua_pushlstring(L, (char *) out, out_size * m0->n);
  free(out);
  return lua_tostring(L, -1);
err_format:
  tk_vec_err(L, raw, 2, "bad format", fmt);
  return NULL;
err_mem:
  tk_vec_err(L, raw, 2, "malloc error");
  return NULL;
err_negative:
  tk_vec_err(L, raw, 2, "unexpected negative");
  goto cleanup;
err_toobig:
  tk_vec_err(L, raw, 2, "number too large");
cleanup:
  free(out);
  return NULL;
}

static inline void tk_vec_pfx(fill_indices) (tk_vec_pfx(t) *v)
{
  for (uint64_t i = 0; i < v->n; i ++)
    v->a[i] = (tk_vec_base) i;
}

#endif

static inline int tk_vec_pfx(create_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peekopt)(L, 1);
  uint64_t m = tk_lua_optunsigned(L, 2, "size", 0);
  if (m0 != NULL) {
    tk_vec_pfx(t) *m1 = tk_vec_pfx(create)(L, m > 0 ? m : m0->n, NULL, NULL);
    tk_vec_pfx(copy)(L, m1, m0, 0, (int64_t) m0->n, 0);
    return 1;
#ifndef tk_vec_limited
  } else if (lua_type(L, 1) == LUA_TTABLE) {
    tk_vec_pfx(t) *m0 = tk_vec_pfx(create)(L, m, NULL, NULL);
    uint64_t l = lua_objlen(L, 1);
    tk_vec_pfx(ensure)(L, m0, l);
    m0->n = l;
    for (uint64_t i = 0; i < l; i ++) {
      lua_pushinteger(L, (int64_t) i + 1); // t i
      lua_gettable(L, 1); // t v
      m0->a[i] = tk_vec_peekbase(L, -1);
      lua_pop(L, 1); // t
    }
    return 1;
#endif
  } else if (lua_type(L, 1) == LUA_TNUMBER || lua_type(L, 1) == LUA_TNIL) {
    uint64_t m = tk_lua_optunsigned(L, 1, "size", 0);
    tk_vec_pfx(create)(L, m, NULL, NULL);
    return 1;
  } else {
    tk_vec_err(L, create, 1, "expected either nil, a table, or a number");
    return 0;
  }
}

#ifndef tk_vec_limited
static inline int tk_vec_pfx(dot_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "a");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "b");
  lua_pushnumber(L, tk_vec_pfx(dot)(L, m0, m1));
  return 1;
}

static inline int tk_vec_pfx(multiply_lua) (lua_State *L)
{
  lua_settop(L, 6);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "a");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "b");
  tk_vec_pfx(t) *m2 = tk_vec_pfx(peek)(L, 3, "c");
  uint64_t k = tk_lua_checkunsigned(L, 4, "middle_dimension");
  bool transpose_a = lua_toboolean(L, 5);
  bool transpose_b = lua_toboolean(L, 6);
  tk_vec_pfx(multiply)(L, m0, m1, m2, k, transpose_a, transpose_b);
  return 0;
}
#endif

static inline int tk_vec_pfx(copy_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "dest");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_vec_pfx(copy)(L, m0, m1, start, end, dest);
  return 0;
}

static inline int tk_vec_pfx(capacity_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  lua_pushinteger(L, (int64_t) tk_vec_pfx(capacity)(m0));
  return 1;
}

static inline int tk_vec_pfx(setn_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t n = tk_lua_checkunsigned(L, 2, "n");
  tk_vec_pfx(setn)(L, m0, n);
  return 0;
}

static inline int tk_vec_pfx(size_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  lua_pushinteger(L, (int64_t) tk_vec_pfx(size)(m0));
  return 1;
}

static inline int tk_vec_pfx(resize_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t m = tk_lua_checkunsigned(L, 2, "size");
  tk_vec_pfx(resize)(L, m0, m, true);
  return 0;
}

static inline int tk_vec_pfx(ensure_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t m = tk_lua_checkunsigned(L, 2, "size");
  tk_vec_pfx(ensure)(L, m0, m);
  return 0;
}

static inline int tk_vec_pfx(clear_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(clear)(m0);
  return 0;
}

static inline int tk_vec_pfx(zero_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(zero)(m0);
  return 0;
}

static inline int tk_vec_pfx(shrink_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(shrink)(L, m0);
  return 0;
}

static inline int tk_vec_pfx(transpose_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "dest");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "source");
  uint64_t cols = tk_lua_checkunsigned(L, 3, "cols");
  tk_vec_pfx(transpose)(L, m0, m1, cols);
  return 0;
}

static inline int tk_vec_pfx(shuffle_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(shuffle)(m0);
  return 1;
}

static inline int tk_vec_pfx(asc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  tk_vec_pfx(asc)(m0, start, end);
  return 1;
}

static inline int tk_vec_pfx(desc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  tk_vec_pfx(desc)(m0, start, end);
  return 1;
}

static inline int tk_vec_pfx(kasc_lua) (lua_State *L)
{
  lua_settop(L, 4);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t k = tk_lua_checkunsigned(L, 2, "k");
  uint64_t start = tk_lua_optunsigned(L, 3, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 4, "end", m0->n);
  tk_vec_pfx(kasc)(m0, k, start, end);
  return 1;
}

static inline int tk_vec_pfx(kdesc_lua) (lua_State *L)
{
  lua_settop(L, 4);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t k = tk_lua_checkunsigned(L, 2, "k");
  uint64_t start = tk_lua_optunsigned(L, 3, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 4, "end", m0->n);
  tk_vec_pfx(kdesc)(m0, k, start, end);
  return 1;
}

#ifdef tk_vec_pushbase
static inline int tk_vec_pfx(get_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  tk_vec_base x = tk_vec_pfx(get)(L, m0, i);
  tk_vec_pushbase(L, x);
  return 1;
}
#endif

#ifdef tk_vec_peekbase
static inline int tk_vec_pfx(set_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  tk_vec_base x = tk_vec_peekbase(L, 3);
  tk_vec_pfx(set)(L, m0, i, x);
  return 0;
}
#endif

#ifdef tk_vec_pushbase

static inline int tk_vec_pfx(table_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start, end;
  if (t == 1) {
    start = 0;
    end =  m0->n;
  } else if (t == 3) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, table, 1, "expected either 1 or 3 arguments (vec and optionally start/end indices)");
    return 0;
  }
  tk_vec_pfx(table)(L, m0, start, end);
  return 1;
}

static inline int tk_vec_pfx(ctable_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  uint64_t start, end;
  if (t == 2) {
    start = 0;
    end = cols;
  } else if (t == 4) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, ctable, 1, "expected either 2 or 4 arguments (vec, cols and optionally start/end rows)");
    return 0;
  }
  tk_vec_pfx(ctable)(L, m0, cols, start, end);
  return 1;
}

static inline int tk_vec_pfx(rtable_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  uint64_t start, end;
  if (t == 2) {
    start = 0;
    end = (m0->n + cols - 1) / cols;
  } else if (t == 4) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, rtable, 1, "expected either 2 or 4 arguments (vec, cols and optionally start/end cols)");
    return 0;
  }
  tk_vec_pfx(rtable)(L, m0, cols, start, end);
  return 1;
}

#endif

#ifndef tk_vec_limited

// TODO: Should be able to support this for complex types (vector of vectors),
// however we'll need to add the pushed subvector as an ephemeron, which likely
// needs to be defined by the user-macros, perhaps tk_vec_link or something
static inline int tk_vec_pfx(push_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base x = tk_vec_peekbase(L, 2);
  tk_vec_pfx(push)(m0, x);
  return 0;
}

static inline int tk_vec_pfx(magnitude_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  lua_pushnumber(L, tk_vec_pfx(magnitude)(m0));
  return 1;
}

static inline int tk_vec_pfx(scale_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base scale;
  uint64_t start, end;
  if (t == 2) {
    scale = tk_lua_checkdouble(L, 2, "scale");
    start = 0;
    end = m0->n;
  } else if (t == 3) {
    scale = tk_lua_checkdouble(L, 2, "scale");
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, rtable, 1, "expected either 2 or 4 arguments (vec, scale, or vec, scale, start, end)");
    return 0;
  }
  tk_vec_pfx(scale)(m0, scale, start, end);
  return 1;
}

static inline int tk_vec_pfx(add_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base add;
  uint64_t start, end;
  if (t == 2) {
    add = tk_lua_checkdouble(L, 2, "add");
    start = 0;
    end = m0->n;
  } else if (t == 3) {
    add = tk_lua_checkdouble(L, 2, "add");
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, rtable, 1, "expected either 2 or 4 arguments (vec, add or vec, add, start, end)");
    return 0;
  }
  tk_vec_pfx(add)(m0, add, start, end);
  return 1;
}

#ifdef tk_vec_abs
static inline int tk_vec_pfx(abs_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start, end;
  if (t == 1) {
    start = 0;
    end = m0->n;
  } else if (t == 2) {
    start = tk_lua_checkunsigned(L, 2, "start");
    end = tk_lua_checkunsigned(L, 3, "end");
  } else {
    tk_vec_err(L, rtable, 1, "expected either 1 or 3 arguments (vec or vec, start, end)");
    return 0;
  }
  tk_vec_pfx(abs)(m0, start, end);
  return 1;
}
#endif

static inline int tk_vec_pfx(exp_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base exp;
  uint64_t start, end;
  if (t == 2) {
    exp = tk_lua_checkdouble(L, 2, "exp");
    start = 0;
    end = m0->n;
  } else if (t == 3) {
    exp = tk_lua_checkdouble(L, 2, "exp");
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, rtable, 1, "expected either 2 or 4 arguments (vec, exp or vec, exp, start, end)");
    return 0;
  }
  tk_vec_pfx(exp)(m0, exp, start, end);
  return 1;
}

static inline int tk_vec_pfx(min_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  if (!m0->n)
    return 0;
  uint64_t minc;
  tk_vec_base minv = tk_vec_pfx(min)(m0, &minc);
  tk_vec_pushbase(L, minv);
  lua_pushinteger(L, (int64_t) minc);
  return 2;
}

static inline int tk_vec_pfx(rmins_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmins)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cmins_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmins)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(max_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  if (!m0->n)
    return 0;
  uint64_t maxc;
  tk_vec_base maxv = tk_vec_pfx(max)(m0, &maxc);
  tk_vec_pushbase(L, maxv);
  lua_pushinteger(L, (int64_t) maxc);
  return 2;
}

static inline int tk_vec_pfx(rmaxs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmaxs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cmaxs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmaxs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(sum_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pushbase(L, tk_vec_pfx(sum)(m0));
  return 1;
}

static inline int tk_vec_pfx(csums_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(csums)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rsums_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rsums)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(each_lua_iter) (lua_State *L)
{
  lua_settop(L, 0);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_vec_base v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  tk_vec_pushbase(L, v);
  return 1;
}

static inline int tk_vec_pfx(ieach_lua_iter) (lua_State *L)
{
  lua_settop(L, 0);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_vec_base v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, (int64_t) n);
  tk_vec_pushbase(L, v);
  return 2;
}

static inline int tk_vec_pfx(each_lua) (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_vec_pfx(each_lua_iter), 2);
  return 1;
}

static inline int tk_vec_pfx(ieach_lua) (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_vec_pfx(ieach_lua_iter), 2);
  return 1;
}

static inline int tk_vec_pfx(fill_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base x = tk_vec_peekbase(L, 2);
  tk_vec_pfx(fill)(v, x);
  return 1;
}

static inline int tk_vec_pfx(fill_indices_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(fill_indices)(v);
  return 1;
}

static inline int tk_vec_pfx(raw_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  const char *fmt = tk_lua_optstring(L, 2, "fmt", NULL);
  tk_vec_pfx(raw)(L, m0, fmt);
  return 1;
}

static inline int tk_vec_pfx(from_raw_lua) (lua_State *L)
{
  lua_settop(L, 2);
  const char *data = tk_lua_checkustring(L, 1, "raw");
  uint64_t s = tk_lua_checkunsigned(L, 2, "size");
  tk_vec_pfx(create)(L, s, (tk_vec_base *) data, NULL);
  return 1;
}

#endif

static luaL_Reg tk_vec_pfx(lua_mt_fns)[] =
{
  // Create, copy, etc.
  { "copy", tk_vec_pfx(copy_lua) }, // to other
  { "size", tk_vec_pfx(size_lua) },
  { "capacity", tk_vec_pfx(capacity_lua) },
  { "resize", tk_vec_pfx(resize_lua) },
  { "setn", tk_vec_pfx(setn_lua) },
  { "ensure", tk_vec_pfx(ensure_lua) },
  { "shrink", tk_vec_pfx(shrink_lua) },
  { "clear", tk_vec_pfx(clear_lua) },
  { "zero", tk_vec_pfx(zero_lua) },
  { "transpose", tk_vec_pfx(transpose_lua) },

#ifdef tk_vec_peekbase
  { "get", tk_vec_pfx(get_lua) },
#endif

#ifdef tk_vec_pushbase
  { "set", tk_vec_pfx(set_lua) },

  // Render to lua table (or table of table of columns/rows)
  { "table", tk_vec_pfx(table_lua) },
  { "ctable", tk_vec_pfx(ctable_lua) },
  { "rtable", tk_vec_pfx(rtable_lua) },
#endif

#ifndef tk_vec_limited

  // Update individual values
  { "push", tk_vec_pfx(push_lua) },

  // Scalar manipulation
  { "add", tk_vec_pfx(add_lua) },
  { "scale", tk_vec_pfx(scale_lua) },
  { "exp", tk_vec_pfx(exp_lua) },

#ifdef tk_vec_abs
  { "abs", tk_vec_pfx(abs_lua) },
#endif

  // Matrix multiplication
  { "multiply", tk_vec_pfx(multiply_lua) },
  { "dot", tk_vec_pfx(dot_lua) },

  // Magnitutes over vector columns or rows
  { "magnitude", tk_vec_pfx(magnitude_lua) },

  // Sums over vector columns or rows
  { "sum", tk_vec_pfx(sum_lua) },
  { "csums", tk_vec_pfx(csums_lua) },
  { "rsums", tk_vec_pfx(rsums_lua) },

  // Maximum over vector columns or rows
  { "max", tk_vec_pfx(max_lua) },
  { "cmaxs", tk_vec_pfx(rmaxs_lua) },
  { "rmaxs", tk_vec_pfx(cmaxs_lua) },

  // Minimum over vector columns or rows
  { "min", tk_vec_pfx(min_lua) },
  { "cmins", tk_vec_pfx(rmins_lua) },
  { "rmins", tk_vec_pfx(cmins_lua) },

  // Sort a vector (full or k)
  { "shuffle", tk_vec_pfx(shuffle_lua) },
  { "asc", tk_vec_pfx(asc_lua) },
  { "desc", tk_vec_pfx(desc_lua) },
  { "kasc", tk_vec_pfx(kasc_lua) },
  { "kdesc", tk_vec_pfx(kdesc_lua) },

  // Iterate
  { "each", tk_vec_pfx(each_lua) },
  { "ieach", tk_vec_pfx(ieach_lua) },

  // Fill with value or index
  { "fill", tk_vec_pfx(fill_lua) },
  { "fill_indices", tk_vec_pfx(fill_indices_lua) },

  // To/from char *
  { "raw", tk_vec_pfx(raw_lua) },

#endif

  { NULL, NULL }
};

static inline void tk_vec_pfx(suppress_unused_lua_mt_fns) (void)
  { (void) tk_vec_pfx(lua_mt_fns); }

static inline tk_vec_pfx(t) *tk_vec_pfx(create) (lua_State *L, size_t n, tk_vec_base *data, tk_vec_pfx(t) *v)
{
  tk_vec_pfx(t) v0;
  if (v == NULL) {
    v = tk_lua_newuserdata(L, tk_vec_pfx(t), tk_vec_mt, tk_vec_pfx(lua_mt_fns), tk_vec_pfx(destroy_lua)); // v (with mt)
    v0 = *v;
    kv_init(v0);
    kv_resize(tk_vec_base, v0, n);
    v0.n = n;
    *v = v0;
    return v;
  } else if (v->a && v->a != data) {
    free(v->a);
  }
  v->n = n;
  v->m = n;
  v->a = data;
  return v;
}

static luaL_Reg tk_vec_pfx(lua_fns)[] =
{
  { "create", tk_vec_pfx(create_lua) },
  { "destroy", tk_vec_pfx(destroy_lua) },
#ifndef tk_vec_limited
  { "from_raw", tk_vec_pfx(from_raw_lua) },
#endif
  { NULL, NULL }
};

static inline void tk_vec_pfx(suppress_unused_lua_fns) (void)
  { (void) tk_vec_pfx(lua_fns); }

#include <santoku/vec/undef.h>
