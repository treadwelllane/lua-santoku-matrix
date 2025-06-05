#define _GNU_SOURCE

#include <santoku/ivec.h>

static inline int tk_ivec_flip_interleave_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
  tk_ivec_flip_interleave(L, m0, n_samples, n_features);
  return 0;
}

static inline int tk_ivec_score_entropy_lua (
  lua_State *L
) {
  lua_settop(L, 4);
  const char *codes = tk_lua_checkustring(L, 1, "codes");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 3, "hidden");
  unsigned int n_threads = tk_threads_getn(L, 4, "threads", NULL);
  tk_ivec_score_entropy(L, (char *) codes, n_samples, n_hidden, n_threads);
  return 1;
}

static inline int tk_ivec_score_chi2_lua (
  lua_State *L
) {
  lua_settop(L, 6);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  unsigned int n_threads = tk_threads_getn(L, 6, "threads", NULL);
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_type(L, 2) == LUA_TSTRING) {
    codes = (char *) luaL_checkstring(L, 2);
  } else if (lua_type(L, 2) == LUA_TLIGHTUSERDATA) {
    codes = (char *) lua_touserdata(L, 2);
  } else {
    tk_ivec_t *m1 = tk_ivec_peek(L, 2, "labels");
    n_samples = m1->n < n_samples ? m1->n : n_samples;
    labels = m1;
  }
  tk_ivec_score_chi2(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  return 1;
}

static inline int tk_ivec_score_mi_lua (
  lua_State *L
) {
  lua_settop(L, 6);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  unsigned int n_threads = tk_threads_getn(L, 6, "threads", NULL);
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_type(L, 2) == LUA_TSTRING) {
    codes = (char *) luaL_checkstring(L, 2);
  } else if (lua_type(L, 2) == LUA_TLIGHTUSERDATA) {
    codes = (char *) lua_touserdata(L, 2);
  } else {
    tk_ivec_t *m1 = tk_ivec_peek(L, 2, "labels");
    n_samples = m1->n < n_samples ? m1->n : n_samples;
    labels = m1;
  }
  tk_ivec_score_mi(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  return 1;
}

static inline int tk_ivec_top_mi_lua (lua_State *L)
{
  lua_settop(L, 7);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 6, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 7, "threads", NULL);
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_type(L, 2) == LUA_TSTRING) {
    codes = (char *) luaL_checkstring(L, 2);
  } else if (lua_type(L, 2) == LUA_TLIGHTUSERDATA) {
    codes = (char *) lua_touserdata(L, 2);
  } else {
    tk_ivec_t *m1 = tk_ivec_peek(L, 2, "labels");
    n_samples = m1->n < n_samples ? m1->n : n_samples;
    labels = m1;
  }
  tk_ivec_top_mi(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, top_k, n_threads);
  return 1;
}

static inline int tk_ivec_top_entropy_lua (lua_State *L)
{
  lua_settop(L, 5);
  const char *codes = tk_lua_checkustring(L, 1, "codes");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 3, "hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 4, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 5, "threads", NULL);
  tk_ivec_top_entropy(L, (char *) codes, n_samples, n_hidden, top_k, n_threads);
  return 1;
}

static inline int tk_ivec_top_chi2_lua (lua_State *L)
{
  lua_settop(L, 7);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 6, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 7, "threads", NULL);
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_type(L, 2) == LUA_TSTRING) {
    codes = (char *) luaL_checkstring(L, 2);
  } else if (lua_type(L, 2) == LUA_TLIGHTUSERDATA) {
    codes = (char *) lua_touserdata(L, 2);
  } else {
    tk_ivec_t *m1 = tk_ivec_peek(L, 2, "labels");
    n_samples = m1->n < n_samples ? m1->n : n_samples;
    labels = m1;
  }
  tk_ivec_top_chi2(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, top_k, n_threads);
  return 1;
}

static inline int tk_ivec_filter_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  tk_ivec_t *top_v = tk_ivec_peek(L, 2, "top_v");
  uint64_t n_visible = tk_lua_checkunsigned(L, 3, "visible");
  tk_ivec_filter(L, set_bits, top_v, n_visible);
  return 0;
}

static inline int tk_ivec_raw_bitmap_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
  size_t l;
  char *r = tk_ivec_raw_bitmap(L, set_bits, n_samples, n_features, &l);
  lua_pushlstring(L, r, l);
  free(r);
  return 1;
}

static inline int tk_ivec_from_bitmap_lua (lua_State *L)
{
  lua_settop(L, 3);
  const char *bm = tk_lua_checkustring(L, 1, "bitmap");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
  tk_ivec_from_bitmap(L, bm, n_samples, n_features);
  return 1;
}

static inline int tk_ivec_extend_bits_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_ivec_t *base = tk_ivec_peek(L, 1, "base_bits");
  tk_ivec_t *ext = tk_ivec_peek(L, 2, "ext_bits");
  uint64_t n_feat = tk_lua_checkunsigned(L, 3, "features");
  uint64_t n_extfeat = tk_lua_checkunsigned(L, 4, "extended");
  tk_ivec_extend_bits(L, base, ext, n_feat, n_extfeat);
  return 0;
}

static luaL_Reg tk_ivec_lua_mt_ext2_fns[] =
{
  { "top_chi2", tk_ivec_top_chi2_lua },
  { "top_mi", tk_ivec_top_mi_lua },
  { "top_entropy", tk_ivec_top_entropy_lua },
  { "score_chi2", tk_ivec_score_chi2_lua },
  { "score_mi", tk_ivec_score_mi_lua },
  { "score_entropy", tk_ivec_score_entropy_lua },
  { "flip_interleave", tk_ivec_flip_interleave_lua },
  { "filter", tk_ivec_filter_lua },
  { "raw_bitmap", tk_ivec_raw_bitmap_lua },
  { "from_bitmap", tk_ivec_from_bitmap_lua },
  { "extend_bits", tk_ivec_extend_bits_lua },
  { NULL, NULL }
};

int luaopen_santoku_ivec (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, tk_ivec_lua_fns); // t
  tk_ivec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_ivec_lua_mt_fns); // t
  luaL_register(L, NULL, tk_ivec_lua_mt_ext_fns); // t
  luaL_register(L, NULL, tk_ivec_lua_mt_ext2_fns); // t
  lua_pop(L, 2);
  return 1;
}
