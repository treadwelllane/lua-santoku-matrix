#include <santoku/iuset.h>
#include <santoku/dvec.h>

static inline int tk_dvec_center_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_dvec_t *M = tk_dvec_peek(L, 1, "dvec");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_dims = tk_lua_checkunsigned(L, 3, "n_dims");
  tk_dvec_center(M->a, n_samples, n_dims);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_dvec_rnorml2_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_dvec_t *M = tk_dvec_peek(L, 1, "dvec");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_dims = tk_lua_checkunsigned(L, 3, "n_dims");
  tk_dvec_rnorml2(M->a, n_samples, n_dims);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_dvec_scores_kaiser_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_kaiser(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_max_curvature_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_curvature(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_lmethod_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_lmethod(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_max_gap_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_gap(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_max_drop_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_drop(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_max_acceleration_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_acceleration(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_tolerance_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double tolerance = luaL_optnumber(L, 2, 1e-3);
  size_t start, end;
  tk_dvec_scores_tolerance(scores->a, scores->n, tolerance, &start, &end);
  double start_val = (start < scores->n) ? scores->a[start] : 0.0;
  double end_val = (end < scores->n) ? scores->a[end] : 0.0;
  lua_pushnumber(L, (lua_Number)start_val);
  lua_pushnumber(L, (lua_Number)start);
  lua_pushnumber(L, (lua_Number)end_val);
  lua_pushnumber(L, (lua_Number)end);
  return 4;
}

static inline int tk_dvec_scores_plateau_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double tolerance = luaL_optnumber(L, 2, 1e-3);
  double val;
  size_t idx = tk_dvec_scores_plateau(scores->a, scores->n, tolerance, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_multiply_blas_lua (lua_State *L)
{
  lua_settop(L, 6);
  tk_dvec_t *a = tk_dvec_peek(L, 1, "a");
  tk_dvec_t *b = tk_dvec_peek(L, 2, "b");
  tk_dvec_t *c = tk_dvec_peek(L, 3, "c");
  uint64_t k = tk_lua_checkunsigned(L, 4, "k");
  bool transpose_a = lua_toboolean(L, 5);
  bool transpose_b = lua_toboolean(L, 6);
  tk_dvec_multiply_override(a, b, c, k, transpose_a, transpose_b);
  lua_pushvalue(L, 3);
  return 1;
}

static inline int tk_dvec_rsums_blas_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *m = tk_dvec_peek(L, 1, "dvec");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_dvec_rsums_override(L, m, cols);
  return 1;
}

static inline int tk_dvec_csums_blas_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *m = tk_dvec_peek(L, 1, "dvec");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_dvec_csums_override(L, m, cols);
  return 1;
}

static inline int tk_dvec_rmags_blas_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *m = tk_dvec_peek(L, 1, "dvec");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_dvec_rmags_override(L, m, cols);
  return 1;
}

static inline int tk_dvec_cmags_blas_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *m = tk_dvec_peek(L, 1, "dvec");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_dvec_cmags_override(L, m, cols);
  return 1;
}

static luaL_Reg tk_dvec_lua_mt_ext2_fns[] =
{
  { "center", tk_dvec_center_lua },
  { "rnorml2", tk_dvec_rnorml2_lua },
  { "multiply", tk_dvec_multiply_blas_lua },
  { "rsums", tk_dvec_rsums_blas_lua },
  { "csums", tk_dvec_csums_blas_lua },
  { "rmags", tk_dvec_rmags_blas_lua },
  { "cmags", tk_dvec_cmags_blas_lua },
  { "scores_kaiser", tk_dvec_scores_kaiser_lua },
  { "scores_max_curvature", tk_dvec_scores_max_curvature_lua },
  { "scores_lmethod", tk_dvec_scores_lmethod_lua },
  { "scores_max_gap", tk_dvec_scores_max_gap_lua },
  { "scores_max_drop", tk_dvec_scores_max_drop_lua },
  { "scores_max_acceleration", tk_dvec_scores_max_acceleration_lua },
  { "scores_tolerance", tk_dvec_scores_tolerance_lua },
  { "scores_plateau", tk_dvec_scores_plateau_lua },
  { NULL, NULL }
};

int luaopen_santoku_dvec (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, tk_dvec_lua_fns); // t
  tk_dvec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_dvec_lua_mt_fns); // t
  luaL_register(L, NULL, tk_dvec_lua_mt_ext_fns); // t
  luaL_register(L, NULL, tk_dvec_lua_mt_ext2_fns); // t
  lua_pop(L, 2);
  return 1;
}
