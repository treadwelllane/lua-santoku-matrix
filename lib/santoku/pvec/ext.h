#ifndef TK_PVEC_EXT_H
#define TK_PVEC_EXT_H

static inline tk_pvec_t *tk_pvec_from_ivec (
  lua_State *L,
  tk_ivec_t *I
) {
  tk_pvec_t *P = tk_pvec_create(L, I->n, 0, 0);
  for (int64_t i = 0; i < (int64_t) I->n; i ++)
    P->a[i] = tk_pair(i, I->a[i]);
  return P;
}

static inline tk_ivec_t *tk_pvec_keys (
  lua_State *L,
  tk_pvec_t *P
) {
  tk_ivec_t *I = tk_ivec_create(L, P->n, 0, 0);
  for (uint64_t i = 0; i < P->n; i ++)
    I->a[i] = P->a[i].i;
  return I;
}

static inline tk_ivec_t *tk_pvec_values (
  lua_State *L,
  tk_pvec_t *P
) {
  tk_ivec_t *I = tk_ivec_create(L, P->n, 0, 0);
  for (uint64_t i = 0; i < P->n; i ++)
    I->a[i] = P->a[i].p;
  return I;
}

static inline int tk_pvec_each0_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_pvec_t *m0 = tk_pvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_pair_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, v.i);
  lua_pushinteger(L, v.p);
  return 2;
}

static inline int tk_pvec_ieach0_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_pvec_t *m0 = tk_pvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_pair_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, (int64_t) n);
  lua_pushinteger(L, v.i);
  lua_pushinteger(L, v.p);
  return 3;
}

static inline int tk_pvec_each0_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_pvec_each0_lua_iter, 2);
  return 1;
}

static inline int tk_pvec_ieach0_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_pvec_ieach0_lua_iter, 2);
  return 1;
}

#endif
