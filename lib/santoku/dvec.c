#define _GNU_SOURCE

#include <santoku/dvec.h>
#include <santoku/threads.h>

int luaopen_santoku_dvec (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, tk_dvec_lua_fns); // t
  tk_dvec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_dvec_lua_mt_fns); // t
  luaL_register(L, NULL, tk_dvec_lua_mt_ext_fns); // t
  lua_pop(L, 2);
  return 1;
}
