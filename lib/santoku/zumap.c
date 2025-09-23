#include <santoku/zumap.h>
#include <santoku/cvec.h>

int luaopen_santoku_zumap (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, tk_zumap_lua_fns); // t
  tk_zumap_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_zumap_lua_mt_fns); // t
  luaL_register(L, NULL, tk_zumap_lua_ext_fns); // t - register extended functions
  lua_pop(L, 2);
  return 1;
}