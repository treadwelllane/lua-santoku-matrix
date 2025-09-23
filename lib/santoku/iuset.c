#include <santoku/iuset.h>
#include <santoku/cvec.h>

int luaopen_santoku_iuset (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, tk_iuset_lua_fns); // t
  tk_iuset_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_iuset_lua_mt_fns); // t
  luaL_register(L, NULL, tk_iuset_lua_ext_fns); // t - register extended functions
  lua_pop(L, 2);
  return 1;
}
