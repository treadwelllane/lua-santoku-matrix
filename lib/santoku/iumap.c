#include <santoku/iuset.h>
#include <santoku/iumap.h>
#include <santoku/cvec.h>

int luaopen_santoku_iumap (lua_State *L)
{
  lua_newtable(L); // t
  luaL_register(L, NULL, tk_iumap_lua_fns); // t
  tk_iumap_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_iumap_lua_mt_fns); // t
  luaL_register(L, NULL, tk_iumap_lua_ext_fns); // t - register extended functions
  lua_pop(L, 2);
  return 1;
}
