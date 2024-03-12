local env = {

  name = "santoku-matrix",
  version = "0.0.5-1",
  variable_prefix = "TK_MATRIX",
  public = true,

  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.202-1",
  },

  test = {
    dependencies = {
      "luacov == 0.15.0-1",
    }
  },

}

env.homepage = "https://github.com/treadwelllane/lua-" .. env.name
env.tarball = env.name .. "-" .. env.version .. ".tar.gz"
env.download = env.homepage .. "/releases/download/" .. env.version .. "/" .. env.tarball

return {
  type = "lib",
  env = env,
}
