local env = {

  name = "santoku-matrix",
  version = "0.0.38-1",
  variable_prefix = "TK_MATRIX",
  license = "MIT",
  public = true,

  cflags = { "-fopenmp", "-I$(shell luarocks show santoku --rock-dir)/include" },
  ldflags = { "-fopenmp" },

  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.261-1",
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
