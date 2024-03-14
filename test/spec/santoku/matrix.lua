local mtx = require("santoku.matrix")

local m0, m1

m0 = mtx.create(1, 10)

for i = 1, 5 do
  mtx.set(m0, 1, i, i)
end

for i = 1, 5 do
  assert(mtx.get(m0, 1, i) == i)
end

mtx.add(m0, -1)

for i = 1, 5 do
  assert(mtx.get(m0, 1, i) == i - 1)
end

m0 = mtx.create({ 1, 2, 3, 4, 5, 6 }, 2, 4)
assert(mtx.rows(m0) == 1)
assert(mtx.columns(m0) == 3)
assert(mtx.get(m0, 1, 1) == 2)
assert(mtx.get(m0, 1, 2) == 3)
assert(mtx.get(m0, 1, 3) == 4)

m0 = mtx.create({ 1, 2, 3, 4, 5, 6 }, 1, 3)
m1 = mtx.create({ 1, 2, 3, 4, 5, 6 }, 4, 6)

assert(mtx.dot(m0, m1) == 32)
