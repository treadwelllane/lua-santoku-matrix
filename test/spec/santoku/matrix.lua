local serialize  = require("santoku.serialize") -- luacheck: ignore
local it  = require("santoku.iter")
local err = require("santoku.error")
local ivec = require("santoku.ivec")
local rvec = require("santoku.rvec")
local dvec = require("santoku.dvec")
local tbl = require("santoku.table")

for _, vec in ipairs({ ivec, dvec }) do

  local m0, m1, m2

  m0 = vec.create(10)

  for i = 0, 4 do
    m0:set(i, i)
  end

  for i = 0, 4 do
    assert(m0:get(i) == i)
  end

  m0:add(-1)

  for i = 0, 4 do
    assert(m0:get(i) == i - 1)
  end

  m0 = vec.create({ 1, 2, 3, 4, 5, 6 })
  assert(m0:size() == 6)
  assert(m0:get(1) == 2)
  assert(m0:get(2) == 3)
  assert(m0:get(3) == 4)

  m0 = vec.create({ 1, 2, 3 })
  m1 = vec.create({ 4, 5, 6 })
  assert(m0:dot(m1) == 32)

  m0 = vec.create({ 1, 2, 3 })
  m1 = vec.create({ 4, 5, 6 })

  assert(m0:size() == 3)
  assert(m1:size() == 3)

  m2 = vec.from_raw(m0:raw(), 3)

  assert(m2:size() == 3)

  m0 = vec.create({ 1, 2, 3, 4, 5, 6 })

  assert(tbl.equals(m0:table(), {
    1, 2, 3, 4, 5, 6
  }))

  assert(tbl.equals(m0:rtable(3), {
    { 1, 2, 3 },
    { 4, 5, 6 }
  }))

  m0 = vec.create({
    10, 2, 3,
    4, 5, 6
  })

  assert(tbl.equals(m0:rdesc(3):rtable(3), {
    { 0, 2, 1 },
    { 2, 1, 0 }
  }))

  assert(tbl.equals(m0:casc(3):rtable(2), {
    { 1, 0 },
    { 0, 1 },
    { 0, 1 }
  }))

  m0 = vec.create({ 1, 6, 3, 2, 5, 2, 3, 6, 1, 4 })
  m0:asc()
  assert(tbl.equals(m0:table(), { 1, 1, 2, 2, 3, 3, 4, 5, 6, 6 }))

  if m0.flip_interleave then
    m0 = vec.create({ 0, 3 }) -- really 1001
    m0:flip_interleave(2, 2) -- should be 10011001 or 0 3 4 7
  end

  m0 = vec.create({ 1, 2, 3, 4 })
  local g = m0:each()
  assert(g() == 1)
  assert(g() == 2)
  assert(g() == 3)
  assert(g() == 4)
  assert(g() == nil)

  g = it.take(2, vec.create({ 1, 2, 3, 4 }):each())
  assert(g() == 1)
  assert(g() == 2)
  assert(g() == nil)

  if m0.flip_interleave then
    m0 = vec.create({ 0 })
    m1 = vec.create({ 1 })
    m1:add(2)
    m0:copy(m1)
    m0:flip_interleave(2, 2)
  end

  m0:persist(".vec.bin")
  m1 = vec.load(".vec.bin")

  assert(tbl.equals(m0:table(), m1:table()))

end

local heap0 = rvec.create(10)
local heap1 = rvec.create()
heap0:setn(0)
for i = 1, 100 do
  local r = math.random(100)
  heap0:hmax(i, r)
  heap1:push(i, r)
end
heap0:asc()
heap1:asc()
heap1:setn(10)

for i = 1, 10 do
  local i0, d0 = heap0:get(i - 1)
  local i1, d1 = heap1:get(i - 1)
  err.assert(i0 == i1)
  err.assert(d0 == d1)
end
