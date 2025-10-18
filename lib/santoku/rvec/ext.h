#ifndef TK_RVEC_EXT_H
#define TK_RVEC_EXT_H

#include <santoku/dumap.h>
#include <santoku/iumap.h>
#include <santoku/pvec.h>
#include <santoku/evec.h>
#include <math.h>

static inline tk_rvec_t *tk_rvec_from_dvec (
  lua_State *L,
  tk_dvec_t *D
) {
  tk_rvec_t *R = tk_rvec_create(L, D->n, 0, 0);
  for (int64_t i = 0; i < (int64_t) D->n; i ++)
    R->a[i] = tk_rank(i, D->a[i]);
  return R;
}

static inline int tk_rvec_split (
  tk_rvec_t *P,
  tk_ivec_t *K,
  tk_dvec_t *V,
  bool append
) {
  if (!append) {
    tk_ivec_clear(K);
    tk_dvec_clear(V);
  }
  if (tk_ivec_ensure(K, K->n + P->n) != 0)
    return -1;
  if (tk_dvec_ensure(V, V->n + P->n) != 0)
    return -1;
  for (uint64_t i = 0; i < P->n; i ++) {
    tk_rank_t r = P->a[i];
    K->a[K->n ++] = r.i;
    V->a[V->n ++] = r.d;
  }
  return 0;
}

static inline tk_ivec_t *tk_rvec_keys (
  lua_State *L,
  tk_rvec_t *P,
  tk_ivec_t *out
) {
  tk_ivec_t *I = out ? out : tk_ivec_create(L, P->n, 0, 0);
  if (out)
    tk_ivec_ensure(I, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    I->a[i] = P->a[i].i;
  if (out)
    I->n = P->n;
  return I;
}

static inline tk_dvec_t *tk_rvec_values (
  lua_State *L,
  tk_rvec_t *P,
  tk_dvec_t *out
) {
  tk_dvec_t *I = out ? out : tk_dvec_create(L, P->n, 0, 0);
  if (out)
    tk_dvec_ensure(I, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    I->a[i] = P->a[i].d;
  if (out)
    I->n = P->n;
  return I;
}

static inline int tk_rvec_each_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_rvec_t *m0 = tk_rvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_rank_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, v.i);
  lua_pushnumber(L, v.d);
  return 2;
}

static inline int tk_rvec_ieach_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_rvec_t *m0 = tk_rvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_rank_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, (int64_t) n);
  lua_pushinteger(L, v.i);
  lua_pushnumber(L, v.d);
  return 3;
}

static inline int tk_rvec_each_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_rvec_each_lua_iter, 2);
  return 1;
}

static inline int tk_rvec_ieach_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_rvec_ieach_lua_iter, 2);
  return 1;
}

static inline void tk_rvec_ranks (
  tk_rvec_t *v,
  tk_dumap_t *r
) {
  tk_dumap_clear(r);
  int kha;
  uint32_t khi;
  for (size_t i = 0; i < v->n; i ++) {
    double rank = (double) i;
    size_t count = 1;
    while (i + 1 < v->n && v->a[i].d == v->a[i + 1].d) {
      count ++;
      i ++;
    }
    double average_rank = (rank + (rank + count - 1)) / 2.0;
    for (size_t j = 0; j < count; j ++) {
      khi = tk_dumap_put(r, v->a[i - j].i, &kha);
      tk_dumap_setval(r, khi, average_rank);
    }
  }
}

static inline double tk_csr_spearman(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_pvec_t *ranks_b,
  tk_dumap_t *rank_buffer_a,
  tk_dumap_t *rank_buffer_b
) {
  uint64_t n_a = (uint64_t)(end_a - start_a);
  if (n_a == 0 || !ranks_b || ranks_b->n == 0)
    return 0.0;
  tk_dumap_clear(rank_buffer_b);
  int kha;
  for (uint64_t i = 0; i < ranks_b->n; i++) {
    double rank = (double)i;
    uint64_t count = 1;
    int64_t hamming = ranks_b->a[i].p;
    while (i + 1 < ranks_b->n && ranks_b->a[i + 1].p == hamming) {
      count++;
      i++;
    }
    double average_rank = (rank + (rank + count - 1)) / 2.0;
    for (uint64_t j = 0; j < count; j++) {
      uint32_t khi = tk_dumap_put(rank_buffer_b, ranks_b->a[i - j].i, &kha);
      tk_dumap_setval(rank_buffer_b, khi, average_rank);
    }
  }
  tk_dumap_clear(rank_buffer_a);
  for (int64_t j = start_a; j < end_a; j++) {
    double rank = (double)(j - start_a);
    uint64_t count = 1;
    double weight = weights_a->a[j];
    while (j + 1 < end_a && weights_a->a[j + 1] == weight) {
      count++;
      j++;
    }
    double average_rank = (rank + (rank + count - 1)) / 2.0;
    for (uint64_t k = 0; k < count; k++) {
      uint32_t khi = tk_dumap_put(rank_buffer_a, neighbors_a->a[(uint64_t) j - k], &kha);
      tk_dumap_setval(rank_buffer_a, khi, average_rank);
    }
  }
  double sum_squared_diff = 0.0;
  uint64_t n = 0;
  int64_t neighbor;
  double ra;
  tk_umap_foreach(rank_buffer_a, neighbor, ra, ({
    uint32_t khi = tk_dumap_get(rank_buffer_b, neighbor);
    if (khi != tk_dumap_end(rank_buffer_b)) {
      double rb = tk_dumap_val(rank_buffer_b, khi);
      double diff = ra - rb;
      sum_squared_diff += diff * diff;
      n++;
    }
  }))
  if (n <= 1)
    return 1.0;
  double n_double = (double)n;
  return 1.0 - (6.0 * sum_squared_diff) / (n_double * (n_double * n_double - 1.0));
}

static inline double tk_csr_biserial(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_iuset_t *group_1
) {
  uint64_t n_a = (uint64_t)(end_a - start_a);
  if (n_a == 0 || !group_1)
    return 0.0;
  double rank_sum_1 = 0.0;
  uint64_t count_0 = 0, count_1 = 0;

  // Handle tied weights using average ranks
  for (int64_t j = start_a; j < end_a; j++) {
    double rank = (double)(j - start_a);
    uint64_t count = 1;
    double weight = weights_a->a[j];

    // Detect ties: consecutive positions with same weight
    while (j + 1 < end_a && weights_a->a[j + 1] == weight) {
      count++;
      j++;
    }

    // Average rank for tied positions
    double average_rank = (rank + (rank + count - 1)) / 2.0 + 1.0;

    // Assign average rank to all tied neighbors
    for (uint64_t k = 0; k < count; k++) {
      int64_t neighbor = neighbors_a->a[(uint64_t)j - k];
      if (tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1)) {
        rank_sum_1 += average_rank;
        count_1++;
      } else {
        count_0++;
      }
    }
  }

  if (count_0 == 0 || count_1 == 0 || n_a == 0)
    return 0.0;
  double U1 = rank_sum_1 - (count_1 * (count_1 + 1)) / 2.0;
  return 1.0 - (2.0 * U1) / (count_0 * count_1);
}

static inline double tk_csr_variance_ratio(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_iuset_t *group_1
) {
  uint64_t n_a = (uint64_t)(end_a - start_a);
  if (n_a == 0 || !group_1)
    return 0.0;

  // Build map of neighbor -> average_rank (handling ties)
  tk_dumap_t *rank_map = tk_dumap_create(NULL, 0);
  if (!rank_map)
    return 0.0;

  int kha;
  for (int64_t j = start_a; j < end_a; j++) {
    double rank = (double)(j - start_a);
    uint64_t count = 1;
    double weight = weights_a->a[j];

    // Detect ties: consecutive positions with same weight
    while (j + 1 < end_a && weights_a->a[j + 1] == weight) {
      count++;
      j++;
    }

    // Average rank for tied positions
    double average_rank = (rank + (rank + count - 1)) / 2.0 + 1.0;

    // Store averaged rank for each tied neighbor
    for (uint64_t k = 0; k < count; k++) {
      int64_t neighbor = neighbors_a->a[(uint64_t)j - k];
      uint32_t khi = tk_dumap_put(rank_map, neighbor, &kha);
      tk_dumap_setval(rank_map, khi, average_rank);
    }
  }

  // Compute group means and overall mean
  double rank_sum_0 = 0.0, rank_sum_1 = 0.0;
  uint64_t count_0 = 0, count_1 = 0;
  int64_t neighbor;
  double rank;
  tk_umap_foreach(rank_map, neighbor, rank, ({
    if (tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1)) {
      rank_sum_1 += rank;
      count_1++;
    } else {
      rank_sum_0 += rank;
      count_0++;
    }
  }))

  if (count_0 == 0 || count_1 == 0) {
    tk_dumap_destroy(rank_map);
    return 0.0;
  }

  double mean_rank_0 = rank_sum_0 / count_0;
  double mean_rank_1 = rank_sum_1 / count_1;
  double overall_mean = (rank_sum_0 + rank_sum_1) / (count_0 + count_1);

  // Compute sum of squares between and within groups
  // SS_between = sum of n_i * (mean_i - overall_mean)^2
  double ss_between = count_0 * (mean_rank_0 - overall_mean) * (mean_rank_0 - overall_mean)
                    + count_1 * (mean_rank_1 - overall_mean) * (mean_rank_1 - overall_mean);

  // SS_within = sum of (x_ij - mean_i)^2 for all observations
  double ss_within = 0.0;
  tk_umap_foreach(rank_map, neighbor, rank, ({
    if (tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1)) {
      double diff = rank - mean_rank_1;
      ss_within += diff * diff;
    } else {
      double diff = rank - mean_rank_0;
      ss_within += diff * diff;
    }
  }))

  tk_dumap_destroy(rank_map);

  // Compute eta-squared (effect size measure from ANOVA)
  // η² = SS_between / SS_total, ranges [0, 1]
  // 0 = no separation (groups have same mean), 1 = perfect separation
  double ss_total = ss_between + ss_within;
  if (ss_total <= 1e-12) {
    // All ranks identical - no variation to explain
    return 0.0;
  }

  return ss_between / ss_total;
}

#endif
