#ifndef TK_IVEC_EXT_H
#define TK_IVEC_EXT_H

#include <santoku/iuset.h>

static inline void tk_ivec_push_selected (
  lua_State *L,
  tk_iuset_t *selected
) {
  // Create matrix
  tk_ivec_t *top_v = tk_ivec_create(L, tk_iuset_size(selected), 0, 0);
  uint64_t write = 0;
  int64_t sel;
  tk_iuset_foreach(selected, sel, ({
    top_v->a[write ++] = sel;
  }));
  tk_iuset_destroy(selected);
}

static inline tk_rvec_t *tk_rvec_rankings (
  lua_State *L,
  tk_dvec_t *scores,
  uint64_t n_visible,
  uint64_t n_hidden
) {
  // Pull scores into ranking structure
  tk_rvec_t *rankings = tk_rvec_create(L, n_hidden * n_visible, NULL, NULL);
  for (uint64_t h = 0; h < n_hidden; h ++)
    for (uint64_t v = 0; v < n_visible; v ++)
      rankings->a[h * n_visible + v] = tk_rank((int64_t) v, scores->a[h * n_visible + v]);

  // Sort best visible by hidden
  for (uint64_t j = 0; j < n_hidden; j ++)
    tk_rvec_desc(rankings);

  // Return rankings
  return rankings;
}

static inline void tk_ivec_select_union (
  lua_State *L,
  tk_iuset_t *selected,
  tk_rvec_t *rankings,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  uint64_t trunc
) {
  // Select top-k union
  int kha;
  uint64_t *offsets = tk_malloc(L, n_hidden * sizeof(uint64_t));
  memset(offsets, 0, n_hidden * sizeof(uint64_t));
  while (tk_iuset_size(selected) < top_k) {
    bool advanced = false;
    for (uint64_t j = 0; j < n_hidden; j++) {
      tk_rank_t *rankings_h = rankings->a + j * n_visible;
      while (offsets[j] < n_visible) {
        tk_rank_t candidate = rankings_h[offsets[j]];
        offsets[j] ++;
        if (candidate.i >= (int64_t) (n_visible - trunc))
          continue;
        tk_iuset_put(selected, (int64_t) candidate.i, &kha);
        advanced = true;
        break;
      }
      if (tk_iuset_size(selected) >= top_k)
        break;
    }
    if (!advanced)
      break;
  }
  // Cleanup
  free(offsets);
}

static inline void tk_ivec_top_generic (
  lua_State *L,
  tk_dvec_t *scores,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  uint64_t trunc
) {
  // Select top-k
  tk_iuset_t *selected = tk_iuset_create();
  tk_rvec_t *rankings = tk_rvec_rankings(L, scores, n_visible, n_hidden);
  tk_ivec_select_union(L, selected, rankings, n_visible, n_hidden, top_k, trunc);
  tk_ivec_destroy(rankings);
  // Push final matrix
  tk_ivec_push_selected(L, selected);
}

#endif
