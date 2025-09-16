#ifndef TK_IUSET_H
#define TK_IUSET_H

#include <santoku/klib.h>
#include <santoku/ivec.h>
#include <santoku/cvec.h>
#include <santoku/cvec/ext.h>
#include <santoku/rvec.h>
#include <santoku/rvec/ext.h>
#include <santoku/dvec.h>
#include <santoku/iumap.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

KHASH_INIT(tk_iuset, int64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal)
typedef khash_t(tk_iuset) tk_iuset_t;

#define tk_iuset_put(...) kh_put(tk_iuset, __VA_ARGS__)
#define tk_iuset_get(...) kh_get(tk_iuset, __VA_ARGS__)
#define tk_iuset_del(...) kh_del(tk_iuset, __VA_ARGS__)
#define tk_iuset_exist(...) kh_exist(__VA_ARGS__)
#define tk_iuset_key(...) kh_key(__VA_ARGS__)
#define tk_iuset_value(...) kh_value(__VA_ARGS__)
#define tk_iuset_begin(...) kh_begin(__VA_ARGS__)
#define tk_iuset_end(...) kh_end(__VA_ARGS__)
#define tk_iuset_size(...) kh_size(__VA_ARGS__)
#define tk_iuset_resize(...) kh_resize(tk_iuset, __VA_ARGS__)
#define tk_iuset_clear(...) kh_clear(tk_iuset, __VA_ARGS__)
#define tk_iuset_destroy(...) kh_destroy(tk_iuset, __VA_ARGS__)
#define tk_iuset_create() kh_init(tk_iuset)
#define tk_iuset_contains(h, v) (tk_iuset_get(h, v) != tk_iuset_end(h))
// Note: this changes the default behavior so that a value var isn't required
#define tk_iuset_foreach(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); __i ++) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

static inline void tk_iuset_union (tk_iuset_t *a, tk_iuset_t *b)
{
  int kha;
  int64_t x;
  tk_iuset_foreach(b, x, ({
    tk_iuset_put(a, x, &kha);
  }));
}

static inline double tk_iuset_jaccard (tk_iuset_t *a, tk_iuset_t *b)
{
  uint64_t intersection = 0;
  uint64_t union_count = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (tk_iuset_contains(b, x))
      intersection ++;
  }));
  union_count = tk_iuset_size(a) + tk_iuset_size(b) - intersection;
  if (union_count == 0)
    return 0.0;
  return (double) intersection / (double) union_count;
}

static inline void tk_iuset_intersect (tk_iuset_t *a, tk_iuset_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (!tk_iuset_contains(b, x))
      kh_del(tk_iuset, a, i);
  }))
}

static inline void tk_iuset_difference (tk_iuset_t *a, tk_iuset_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (tk_iuset_contains(b, x))
      kh_del(tk_iuset, a, i);
  }))
}

// Operations with iumap (using map keys as the secondary set)
#include <santoku/iumap.h>

static inline void tk_iuset_union_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  int kha;
  for (khint_t i = tk_iumap_begin(b); i < tk_iumap_end(b); i ++) {
    if (tk_iumap_exist(b, i)) {
      int64_t key = tk_iumap_key(b, i);
      tk_iuset_put(a, key, &kha);
    }
  }
}

static inline void tk_iuset_intersect_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    khint_t khi = tk_iumap_get(b, x);
    if (khi == tk_iumap_end(b))
      kh_del(tk_iuset, a, i);
  }))
}

static inline void tk_iuset_difference_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    khint_t khi = tk_iumap_get(b, x);
    if (khi != tk_iumap_end(b))
      kh_del(tk_iuset, a, i);
  }))
}

static inline tk_iuset_t *tk_iuset_from_ivec (tk_ivec_t *v)
{
  int kha;
  tk_iuset_t *s = tk_iuset_create();
  for (uint64_t i = 0; i < v->n; i ++)
    tk_iuset_put(s, v->a[i], &kha);
  return s;
}

static inline int tk_ivec_bits_filter (
  tk_ivec_t *set_bits,
  tk_ivec_t *top_v,
  tk_ivec_t *sample_ids,
  uint64_t n_visible
) {
  // Return immediately if neither filter list is provided
  if ((top_v == NULL || top_v->n == 0) && (sample_ids == NULL || sample_ids->n == 0))
    return 0;

  // Create feature map if filtering by features
  int64_t *feature_map = NULL;
  uint64_t n_new_features = n_visible;
  if (top_v != NULL && top_v->n > 0) {
    feature_map = malloc(n_visible * sizeof(int64_t));
    if (!feature_map) return -1;
    for (uint64_t i = 0; i < n_visible; i ++)
      feature_map[i] = -1;
    for (uint64_t i = 0; i < top_v->n; i ++) {
      int64_t feat = top_v->a[i];
      if (feat >= 0 && (uint64_t) feat < n_visible)
        feature_map[feat] = (int64_t) i;
    }
    n_new_features = top_v->n;
  }

  // Create sample map if filtering by samples
  tk_iuset_t *sample_set = NULL;
  int64_t *sample_map = NULL;
  uint64_t max_sample = 0;
  if (sample_ids != NULL && sample_ids->n > 0) {
    // First find max sample to allocate map
    for (uint64_t i = 0; i < set_bits->n; i ++) {
      if (set_bits->a[i] >= 0) {
        uint64_t s = (uint64_t) set_bits->a[i] / n_visible;
        if (s > max_sample) max_sample = s;
      }
    }

    // Create sample set for quick lookup
    sample_set = tk_iuset_from_ivec(sample_ids);

    // Create sample remapping
    sample_map = malloc((max_sample + 1) * sizeof(int64_t));
    if (!sample_map) {
      if (feature_map) free(feature_map);
      tk_iuset_destroy(sample_set);
      return -1;
    }
    for (uint64_t i = 0; i <= max_sample; i ++)
      sample_map[i] = -1;

    // Map old sample indices to new compacted indices
    uint64_t new_idx = 0;
    for (uint64_t i = 0; i < sample_ids->n; i ++) {
      int64_t sid = sample_ids->a[i];
      if (sid >= 0 && (uint64_t) sid <= max_sample)
        sample_map[sid] = (int64_t) new_idx ++;
    }
  }

  size_t write = 0;
  for (size_t i = 0; i < set_bits->n; i ++) {
    int64_t val = set_bits->a[i];
    if (val < 0)
      continue;
    uint64_t sample = (uint64_t) val / n_visible;
    uint64_t feature = (uint64_t) val % n_visible;

    // Check and remap sample if filtering
    int64_t new_sample = (int64_t) sample;
    if (sample_set != NULL) {
      if (!tk_iuset_contains(sample_set, (int64_t) sample))
        continue;
      if (sample <= max_sample)
        new_sample = sample_map[sample];
      if (new_sample < 0)
        continue;
    }

    // Check and remap feature if filtering
    int64_t new_feature = (int64_t) feature;
    if (feature_map != NULL) {
      new_feature = feature_map[feature];
      if (new_feature < 0)
        continue;
    }

    // Store remapped index
    set_bits->a[write ++] = new_sample * (int64_t) n_new_features + new_feature;
  }
  set_bits->n = write;

  if (feature_map != NULL)
    free(feature_map);
  if (sample_map != NULL)
    free(sample_map);
  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);
  return 0;
}

static inline void tk_cvec_bits_filter (
  tk_cvec_t *bitmap,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_features
) {
  uint64_t n_samples = bitmap->n / TK_CVEC_BITS_BYTES(n_features);
  // Return immediately if neither filter list is provided
  if ((selected_features == NULL || selected_features->n == 0) &&
      (sample_ids == NULL || sample_ids->n == 0))
    return;

  // Create sample set if filtering by samples
  tk_iuset_t *sample_set = NULL;
  if (sample_ids != NULL && sample_ids->n > 0) {
    sample_set = tk_iuset_from_ivec(sample_ids);
  }

  uint64_t n_selected_features = (selected_features != NULL && selected_features->n > 0)
    ? selected_features->n : n_features;
  uint64_t in_bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t out_bytes_per_sample = TK_CVEC_BITS_BYTES(n_selected_features);
  uint8_t *data = (uint8_t *)bitmap->a;

  // Count output samples
  uint64_t n_output_samples = 0;
  if (sample_set != NULL) {
    for (uint64_t s = 0; s < n_samples; s ++) {
      if (tk_iuset_contains(sample_set, (int64_t) s))
        n_output_samples ++;
    }
  } else {
    n_output_samples = n_samples;
  }

  // Process samples
  uint64_t write_sample = 0;
  for (uint64_t s = 0; s < n_samples; s ++) {
    // Skip if sample not in filter set
    if (sample_set != NULL && !tk_iuset_contains(sample_set, (int64_t) s))
      continue;

    uint8_t *temp = malloc(out_bytes_per_sample);
    if (!temp) {
      if (sample_set != NULL) tk_iuset_destroy(sample_set);
      return;
    }
    memset(temp, 0, out_bytes_per_sample);

    uint64_t in_offset = s * in_bytes_per_sample;

    if (selected_features != NULL && selected_features->n > 0) {
      // Copy selected feature bits
      for (uint64_t i = 0; i < selected_features->n; i ++) {
        int64_t src_bit = selected_features->a[i];
        if (src_bit >= 0 && (uint64_t) src_bit < n_features) {
          uint64_t src_byte = (uint64_t) src_bit / CHAR_BIT;
          uint8_t src_bit_pos = (uint64_t) src_bit % CHAR_BIT;

          if (data[in_offset + src_byte] & (1u << src_bit_pos)) {
            uint64_t dst_byte = i / CHAR_BIT;
            uint8_t dst_bit_pos = i % CHAR_BIT;
            temp[dst_byte] |= (1u << dst_bit_pos);
          }
        }
      }
    } else {
      // No feature filtering, copy all features
      memcpy(temp, data + in_offset, in_bytes_per_sample);
    }

    // Copy result to write position
    uint64_t out_offset = write_sample * out_bytes_per_sample;
    memcpy(data + out_offset, temp, out_bytes_per_sample);
    free(temp);
    write_sample ++;
  }

  bitmap->n = n_output_samples * out_bytes_per_sample;

  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);
}

static inline void tk_ivec_bits_copy (
  tk_ivec_t *dest,
  tk_ivec_t *src_bits,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_visible,
  uint64_t dest_sample  // 0 = overwrite, >0 = append after this many samples
) {
  // Clear destination if dest_sample == 0 (overwrite mode)
  if (dest_sample == 0) {
    tk_ivec_clear(dest);
  }

  // Return if no source bits
  if (src_bits == NULL || src_bits->n == 0)
    return;

  // Create feature map if filtering by features
  int64_t *feature_map = NULL;
  uint64_t n_new_features = n_visible;
  if (selected_features != NULL && selected_features->n > 0) {
    feature_map = malloc(n_visible * sizeof(int64_t));
    if (!feature_map)
      return;
    for (uint64_t i = 0; i < n_visible; i ++)
      feature_map[i] = -1;
    for (uint64_t i = 0; i < selected_features->n; i ++) {
      int64_t feat = selected_features->a[i];
      if (feat >= 0 && (uint64_t) feat < n_visible)
        feature_map[feat] = (int64_t) i;
    }
    n_new_features = selected_features->n;
  }

  // Create sample map if filtering by samples
  tk_iuset_t *sample_set = NULL;
  int64_t *sample_map = NULL;
  uint64_t max_sample = 0;
  if (sample_ids != NULL && sample_ids->n > 0) {
    // Create sample set for quick lookup
    sample_set = tk_iuset_from_ivec(sample_ids);

    // Always create sample remapping when filtering (ivec is always "packed")
    // First find max sample to allocate map
    for (uint64_t i = 0; i < src_bits->n; i ++) {
      if (src_bits->a[i] >= 0) {
        uint64_t s = (uint64_t) src_bits->a[i] / n_visible;
        if (s > max_sample) max_sample = s;
      }
    }

    // Create sample remapping for compacted output
    sample_map = malloc((max_sample + 1) * sizeof(int64_t));
    if (!sample_map) {
      if (feature_map) free(feature_map);
      tk_iuset_destroy(sample_set);
      return;
    }
    for (uint64_t i = 0; i <= max_sample; i ++)
      sample_map[i] = -1;

    // Map old sample indices to new compacted indices
    uint64_t new_idx = 0;
    for (uint64_t i = 0; i < sample_ids->n; i ++) {
      int64_t sid = sample_ids->a[i];
      if (sid >= 0 && (uint64_t) sid <= max_sample)
        sample_map[sid] = (int64_t) new_idx ++;
    }
  }

  // Process bits
  for (size_t i = 0; i < src_bits->n; i ++) {
    int64_t val = src_bits->a[i];
    if (val < 0)
      continue;
    uint64_t sample = (uint64_t) val / n_visible;
    uint64_t feature = (uint64_t) val % n_visible;

    // Check and remap sample if filtering
    int64_t new_sample = (int64_t) sample;
    if (sample_set != NULL) {
      if (!tk_iuset_contains(sample_set, (int64_t) sample))
        continue;
      // Always use sample_map for remapping when filtering
      if (sample <= max_sample) {
        new_sample = sample_map[sample];
        if (new_sample < 0)
          continue;
      } else {
        continue;  // Sample index out of range
      }
    }

    // Check and remap feature if filtering
    int64_t new_feature = (int64_t) feature;
    if (feature_map != NULL) {
      new_feature = feature_map[feature];
      if (new_feature < 0)
        continue;
    }

    // Add dest_sample offset to the final sample index
    new_sample += (int64_t) dest_sample;

    // Add remapped index to destination
    int64_t new_val = new_sample * (int64_t) n_new_features + new_feature;
    tk_ivec_push(dest, new_val);
  }

  if (feature_map != NULL)
    free(feature_map);
  if (sample_map != NULL)
    free(sample_map);
  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);

  tk_ivec_shrink(dest);
}

static inline void tk_cvec_bits_copy (
  tk_cvec_t *dest,
  tk_cvec_t *src_bitmap,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_features,
  uint64_t dest_sample,  // Sample index in destination
  uint64_t dest_stride  // 0 = byte-aligned, >0 = fixed bit stride between rows
) {
  uint64_t n_samples = src_bitmap->n / TK_CVEC_BITS_BYTES(n_features);
  // Create sample set if filtering by samples
  tk_iuset_t *sample_set = NULL;
  if (sample_ids != NULL && sample_ids->n > 0) {
    sample_set = tk_iuset_from_ivec(sample_ids);
  }

  uint64_t n_selected_features = (selected_features != NULL && selected_features->n > 0)
    ? selected_features->n : n_features;
  uint64_t in_bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t out_bytes_per_sample = TK_CVEC_BITS_BYTES(n_selected_features);
  uint8_t *src_data = (uint8_t *)src_bitmap->a;

  // Count output samples
  uint64_t n_output_samples = 0;
  if (sample_set != NULL) {
    for (uint64_t s = 0; s < n_samples; s ++) {
      if (tk_iuset_contains(sample_set, (int64_t) s))
        n_output_samples ++;
    }
  } else {
    n_output_samples = n_samples;
  }

  // Determine stride
  uint64_t row_stride_bits;
  bool use_packed;

  if (dest_stride > 0) {
    // Explicit stride: packed mode with fixed stride
    // Round up to byte boundary for proper alignment
    use_packed = true;
    row_stride_bits = TK_CVEC_BITS_BYTES(dest_stride) * CHAR_BIT;
  } else {
    // Default: byte-aligned mode
    use_packed = false;
    row_stride_bits = out_bytes_per_sample * CHAR_BIT;
  }

  // Calculate total bytes needed for the destination
  uint64_t final_samples = dest_sample + n_output_samples;
  uint64_t total_bytes;

  if (use_packed) {
    // Packed mode: use bit-level stride
    uint64_t total_bits = dest_sample * row_stride_bits + n_output_samples * n_selected_features;
    // If using fixed stride, need space for full final row
    if (dest_stride > 0) {
      total_bits = (dest_sample + 1) * row_stride_bits;
    }
    total_bytes = TK_CVEC_BITS_BYTES(total_bits);
  } else {
    // Unpacked mode: byte-aligned samples
    total_bytes = final_samples * out_bytes_per_sample;
  }

  // Ensure destination capacity
  tk_cvec_ensure(dest, total_bytes);
  uint8_t *dest_data = (uint8_t *)dest->a;

  if (dest_sample == 0) {
    // Overwrite mode: clear everything and set new size
    dest->n = total_bytes;
    memset(dest_data, 0, total_bytes);
  } else {
    // Append mode: extend size and zero only new portion
    uint64_t old_size = dest->n;
    dest->n = total_bytes;
    if (total_bytes > old_size) {
      memset(dest_data + old_size, 0, total_bytes - old_size);
    }
  }

  // Process samples
  uint64_t write_sample = dest_sample;  // Start writing at dest_sample
  uint64_t write_bit_offset = 0;

  if (use_packed) {
    // Calculate starting bit position using stride
    write_bit_offset = dest_sample * row_stride_bits;
  }

  for (uint64_t s = 0; s < n_samples; s ++) {
    // Skip if sample not in filter set
    if (sample_set != NULL && !tk_iuset_contains(sample_set, (int64_t) s))
      continue;

    uint8_t *temp = malloc(out_bytes_per_sample);
    if (!temp) {
      if (sample_set != NULL) tk_iuset_destroy(sample_set);
      return;
    }
    memset(temp, 0, out_bytes_per_sample);

    uint64_t in_offset = s * in_bytes_per_sample;

    if (selected_features != NULL && selected_features->n > 0) {
      // Copy selected feature bits
      for (uint64_t i = 0; i < selected_features->n; i ++) {
        int64_t src_bit = selected_features->a[i];
        if (src_bit >= 0 && (uint64_t) src_bit < n_features) {
          uint64_t src_byte = (uint64_t) src_bit / CHAR_BIT;
          uint8_t src_bit_pos = (uint64_t) src_bit % CHAR_BIT;

          if (src_data[in_offset + src_byte] & (1u << src_bit_pos)) {
            uint64_t dst_byte = i / CHAR_BIT;
            uint8_t dst_bit_pos = i % CHAR_BIT;
            temp[dst_byte] |= (1u << dst_bit_pos);
          }
        }
      }
    } else {
      // No feature filtering, copy all features
      memcpy(temp, src_data + in_offset, in_bytes_per_sample);
    }

    // Copy result to destination
    if (use_packed) {
      // Packed mode: write bits at current offset
      for (uint64_t i = 0; i < n_selected_features; i++) {
        uint64_t src_byte = i / CHAR_BIT;
        uint8_t src_bit_pos = i % CHAR_BIT;

        if (temp[src_byte] & (1u << src_bit_pos)) {
          uint64_t dst_byte = write_bit_offset / CHAR_BIT;
          uint8_t dst_bit_pos = write_bit_offset % CHAR_BIT;
          dest_data[dst_byte] |= (1u << dst_bit_pos);
        }
        write_bit_offset++;
      }
      // write_bit_offset already advanced by n_selected_features bits for this sample
      // Continue with next sample in the same row (for concatenation)
    } else {
      // In unpacked mode, maintain byte alignment per sample
      uint64_t out_offset = write_sample * out_bytes_per_sample;
      memcpy(dest_data + out_offset, temp, out_bytes_per_sample);
      write_sample++;
    }

    free(temp);
  }

  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);
}

// Additional ivec functions that depend on iuset

static inline tk_ivec_t *tk_iuset_keys (lua_State *L, tk_iuset_t *S)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iuset_size(S), 0, 0);
  int64_t k;
  out->n = 0;
  tk_iuset_foreach(S, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

static inline tk_ivec_t *tk_ivec_from_iuset (lua_State *L, tk_iuset_t *s)
{
  tk_ivec_t *v = tk_ivec_create(L, tk_iuset_size(s), 0, 0);
  int64_t x;
  v->n = 0;
  tk_iuset_foreach(s, x, ({
    v->a[v->n ++] = x;
  }))
  return v;
}

static inline tk_ivec_t *tk_ivec_top_select (
  lua_State *L,
  tk_iuset_t *selected
) {
  // Create matrix
  tk_ivec_t *top_v = tk_ivec_create(L, 0, 0, 0);
  int64_t sel;
  tk_iuset_foreach(selected, sel, ({
    tk_ivec_push(top_v, sel);
  }));
  tk_ivec_shrink(top_v);
  tk_iuset_destroy(selected);
  return top_v;
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
    for (uint64_t j = 0; j < n_hidden; j ++) {
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

static inline tk_ivec_t *tk_ivec_bits_top_mi (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k
) {
  // Sort set_bits for efficient processing
  tk_ivec_asc(set_bits, 0, set_bits->n);

  // Stack: (empty)

  // Branch based on dense codes vs sparse labels
  if (codes) {
    // DENSE PATH: Use existing implementation with dense counting array
    tk_ivec_t *counts = tk_ivec_create(L, n_visible * n_hidden * 4, 0, 0);
    tk_ivec_zero(counts);
    // Stack: counts

    // Process sparse bits to count all 4 combinations
    uint64_t prev_sample = UINT64_MAX;
    uint8_t *sample_codes = NULL;

    // First pass: count features
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit_idx = set_bits->a[i];
      if (bit_idx < 0) continue;

      uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
      uint64_t feature_idx = (uint64_t)bit_idx % n_visible;

      if (sample_idx >= n_samples || feature_idx >= n_visible) continue;

      // Get codes for this sample if it's new
      if (sample_idx != prev_sample) {
        prev_sample = sample_idx;
        sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
      }

      // Count combinations for MI calculation
      for (uint64_t j = 0; j < n_hidden; j++) {
        uint64_t byte_idx = j / CHAR_BIT;
        uint64_t bit_pos = j % CHAR_BIT;
        bool hidden = (sample_codes[byte_idx] & (1u << bit_pos)) != 0;
        // visible=1 (feature is present), hidden=0/1
        counts->a[feature_idx * n_hidden * 4 + j * 4 + 2 + (hidden ? 1 : 0)]++;
      }
    }

    // Second pass: count non-visible cases
    for (uint64_t s = 0; s < n_samples; s++) {
      // Build set of features present in this sample
      tk_iuset_t *sample_features = tk_iuset_create();
      for (uint64_t i = 0; i < set_bits->n; i++) {
        int64_t bit_idx = set_bits->a[i];
        if (bit_idx < 0) continue;
        uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
        uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
        if (sample_idx == s) {
          int absent;
          tk_iuset_put(sample_features, (int64_t)feature_idx, &absent);
        }
      }

      // Count non-visible cases
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t f = 0; f < n_visible; f++) {
        if (!tk_iuset_contains(sample_features, (int64_t)f)) {
          for (uint64_t j = 0; j < n_hidden; j++) {
            uint64_t byte_idx = j / CHAR_BIT;
            uint64_t bit_pos = j % CHAR_BIT;
            bool hidden = (s_codes[byte_idx] & (1u << bit_pos)) != 0;
            // visible=0 (feature not present), hidden=0/1
            counts->a[f * n_hidden * 4 + j * 4 + (hidden ? 1 : 0)]++;
          }
        }
      }

      tk_iuset_destroy(sample_features);
    }

    // Calculate MI scores using fixed-size heap for top-k
    // Stack: counts
    tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
    // Stack: counts, top_heap

    for (uint64_t f = 0; f < n_visible; f++) {
      double max_mi = 0.0;

      for (uint64_t j = 0; j < n_hidden; j++) {
        int64_t c[4];
        int64_t *counts_ptr = counts->a + f * n_hidden * 4 + j * 4;
        for (int k = 0; k < 4; k++) {
          c[k] = counts_ptr[k] + 1; // Add 1 for smoothing
        }

        double total = c[0] + c[1] + c[2] + c[3];
        double mi = 0.0;

        if (total > 0.0) {
          for (unsigned int o = 0; o < 4; o++) {
            if (c[o] == 0) continue;
            double p_fb = c[o] / total;
            unsigned int feat = o >> 1;
            unsigned int hid = o & 1;
            double pf = (c[2] + c[3]) / total;
            if (feat == 0) pf = 1.0 - pf;
            double ph = (c[1] + c[3]) / total;
            if (hid == 0) ph = 1.0 - ph;
            double d = pf * ph;
            if (d > 0) mi += p_fb * log2(p_fb / d);
          }
        }

        if (mi > max_mi) max_mi = mi;
      }

      // Use heap to maintain top-k features
      tk_rank_t r = { (int64_t)f, max_mi };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    // Sort the final heap in descending order (highest scores first)
    tk_rvec_desc(top_heap, 0, top_heap->n);

    // Extract results (scores are already positive)
    // Stack: counts, top_heap
    tk_ivec_t *out = tk_ivec_create(L, top_heap->n, 0, 0);
    // Stack: counts, top_heap, out

    tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
    // Stack: counts, top_heap, out, weights

    for (uint64_t i = 0; i < top_heap->n; i++) {
      out->a[i] = top_heap->a[i].i;
      weights->a[i] = top_heap->a[i].d;
    }

    // Rearrange stack to return: out, weights
    lua_remove(L, -4); // Remove counts
    lua_remove(L, -3); // Remove top_heap

    return out;

  } else if (labels) {
    // SPARSE PATH: Interpret labels as sparse multi-label bits
    // Labels format: sample_idx * n_hidden + label_idx
    tk_ivec_asc(labels, 0, labels->n);

    // Use sparse counting with hashmap
    tk_iumap_t *counts = tk_iumap_create();

    // Track marginals for total-based calculation
    tk_ivec_t *feat_counts = tk_ivec_create(L, n_visible, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_t *label_counts = tk_ivec_create(L, n_hidden, 0, 0);
    tk_ivec_zero(label_counts);
    // Stack: feat_counts, label_counts

    // Count marginals from the sparse vectors
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0) {
        uint64_t f = (uint64_t)bit % n_visible;
        feat_counts->a[f]++;
      }
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        label_counts->a[h]++;
      }
    }

    // Two-pointer traversal for counting (v=1, h=1) combinations
    size_t si = 0, li = 0;
    while (si < set_bits->n) {
      if (set_bits->a[si] < 0) {
        si++;
        continue;
      }

      uint64_t s_sample = (uint64_t)set_bits->a[si] / n_visible;
      uint64_t f = (uint64_t)set_bits->a[si] % n_visible;

      // Skip labels until we reach this sample
      while (li < labels->n && labels->a[li] >= 0 &&
             (uint64_t)labels->a[li] / n_hidden < s_sample) {
        li++;
      }

      // If no labels for this sample, move to next feature
      if (li >= labels->n || labels->a[li] < 0 ||
          (uint64_t)labels->a[li] / n_hidden > s_sample) {
        si++;
        continue;
      }

      // Remember where this sample's labels start
      size_t li_start = li;

      // Count this feature with all labels of this sample
      while (li < labels->n && labels->a[li] >= 0 &&
             (uint64_t)labels->a[li] / n_hidden == s_sample) {
        uint64_t h = (uint64_t)labels->a[li] % n_hidden;
        int64_t key = (int64_t)(f * n_hidden * 4 + h * 4 + 3); // 11 combination
        tk_iumap_inc(counts, key);
        li++;
      }

      // Move to next feature
      si++;

      // If next feature is in same sample, reset label pointer
      if (si < set_bits->n && set_bits->a[si] >= 0 &&
          (uint64_t)set_bits->a[si] / n_visible == s_sample) {
        li = li_start;
      }
    }

    // Calculate MI scores using total-based calculation for absent combinations
    tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
    // Stack: feat_counts, label_counts, top_heap

    for (uint64_t f = 0; f < n_visible; f++) {
      double max_mi = 0.0;
      int64_t feat_total = feat_counts->a[f];

      for (uint64_t h = 0; h < n_hidden; h++) {
        int64_t label_total = label_counts->a[h];

        // Get count_11 from sparse map
        int64_t key_11 = (int64_t)(f * n_hidden * 4 + h * 4 + 3);
        int64_t count_11 = tk_iumap_get_or(counts, key_11, 0);

        // Calculate other combinations using totals
        int64_t count_10 = feat_total - count_11;      // f=1, h=0
        int64_t count_01 = label_total - count_11;     // f=0, h=1
        int64_t count_00 = (int64_t)n_samples - count_11 - count_10 - count_01; // f=0, h=0

        // Add smoothing
        int64_t c[4];
        c[0] = count_00 + 1;
        c[1] = count_01 + 1;
        c[2] = count_10 + 1;
        c[3] = count_11 + 1;

        double total = c[0] + c[1] + c[2] + c[3];
        double mi = 0.0;

        if (total > 0.0) {
          for (unsigned int o = 0; o < 4; o++) {
            if (c[o] == 0) continue;
            double p_fb = c[o] / total;
            unsigned int feat = o >> 1;
            unsigned int hid = o & 1;
            double pf = (c[2] + c[3]) / total;
            if (feat == 0) pf = 1.0 - pf;
            double ph = (c[1] + c[3]) / total;
            if (hid == 0) ph = 1.0 - ph;
            double d = pf * ph;
            if (d > 0) mi += p_fb * log2(p_fb / d);
          }
        }

        if (mi > max_mi) max_mi = mi;
      }

      // Use heap to maintain top-k features
      tk_rank_t r = { (int64_t)f, max_mi };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    // Clean up hashmap
    tk_iumap_destroy(counts);

    // Sort the final heap in descending order (highest scores first)
    tk_rvec_desc(top_heap, 0, top_heap->n);

    // Extract results
    // Stack: feat_counts, label_counts, top_heap
    tk_ivec_t *out = tk_ivec_create(L, top_heap->n, 0, 0);
    // Stack: feat_counts, label_counts, top_heap, out

    tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
    // Stack: feat_counts, label_counts, top_heap, out, weights

    for (uint64_t i = 0; i < top_heap->n; i++) {
      out->a[i] = top_heap->a[i].i;
      weights->a[i] = top_heap->a[i].d;
    }

    // Rearrange stack to return: out, weights
    lua_remove(L, -5); // Remove feat_counts
    lua_remove(L, -4); // Remove label_counts
    lua_remove(L, -3); // Remove top_heap

    return out;

  } else {
    // No codes or labels provided - return empty
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }
}

static inline tk_ivec_t *tk_ivec_bits_top_chi2 (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k
) {
  // Sort set_bits for efficient processing
  tk_ivec_asc(set_bits, 0, set_bits->n);

  // Stack: (empty)

  // Branch based on dense codes vs sparse labels
  if (codes) {
    // DENSE PATH: Use existing implementation with dense counting array
    tk_ivec_t *active_counts = tk_ivec_create(L, n_visible * n_hidden, 0, 0);
    tk_ivec_zero(active_counts);
    // Stack: active_counts

    tk_ivec_t *global_counts = tk_ivec_create(L, n_hidden, 0, 0);
    tk_ivec_zero(global_counts);
    // Stack: active_counts, global_counts

    tk_ivec_t *feat_counts = tk_ivec_create(L, n_visible, 0, 0);
    tk_ivec_zero(feat_counts);
    // Stack: active_counts, global_counts, feat_counts

    // Count features from sparse representation
    uint64_t prev_sample = UINT64_MAX;
    uint8_t *sample_codes = NULL;

    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit_idx = set_bits->a[i];
      if (bit_idx < 0) continue;

      uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
      uint64_t feature_idx = (uint64_t)bit_idx % n_visible;

      if (sample_idx >= n_samples || feature_idx >= n_visible) continue;

      // Track feature counts
      feat_counts->a[feature_idx]++;

      // Get codes for this sample if it's new
      if (sample_idx != prev_sample) {
        prev_sample = sample_idx;
        sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
      }

      // Count actives for chi2 calculation
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (sample_codes[byte_idx] & (1u << bit_pos)) {
          active_counts->a[feature_idx * n_hidden + b]++;
        }
      }
    }

    // Count global occurrences
    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos)) {
          global_counts->a[b]++;
        }
      }
    }

    // Calculate chi2 scores using fixed-size heap for top-k
    // Stack: active_counts, global_counts, feat_counts
    tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
    // Stack: active_counts, global_counts, feat_counts, top_heap

    for (uint64_t f = 0; f < n_visible; f++) {
      double max_chi2 = 0.0;

      for (uint64_t b = 0; b < n_hidden; b++) {
        int64_t A = active_counts->a[f * n_hidden + b];
        int64_t G = global_counts->a[b];
        int64_t C = feat_counts->a[f];

        if (C == 0 || G == 0 || C == (int64_t)n_samples || G == (int64_t)n_samples) {
          continue;
        }

        int64_t B = G - A; // f=0, b=1
        int64_t C_ = C - A; // f=1, b=0
        int64_t D = (int64_t)n_samples - C - B; // f=0, b=0

        double n = (double)n_samples;
        double E_A = ((double)C * (double)G) / n;
        double E_B = ((double)(n - C) * (double)G) / n;
        double E_C = ((double)C * (double)(n - G)) / n;
        double E_D = ((double)(n - C) * (double)(n - G)) / n;

        double chi2 = 0.0;
        if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
        if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
        if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
        if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;

        if (chi2 > max_chi2) max_chi2 = chi2;
      }

      // Use heap to maintain top-k features
      tk_rank_t r = { (int64_t)f, max_chi2 };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    // Sort the final heap in descending order (highest scores first)
    tk_rvec_desc(top_heap, 0, top_heap->n);

    // Extract results
    // Stack: active_counts, global_counts, feat_counts, top_heap
    tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);
    // Stack: active_counts, global_counts, feat_counts, top_heap, out

    tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
    // Stack: active_counts, global_counts, feat_counts, top_heap, out, weights

    for (uint64_t i = 0; i < top_heap->n; i++) {
      weights->a[i] = top_heap->a[i].d;
    }

    // Rearrange stack to return: out, weights
    lua_remove(L, -6); // Remove active_counts
    lua_remove(L, -5); // Remove global_counts
    lua_remove(L, -4); // Remove feat_counts
    lua_remove(L, -3); // Remove top_heap

    return out;

  } else if (labels) {
    // SPARSE PATH: Interpret labels as sparse multi-label bits
    // Labels format: sample_idx * n_hidden + label_idx
    tk_ivec_asc(labels, 0, labels->n);

    // Use sparse counting with hashmap for active counts
    tk_iumap_t *active_counts = tk_iumap_create();

    // Track marginals
    tk_ivec_t *feat_counts = tk_ivec_create(L, n_visible, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_t *global_counts = tk_ivec_create(L, n_hidden, 0, 0);
    tk_ivec_zero(global_counts);
    // Stack: feat_counts, global_counts

    // Count marginals from the sparse vectors
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0) {
        uint64_t f = (uint64_t)bit % n_visible;
        feat_counts->a[f]++;
      }
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        global_counts->a[h]++;
      }
    }

    // Two-pointer traversal for counting (f=1, h=1) combinations
    size_t si = 0, li = 0;
    while (si < set_bits->n) {
      if (set_bits->a[si] < 0) {
        si++;
        continue;
      }

      uint64_t s_sample = (uint64_t)set_bits->a[si] / n_visible;
      uint64_t f = (uint64_t)set_bits->a[si] % n_visible;

      // Skip labels until we reach this sample
      while (li < labels->n && labels->a[li] >= 0 &&
             (uint64_t)labels->a[li] / n_hidden < s_sample) {
        li++;
      }

      // If no labels for this sample, move to next feature
      if (li >= labels->n || labels->a[li] < 0 ||
          (uint64_t)labels->a[li] / n_hidden > s_sample) {
        si++;
        continue;
      }

      // Remember where this sample's labels start
      size_t li_start = li;

      // Count this feature with all labels of this sample
      while (li < labels->n && labels->a[li] >= 0 &&
             (uint64_t)labels->a[li] / n_hidden == s_sample) {
        uint64_t h = (uint64_t)labels->a[li] % n_hidden;
        int64_t key = (int64_t)(f * n_hidden + h);
        tk_iumap_inc(active_counts, key);
        li++;
      }

      // Move to next feature
      si++;

      // If next feature is in same sample, reset label pointer
      if (si < set_bits->n && set_bits->a[si] >= 0 &&
          (uint64_t)set_bits->a[si] / n_visible == s_sample) {
        li = li_start;
      }
    }

    // Calculate chi2 scores using sparse data
    tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
    // Stack: feat_counts, global_counts, top_heap

    for (uint64_t f = 0; f < n_visible; f++) {
      double max_chi2 = 0.0;
      int64_t C = feat_counts->a[f];

      for (uint64_t b = 0; b < n_hidden; b++) {
        int64_t G = global_counts->a[b];

        if (C == 0 || G == 0 || C == (int64_t)n_samples || G == (int64_t)n_samples) {
          continue;
        }

        // Get active count from sparse map
        int64_t key = (int64_t)(f * n_hidden + b);
        int64_t A = tk_iumap_get_or(active_counts, key, 0);

        int64_t B = G - A; // f=0, b=1
        int64_t C_ = C - A; // f=1, b=0
        int64_t D = (int64_t)n_samples - C - B; // f=0, b=0

        double n = (double)n_samples;
        double E_A = ((double)C * (double)G) / n;
        double E_B = ((double)(n - C) * (double)G) / n;
        double E_C = ((double)C * (double)(n - G)) / n;
        double E_D = ((double)(n - C) * (double)(n - G)) / n;

        double chi2 = 0.0;
        if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
        if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
        if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
        if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;

        if (chi2 > max_chi2) max_chi2 = chi2;
      }

      // Use heap to maintain top-k features
      tk_rank_t r = { (int64_t)f, max_chi2 };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    // Clean up hashmap
    tk_iumap_destroy(active_counts);

    // Sort the final heap in descending order (highest scores first)
    tk_rvec_desc(top_heap, 0, top_heap->n);

    // Extract results
    // Stack: feat_counts, global_counts, top_heap
    tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);
    // Stack: feat_counts, global_counts, top_heap, out

    tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
    // Stack: feat_counts, global_counts, top_heap, out, weights

    for (uint64_t i = 0; i < top_heap->n; i++) {
      weights->a[i] = top_heap->a[i].d;
    }

    // Rearrange stack to return: out, weights
    lua_remove(L, -5); // Remove feat_counts
    lua_remove(L, -4); // Remove global_counts
    lua_remove(L, -3); // Remove top_heap

    return out;

  } else {
    // No codes or labels provided - return empty
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }
}

static inline tk_ivec_t *tk_ivec_bits_top_entropy (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k
) {
  // Sort set_bits for efficient processing
  tk_ivec_asc(set_bits, 0, set_bits->n);

  // Stack: (empty)

  // Count bits per hidden unit
  tk_ivec_t *bit_counts = tk_ivec_create(L, n_hidden, 0, 0);
  tk_ivec_zero(bit_counts);
  // Stack: bit_counts

  // Process sparse bits to count occurrences
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    if (bit_idx < 0) continue;

    uint64_t sample_idx = (uint64_t)bit_idx / n_hidden;
    uint64_t hidden_idx = (uint64_t)bit_idx % n_hidden;

    if (sample_idx >= n_samples || hidden_idx >= n_hidden) continue;

    bit_counts->a[hidden_idx]++;
  }

  // Calculate entropy scores using fixed-size heap for top-k
  // Stack: bit_counts
  tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
  // Stack: bit_counts, top_heap

  for (uint64_t h = 0; h < n_hidden; h++) {
    double p = (double)bit_counts->a[h] / (double)n_samples;
    double entropy = 0.0;

    if (p > 0.0 && p < 1.0) {
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    }

    // Use heap to maintain top-k hidden units (highest entropy)
    tk_rank_t r = { (int64_t)h, entropy };
    tk_rvec_hmin(top_heap, top_k, r);
  }

  // Sort the final heap in descending order (highest scores first)
  tk_rvec_desc(top_heap, 0, top_heap->n);

  // Extract results
  // Stack: bit_counts, top_heap
  tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);
  // Stack: bit_counts, top_heap, out

  tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
  // Stack: bit_counts, top_heap, out, weights

  for (uint64_t i = 0; i < top_heap->n; i++) {
    weights->a[i] = top_heap->a[i].d; // Scores are already positive
  }

  // Rearrange stack to return: out, weights
  lua_remove(L, -4); // Remove bit_counts
  lua_remove(L, -3); // Remove top_heap

  return out;
}

// cvec top functions for dense bitmaps
static inline tk_ivec_t *tk_cvec_bits_top_mi (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t top_k
) {
  // Stack: (empty)

  // Branch based on dense codes vs sparse labels
  if (codes) {
    // DENSE PATH: Use existing implementation with dense counting array
    tk_ivec_t *counts = tk_ivec_create(L, n_features * n_hidden * 4, 0, 0);
    tk_ivec_zero(counts);
    // Stack: counts

    // Process dense bitmap
    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

    // Count all 4 combinations
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));

      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        bool visible = (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) != 0;

        for (uint64_t j = 0; j < n_hidden; j++) {
          uint64_t j_byte = j / CHAR_BIT;
          uint8_t j_bit = j % CHAR_BIT;
          bool hidden = (sample_codes[j_byte] & (1u << j_bit)) != 0;
          counts->a[f * n_hidden * 4 + j * 4 + (visible ? 2 : 0) + (hidden ? 1 : 0)]++;
        }
      }
    }

    // Calculate MI scores and build heap
    // Stack: counts
    tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
    // Stack: counts, top_heap

    for (uint64_t f = 0; f < n_features; f++) {
      double max_mi = 0.0;

      for (uint64_t j = 0; j < n_hidden; j++) {
        int64_t c[4];
        int64_t *counts_ptr = counts->a + f * n_hidden * 4 + j * 4;
        for (int k = 0; k < 4; k++) {
          c[k] = counts_ptr[k] + 1; // Add 1 for smoothing
        }

        double total = c[0] + c[1] + c[2] + c[3];
        double mi = 0.0;

        if (total > 0.0) {
          for (unsigned int o = 0; o < 4; o++) {
            if (c[o] == 0) continue;
            double p_fb = c[o] / total;
            unsigned int feat = o >> 1;
            unsigned int hid = o & 1;
            double pf = (c[2] + c[3]) / total;
            if (feat == 0) pf = 1.0 - pf;
            double ph = (c[1] + c[3]) / total;
            if (hid == 0) ph = 1.0 - ph;
            double d = pf * ph;
            if (d > 0) mi += p_fb * log2(p_fb / d);
          }
        }

        if (mi > max_mi) max_mi = mi;
      }

      // Use heap to maintain top-k features (highest MI)
      tk_rank_t r = { (int64_t)f, max_mi };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    // Sort the final heap in descending order (highest scores first)
    tk_rvec_desc(top_heap, 0, top_heap->n);

    // Extract results
    // Stack: counts, top_heap
    tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);
    // Stack: counts, top_heap, out

    tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
    // Stack: counts, top_heap, out, weights

    for (uint64_t i = 0; i < top_heap->n; i++) {
      weights->a[i] = top_heap->a[i].d;
    }

    // Rearrange stack to return: out, weights
    lua_remove(L, -4); // Remove counts
    lua_remove(L, -3); // Remove top_heap

    return out;

  } else if (labels) {
    // SPARSE PATH: Interpret labels as sparse multi-label bits
    // Labels format: sample_idx * n_hidden + label_idx
    tk_ivec_asc(labels, 0, labels->n);

    // Use sparse counting with hashmap
    tk_iumap_t *counts = tk_iumap_create();

    // Track marginals for total-based calculation
    tk_ivec_t *feat_counts = tk_ivec_create(L, n_features, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_t *label_counts = tk_ivec_create(L, n_hidden, 0, 0);
    tk_ivec_zero(label_counts);
    // Stack: feat_counts, label_counts

    // Count feature marginals from dense bitmap
    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          feat_counts->a[f]++;
        }
      }
    }

    // Count label marginals from sparse labels
    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        label_counts->a[h]++;
      }
    }

    // Count (f=1, h=1) combinations using binary search
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;

      // Find labels for this sample using binary search
      int64_t label_start = tk_ivec_set_find(
        labels->a, 0, (int64_t) labels->n, (int64_t)(s * n_hidden)
      );

      if (label_start < 0) {
        label_start = -(label_start + 1);
      }

      // Process each set feature in this sample
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;

        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx))) {
          continue; // Feature not present
        }

        // Count combinations with all labels of this sample
        for (int64_t li = label_start;
             li < (int64_t) labels->n && labels->a[li] >= 0 &&
             (uint64_t)labels->a[li] / n_hidden == s;
             li++) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden * 4 + h * 4 + 3); // 11 combination
          tk_iumap_inc(counts, key);
        }
      }
    }

    // Calculate MI scores using total-based calculation for absent combinations
    tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
    // Stack: feat_counts, label_counts, top_heap

    for (uint64_t f = 0; f < n_features; f++) {
      double max_mi = 0.0;
      int64_t feat_total = feat_counts->a[f];

      for (uint64_t h = 0; h < n_hidden; h++) {
        int64_t label_total = label_counts->a[h];

        // Get count_11 from sparse map
        int64_t key_11 = (int64_t)(f * n_hidden * 4 + h * 4 + 3);
        int64_t count_11 = tk_iumap_get_or(counts, key_11, 0);

        // Calculate other combinations using totals
        int64_t count_10 = feat_total - count_11;      // f=1, h=0
        int64_t count_01 = label_total - count_11;     // f=0, h=1
        int64_t count_00 = (int64_t)n_samples - count_11 - count_10 - count_01; // f=0, h=0

        // Add smoothing
        int64_t c[4];
        c[0] = count_00 + 1;
        c[1] = count_01 + 1;
        c[2] = count_10 + 1;
        c[3] = count_11 + 1;

        double total = c[0] + c[1] + c[2] + c[3];
        double mi = 0.0;

        if (total > 0.0) {
          for (unsigned int o = 0; o < 4; o++) {
            if (c[o] == 0) continue;
            double p_fb = c[o] / total;
            unsigned int feat = o >> 1;
            unsigned int hid = o & 1;
            double pf = (c[2] + c[3]) / total;
            if (feat == 0) pf = 1.0 - pf;
            double ph = (c[1] + c[3]) / total;
            if (hid == 0) ph = 1.0 - ph;
            double d = pf * ph;
            if (d > 0) mi += p_fb * log2(p_fb / d);
          }
        }

        if (mi > max_mi) max_mi = mi;
      }

      // Use heap to maintain top-k features
      tk_rank_t r = { (int64_t)f, max_mi };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    // Clean up hashmap
    tk_iumap_destroy(counts);

    // Sort the final heap in descending order (highest scores first)
    tk_rvec_desc(top_heap, 0, top_heap->n);

    // Extract results
    // Stack: feat_counts, label_counts, top_heap
    tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);
    // Stack: feat_counts, label_counts, top_heap, out

    tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
    // Stack: feat_counts, label_counts, top_heap, out, weights

    for (uint64_t i = 0; i < top_heap->n; i++) {
      out->a[i] = top_heap->a[i].i;
      weights->a[i] = top_heap->a[i].d;
    }

    // Rearrange stack to return: out, weights
    lua_remove(L, -5); // Remove feat_counts
    lua_remove(L, -4); // Remove label_counts
    lua_remove(L, -3); // Remove top_heap

    return out;

  } else {
    // No codes or labels provided - return empty
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }
}

static inline tk_ivec_t *tk_cvec_bits_top_chi2 (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t top_k
) {
  // Stack: (empty)

  // Branch based on dense codes vs sparse labels
  if (codes) {
    // DENSE PATH: Use existing implementation with dense counting array
    tk_ivec_t *active_counts = tk_ivec_create(L, n_features * n_hidden, 0, 0);
    tk_ivec_zero(active_counts);
    // Stack: active_counts

    tk_ivec_t *global_counts = tk_ivec_create(L, n_hidden, 0, 0);
    tk_ivec_zero(global_counts);
    // Stack: active_counts, global_counts

    tk_ivec_t *feat_counts = tk_ivec_create(L, n_features, 0, 0);
    tk_ivec_zero(feat_counts);
    // Stack: active_counts, global_counts, feat_counts

    // Process dense bitmap
    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

    // Count features and build mappings
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));

      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;

        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          feat_counts->a[f]++;

          // Count actives for chi2 calculation
          for (uint64_t b = 0; b < n_hidden; b++) {
            uint64_t b_byte = b / CHAR_BIT;
            uint8_t b_bit = b % CHAR_BIT;
            if (sample_codes[b_byte] & (1u << b_bit)) {
              active_counts->a[f * n_hidden + b]++;
            }
          }
        }
      }
    }

    // Count global occurrences
    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint8_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos)) {
          global_counts->a[b]++;
        }
      }
    }

    // Calculate chi2 scores and build heap
    // Stack: active_counts, global_counts, feat_counts
    tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
    // Stack: active_counts, global_counts, feat_counts, top_heap

    for (uint64_t f = 0; f < n_features; f++) {
      double max_chi2 = 0.0;

      for (uint64_t b = 0; b < n_hidden; b++) {
        int64_t A = active_counts->a[f * n_hidden + b];
        int64_t G = global_counts->a[b];
        int64_t C = feat_counts->a[f];

        if (C == 0 || G == 0 || C == (int64_t)n_samples || G == (int64_t)n_samples) {
          continue;
        }

        int64_t B = G - A; // f=0, b=1
        int64_t C_ = C - A; // f=1, b=0
        int64_t D = (int64_t)n_samples - C - B; // f=0, b=0

        double n = (double)n_samples;
        double E_A = ((double)C * (double)G) / n;
        double E_B = ((double)(n - C) * (double)G) / n;
        double E_C = ((double)C * (double)(n - G)) / n;
        double E_D = ((double)(n - C) * (double)(n - G)) / n;

        double chi2 = 0.0;
        if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
        if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
        if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
        if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;

        if (chi2 > max_chi2) max_chi2 = chi2;
      }

      // Use heap to maintain top-k features
      tk_rank_t r = { (int64_t)f, max_chi2 };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    // Sort the final heap in descending order (highest scores first)
    tk_rvec_desc(top_heap, 0, top_heap->n);

    // Extract results
    // Stack: active_counts, global_counts, feat_counts, top_heap
    tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);
    // Stack: active_counts, global_counts, feat_counts, top_heap, out

    tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
    // Stack: active_counts, global_counts, feat_counts, top_heap, out, weights

    for (uint64_t i = 0; i < top_heap->n; i++) {
      weights->a[i] = top_heap->a[i].d;
    }

    // Rearrange stack to return: out, weights
    lua_remove(L, -6); // Remove active_counts
    lua_remove(L, -5); // Remove global_counts
    lua_remove(L, -4); // Remove feat_counts
    lua_remove(L, -3); // Remove top_heap

    return out;

  } else if (labels) {
    // SPARSE PATH: Interpret labels as sparse multi-label bits
    // Labels format: sample_idx * n_hidden + label_idx
    tk_ivec_asc(labels, 0, labels->n);

    // Use sparse counting with hashmap for active counts
    tk_iumap_t *active_counts = tk_iumap_create();

    // Track marginals
    tk_ivec_t *feat_counts = tk_ivec_create(L, n_features, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_t *global_counts = tk_ivec_create(L, n_hidden, 0, 0);
    tk_ivec_zero(global_counts);
    // Stack: feat_counts, global_counts

    // Count feature marginals from dense bitmap
    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          feat_counts->a[f]++;
        }
      }
    }

    // Count label marginals from sparse labels
    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        global_counts->a[h]++;
      }
    }

    // Count (f=1, h=1) combinations using binary search
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;

      // Find labels for this sample using binary search
      int64_t label_start = tk_ivec_set_find(
        labels->a, 0, (int64_t) labels->n, (int64_t)(s * n_hidden)
      );

      if (label_start < 0) {
        label_start = -(label_start + 1);
      }

      // Process each set feature in this sample
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;

        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx))) {
          continue; // Feature not present
        }

        // Count combinations with all labels of this sample
        for (int64_t li = label_start;
             li < (int64_t) labels->n && labels->a[li] >= 0 &&
             (uint64_t)labels->a[li] / n_hidden == s;
             li++) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden + h);
          tk_iumap_inc(active_counts, key);
        }
      }
    }

    // Calculate chi2 scores using sparse data
    tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
    // Stack: feat_counts, global_counts, top_heap

    for (uint64_t f = 0; f < n_features; f++) {
      double max_chi2 = 0.0;
      int64_t C = feat_counts->a[f];

      for (uint64_t b = 0; b < n_hidden; b++) {
        int64_t G = global_counts->a[b];

        if (C == 0 || G == 0 || C == (int64_t)n_samples || G == (int64_t)n_samples) {
          continue;
        }

        // Get active count from sparse map
        int64_t key = (int64_t)(f * n_hidden + b);
        int64_t A = tk_iumap_get_or(active_counts, key, 0);

        int64_t B = G - A; // f=0, b=1
        int64_t C_ = C - A; // f=1, b=0
        int64_t D = (int64_t)n_samples - C - B; // f=0, b=0

        double n = (double)n_samples;
        double E_A = ((double)C * (double)G) / n;
        double E_B = ((double)(n - C) * (double)G) / n;
        double E_C = ((double)C * (double)(n - G)) / n;
        double E_D = ((double)(n - C) * (double)(n - G)) / n;

        double chi2 = 0.0;
        if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
        if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
        if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
        if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;

        if (chi2 > max_chi2) max_chi2 = chi2;
      }

      // Use heap to maintain top-k features
      tk_rank_t r = { (int64_t)f, max_chi2 };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    // Clean up hashmap
    tk_iumap_destroy(active_counts);

    // Sort the final heap in descending order (highest scores first)
    tk_rvec_desc(top_heap, 0, top_heap->n);

    // Extract results
    // Stack: feat_counts, global_counts, top_heap
    tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);
    // Stack: feat_counts, global_counts, top_heap, out

    tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
    // Stack: feat_counts, global_counts, top_heap, out, weights

    for (uint64_t i = 0; i < top_heap->n; i++) {
      weights->a[i] = top_heap->a[i].d;
    }

    // Rearrange stack to return: out, weights
    lua_remove(L, -5); // Remove feat_counts
    lua_remove(L, -4); // Remove global_counts
    lua_remove(L, -3); // Remove top_heap

    return out;

  } else {
    // No codes or labels provided - return empty
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }
}

static inline tk_ivec_t *tk_cvec_bits_top_entropy (
  lua_State *L,
  tk_cvec_t *codes,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k
) {
  // Stack: (empty)

  // Count bits per hidden unit
  tk_ivec_t *bit_counts = tk_ivec_create(L, n_hidden, 0, 0);
  tk_ivec_zero(bit_counts);
  // Stack: bit_counts

  // Process codes to count occurrences
  for (uint64_t s = 0; s < n_samples; s++) {
    uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));

    for (uint64_t h = 0; h < n_hidden; h++) {
      uint64_t byte_idx = h / CHAR_BIT;
      uint8_t bit_idx = h % CHAR_BIT;

      if (sample_codes[byte_idx] & (1u << bit_idx)) {
        bit_counts->a[h]++;
      }
    }
  }

  // Calculate entropy scores using fixed-size heap for top-k
  // Stack: bit_counts
  tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
  // Stack: bit_counts, top_heap

  for (uint64_t h = 0; h < n_hidden; h++) {
    double p = (double)bit_counts->a[h] / (double)n_samples;
    double entropy = 0.0;

    if (p > 0.0 && p < 1.0) {
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    }

    // Use heap to maintain top-k hidden units (highest entropy)
    tk_rank_t r = { (int64_t)h, entropy };
    tk_rvec_hmin(top_heap, top_k, r);
  }

  // Sort the final heap in descending order (highest scores first)
  tk_rvec_desc(top_heap, 0, top_heap->n);

  // Extract results
  // Stack: bit_counts, top_heap
  tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);
  // Stack: bit_counts, top_heap, out

  tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);
  // Stack: bit_counts, top_heap, out, weights

  for (uint64_t i = 0; i < top_heap->n; i++) {
    weights->a[i] = top_heap->a[i].d; // Scores are already positive
  }

  // Rearrange stack to return: out, weights
  lua_remove(L, -4); // Remove bit_counts
  lua_remove(L, -3); // Remove top_heap

  return out;
}

// Document frequency based feature selection with minimum coverage guarantee
static inline tk_ivec_t *tk_ivec_bits_top_df (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_features,
  double min_df,
  double max_df,
  uint64_t top_k
) {
  // Sort set_bits for efficient processing
  tk_ivec_asc(set_bits, 0, set_bits->n);

  // Stack: (empty)
  // Create IDF scores vector
  tk_dvec_t *idf_scores = tk_dvec_create(L, n_features, 0, 0);
  // Stack: idf_scores
  tk_dvec_zero(idf_scores);
  int i_scores = tk_lua_absindex(L, -1);

  // Count documents per feature using sets
  tk_iuset_t **feature_docs = (tk_iuset_t **)calloc(n_features, sizeof(tk_iuset_t *));
  for (uint64_t i = 0; i < n_features; i++) {
    feature_docs[i] = tk_iuset_create();
  }

  // Process set_bits to build document-feature mapping
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    if (bit_idx < 0) continue;

    uint64_t sample_idx = (uint64_t)bit_idx / n_features;
    uint64_t feature_idx = (uint64_t)bit_idx % n_features;

    if (sample_idx < n_samples && feature_idx < n_features) {
      int absent;
      tk_iuset_put(feature_docs[feature_idx], (int64_t)sample_idx, &absent);
    }
  }

  // Build heap for top-k features with IDF weighting
  // Stack: idf_scores
  tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
  // Stack: idf_scores, top_heap

  // Convert thresholds to absolute counts
  // Negative values = absolute counts, positive = percentages
  double min_df_abs = min_df < 0 ? -min_df : min_df * n_samples;
  double max_df_abs = max_df < 0 ? -max_df : max_df * n_samples;

  for (uint64_t i = 0; i < n_features; i++) {
    double df_count = (double)tk_iuset_size(feature_docs[i]);

    // Calculate IDF score: log((N + 1) / (df + 1))
    double idf = log((double)(n_samples + 1) / (df_count + 1));
    idf_scores->a[i] = idf;

    // Only include features within bounds
    if (df_count >= min_df_abs && df_count <= max_df_abs) {
      // Use heap to maintain top-k features (highest IDF)
      tk_rank_t r = { (int64_t)i, idf };
      tk_rvec_hmin(top_heap, top_k, r);
    }
  }

  // Sort the final heap in descending order (highest scores first)
  tk_rvec_desc(top_heap, 0, top_heap->n);

  // Extract results
  // Stack: idf_scores, top_heap
  tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);  // Get feature indices
  // Stack: idf_scores, top_heap, out
  tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);  // Create weights vector
  // Stack: idf_scores, top_heap, out, weights

  // Fill weights with IDF scores
  for (uint64_t i = 0; i < top_heap->n; i++) {
    weights->a[i] = top_heap->a[i].d;  // Scores are already positive
  }

  // Cleanup
  for (uint64_t i = 0; i < n_features; i++) {
    tk_iuset_destroy(feature_docs[i]);
  }
  free(feature_docs);

  // Rearrange stack to return: out, weights
  // Current: idf_scores, top_heap, out, weights
  lua_remove(L, i_scores);       // top_heap, out, weights
  lua_remove(L, -3);             // out, weights

  return out;
}

// Document frequency based feature selection for dense cvec bitmaps
static inline tk_ivec_t *tk_cvec_bits_top_df (
  lua_State *L,
  tk_cvec_t *bitmap,
  uint64_t n_samples,
  uint64_t n_features,
  double min_df,
  double max_df,
  uint64_t top_k
) {
  // Stack: (empty)

  // Create IDF scores vector
  tk_dvec_t *idf_scores = tk_dvec_create(L, n_features, 0, 0);
  // Stack: idf_scores
  tk_dvec_zero(idf_scores);
  int i_scores = tk_lua_absindex(L, -1);

  // Count documents per feature using sets
  tk_iuset_t **feature_docs = (tk_iuset_t **)calloc(n_features, sizeof(tk_iuset_t *));
  for (uint64_t i = 0; i < n_features; i++) {
    feature_docs[i] = tk_iuset_create();
  }

  // Process bitmap to build document-feature mapping
  uint8_t *data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (data[sample_offset + byte_idx] & (1u << bit_idx)) {
        int absent;
        tk_iuset_put(feature_docs[f], (int64_t)s, &absent);
      }
    }
  }

  // Build heap for top-k features with IDF weighting
  // Stack: idf_scores
  tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
  // Stack: idf_scores, top_heap

  // Convert thresholds to absolute counts
  // Negative values = absolute counts, positive = percentages
  double min_df_abs = min_df < 0 ? -min_df : min_df * n_samples;
  double max_df_abs = max_df < 0 ? -max_df : max_df * n_samples;

  for (uint64_t i = 0; i < n_features; i++) {
    double df_count = (double)tk_iuset_size(feature_docs[i]);

    // Calculate IDF score: log((N + 1) / (df + 1))
    double idf = log((double)(n_samples + 1) / (df_count + 1));
    idf_scores->a[i] = idf;

    // Only include features within bounds
    if (df_count >= min_df_abs && df_count <= max_df_abs) {
      // Use heap to maintain top-k features (highest IDF)
      tk_rank_t r = { (int64_t)i, idf };
      tk_rvec_hmin(top_heap, top_k, r);
    }
  }

  // Sort the final heap in descending order (highest scores first)
  tk_rvec_desc(top_heap, 0, top_heap->n);

  // Extract results
  // Stack: idf_scores, top_heap
  tk_ivec_t *out = tk_ivec_from_rvec(L, top_heap);  // Get feature indices
  // Stack: idf_scores, top_heap, out
  tk_dvec_t *weights = tk_dvec_create(L, top_heap->n, 0, 0);  // Create weights vector
  // Stack: idf_scores, top_heap, out, weights

  // Fill weights with IDF scores
  for (uint64_t i = 0; i < top_heap->n; i++) {
    weights->a[i] = top_heap->a[i].d;  // Scores are already positive
  }

  // Cleanup
  for (uint64_t i = 0; i < n_features; i++) {
    tk_iuset_destroy(feature_docs[i]);
  }
  free(feature_docs);

  // Rearrange stack to return: out, weights
  // Current: idf_scores, top_heap, out, weights
  lua_remove(L, i_scores);       // top_heap, out, weights
  lua_remove(L, -3);             // out, weights

  return out;
}

#endif
