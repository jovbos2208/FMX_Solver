#include "geom/BVH.hpp"
#include <algorithm>
#include <stack>

namespace fmx::geom {

static inline double max3(double a, double b, double c) { return std::max(a, std::max(b, c)); }

Aabb BVHOccluder::tri_bounds(const Triangle& t) {
  Aabb b; b.expand(t.v0); b.expand(t.v1); b.expand(t.v2); return b;
}

BVHOccluder::BVHOccluder(const std::vector<Triangle>& tris) : m_tris(tris) {
  m_indices.resize(m_tris.size());
  for (size_t i = 0; i < m_tris.size(); ++i) m_indices[i] = static_cast<int>(i);
  m_nodes.reserve(2 * m_tris.size());
  if (!m_tris.empty()) build_node(0, static_cast<int>(m_tris.size()));
}

int BVHOccluder::build_node(int start, int count) {
  BVHNode node; node.start = start; node.count = count; node.leaf = (count <= 8);
  // bounds
  Aabb box; for (int i = 0; i < count; ++i) box.expand(tri_bounds(m_tris[m_indices[start+i]]));
  node.box = box;
  int idx = static_cast<int>(m_nodes.size());
  m_nodes.push_back(node);
  if (node.leaf) return idx;

  // choose split axis by largest extent of centroids
  Aabb cb; for (int i = 0; i < count; ++i) {
    const auto& t = m_tris[m_indices[start+i]];
    fmx::Vec3 c = (t.v0 + t.v1 + t.v2) / 3.0; cb.expand(c);
  }
  fmx::Vec3 e = cb.extent(); int axis = 0;
  if (e.y > e.x && e.y >= e.z) axis = 1; else if (e.z > e.x && e.z >= e.y) axis = 2;

  int mid = start + count/2;
  std::nth_element(m_indices.begin()+start, m_indices.begin()+mid, m_indices.begin()+start+count,
    [&](int a, int b){
      auto ca = (m_tris[a].v0 + m_tris[a].v1 + m_tris[a].v2) / 3.0;
      auto cb2 = (m_tris[b].v0 + m_tris[b].v1 + m_tris[b].v2) / 3.0;
      double va = (axis==0?ca.x:(axis==1?ca.y:ca.z));
      double vb = (axis==0?cb2.x:(axis==1?cb2.y:cb2.z));
      return va < vb;
    });

  int left = build_node(start, mid - start);
  int right = build_node(mid, start + count - mid);
  m_nodes[idx].left = left;
  m_nodes[idx].right = right;
  return idx;
}

bool BVHOccluder::ray_triangle(const fmx::Vec3& ro, const fmx::Vec3& rd, const Triangle& t, double t_max) {
  // Möller–Trumbore
  const double eps = 1e-7;
  fmx::Vec3 v0v1 = t.v1 - t.v0;
  fmx::Vec3 v0v2 = t.v2 - t.v0;
  fmx::Vec3 pvec = fmx::Vec3::cross(rd, v0v2);
  double det = fmx::Vec3::dot(v0v1, pvec);
  if (std::abs(det) < eps) return false;
  double invDet = 1.0 / det;
  fmx::Vec3 tvec = ro - t.v0;
  double u = fmx::Vec3::dot(tvec, pvec) * invDet;
  if (u < 0.0 || u > 1.0) return false;
  fmx::Vec3 qvec = fmx::Vec3::cross(tvec, v0v1);
  double v = fmx::Vec3::dot(rd, qvec) * invDet;
  if (v < 0.0 || u + v > 1.0) return false;
  double tparam = fmx::Vec3::dot(v0v2, qvec) * invDet;
  return (tparam > 1e-5 && tparam < t_max);
}

bool BVHOccluder::traverse_any(const fmx::Vec3& ro, const fmx::Vec3& rd, double t_max) const {
  if (m_nodes.empty()) return false;
  std::stack<int> st; st.push(0);
  while (!st.empty()) {
    int ni = st.top(); st.pop();
    const BVHNode& n = m_nodes[ni];
    if (!n.box.intersect(ro, rd, t_max)) continue;
    if (n.leaf) {
      for (int i = 0; i < n.count; ++i) {
        int triIdx = m_indices[n.start + i];
        if (ray_triangle(ro, rd, m_tris[triIdx], t_max)) return true;
      }
    } else {
      st.push(n.left);
      st.push(n.right);
    }
  }
  return false;
}

bool BVHOccluder::any_hit(const Ray& r, double t_max) const {
  // Offset origin by small step along ray to avoid self-intersection
  fmx::Vec3 ro = r.o + r.d * 1e-6;
  return traverse_any(ro, r.d, t_max);
}

} // namespace fmx::geom

