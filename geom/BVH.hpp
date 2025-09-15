// Simple median-split BVH occluder with ray-triangle any-hit
#pragma once

#include <vector>
#include <limits>
#include "core/types.hpp"
#include "geom/Mesh.hpp"
#include "geom/Occluder.hpp"

namespace fmx::geom {

struct Aabb {
  fmx::Vec3 lo{ std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
  fmx::Vec3 hi{ -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
  void expand(const fmx::Vec3& p) {
    lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y); lo.z = std::min(lo.z, p.z);
    hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y); hi.z = std::max(hi.z, p.z);
  }
  void expand(const Aabb& b) { expand(b.lo); expand(b.hi); }
  fmx::Vec3 extent() const { return {hi.x - lo.x, hi.y - lo.y, hi.z - lo.z}; }
  bool intersect(const fmx::Vec3& ro, const fmx::Vec3& rd, double t_max) const {
    // Slab test
    double t0 = 0.0, t1 = t_max;
    for (int i = 0; i < 3; ++i) {
      double o = (i==0?ro.x:(i==1?ro.y:ro.z));
      double d = (i==0?rd.x:(i==1?rd.y:rd.z));
      double invd = 1.0 / d;
      double tNear = (((i==0)?lo.x:(i==1?lo.y:lo.z)) - o) * invd;
      double tFar  = (((i==0)?hi.x:(i==1?hi.y:hi.z)) - o) * invd;
      if (invd < 0.0) std::swap(tNear, tFar);
      t0 = tNear > t0 ? tNear : t0;
      t1 = tFar  < t1 ? tFar  : t1;
      if (t1 < t0) return false;
    }
    return t1 >= t0;
  }
};

struct BVHNode { Aabb box; int left{-1}, right{-1}; int start{0}, count{0}; bool leaf{false}; };

class BVHOccluder : public Occluder {
public:
  explicit BVHOccluder(const std::vector<Triangle>& tris);
  bool any_hit(const Ray& r, double t_max) const override;

private:
  std::vector<Triangle> m_tris;
  std::vector<int> m_indices;
  std::vector<BVHNode> m_nodes;

  int build_node(int start, int count);
  static Aabb tri_bounds(const Triangle& t);
  bool traverse_any(const fmx::Vec3& ro, const fmx::Vec3& rd, double t_max) const;
  static bool ray_triangle(const fmx::Vec3& ro, const fmx::Vec3& rd, const Triangle& t, double t_max);
};

} // namespace fmx::geom

