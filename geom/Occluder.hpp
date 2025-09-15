// Occlusion interface (stub). Real BVH/Embree backends can implement this.
#pragma once

#include <vector>
#include "core/types.hpp"
#include "geom/Mesh.hpp"

namespace fmx::geom {

struct Ray { fmx::Vec3 o, d; }; // origin, direction (normalized)

class Occluder {
public:
  virtual ~Occluder() = default;
  virtual bool any_hit(const Ray& r, double t_max) const = 0;
};

// No-occlusion implementation (always returns false)
class NoneOccluder : public Occluder {
public:
  bool any_hit(const Ray&, double) const override { return false; }
};

} // namespace fmx::geom

