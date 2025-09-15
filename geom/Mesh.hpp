// Minimal OBJ/ASCII STL mesh loader and facet extraction
#pragma once

#include <string>
#include <vector>
#include <optional>
#include "core/types.hpp"

namespace fmx::geom {

struct Triangle { fmx::Vec3 v0, v1, v2; };

struct Mesh {
  std::vector<Triangle> tris;

  static std::optional<Mesh> load(const std::string& path, std::string* err = nullptr);
  static std::optional<Mesh> loadOBJ(const std::string& path, std::string* err = nullptr);
  static std::optional<Mesh> loadSTL(const std::string& path, std::string* err = nullptr);

  // Convert triangles to solver facets with centers, normals, and areas
  std::vector<fmx::Facet> to_facets(std::size_t material_id = 0) const;
};

} // namespace fmx::geom

