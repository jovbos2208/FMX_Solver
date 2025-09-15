#include "geom/Mesh.hpp"

#include <fstream>
#include <sstream>
#include <cctype>

namespace fmx::geom {

using fmx::Vec3;

static inline std::string to_lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::optional<Mesh> Mesh::load(const std::string& path, std::string* err) {
  auto lower = to_lower(path);
  if (lower.size() >= 4 && lower.substr(lower.size()-4) == ".obj")
    return loadOBJ(path, err);
  if (lower.size() >= 4 && lower.substr(lower.size()-4) == ".stl")
    return loadSTL(path, err);
  if (err) *err = "Unsupported mesh extension: " + path;
  return std::nullopt;
}

std::optional<Mesh> Mesh::loadOBJ(const std::string& path, std::string* err) {
  std::ifstream in(path);
  if (!in) { if (err) *err = "Failed to open OBJ: " + path; return std::nullopt; }

  std::vector<Vec3> verts;
  Mesh m;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string tok; ss >> tok;
    if (tok == "v") {
      double x,y,z; if (!(ss >> x >> y >> z)) continue; verts.emplace_back(x,y,z);
    } else if (tok == "f") {
      // Support simple faces: f i j k (1-based indices, no texture/normals)
      int i1=0,i2=0,i3=0;
      auto read_index = [&](const std::string& s)->int{
        // Handle forms like "3", "3/..", "3/../.."
        size_t pos = s.find('/');
        return std::stoi(pos==std::string::npos ? s : s.substr(0,pos));
      };
      std::string s1,s2,s3; if (!(ss >> s1 >> s2 >> s3)) continue;
      i1 = read_index(s1); i2 = read_index(s2); i3 = read_index(s3);
      if (i1==0||i2==0||i3==0) continue;
      auto v0 = verts[static_cast<size_t>(i1-1)];
      auto v1 = verts[static_cast<size_t>(i2-1)];
      auto v2 = verts[static_cast<size_t>(i3-1)];
      m.tris.push_back({v0,v1,v2});
    }
  }
  return m;
}

std::optional<Mesh> Mesh::loadSTL(const std::string& path, std::string* err) {
  std::ifstream in(path);
  if (!in) { if (err) *err = "Failed to open STL: " + path; return std::nullopt; }
  // Only ASCII STL supported in this minimal loader.
  Mesh m; std::string tok;
  while (in >> tok) {
    if (tok == "facet") {
      // skip normal header
      std::string normal_label; in >> normal_label; double nx,ny,nz; in >> nx >> ny >> nz;
      std::string outer; in >> outer; // "outer"
      std::string loop; in >> loop;   // "loop"
      std::string vertex;
      Vec3 v[3];
      for (int k=0;k<3;k++) {
        in >> vertex; // "vertex"
        in >> v[k].x >> v[k].y >> v[k].z;
      }
      std::string endloop, endfacet; in >> endloop; in >> endfacet;
      m.tris.push_back({v[0],v[1],v[2]});
    }
  }
  if (m.tris.empty()) { if (err) *err = "No triangles parsed (ASCII STL expected): " + path; }
  return m;
}

std::vector<fmx::Facet> Mesh::to_facets(std::size_t material_id) const {
  std::vector<fmx::Facet> facets;
  facets.reserve(tris.size());
  for (const auto& t : tris) {
    Vec3 e1 = t.v1 - t.v0;
    Vec3 e2 = t.v2 - t.v0;
    Vec3 n = Vec3::cross(e1, e2);
    double area2 = n.norm();
    double area = 0.5 * area2;
    Vec3 nn = (area2 > 0.0) ? (n / area2) : Vec3{0,0,1};
    Vec3 rc = (t.v0 + t.v1 + t.v2) / 3.0;
    facets.push_back({area, nn, rc, material_id});
  }
  return facets;
}

} // namespace fmx::geom

