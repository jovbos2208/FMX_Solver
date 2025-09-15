// Minimal core types for FMX
#pragma once

#include <cmath>
#include <cstddef>

namespace fmx {

struct Vec3 {
  double x{0}, y{0}, z{0};

  Vec3() = default;
  constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  constexpr Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
  constexpr Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
  Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
  Vec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

  friend constexpr Vec3 operator*(double s, const Vec3& v) { return v * s; }

  static constexpr double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
  }

  static constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
      a.y*b.z - a.z*b.y,
      a.z*b.x - a.x*b.z,
      a.x*b.y - a.y*b.x
    };
  }

  double norm() const { return std::sqrt(x*x + y*y + z*z); }
  double norm2() const { return x*x + y*y + z*z; }

  Vec3 normalized() const {
    double n = norm();
    if (n == 0.0) return *this;
    return *this / n;
  }
};

struct Facet {
  // Surface area [m^2]
  double area{0.0};
  // Unit normal (outward)
  Vec3 n{0,0,1};
  // Center position [m]
  Vec3 r_center{0,0,0};
  // Material id (index into materials table)
  std::size_t material_id{0};
};

} // namespace fmx
