// Fundamental constants and species masses (SI units)
#pragma once

namespace fmx::units {

// Mathematical constants
inline constexpr double pi = 3.141592653589793238462643383279502884;

// Boltzmann constant [J/K]
inline constexpr double k_B = 1.380649e-23;

// Atomic mass unit [kg]
inline constexpr double m_u = 1.66053906660e-27;

// Species masses [kg]
inline constexpr double m_H  = 1.00784 * m_u;      // Hydrogen atom
inline constexpr double m_He = 4.002602 * m_u;     // Helium atom
inline constexpr double m_O  = 15.999 * m_u;       // Oxygen atom
inline constexpr double m_N2 = 28.0134 * m_u;      // Nitrogen molecule
inline constexpr double m_O2 = 31.9988 * m_u;      // Oxygen molecule
// Additional species [kg]
inline constexpr double m_Ar = 39.948 * m_u;       // Argon atom
inline constexpr double m_N  = 14.0067 * m_u;      // Nitrogen atom
inline constexpr double m_NO = (14.0067 + 15.999) * m_u; // Nitric oxide molecule

} // namespace fmx::units
