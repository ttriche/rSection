// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup constants
//! @{


// the long lengths of the constants are for future support of "long double"
// and any smart compiler that does high-precision computation at compile-time

template<typename eT>
class Math
  {
  public:
  
  //! ratio of any circle's circumference to its diameter
  static const eT pi()        { return eT(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679); }
  
  //! base of the natural logarithm
  static const eT e()         { return eT(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274); }
  
  //! Euler's constant, aka Euler-Mascheroni constant
  static const eT euler()     { return eT(0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495); }
  
  //! golden ratio
  static const eT gratio()    { return eT(1.6180339887498948482045868343656381177203091798057628621354486227052604628189024497072072041893911374); }
  
  //! square root of 2
  static const eT sqrt2()     { return eT(1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727); }
  
  //! the difference between 1 and the least value greater than 1 that is representable
  static const eT eps()       { return std::numeric_limits<eT>::epsilon(); }
  
  //! log of the minimum representable value
  static const eT log_min()   { static const eT out = std::log(std::numeric_limits<eT>::min()); return out; }
    
  //! log of the maximum representable value
  static const eT log_max()   { static const eT out = std::log(std::numeric_limits<eT>::max()); return out; }
  };



//! Physical constants taken from NIST and WolframAlpha on 2009-06-23
//! http://physics.nist.gov/cuu/Constants
//! http://www.wolframalpha.com
//! See also http://en.wikipedia.org/wiki/Physical_constant
template<typename eT>
class Phy
  {
  public:

  //! atomic mass constant (in kg)
  static const eT m_u()       {  return eT(1.660538782e-27); }

  //! Avogadro constant
  static const eT N_A()       {  return eT(6.02214179e23); }

  //! Boltzmann constant (in joules per kelvin)
  static const eT k()         {  return eT(1.3806504e-23); }

  //! Boltzmann constant (in eV/K)
  static const eT k_evk()     {  return eT(8.617343e-5); }

  //! Bohr radius (in meters)
  static const eT a_0()       { return eT(0.52917720859e-10); }

  //! Bohr magneton
  static const eT mu_B()      { return(927.400915e-26); }

  //! characteristic impedance of vacuum (in ohms)
  static const eT Z_0()       { return eT(3.76730313461771e-2); }

  //! conductance quantum (in siemens)
  static const eT G_0()       { return eT(7.7480917004e-5); }

  //! Coulomb's constant (in meters per farad)
  static const eT k_e()       { return eT(8.9875517873681764e9); }

  //! electric constant (in farads per meter)
  static const eT eps_0()     { return eT(8.85418781762039e-12); }

  //! electron mass (in kg)
  static const eT m_e()       { return eT(9.10938215e-31); }

  //! electron volt (in joules)
  static const eT eV()        { return eT(1.602176487e-19); }

  //! elementary charge (in coulombs)
  static const eT e()         { return eT(1.602176487e-19); }
  
  //! Faraday constant (in coulombs)
  static const eT F()         { return eT(96485.3399); }

  //! fine-structure constant
  static const eT alpha()     { return eT(7.2973525376e-3); }

  //! inverse fine-structure constant
  static const eT alpha_inv() { return eT(137.035999679); }

  //! Josephson constant
  static const eT K_J()       { return eT(483597.891e9); }

  //! magnetic constant (in henries per meter)
  static const eT mu_0()      { return eT(1.25663706143592e-06); }

  //! magnetic flux quantum (in webers)
  static const eT phi_0()     { return eT(2.067833667e-15); }

  //! molar gas constant (in joules per mole kelvin)
  static const eT R()         { return eT(8.314472); }

  //! Newtonian constant of gravitation (in newton square meters per kilogram squared)
  static const eT G()         { return eT(6.67428e-11); }

  //! Planck constant (in joule seconds)
  static const eT h()         { return eT(6.62606896e-34); }

  //! Planck constant over 2 pi, aka reduced Planck constant (in joule seconds)
  static const eT h_bar()     { return eT(1.054571628e-34); }

  //! proton mass (in kg)
  static const eT m_p()       { return eT(1.672621637e-27); }

  //! Rydberg constant (in reciprocal meters)
  static const eT R_inf()     { return eT(10973731.568527); }

  //! speed of light in vacuum (in meters per second)
  static const eT c_0()       { return eT(299792458.0); }

  //! Stefan-Boltzmann constant
  static const eT sigma()     { return eT(5.670400e-8); }

  //! von Klitzing constant (in ohms)
  static const eT R_k()       { return eT(25812.807557); }

  //! Wien wavelength displacement law constant
  static const eT b()         { return eT(2.8977685e-3); }
  };



typedef Math<float>  fmath;
typedef Math<double> math;

typedef Phy<float>   fphy;
typedef Phy<double>  phy;



//! @}
