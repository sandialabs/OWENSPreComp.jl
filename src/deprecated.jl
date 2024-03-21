"""
Struct for holding inputs to OWENSPreComp.properties()
"""
struct input{I<:Integer,R<:Real}
  chord::R
  tw_aero_d::R
  tw_prime_d::R
  le_loc::R
  xnode::Array{R,1}
  ynode::Array{R,1}
  e1::Array{R,1}
  e2::Array{R,1}
  g12::Array{R,1}
  anu12::Array{R,1}
  density::Array{R,1}
  xsec_nodeU::Array{R,1}
  n_laminaU::Array{I,1}
  n_pliesU::Array{I,1}
  t_lamU::Array{R,1}
  tht_lamU::Array{R,1}
  mat_lamU::Array{I,1}
  xsec_nodeL::Array{R,1}
  n_laminaL::Array{I,1}
  n_pliesL::Array{I,1}
  t_lamL::Array{R,1}
  tht_lamL::Array{R,1}
  mat_lamL::Array{I,1}
  loc_web::Array{R,1}
  n_laminaW::Array{I,1}
  n_pliesW::Array{I,1}
  t_lamW::Array{R,1}
  tht_lamW::Array{R,1}
  mat_lamW::Array{I,1}
end

output = Output
