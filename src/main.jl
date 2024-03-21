#parameters
eps=1e-10  #from OWENSPreCompPy.f90

"""
Struct for holding inputs to OWENSPreComp.properties()
"""
struct Input{R<:Real}
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
  n_laminaU::Array{Int,1}
  n_pliesU::Array{R,1}
  t_lamU::Array{R,1}
  tht_lamU::Array{R,1}
  mat_lamU::Array{Int,1}
  xsec_nodeL::Array{R,1}
  n_laminaL::Array{Int,1}
  n_pliesL::Array{R,1}
  t_lamL::Array{R,1}
  tht_lamL::Array{R,1}
  mat_lamL::Array{Int,1}
  loc_web::Array{R,1}
  n_laminaW::Array{Int,1}
  n_pliesW::Array{R,1}
  t_lamW::Array{R,1}
  tht_lamW::Array{R,1}
  mat_lamW::Array{Int,1}
end

"""
Struct type for holding outputs of OWENSPreComp.properties()
"""
struct Output{R<:Real}
  ei_flap::R
  ei_lag::R
  gj::R
  ea::R
  s_fl::R
  s_af::R
  s_al::R
  s_ft::R
  s_lt::R
  s_at::R
  x_sc::R
  y_sc::R
  x_tc::R
  y_tc::R
  mass::R
  flap_iner::R
  lag_iner::R
  tw_iner_d::R
  x_cm::R
  y_cm::R
end

"""
    properties(pc_input::Input)
Calculates span-variant structural properties for composite blades. Holds all
inputs and outputs to properties function in structs
"""
function properties(pc_input::Input{R}) where R<:Real
  pc_output = Output{R}(OWENSPreComp.properties(pc_input.chord,
    pc_input.tw_aero_d, pc_input.tw_prime_d, pc_input.le_loc, pc_input.xnode,
    pc_input.ynode, pc_input.e1, pc_input.e2,pc_input.g12, pc_input.anu12,
    pc_input.density,pc_input.xsec_nodeU, pc_input.n_laminaU, pc_input.n_pliesU,
    pc_input.t_lamU, pc_input.tht_lamU, pc_input.mat_lamU, pc_input.xsec_nodeL,
    pc_input.n_laminaL, pc_input.n_pliesL, pc_input.t_lamL, pc_input.tht_lamL,
    pc_input.mat_lamL, pc_input.loc_web, pc_input.n_laminaW, pc_input.n_pliesW,
    pc_input.t_lamW, pc_input.tht_lamW, pc_input.mat_lamW)[1:20]...)
  return pc_output
end

"""
    properties(chord::Array{<:Real,1}, tw_aero_d::Array{<:Real,1},
    tw_prime_d::Array{<:Real,1}, le_loc::Real, xnode::Array{<:Real,1},
    ynode::Array{<:Real,1}, e1::Array{<:Real,1}, e2::Array{<:Real,1},
    g12::Array{<:Real,1}, anu12::Array{<:Real,1}, density::Array{<:Real,1},
    xsec_nodeU::Array{<:Real,1}, n_laminaU::Array{Int64,1},
    n_pliesU::Array{Int64,1}, t_lamU::Array{<:Real,1},
    tht_lamU::Array{<:Real,1}, mat_lamU::Array{Int64,1},
    xsec_nodeL::Array{<:Real,1}, n_laminaL::Array{Int64,1},
    n_pliesL::Array{Int64,1}, t_lamL::Array{<:Real,1},
    tht_lamL::Array{<:Real,1}, mat_lamL::Array{Int64,1},
    loc_web::Array{<:Real,1}, n_laminaW::Array{Int64,1},
    n_pliesW::Array{Int64,1}, t_lamW::Array{<:Real,1},
    tht_lamW::Array{<:Real,1}, mat_lamW::Array{Int64,1})

Calculates span-variant structural properties for composite blades

# Inputs
* `chord::Real`: section chord length (m)
* `tw_aero_d::Real`: section twist angle (deg)
* `tw_prime_d::Real`: derivative of section twist angle w.r.t. span location (deg/m)
* `le_loc::Real`: leading edge location relative to reference axis (normalized by chord)
* `xnode::Array{<:Real,1}`: x airfoil coordinates starting at leading edge traversing upper surface and back around lower surface
* `ynode::Array{<:Real,1}`: y airfoil coordinates starting at leading edge traversing upper surface and back around lower surface
* `e1::Array{<:Real,1}`: E1
* `e2::Array{<:Real,1}`: E2
* `g12::Array{<:Real,1}`: G12
* `anu12::Array{<:Real,1}`: Nu12
* `density::Array{<:Real,1}`: density
* `xsec_nodeU::Array{<:Real,1}`: upper surface normalized chord location of sector boundaries
* `n_laminaU::Array{Int64,1}`: upper surface number of lamina in each sector
* `n_pliesU::Array{Int64,1}`: upper surface number of plies
* `t_lamU::Array{<:Real,1}`: upper surface ply thickness (m) for the lamina
* `tht_lamU::Array{<:Real,1}`: upper surface orientation (deg) for the lamina
* `mat_lamU::Array{Int64,1}`: upper surface material id for the lamina
* `xsec_nodeL::Array{<:Real,1}`: lower surface normalized chord location of sector boundaries
* `n_laminaL::Array{Int64,1}`: lower surface number of lamina in each sector
* `n_pliesL::Array{Int64,1}`: lower surface number of plies
* `t_lamL::Array{<:Real,1}`: lower surface ply thickness (m) for the lamina
* `tht_lamL::Array{<:Real,1}`: lower surface orientation (deg) for the lamina
* `mat_lamL::Array{Int64,1}`: lower surface material id for the lamina
* `loc_web::Array{<:Real,1}`: web normalized chord location of sector boundaries
* `n_laminaW::Array{Int64,1}`: web number of lamina in each sector
* `n_pliesW::Array{Int64,1}`: web number of plies
* `t_lamW::Array{<:Real,1}`: web ply thickness (m) for the lamina
* `tht_lamW::Array{<:Real,1}`: web orientation (deg) for the lamina
* `mat_lamW::Array{Int64,1}`: web material id for the lamina

# Outputs:
* `eifbar`: ei_flap, Section flap bending stiffness about the YE axis (Nm2)
* `eilbar`: ei_lag, Section lag (edgewise) bending stiffness about the XE axis (Nm2)
* `gjbar`: gj, Section torsion stiffness (Nm2)
* `eabar`: ea, Section axial stiffness (N)
* `eiflbar`: s_fl, Coupled flap-lag stiffness with respect to the XE-YE frame (Nm2)
* `sfbar`: s_af, Coupled axial-flap stiffness with respect to the XE-YE frame (Nm)
* `slbar`: s_al, Coupled axial-lag stiffness with respect to the XE-YE frame (Nm.)
* `sftbar`: s_ft, Coupled flap-torsion stiffness with respect to the XE-YE frame (Nm2)
* `sltbar`: s_lt, Coupled lag-torsion stiffness with respect to the XE-YE frame (Nm2)
* `satbar`: s_at, Coupled axial-torsion stiffness (Nm)
* `z_sc`: x_sc, X-coordinate of the shear-center offset with respect to the XR-YR axes (m)
* `y_sc`: y_sc, Chordwise offset of the section shear-center with respect to the reference frame, XR-YR (m)
* `ztc_ref`: x_tc, X-coordinate of the tension-center offset with respect to the XR-YR axes (m)
* `ytc_ref`: y_tc, Chordwise offset of the section tension-center with respect to the XR-YR axes (m)
* `mass`: mass, Section mass per unit length (Kg/m)
* `iflap_eta`: flap_iner, Section flap inertia about the YG axis per unit length (Kg-m)
* `ilag_zeta`: lag_iner, Section lag inertia about the XG axis per unit length (Kg-m)
* `tw_iner_d`: tw_iner_d, Orientation of the section principal inertia axes with respect the blade reference plane, θ (deg)
* `zcm_ref`: x_cm, X-coordinate of the center-of-mass offset with respect to the XR-YR axes (m)
* `ycm_ref`: y_cm, Chordwise offset of the section center of mass with respect to the XR-YR axes (m)
* `n_af_nodes`: number of airfoil nodes
* `n_materials`: number of materials
* `n_sctU`: number of sectors on upper
* `n_sctL`: number of sectors on lower
* `nwebin`: number of webs
* `n_laminaTotalU`: total number of lamina on upper
* `n_laminaTotalL`: total number of lamina on lower
* `n_laminaTotalW`: total number of lamina on webs
"""
function properties(chord::Real, tw_aero_d::Real,
  tw_prime_d::Real, le_loc::Real, xnode::AbstractArray{<:Real,1},
  ynode::AbstractArray{<:Real,1}, e1::AbstractArray{<:Real,1}, e2::AbstractArray{<:Real,1},
  g12::AbstractArray{<:Real,1}, anu12::AbstractArray{<:Real,1}, density::AbstractArray{<:Real,1},
  xsec_nodeU::AbstractArray{<:Real,1}, n_laminaU::AbstractArray{<:Real,1},
  n_pliesU::AbstractArray{<:Real,1}, t_lamU::AbstractArray{<:Real,1},
  tht_lamU::AbstractArray{<:Real,1}, mat_lamU::AbstractArray{<:Real,1},
  xsec_nodeL::AbstractArray{<:Real,1}, n_laminaL::AbstractArray{<:Real,1},
  n_pliesL::AbstractArray{<:Real,1}, t_lamL::AbstractArray{<:Real,1},
  tht_lamL::AbstractArray{<:Real,1}, mat_lamL::AbstractArray{<:Real,1},
  loc_web::AbstractArray{<:Real,1}, n_laminaW::AbstractArray{<:Real,1},
  n_pliesW::AbstractArray{<:Real,1}, t_lamW::AbstractArray{<:Real,1},
  tht_lamW::AbstractArray{<:Real,1}, mat_lamW::AbstractArray{<:Real,1})

# Implicitly assigned inputs
if length(xnode) != length(ynode)
  error("x and y node lengths do not match")
end
n_af_nodes = length(xnode)

if !(length(e1) == length(e2) == length(g12) == length(anu12) == length(density))
  error("lengths of specified material properties do not match")
end
n_materials = length(e1)

if ((length(xsec_nodeU)-1) != length(n_laminaU))
  error("length of `xsec_nodeU` must be one greater than `n_laminaU`")
end
n_sctU = length(n_laminaU)

if ((length(xsec_nodeL)-1) != length(n_laminaL))
  error("length of `xsec_nodeL` must be one greater than `n_laminaL`")
end
n_sctL = length(n_laminaL)

if (length(loc_web) != length(n_laminaW))
  error("length of `loc_web` must equal length of `n_laminaW`")
end
nwebin = length(loc_web)

if !(length(n_pliesU) == length(mat_lamU) == length(t_lamU) == length(tht_lamU))
  error("lengths of upper surface lamina properties do not match")
end
n_laminaTotalU = length(n_pliesU)

if !(length(n_pliesL) == length(mat_lamL) == length(t_lamL) == length(tht_lamL))
  error("lengths of lower surface lamina properties do not match")
end
n_laminaTotalL = length(n_pliesL)

if !(length(n_pliesW) == length(mat_lamW) == length(t_lamW) == length(tht_lamW))
  error("lengths of web lamina properties do not match")
end
n_laminaTotalW = length(n_pliesW)

# Declare all scalars as local
local tw_aero, tw_prime, tphip
local max_sectors, max_laminates
local webs_exist, tenode_u, tenode_l
local ieta1, izeta1, iepz, iemz, ipp,
    iqq, ipq, iflap_sc, ilag_sc, ifl_sc, iflap_cm, ilag_cm, ifl_cm,
    m_inertia, r_inertia
local x, xl, xr, y, yl, yr
local sths, s2ths, cths, c2ths, em_stiff, er_stiff, y_tc, z_tc,
    sigm2, ycm_sc, zcm_sc, sigma, tbar, q11t, q11yt, q11zt, dtbar,
    q2bar, zbart, ybart, tbart, q11ysqt, q11zsqt, q11yzt, rhot, rhoyt,
    rhozt, rhoysqt, rhozsqt, rhoyzt, pflap_stff, plag_stff,
    q11ya, q11za, ap, bp, cp, dp, ep, q11ysqa, q11zsqa, q11yza,
    rhoya, rhoza, rhoysqa, rhozsqa, rhoyza, q11yt_u, q11zt_u,
    q11yt_l, q11zt_l, qtil12t, qtil11t, qtil22t, rot, t, the_pa, th_pa,
    wdq2bar, w, xu1, xu2, xl1, xl2, yu1, yu2, yl1, yl2, y0,
    y0sq, ysg, zsg, z0, z0sq, ynd, str_tw
local idsec, ilam, is, iseg, iweb, nlam, ncounter, ndl2,
    ndu2, nsects, wreq, id_form
local newnode, nodes_u, nodes_l
local nseg, nseg_l, nseg_u, nseg_p, ndl1, ndu1
local rho_m, thp, qbar11, qbar22, qbar12, qbar16, qbar26, qbar66
local mat

# Initialize Arrays
tht_wlam = zeros(Real,nwebin,15)
twlam = zeros(Real,nwebin,15)
wmat_id = zeros(Int64,nwebin,15)
n_weblams = zeros(Int64,nwebin)

location = zeros(Int64,1)

xnode_u = zeros(Real,n_af_nodes+20)
xnode_l = zeros(Real,n_af_nodes+20)
ynode_u = zeros(Real,n_af_nodes+20)
ynode_l = zeros(Real,n_af_nodes+20)

weby_u = zeros(Real,nwebin)
weby_l = zeros(Real,nwebin)

yseg = zeros(Real,n_af_nodes+20)
zseg = zeros(Real,n_af_nodes+20)
wseg = zeros(Real,n_af_nodes+20)
sthseg = zeros(Real,n_af_nodes+20)
cthseg = zeros(Real,n_af_nodes+20)
s2thseg = zeros(Real,n_af_nodes+20)
c2thseg = zeros(Real,n_af_nodes+20)

n_scts = zeros(Int64,2)
isur = zeros(Int64,n_af_nodes+20)
idsect = zeros(Int64,n_af_nodes+20)

q11 = zeros(Real,n_materials)
q22 = zeros(Real,n_materials)
q12 = zeros(Real,n_materials)
q66 = zeros(Real,n_materials)
anud = zeros(Real,n_materials)
qtil = zeros(Real,2,2)

max_sectors = max(n_sctU,n_sctL,nwebin)
if length(n_laminaW) == 0
  max_laminates = max(maximum(n_laminaU),maximum(n_laminaL))
else
  max_laminates = max(maximum(n_laminaU),maximum(n_laminaL),maximum(n_laminaW))
end

n_laminas = zeros(Int64,2,max_sectors)
tht_lam = zeros(Real,2,max_sectors,max_laminates)
tlam = zeros(Real,2,max_sectors,max_laminates)
mat_id = zeros(Int64,2,max_sectors,max_laminates)
xsec_node = zeros(Real,2,max_sectors+1)


# convert twist angle to radians
tw_aero = tw_aero_d * (pi/180)
tw_prime = tw_prime_d * (pi/180)

# webs?

webs_exist = 1

#     check number of webs
if (nwebin == 0)
#             println(" ** no webs in this blade **")
  webs_exist = 0
end


# ---- checks --------------
#  check leading edge location
if (le_loc < 0)
    @warn("leading edge aft of reference axis")
end


# check materials
for i = 1:Int(n_materials)

    if (anu12[i] > sqrt(e1[i]/e2[i]))
      error("material ", i, " properties not consistent")
    end

end

# check airfoil nodes
if (n_af_nodes <= 2)
    error("min 3 nodes reqd to define airfoil geom")
end


#   check if the first airfoil node is a leading-edge node and at (0,0)
location = argmin(xnode)
if (location[1] != 1)
    error("the first airfoil node not a leading node")
end


if (abs(xnode[1]) > eps || abs(ynode[1]) > eps)
    error("leading-edge node not located at (0,0)")
end

#   identify trailing-edge end nodes on upper and lower surfaces
location = argmax(xnode)
if(abs(xnode[location[1]]) > 1)
    error("trailing-edge node exceeds chord boundary")
end



# ----------------


#   get th_prime and phi_prime
#tw_rate(naf, sloc, tw_aero, th_prime)

#for i in 1:naf
#    phi_prime[i] = 0.  # later: get it from aeroelastic code
#    tphip[i] = th_prime[i] + 0.5*phi_prime[i]
#end

tphip = tw_prime


# material properties
anud = 1.0 .- anu12.*anu12.*e2./e1
q11 = e1 ./ anud
q22 = e2 ./ anud
q12 = anu12.*e2 ./ anud
q66 = g12



# begin blade sections loop sec-sec-sec-sec-sec-sec-sec-sec-sec-sec--------


# ----------- airfoil data -------------------

#   identify trailing-edge end nodes on upper and lower surfaces
location = argmax(xnode)
tenode_u = location[1]

for i in location[1]:(n_af_nodes)
    if (abs(xnode[i]-xnode[location[1]])< eps)
        ncounter = i - location[1]
    end
end

tenode_l = tenode_u + ncounter


#   renumber airfoil nodes
#   (modify later using equivalence or pointers)
nodes_u = tenode_u
nodes_l = n_af_nodes - tenode_l + 2

for i = 1:nodes_u
    xnode_u[i] = xnode[i]
    ynode_u[i] = ynode[i]
end

xnode_l[1] = xnode[1]
ynode_l[1] = ynode[1]

for i = 2:tenode_l
    xnode_l[i] = xnode[n_af_nodes+2-i]
    ynode_l[i] = ynode[n_af_nodes+2-i]
end

# ----------------------------------------------


# ------ more checks -------------
#   ensure surfaces are single-valued functions

for i = 2:nodes_u
    if ((xnode_u[i] - xnode_u[i-1]) <= eps )
        error("upper surface not single-valued")
    end
end

for i = 2:nodes_l
    if ((xnode_l[i] - xnode_l[i-1]) <= eps )
        error("lower surface not single-valued")
    end
end

#   check clockwise node numbering at the le

if (ynode_u[2]/xnode_u[2] <= ynode_l[2]/xnode_l[2])
    error("airfoil node numbering not clockwise")
end

#   check for single-connectivity of the airfoil shape
#   (modify later using binary search)


for j = 2:(nodes_l - 1)   # loop over lower-surf nodes
    x = xnode_l[j]

    for i = 1:(nodes_u - 1)  # loop over upper-surf nodes

        xl = xnode_u[i]
        xr = xnode_u[i+1]

        if(x >= xl && x <= xr)
            yl = ynode_u[i]
            yr = ynode_u[i+1]
            y = yl + (yr-yl)*(x-xl)/(xr-xl)

            if(ynode_l[j] >= y)
                error("airfoil shape self-crossing")
            end
        end

    end    # end loop over upper-surf nodes
end   # end loop over lower-surf nodes

# ---------- end checks ---------------------


# -------------- webs ------------------


#   embed airfoil nodes at web-to-airfoil intersections


if (webs_exist == 1)
    for i = 1:nwebin
        ynd,nodes_u,xnode_u,ynode_u,newnode = embed_us(loc_web[i], nodes_u,
            xnode_u, ynode_u)
        weby_u[i] = ynd
        nodes_l,xnode_l,ynode_l,ynd,newnode = embed_ls(loc_web[i], nodes_l,
            xnode_l, ynode_l)
        weby_l[i] = ynd
    end
end


# ----------------------------------------------


# ------ internal structure data ------------
n_scts[1] = n_sctU
n_scts[2] = n_sctL
xsec_node[1, :] = xsec_nodeU
xsec_node[2, :] = xsec_nodeL

# unpack data
k = 1
for i = 1:n_sctU
    n_laminas[1, i] = n_laminaU[i]

    for j = 1:n_laminaU[i]
        tlam[1, i, j] = n_pliesU[k] * t_lamU[k]
        tht_lam[1, i, j] = tht_lamU[k] *(pi/180)
        mat_id[1, i, j] = mat_lamU[k]

        k = k + 1
    end
end

k = 1
for i = 1:n_sctL
    n_laminas[2, i] = n_laminaL[i]

    for j = 1:n_laminaL[i]
        tlam[2, i, j] = n_pliesL[k] * t_lamL[k]
        tht_lam[2, i, j] = tht_lamL[k] *(pi/180)
        mat_id[2, i, j] = mat_lamL[k]

        k = k + 1
    end
end

k = 1
for i = 1:nwebin
    n_weblams[i] = n_laminaW[i]

    for j = 1:n_laminaW[i]
        twlam[i, j] = n_pliesW[k] * t_lamW[k]
        tht_wlam[i, j] = tht_lamW[k] *(pi/180)
        wmat_id[i, j] = mat_lamW[k]

        k = k + 1
    end
end

for is = 1:2  # begin loop for blade surfaces

    nsects = n_scts[is]

    if (nsects <= 0)
        error("no of sectors not positive")
    end

    if (xsec_node[is,1] < 0)
        error("sector node x-location not positive")
    end

    if (is == 1)
        xu1 = xsec_node[is,1]
        xu2 = xsec_node[is,nsects+1]
        if (xu2 > xnode_u[nodes_u])
            error("upper-surf last sector node out of bounds")
        end
    else
        xl1 = xsec_node[is,1]
        xl2 = xsec_node[is,nsects+1]
        if (xl2 > xnode_l[nodes_l])
            error("lower-surf last sector node out of bounds") #this error is triggering
        end
    end

    for i = 1:nsects
        if (xsec_node[is,i+1] <= xsec_node[is,i])
            error("sector nodal x-locations not in ascending order")
        end
    end



    #embed airfoil nodes representing sectors bounds
    for i = 1:(nsects + 1)

      if(is == 1)
        ynd,nodes_u,xnode_u,ynode_u,newnode = embed_us(xsec_node[is,i], nodes_u,
            xnode_u, ynode_u)
        if(i == 1)
          yu1 = ynd
          ndu1 = newnode
        end
        if(i == nsects+1)
          yu2 = ynd
          ndu2 = newnode
        end

      end

      if(is == 2)

        nodes_l,xnode_l,ynode_l,ynd,newnode = embed_ls(xsec_node[is,i], nodes_l,
          xnode_l, ynode_l)

        if(i == 1)
          yl1 = ynd
          ndl1 = newnode
        end
        if(i == nsects+1)
          yl2 = ynd
          ndl2 = newnode
        end

      end

    end

end      # end blade surfaces loop


        #.... check for le and te non-closures and issue warning ....

        if (abs(xu1-xl1) > eps)

            @warn("the leading edge may be open; check closure")

        else

            if ((yu1-yl1) > eps)
                wreq = 1

                if (webs_exist != 0)
                    if (abs(xu1-loc_web[1]) < eps)
                        wreq = 0
                    end
                end


                if (wreq == 1)
                    @warn("open leading edge; check web requirement")
                end

            end

        end

    #
        if (abs(xu2-xl2) > eps)

            @warn("the trailing edge may be open; check closure")

        else

            if ((yu2-yl2) > eps)
                wreq = 1

                if (webs_exist != 0)
                    if (abs(xu2-loc_web[nwebin]) < eps)
                        wreq = 0
                    end
                end


                if (wreq == 1)
                    @warn("open trailing edge; check web requirement")
                end

            end

        end
        #................

        if (webs_exist == 1)

            if(loc_web[1] < xu1 || loc_web[1] < xl1)
                error("first web out of sectors-bounded airfoil")
            end

            if(loc_web[nwebin] > xu2 || loc_web[nwebin] > xl2)
                error("last web out of sectors-bounded airfoil")
            end

        end


        # ------------- Done Processing Inputs ----------------------


# ----------- Start Computations ------------------


#   identify segments groupings and get segs attributes
nseg_u = ndu2 - ndu1
nseg_l = ndl2 - ndl1
nseg_p = nseg_u + nseg_l    # no of peripheral segments

if(webs_exist == 1)
    nseg = nseg_p + nwebin  # total no of segments (webs in section)
else
    nseg = nseg_p      # total no of segments (no webs in section)
end

isur,idsect,yseg,zseg,wseg,sthseg,cthseg,s2thseg,c2thseg =
seg_info(chord, le_loc, nseg, nseg_u, nseg_p, xnode_u, ynode_u, xnode_l,
ynode_l, ndl1, ndu1, loc_web, weby_u, weby_l, n_scts, xsec_node)


#------------------------------------------

#   initialize for section (sc)
sigma = 0.0
eabar = 0.0
q11ya = 0.0
q11za = 0.0

#   segments loop for sc


for iseg = 1:nseg_p #begin paeripheral segments loop (sc)

#     retrieve seg attributes
    is = isur[iseg]
    idsec = idsect[iseg]
    ysg = yseg[iseg]
    zsg = zseg[iseg]
    w = wseg[iseg]
    sths = sthseg[iseg]
    cths = cthseg[iseg]
#       s2ths = s2thseg[iseg]
#       c2ths = c2thseg[iseg]  # commented out

    nlam = n_laminas[is,idsec]    # for sector seg

    #     initialization for seg (sc)
    tbar = 0.0
    q11t = 0.0
    q11yt_u = 0.0
    q11zt_u = 0.0
    q11yt_l = 0.0
    q11zt_l = 0.0

    for ilam = 1:nlam #laminas loop (sc)

        t = tlam[is,idsec,ilam]          #thickness
        thp = tht_lam[is,idsec,ilam]  # ply angle
        mat = mat_id[is,idsec,ilam]    #material

        tbar = tbar + t/2.0
        y0 = ysg - ((-1.0)^is*tbar*sths )
        z0 = zsg + ((-1.0)^is*tbar*cths )


        # obtain qtil for specified mat
        qbar11,qbar22,qbar12,qbar16,qbar26,qbar66,rho_m = q_bars(mat, thp,
          density, q11, q22, q12, q66)
        qtil=q_tildas(qbar11, qbar22, qbar12, qbar16, qbar26, qbar66, mat)

        # add seg-laminas contributions (sc)
        qtil11t = qtil[1,1]*t
        q11t = q11t + qtil11t
        if iseg <= nseg_u
            q11yt_u = q11yt_u + qtil11t*y0
            q11zt_u = q11zt_u + qtil11t*z0
        else
            q11yt_l = q11yt_l + qtil11t*y0
            q11zt_l = q11zt_l + qtil11t*z0
        end

        tbar = tbar + t/2.0

    end         # end laminas loop

    # add seg contributions (sc)
    sigma = sigma + w*abs(zsg + ((-1.0)^is)*0.5*tbar*cths)*cths
    eabar = eabar + q11t*w
    q11ya = q11ya + (q11yt_u+q11yt_l)*w
    q11za = q11za + (q11zt_u+q11zt_l)*w


end           #end af_periph segment loop (sc)

# get section sc
y_sc = q11ya/eabar     #wrt r-frame
z_sc = q11za/eabar     #wrt r-frame


#---------------- end section sc -----------

#   initializations for section (properties)

eabar = 0.0
q11ya = 0.0
q11za = 0.0
ap = 0.0
bp = 0.0
cp = 0.0
dp = 0.0
ep = 0.0
q11ysqa = 0.0
q11zsqa = 0.0
q11yza = 0.0

mass = 0.0
rhoya = 0.0
rhoza = 0.0
rhoysqa = 0.0
rhozsqa = 0.0
rhoyza = 0.0

#   segments loop (for properties)

for iseg = 1:nseg   #begin segment loop (properties)

    # retrieve seg attributes
    is = isur[iseg]
    idsec = idsect[iseg]
    ysg = yseg[iseg]
    zsg = zseg[iseg]
    w = wseg[iseg]
    sths = sthseg[iseg]
    cths = cthseg[iseg]
    s2ths = s2thseg[iseg]
    c2ths = c2thseg[iseg]

    if is > 0
        nlam = n_laminas[is,idsec]  # for sector seg
    else
        iweb = idsec
        nlam = n_weblams[iweb]      # for web seg
    end

    # initialization for seg (properties)
    tbar = 0.0
    q11t = 0.0
    q11yt = 0.0
    q11zt = 0.0
    dtbar = 0.0
    q2bar = 0.0
    zbart = 0.0
    ybart = 0.0
    tbart = 0.0
    q11ysqt = 0.0
    q11zsqt = 0.0
    q11yzt = 0.0

    rhot = 0.0
    rhoyt = 0.0
    rhozt = 0.0
    rhoysqt = 0.0
    rhozsqt = 0.0
    rhoyzt = 0.0

    for ilam = 1:nlam #laminas loop (properties)

        if (is > 0)
            t = tlam[is,idsec,ilam]          #thickness
            thp = tht_lam[is,idsec,ilam]  # ply angle
            mat = mat_id[is,idsec,ilam]      # material
            tbar = tbar + t/2.0
            y0 = ysg - ((-1.0)^is)*tbar*sths - y_sc
            z0 = zsg + ((-1.0)^is)*tbar*cths - z_sc
        else
            t = twlam[iweb,ilam]
            thp = tht_wlam[iweb,ilam]
            mat = wmat_id[iweb,ilam]
            tbar = tbar + t/2.0
            y0 = ysg - tbar/2.0 - y_sc
            z0 = zsg - z_sc
        end

        y0sq = y0*y0
        z0sq = z0*z0

        # obtain qtil and rho for specified mat
        qbar11,qbar22,qbar12,qbar16,qbar26,qbar66,rho_m = q_bars(mat, thp,
          density, q11, q22, q12, q66)
        qtil=q_tildas(qbar11, qbar22, qbar12, qbar16, qbar26, qbar66, mat)


        ieta1 = (t^2)/12.0
        izeta1 = (w^2)/12.0
        iepz = 0.5*(ieta1+izeta1)
        iemz = 0.5*(ieta1-izeta1)
        ipp = iepz+(iemz*c2ths)   # check this block later
        iqq = iepz-(iemz*c2ths)
        ipq = iemz*s2ths

        qtil11t = qtil[1,1]*t
        rot = rho_m*t

        #add laminas contributions (properties) at current segment

        if (is > 0) # peripheral segs contribution

            qtil12t = qtil[1,2]*t
            qtil22t = qtil[2,2]*t

            q11t = q11t + qtil11t
            q11yt = q11yt + qtil11t*y0
            q11zt = q11zt + qtil11t*z0

            dtbar = dtbar + qtil12t*(y0sq+z0sq)*tphip*t
            q2bar = q2bar + qtil22t    # later: retain only this block
            zbart = zbart + z0*qtil12t
            ybart = ybart + y0*qtil12t
            tbart = tbart + qtil12t

            q11ysqt = q11ysqt + qtil11t*(y0sq+iqq)
            q11zsqt = q11zsqt + qtil11t*(z0sq+ipp)
            q11yzt = q11yzt + qtil11t*(y0*z0+ipq)

            rhot = rhot + rot
            rhoyt = rhoyt + rot*y0
            rhozt = rhozt + rot*z0
            rhoysqt = rhoysqt + rot*(y0sq+iqq)
            rhozsqt = rhozsqt + rot*(z0sq+ipp)
            rhoyzt = rhoyzt + rot*(y0*z0+ipq)

        else            #web segs contribution

            q11t = q11t + qtil11t
            q11yt = q11yt + qtil11t*y0
            q11zt = q11zt + qtil11t*z0
            q11ysqt = q11ysqt + qtil11t*(y0sq+iqq)
            q11zsqt = q11zsqt + qtil11t*(z0sq+ipp)
            q11yzt = q11yzt + qtil11t*(y0*z0+ipq)

            rhot = rhot + rot
            rhoyt = rhoyt + rot*y0
            rhozt = rhozt + rot*z0
            rhoysqt = rhoysqt + rot*(y0sq+iqq)
            rhozsqt = rhozsqt + rot*(z0sq+ipp)
            rhoyzt = rhoyzt + rot*(y0*z0+ipq)

        end

        tbar = tbar + t/2.0

    end         # end laminas loop


    # add seg contributions to obtain sec props about ref_parallel axes at sc
    eabar = eabar + q11t*w
    q11ya = q11ya + q11yt*w
    q11za = q11za + q11zt*w
    q11ysqa = q11ysqa + q11ysqt*w
    q11zsqa = q11zsqa + q11zsqt*w
    q11yza = q11yza + q11yzt*w

    if (is > 0)
        wdq2bar = w/q2bar
        ap = ap + wdq2bar
        bp = bp + wdq2bar*tbart
        cp = cp + wdq2bar*dtbar
        dp = dp + wdq2bar*zbart
        ep = ep + wdq2bar*ybart
    end

    mass = mass + rhot*w
    rhoya = rhoya + rhoyt*w
    rhoza = rhoza + rhozt*w
    rhoysqa = rhoysqa + rhoysqt*w
    rhozsqa = rhozsqa + rhozsqt*w
    rhoyza = rhoyza + rhoyzt*w

end       #end af_periph segment loop (properties)

#  get more section properties # about ref_parallel axes at sc

y_tc = q11ya/eabar
z_tc = q11za/eabar

sfbar = q11za
slbar = q11ya
eifbar = q11zsqa
eilbar = q11ysqa
eiflbar = q11yza

sigm2 = sigma*2.0
gjbar = sigm2*(sigm2+cp)/ap
sftbar = -sigm2*dp/ap
sltbar = -sigm2*ep/ap
satbar = sigm2*bp/ap

ycm_sc =   rhoya/mass #wrt sc
zcm_sc =   rhoza/mass #wrt sc

iflap_sc = rhozsqa #wrt sc
ilag_sc = rhoysqa   #wrt sc
ifl_sc = rhoyza     #wrt sc

# get section tc and cm

ytc_ref =  y_tc + y_sc  #wrt the ref axes
ztc_ref =  z_tc + z_sc  #wrt the ref axes

ycm_ref =  ycm_sc + y_sc    #wrt the ref axes
zcm_ref =  zcm_sc + z_sc    #wrt the ref axes

# moments of inertia # about ref_parallel axes at cm

iflap_cm = iflap_sc - mass*zcm_sc^2
ilag_cm = ilag_sc - mass*ycm_sc^2
ifl_cm = ifl_sc - mass*ycm_sc*zcm_sc

# inertia principal axes orientation and moments of inertia

m_inertia = 0.5*(ilag_cm + iflap_cm)
r_inertia = sqrt(0.25*((ilag_cm-iflap_cm)^2)+ifl_cm^2)

if (iflap_cm <= ilag_cm)
    iflap_eta = m_inertia - r_inertia
    ilag_zeta = m_inertia + r_inertia
else
    iflap_eta = m_inertia + r_inertia
    ilag_zeta = m_inertia - r_inertia
end

if (ilag_cm == iflap_cm)
    th_pa = pi/4.0
    if (abs(ifl_cm/iflap_cm) < 1e-6)
      th_pa = 0.0
    end
else
    th_pa = 0.5*abs(atan(2.0*ifl_cm/(ilag_cm-iflap_cm)))
end

if (abs(ifl_cm) < eps)
    th_pa = 0.0
else          # check this block later
    if (iflap_cm >= ilag_cm && ifl_cm > 0.0)
      th_pa = -th_pa
    end
    if (iflap_cm >= ilag_cm && ifl_cm < 0.0)
      th_pa = th_pa
    end
    if (iflap_cm < ilag_cm && ifl_cm > 0.0)
      th_pa = th_pa
    end
    if (iflap_cm < ilag_cm && ifl_cm < 0.0)
      th_pa = -th_pa
    end
end

# elastic principal axes orientation and principal bending stiffneses

em_stiff = 0.5*(eilbar+eifbar)
er_stiff = sqrt(0.25*((eilbar-eifbar)^2)+eiflbar^2)

if (eifbar <= eilbar)
    pflap_stff = em_stiff - er_stiff
    plag_stff = em_stiff + er_stiff
else
    pflap_stff = em_stiff + er_stiff
    plag_stff = em_stiff - er_stiff
end

if (eilbar == eifbar)
    the_pa = pi/4.0
else
    the_pa = 0.5*abs(atan(2.0*eiflbar/(eilbar-eifbar)))
end

if (abs(eiflbar) < eps)
    the_pa = 0.0
# check this block later
elseif (eifbar >= eilbar && eiflbar > 0)
      the_pa = -the_pa
elseif (eifbar >= eilbar && eiflbar < 0)
      the_pa = the_pa
elseif (eifbar < eilbar && eiflbar > 0)
      the_pa = the_pa
elseif (eifbar < eilbar && eiflbar < 0)
      the_pa = -the_pa
end

#---------------- end properties computation -----------


# ---------- prepare outputs --------------



id_form = 1  # hardwired for wt's

if (id_form == 1)
    tw_iner_d = tw_aero - th_pa
    str_tw =  tw_aero - the_pa
    y_sc = -y_sc
    ytc_ref = -ytc_ref
    ycm_ref = -ycm_ref
else         # for h/c
    #       note: for h/c, th_aero input is +ve acc to h/c convention
    tw_iner_d = tw_aero + th_pa
    str_tw =  tw_aero + the_pa
end



# conversions
eiflbar = -eiflbar
sfbar = -sfbar
sltbar = -sltbar
tw_iner_d = tw_iner_d * (180/pi)

return (eifbar, eilbar, gjbar, eabar, eiflbar,
sfbar, slbar, sftbar, sltbar, satbar,
z_sc, y_sc, ztc_ref, ytc_ref, mass, iflap_eta,
ilag_zeta, tw_iner_d, zcm_ref, ycm_ref,
n_af_nodes, n_materials, n_sctU, n_sctL, nwebin,
n_laminaTotalU, n_laminaTotalL, n_laminaTotalW)

end #properties

"""
    embed_us(x, nodes_u, xnode_u, ynode_u)

purpose: embed a node in the upper-surface airfoil section nodes
NOTE: nodal x coordinates must be in ascending order

# Arguments:
* `x::Real`: x-coordinate of node to be embedded in the u-surf
* `nodes_u::Int`: no of current nodes on the upper surface
* `xnode_u::Array{<:Real,1}`: x-nodes on upper surface
* `ynode_u::Array{<:Real,1}`: y-nodes on upper surface
# Outputs:
* `nodes_u`: revised no of current nodes on upper surface
* `xnode_u`: x-nodes on upper surface
* `ynode_u`: y-nodes on upper surface
* `y`: y-coordinate of node embedded in the u-surf
* `newnode` : number of the embedded node
"""
function embed_us(x::Real, nodes_u::Int, xnode_u, ynode_u)

newnode = -1
isave = 0

for i = 1:(nodes_u - 1)   # loop over upper-surf nodes
    xl = xnode_u[i]
    xr = xnode_u[i + 1]
    yl = ynode_u[i]

    if (abs(x - xl) <= eps)
        newnode = 0
        isave = i
        y = yl
        break
    elseif (x < (xr - eps))
        yr = ynode_u[i + 1]
        y = yl + ((yr - yl) * (x - xl) / (xr-xl))
        newnode = i + 1
        break
    end
end       # end loop over upper-surf nodes

if (newnode == -1)
    if (abs(x - xnode_u[nodes_u]) <= eps)
        newnode = 0
        isave = nodes_u
        y =  ynode_u[nodes_u]
    else
        error(" ERROR unknown, consult NWTC")  #this error is triggering. what's going on?
    end
end

if (newnode > 0)
    nodes_u = nodes_u + 1

    for i in nodes_u:-1:(newnode + 1)
        xnode_u[i] = xnode_u[i-1]
        ynode_u[i] = ynode_u[i-1]
    end

    xnode_u[newnode] = x
    ynode_u[newnode] = y

else
    newnode = isave
end
return y,nodes_u,xnode_u,ynode_u,newnode

end # end of function embed_us

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

"""
    embed_ls(x, nodes_l, xnode_l, ynode_l)

purpose: embed a node in the lower-surface airfoil section nodes
NOTE: nodal x coordinates must be in ascending order

# Arguments:
* `x::Real`: x-coordinate of node to be embedded in the l-surf
* `nodes_l::Int`: no of current nodes on the lower surface
* `xnode_l::Array{<:Real,1}`: x-nodes on lower surface
* `ynode_l::Array{<:Real,1}`: y-nodes on lower surface
# Outputs:
* `nodes_l`: revised no of current nodes on lower surface
* `xnode_l`: x-nodes on lower surface
* `ynode_l`: y-nodes on lower surface
* `y`: y-coordinate of node embedded in the l-surf
* `newnode` : number of the embedded node
"""
function embed_ls(x, nodes_l, xnode_l, ynode_l)

newnode = -1

for i = 1:(nodes_l - 1)   # loop over lower-surf nodes
    xl = xnode_l[i]
    xr = xnode_l[i + 1]
    yl = ynode_l[i]

    if (abs(x-xl) <= eps)
        newnode = 0
        isave = i
        y = yl
        break
    elseif (x < (xr - eps))
        yr = ynode_l[i + 1]
        y = yl + ((yr - yl) * (x - xl) / (xr - xl))
        newnode = i + 1
        break

    end

end       # end loop over lower-surf nodes

if newnode == -1
    if (abs(x - xnode_l[nodes_l]) <= eps)
        newnode = 0
        isave = nodes_l
        y =  ynode_l[nodes_l]
    else
        error("unknown, consult NWTC")
    end
end

if newnode > 0
    nodes_l = nodes_l + 1

    for i in nodes_l:-1:(newnode + 1)
        xnode_l[i] = xnode_l[i - 1]
        ynode_l[i] = ynode_l[i - 1]
    end

    xnode_l[newnode] = x
    ynode_l[newnode] = y

else
    newnode = isave
end
return (nodes_l,xnode_l,ynode_l,y,newnode)

end #embed_ls

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

"""
    seg_info(ch, rle, nseg, nseg_u, nseg_p, xnode_u, ynode_u,
    xnode_l, ynode_l, ndl1, ndu1, loc_web, weby_u, weby_l, n_scts,
    xsec_node)

NOTE: coord transformation from xaf-yaf to yre-zref and seg info
"""
function seg_info(ch, rle, nseg, nseg_u, nseg_p, xnode_u, ynode_u,
  xnode_l, ynode_l, ndl1, ndu1, loc_web, weby_u, weby_l, n_scts,
  xsec_node)

# Initialize outputs
isur=zeros(Int64,nseg)
idsect=zeros(Int64,nseg)
yseg=zeros(Real,nseg)
zseg=zeros(Real,nseg)
wseg=zeros(Real,nseg)
sthseg=zeros(Real,nseg)
cthseg=zeros(Real,nseg)
s2thseg=zeros(Real,nseg)
c2thseg=zeros(Real,nseg)

local iseg, icheck, iweb, nd_a
local xa, ya, xb, yb, xba, yba, thseg

for iseg = 1:nseg   # seg numbering from le clockwise
    is = -1
    if iseg <= nseg_u  # upper surface segs

        nd_a = ndu1 + iseg - 1
        xa = xnode_u[nd_a]
        ya = ynode_u[nd_a]
        xb = xnode_u[nd_a + 1]
        yb = ynode_u[nd_a + 1]
        is = 1
    else
        if iseg <= nseg_p   # lower surface segs
            nd_a = ndl1+iseg-nseg_u-1
            xa = xnode_l[nd_a]        #xref of node toward le (in a/f ref frame)
            ya = ynode_l[nd_a]        #yref of node toward le (in new ref frame)
            xb = xnode_l[nd_a + 1]    #xref of node toward te (in a/f ref frame)
            yb = ynode_l[nd_a + 1]    #yref of node toward te (in new ref frame)
            is = 2
        end

        if iseg > nseg_p  # web segs
            iweb = iseg - nseg_p
            xa = loc_web[iweb]
            xb = xa
            ya = weby_u[iweb]
            yb = weby_l[iweb]
            is = 0
        end

    end  # end seg group identification


    if is == -1
        println("iseg=", iseg)
        error("unknown, contact NREL")
    end

    isur[iseg] = is


    if is > 0 #id assocaited sect number
        icheck = 0
        for i = 1:n_scts[is]
            if (xa > (xsec_node[is,i]-eps)) &&
               (xb < (xsec_node[is,i+1]+eps))
                idsect[iseg] = i
                icheck = 1
                break
            end
        end
    end


    if (icheck == 0)
        error("unknown, contact NREL")
    end

    if is == 0
      idsect[iseg] = iweb   #id assocaited web number
    end

    xba = xb - xa
    yba = ya - yb
    yseg[iseg] = ch*(2.0*rle-xa-xb)/2.0 #yref coord of mid-seg pt (in r-frame)
    zseg[iseg] = ch*(ya+yb)/2.0    #zref coord of mid-seg pt (in r-frame)
    wseg[iseg] = ch*sqrt(xba^2 + yba^2)

    if (is == 0)
        thseg = -pi/2.0
    else
        thseg = atan(yba/xba) # thseg +ve in new y-z ref frame
    end

    sthseg[iseg] = sin(thseg)
    cthseg[iseg] = cos(thseg)
    s2thseg[iseg] = sin(2.0*thseg)
    c2thseg[iseg] = cos(2.0*thseg)


end   # end seg loop

return isur,idsect,yseg,zseg,wseg,sthseg,cthseg,s2thseg,c2thseg

end #seg_info

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

"""
    tw_rate(naf, sloc, tw_aero)

# Arguments
* `naf`: no of blade stations
* `sloc`: vector of station locations
* `tw_aero_d`: vector of twist distribution in degrees
# Outputs
* `th_prime_d`: vector of twist rates in degrees
"""
function tw_rate(naf, sloc, tw_aero_d)

tw_aero = tw_aero_d*(pi/180)

th_prime = Array{Real, 1}(undef, naf)
for i = 2:(naf-1)
    f0 = tw_aero[i]
    f1 = tw_aero[i-1]
    f2 = tw_aero[i+1]
    h1 = sloc[i] - sloc[i-1]
    h2 = sloc[i+1] - sloc[i]
    th_prime[i] = (h1*(f2-f0) + h2*(f0-f1))/(2.0*h1*h2)

end

th_prime[1] = (tw_aero[2]-tw_aero[1])/(sloc[2]-sloc[1])
th_prime[naf]=(tw_aero[naf]-tw_aero[naf-1])/(sloc[naf]-sloc[naf-1])

return th_prime*(180/pi)
end #tw_rate

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

"""
    q_bars(mat, thp, density, q11, q22, q12, q66)

# Arguments
* `mat`: material id
* `thp`: ply orientation
* `density`:
* `q11`:
* `q22`:
* `q12`:
* `q66`:
# Outputs
* `qbar11`:
* `qbar22`:
* `qbar12`:
* `qbar16`:
* `qbar26`:
* `qbar66`:
* `rho_m`:
"""
function q_bars(mat, thp, density, q11, q22, q12, q66)

ct = cos(thp)
st = sin(thp)

c2t = ct*ct
c3t = c2t*ct
c4t = c3t*ct
s2t = st*st
s3t = s2t*st
s4t = s3t*st
s2thsq = 4.0*s2t*c2t

k11 = q11[mat]
k22 = q22[mat]
k12 = q12[mat]
k66 = q66[mat]
kmm = k11 -k12 -2.0*k66
kmp = k12 -k22 +2.0*k66

qbar11 = k11*c4t + 0.5*(k12+2.0*k66)*s2thsq + k22*s4t
qbar22 = k11*s4t + 0.5*(k12+2.0*k66)*s2thsq + k22*c4t
qbar12 = 0.25*(k11+k22-4.0*k66)*s2thsq + k12*(s4t+c4t)
qbar16 = kmm*st*c3t + kmp*s3t*ct
qbar26 = kmm*s3t*ct + kmp*st*c3t
qbar66 = 0.25*(kmm+k22-k12)*s2thsq  + k66*(s4t+c4t)

rho_m = density[mat]

return qbar11,qbar22,qbar12,qbar16,qbar26,qbar66,rho_m

end #q_bars

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function q_tildas(qbar11, qbar22, qbar12, qbar16, qbar26, qbar66, mat)

# Initialize Output
qtil=zeros(Real,2,2)

qtil[1,1] = qbar11 - qbar12*qbar12/qbar22
if (qtil[1,1] < 0.0)
    error(" check material no ", mat, " properties; these are not physically realizable.")

end

qtil[1,2] = qbar16 - (qbar12*qbar26/qbar22)
qtil[2,2] = qbar66 - (qbar26*qbar26/qbar22)

return qtil

end #q_tildas

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
