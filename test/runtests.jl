import OWENSPreComp
using Test

function issame(x1,x2)
    return abs(x1-x2) < 0.001*10.0^(Int(floor(log10(abs(x1)))))
end

# Read main input file
sloc, le_loc, chord, tw_aero_d, af_shape_file, int_str_file, ib_sp_stn,
    ob_sp_stn, ib_ch_loc, ob_ch_loc = OWENSPreComp.readmain("test01_composite_blade.pci")

naf = length(sloc)
input = Array{OWENSPreComp.Input,1}(undef, naf)
output = Array{OWENSPreComp.Output,1}(undef, naf)
tw_prime_d = OWENSPreComp.tw_rate(naf, sloc, tw_aero_d)
for i = 1:naf
    xu, yu, xl, yl = OWENSPreComp.readprecompprofile(af_shape_file[i])
    e1, e2, g12, nu12, rho, name = OWENSPreComp.readmaterials("materials.inp")

    nweb = length(ib_ch_loc)
    locW = Array{Float64,1}(undef, nweb)
    for j = 1:nweb
        if i < ib_sp_stn || i > ob_sp_stn
            nweb = 0
            locW = Float64[]
        else
            rle = le_loc[i]
            r1w = le_loc[ib_sp_stn]
            r2w = le_loc[ob_sp_stn]
            p1w = ib_ch_loc[j]
            p2w = ob_ch_loc[j]
            ch1 = chord[ib_sp_stn]
            ch2 = chord[ob_sp_stn]
            x1w = sloc[ib_sp_stn]
            l_web = sloc[ob_sp_stn] - x1w
            xlocn = (sloc[i]-x1w)/l_web
            ch = chord[i]
            locW[j] = rle - (r1w-p1w)*ch1*(1-xlocn)/ch-(r2w-p2w)*ch2*xlocn/ch
        end
    end

    locU, n_laminaU, n_pliesU, tU, thetaU, mat_idxU, locL, n_laminaL, n_pliesL,
    tL, thetaL, mat_idxL, locW, n_laminaW, n_pliesW, tW, thetaW, mat_idxW =
    OWENSPreComp.readcompositesection(int_str_file[i], locW)

    input[i] = OWENSPreComp.Input(chord[i], tw_aero_d[i], tw_prime_d[i], le_loc[i],
    vcat(xu, reverse(xl)[2:(end-1)]), vcat(yu, reverse(yl)[2:(end-1)]),
    e1, e2, g12, nu12, rho, locU, n_laminaU, float.(n_pliesU), tU, thetaU, mat_idxU,
    locL, n_laminaL, float.(n_pliesL), tL, thetaL, mat_idxL, locW, n_laminaW,
    float.(n_pliesW), tW, thetaW, mat_idxW)

    output[i] = OWENSPreComp.properties(input[i])
end

# Check fields vs. precomp generated file
open("./test01_composite_blade.out_gen") do f
    readline(f)
    readline(f)
    readline(f)
    readline(f)
    readline(f)
    readline(f)
    readline(f)
    readline(f)
    readline(f)
    istat = 0
    while !eof(f)
        istat += 1
        @testset "Station $istat" begin
        substr = split(chomp(readline(f)))
        #Check span location
        @test abs(round(sloc[istat]/sloc[end], digits = 3) -
            Meta.parse(substr[1])) < 0.001
        #Check chord length
        @test abs(round(input[istat].chord, digits = 3) -
            Meta.parse(substr[2])) < 0.001
        #Check aerodynamic twist
        @test abs(round(input[istat].tw_aero_d, digits = 3) -
            Meta.parse(substr[3])) < 0.001
        #Check ei_flap
        @test issame(output[istat].ei_flap,Meta.parse(substr[4]))
        #Check ei_lag
        @test issame(output[istat].ei_lag,Meta.parse(substr[5]))
        #Check gj
        @test issame(output[istat].gj,Meta.parse(substr[6]))
        #Check ea
        @test issame(output[istat].ea,Meta.parse(substr[7]))
        #Check s_fl
        @test issame(output[istat].s_fl,Meta.parse(substr[8]))
        #Check s_a
        @test issame(output[istat].s_af,Meta.parse(substr[9]))
        #Check s_al
        @test issame(output[istat].s_al,Meta.parse(substr[10]))
        #Check s_ft
        @test issame(output[istat].s_ft,Meta.parse(substr[11]))
        #Check s_lt
        @test issame(output[istat].s_lt,Meta.parse(substr[12]))
        #Check s_at
        @test issame(output[istat].s_at,Meta.parse(substr[13]))
        #Check x_sc
        @test abs(round(output[istat].x_sc, digits = 3)-Meta.parse(substr[14])) < 0.001
        #Check y_sc
        @test abs(round(output[istat].y_sc, digits = 3)-Meta.parse(substr[15])) < 0.001
        #Check x_sc
        @test abs(round(output[istat].x_tc, digits = 3)-Meta.parse(substr[16])) < 0.001
        #Check y_sc
        @test abs(round(output[istat].y_tc, digits = 3)-Meta.parse(substr[17])) < 0.001
        #Check mass
        @test issame(output[istat].mass,Meta.parse(substr[18]))
        #Check flap_iner
        @test issame(output[istat].flap_iner,Meta.parse(substr[19]))
        #Check lag_iner
        @test issame(output[istat].lag_iner,Meta.parse(substr[20]))
        #Check tw_iner
        @test abs(round(output[istat].tw_iner_d, digits = 3)-Meta.parse(substr[21])) < 0.001
        #Check x_cm
        @test abs(round(output[istat].x_cm, digits = 3)-Meta.parse(substr[22])) < 0.001
        #Check y_cm
        @test abs(round(output[istat].y_cm, digits = 3)-Meta.parse(substr[23])) < 0.001
    end
end
end
