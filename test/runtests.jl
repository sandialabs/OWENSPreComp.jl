import OWENSPreComp
using Test

const FIXTURE_DIR = @__DIR__
const OUT_GEN_FIELDS = (
    :span_loc,
    :chord,
    :tw_aero,
    :ei_flap,
    :ei_lag,
    :gj,
    :ea,
    :s_fl,
    :s_af,
    :s_al,
    :s_ft,
    :s_lt,
    :s_at,
    :x_sc,
    :y_sc,
    :x_tc,
    :y_tc,
    :mass,
    :flap_iner,
    :lag_iner,
    :tw_iner,
    :x_cm,
    :y_cm,
)
const OutGenRow = NamedTuple{OUT_GEN_FIELDS,NTuple{length(OUT_GEN_FIELDS),Float64}}
const THREE_DECIMAL_ATOL = 5e-4
const SCIENTIFIC_RTOL = 5e-4
const SCIENTIFIC_ATOL = 1e-12

fixture_path(name) = joinpath(FIXTURE_DIR, name)

function read_out_gen(path)
    rows = OutGenRow[]
    open(path) do f
        for _ = 1:9
            readline(f)
        end
        while !eof(f)
            tokens = split(chomp(readline(f)))
            isempty(tokens) && continue
            length(tokens) == length(OUT_GEN_FIELDS) ||
                error("Expected $(length(OUT_GEN_FIELDS)) columns in $path, got $(length(tokens))")
            values = Tuple(parse(Float64, token) for token in tokens)
            push!(rows, OutGenRow(values))
        end
    end
    return rows
end

# Read main input file
sloc, le_loc, chord, tw_aero_d, af_shape_file, int_str_file, ib_sp_stn,
    ob_sp_stn, ib_ch_loc, ob_ch_loc = OWENSPreComp.readmain(fixture_path("test01_composite_blade.pci"))

naf = length(sloc)
input = Array{OWENSPreComp.Input,1}(undef, naf)
output = Array{OWENSPreComp.Output,1}(undef, naf)
tw_prime_d = OWENSPreComp.tw_rate(naf, sloc, tw_aero_d)
for i = 1:naf
    xu, yu, xl, yl = OWENSPreComp.readprecompprofile(fixture_path(af_shape_file[i]))
    e1, e2, g12, nu12, rho, name = OWENSPreComp.readmaterials(fixture_path("materials.inp"))

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
    OWENSPreComp.readcompositesection(fixture_path(int_str_file[i]), locW)

    input[i] = OWENSPreComp.Input(chord[i], tw_aero_d[i], tw_prime_d[i], le_loc[i],
    vcat(xu, reverse(xl)[2:(end-1)]), vcat(yu, reverse(yl)[2:(end-1)]),
    e1, e2, g12, nu12, rho, locU, n_laminaU, float.(n_pliesU), tU, thetaU, mat_idxU,
    locL, n_laminaL, float.(n_pliesL), tL, thetaL, mat_idxL, locW, n_laminaW,
    float.(n_pliesW), tW, thetaW, mat_idxW)

    output[i] = OWENSPreComp.properties(input[i])
end

expected_rows = read_out_gen(fixture_path("test01_composite_blade.out_gen"))

function modified_input(base; kwargs...)
    values = map(fieldnames(typeof(base))) do field
        haskey(kwargs, field) ? kwargs[field] : getfield(base, field)
    end
    return OWENSPreComp.Input(values...)
end

function with_value(values, index, value)
    copy_values = copy(values)
    copy_values[index] = value
    return copy_values
end

function thrown_message(f)
    try
        f()
    catch err
        return sprint(showerror, err)
    end
    return "NO_ERROR"
end

properties_error_message(pc_input) = thrown_message(() -> OWENSPreComp.properties(pc_input))

@testset "tw_rate linear twist" begin
    @test OWENSPreComp.tw_rate(3, [0, 1, 2], [0, 10, 20]) == [10, 10, 10]
end

@testset "Profile and materials readers" begin
    te_x = [1.0, 0.5, 0.0, 0.5, 1.0]
    te_y = [0.02, 0.1, 0.0, -0.1, -0.02]
    @test OWENSPreComp.TEtoTEdata(te_x, te_y) == (
        [0.0, 0.5, 1.0],
        [0.0, 0.1, 0.02],
        [0.0, 0.5, 1.0],
        [0.0, -0.1, -0.02],
    )

    reversed_te_y = [-0.02, -0.1, 0.0, 0.1, 0.02]
    @test OWENSPreComp.TEtoTEdata(te_x, reversed_te_y) == (
        [0.0, 0.5, 1.0],
        [0.0, 0.1, 0.02],
        [0.0, 0.5, 1.0],
        [0.0, -0.1, -0.02],
    )

    le_x = [0.0, 0.5, 1.0, 1.0, 0.5, 0.0]
    le_y = [0.0, -0.1, -0.02, 0.02, 0.1, 0.0]
    @test OWENSPreComp.LEtoLEdata(le_x, le_y) == (
        [0.0, 0.5, 1.0],
        [0.0, 0.1, 0.02],
        [0.0, 0.5, 1.0],
        [0.0, -0.1, -0.02],
    )

    mktempdir() do tmp
        profile_path = joinpath(tmp, "te_profile.dat")
        open(profile_path, "w") do io
            println(io, "header one")
            println(io, "header two")
            for (x, y) in zip(te_x, te_y)
                println(io, x, " ", y)
            end
            println(io)
            println(io, "ignored trailing row")
        end
        @test OWENSPreComp.readprofile(profile_path, 2, false) == (
            [0.0, 0.5, 1.0],
            [0.0, 0.1, 0.02],
            [0.0, 0.5, 1.0],
            [0.0, -0.1, -0.02],
        )

        materials_path = joinpath(tmp, "materials.inp")
        open(materials_path, "w") do io
            println(io, "header one")
            println(io, "header two")
            println(io, "header three")
            println(io, "2 20.0 10.0 5.0 0.2 1000.0 (Named Material)")
            println(io, "1 10.0 5.0 2.0 0.1 500.0")
        end
        e1, e2, g12, nu12, rho, name = OWENSPreComp.readmaterials(materials_path)
        @test e1 == [10.0, 20.0]
        @test e2 == [5.0, 10.0]
        @test g12 == [2.0, 5.0]
        @test nu12 == [0.1, 0.2]
        @test rho == [500.0, 1000.0]
        @test name == ["", "Named Material"]
        @test eltype(e1) === Float64
        @test eltype(name) === String
    end
end

@testset "properties input validation" begin
    base = input[1]
    webbed = input[ib_sp_stn]

    @test properties_error_message(modified_input(base; xnode = base.xnode[1:(end-1)])) ==
        "x and y node lengths do not match"
    @test properties_error_message(modified_input(base; e2 = base.e2[1:(end-1)])) ==
        "lengths of specified material properties do not match"
    @test properties_error_message(modified_input(base; xsec_nodeU = base.xsec_nodeU[1:(end-1)])) ==
        "length of `xsec_nodeU` must be one greater than `n_laminaU`"
    @test properties_error_message(modified_input(base; xsec_nodeL = base.xsec_nodeL[1:(end-1)])) ==
        "length of `xsec_nodeL` must be one greater than `n_laminaL`"
    @test properties_error_message(modified_input(webbed; loc_web = webbed.loc_web[1:(end-1)])) ==
        "length of `loc_web` must equal length of `n_laminaW`"
    @test properties_error_message(modified_input(base; t_lamU = base.t_lamU[1:(end-1)])) ==
        "lengths of upper surface lamina properties do not match"
    @test properties_error_message(modified_input(base; t_lamL = base.t_lamL[1:(end-1)])) ==
        "lengths of lower surface lamina properties do not match"
    @test properties_error_message(modified_input(webbed; t_lamW = webbed.t_lamW[1:(end-1)])) ==
        "lengths of web lamina properties do not match"

    inconsistent_nu = with_value(base.anu12, 1, sqrt(base.e1[1] / base.e2[1]) + 0.1)
    @test properties_error_message(modified_input(base; anu12 = inconsistent_nu)) ==
        "material 1 properties not consistent"
    @test properties_error_message(modified_input(base; xnode = [0.0, 1.0], ynode = [0.0, 0.0])) ==
        "min 3 nodes reqd to define airfoil geom"
    @test properties_error_message(modified_input(base; xnode = circshift(base.xnode, -1),
        ynode = circshift(base.ynode, -1))) == "the first airfoil node not a leading node"
    @test properties_error_message(modified_input(base; ynode = with_value(base.ynode, 1, 0.01))) ==
        "leading-edge node not located at (0,0)"
    @test properties_error_message(modified_input(base; xnode = with_value(base.xnode,
        argmax(base.xnode), 1.1))) == "trailing-edge node exceeds chord boundary"
    @test properties_error_message(modified_input(base; xnode = with_value(base.xnode, 2, base.xnode[1]))) ==
        "upper surface not single-valued"
    @test properties_error_message(modified_input(base; xnode = with_value(base.xnode,
        length(base.xnode) - 1, base.xnode[end]))) == "lower surface not single-valued"
    @test properties_error_message(modified_input(base; xsec_nodeU = with_value(base.xsec_nodeU, 1, -0.01))) ==
        "sector node x-location not positive"
    @test properties_error_message(modified_input(base; xsec_nodeU = with_value(base.xsec_nodeU,
        length(base.xsec_nodeU), 1.1))) == "upper-surf last sector node out of bounds"
    @test properties_error_message(modified_input(base; xsec_nodeL = with_value(base.xsec_nodeL,
        length(base.xsec_nodeL), 1.1))) == "lower-surf last sector node out of bounds"
    @test properties_error_message(modified_input(base; xsec_nodeU = with_value(base.xsec_nodeU, 2,
        base.xsec_nodeU[1]))) == "sector nodal x-locations not in ascending order"

    warned = modified_input(base; le_loc = -0.1)
    warned_output = @test_logs (:warn, "leading edge aft of reference axis") OWENSPreComp.properties(warned)
    @test warned_output.mass == output[1].mass
end

@testset "material transformed stiffness validation" begin
    @test OWENSPreComp.q_tildas(10.0, 5.0, 1.0, 0.2, 0.3, 2.0, 1) == [9.8 0.14; 0 1.982]
    @test thrown_message(() -> OWENSPreComp.q_tildas(1.0, 1.0, 2.0, 0.0, 0.0, 1.0, 3)) ==
        " check material no 3 properties; these are not physically realizable."
end

@testset "geometry helper bounds validation" begin
    @test thrown_message(() -> OWENSPreComp.embed_us(2.0, 2, [0.0, 1.0], [0.0, 0.0])) ==
        " ERROR unknown, consult NWTC"
    @test thrown_message(() -> OWENSPreComp.embed_ls(2.0, 2, [0.0, 1.0], [0.0, 0.0])) ==
        "unknown, consult NWTC"
end

@testset "Pinned first station reference" begin
    first_station = expected_rows[1]

    @test first_station.span_loc == 0.0
    @test first_station.chord == 0.664
    @test first_station.tw_aero == 0.0
    @test first_station.ei_flap == 0.5481e8
    @test first_station.ei_lag == 0.2739e8
    @test first_station.gj == 0.3003e8
    @test first_station.ea == 0.8840e9
    @test first_station.mass == 0.1080e3
    @test first_station.flap_iner == 0.6343e1
    @test first_station.lag_iner == 0.4066e1
end

# Check fields vs. precomp generated file
@test length(expected_rows) == naf
for (istat, expected) in enumerate(expected_rows)
    @testset "Station $istat" begin
        three_decimal_values = (
            span_loc = sloc[istat]/sloc[end],
            chord = input[istat].chord,
            tw_aero = input[istat].tw_aero_d,
            x_sc = output[istat].x_sc,
            y_sc = output[istat].y_sc,
            x_tc = output[istat].x_tc,
            y_tc = output[istat].y_tc,
            tw_iner = output[istat].tw_iner_d,
            x_cm = output[istat].x_cm,
            y_cm = output[istat].y_cm,
        )
        scientific_values = (
            ei_flap = output[istat].ei_flap,
            ei_lag = output[istat].ei_lag,
            gj = output[istat].gj,
            ea = output[istat].ea,
            s_fl = output[istat].s_fl,
            s_af = output[istat].s_af,
            s_al = output[istat].s_al,
            s_ft = output[istat].s_ft,
            s_lt = output[istat].s_lt,
            s_at = output[istat].s_at,
            mass = output[istat].mass,
            flap_iner = output[istat].flap_iner,
            lag_iner = output[istat].lag_iner,
        )

        for (name, actual) in pairs(three_decimal_values)
            @test isapprox(actual, getproperty(expected, name); rtol = 0.0, atol = THREE_DECIMAL_ATOL)
        end
        for (name, actual) in pairs(scientific_values)
            @test isapprox(actual, getproperty(expected, name); rtol = SCIENTIFIC_RTOL, atol = SCIENTIFIC_ATOL)
        end
    end
end
