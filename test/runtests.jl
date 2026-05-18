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

@testset "tw_rate linear twist" begin
    @test OWENSPreComp.tw_rate(3, [0, 1, 2], [0, 10, 20]) == [10, 10, 10]
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
