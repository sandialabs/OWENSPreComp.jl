"""
readmain("fname::String")

Reads a main OWENSPreComp input file.
Returns sloc,le_loc,chord,tw_aero,af_shape_file,int_str_file,
ib_sp_stn,ob_sp_stn,ib_ch_loc,ob_ch_loc
"""
function readmain(fname::String)
    local sloc,le_loc,chord,tw_aero,af_shape_file,int_str_file,ib_sp_stn,
    ob_sp_stn,ib_ch_loc,ob_ch_loc

    open(fname) do f
        readline(f)
        title = chomp(readline(f))
        readline(f)
        readline(f)
        bl_length = Meta.parse(split(chomp(readline(f)))[1])
        naf = Meta.parse(split(chomp(readline(f)))[1])
        n_materials = Meta.parse(split(chomp(readline(f)))[1])
        out_format = Meta.parse(split(chomp(readline(f)))[1])
        TabDelim = lowercase(split(chomp(readline(f)))[1]) == "t"
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)

        sloc = Array{Float64,1}(undef, naf)
        le_loc = Array{Float64,1}(undef, naf)
        chord = Array{Float64,1}(undef, naf)
        tw_aero = Array{Float64,1}(undef, naf)
        af_shape_file = Array{String,1}(undef, naf)
        int_str_file = Array{String,1}(undef, naf)
        for i = 1:naf
            substr = split(chomp(readline(f)))
            sloc[i] = Meta.parse(substr[1])
            le_loc[i] = Meta.parse(substr[2])
            chord[i] = Meta.parse(substr[3])
            tw_aero[i] = Meta.parse(substr[4])
            af_shape_file[i] = replace(substr[5],"'"=>"")
            int_str_file[i] = replace(substr[6],"'"=>"")
        end
        readline(f)
        readline(f)
        readline(f)

        nweb = Meta.parse(split(chomp(readline(f)))[1])
        ib_ch_loc = Array{Float64,1}(undef, nweb)
        ob_ch_loc = Array{Float64,1}(undef, nweb)
        if nweb >= 1
            ib_sp_stn = Meta.parse(split(chomp(readline(f)))[1])
            ob_sp_stn = Meta.parse(split(chomp(readline(f)))[1])

            readline(f)
            readline(f)

            for iw = 1:nweb
                substr = split(chomp(readline(f)))
                ib_ch_loc[iw] = Meta.parse(substr[2])
                ob_ch_loc[iw] = Meta.parse(substr[3])
            end
        end
        # Dimensionalize station locations
        sloc = sloc*bl_length
    end

    return (sloc,le_loc,chord,tw_aero,af_shape_file,int_str_file,
    ib_sp_stn,ob_sp_stn,ib_ch_loc,ob_ch_loc)
end

"""
readcompositesection(fname::String,locw::Array{Float64,1})

Reads a composite section input file.
Returns locU, n_laminaU, n_pliesU, tU, thetaU, mat_idxU, locL, n_laminaL,
n_pliesL, tL, thetaL, mat_idxL, locW, n_laminaW, n_pliesW, tW, thetaW, mat_idxW
"""
function readcompositesection(fname::String,locW::Array{Float64,1})
    open(fname) do f
        readline(f)
        readline(f)
        readline(f)

        # number of sectors
        n_sector = Int64(Meta.parse(split(chomp(readline(f)))[1]))

        readline(f)
        readline(f)

        # read normalized chord locations
        locU = [Float64(Meta.parse(x)) for x in split(chomp(readline(f)))]

        n_laminaU, n_pliesU, tU, thetaU, mat_idxU = readsectorsfromfile(f, n_sector)

        readline(f)
        readline(f)
        readline(f)

        # number of sectors
        n_sector = Int64(Meta.parse(split(chomp(readline(f)))[1]))

        readline(f)
        readline(f)

        locL = [Float64(Meta.parse(x)) for x in split(chomp(readline(f)))]

        n_laminaL, n_pliesL, tL, thetaL, mat_idxL = readsectorsfromfile(f, n_sector)

        readline(f)
        readline(f)
        readline(f)
        readline(f)

        n_sector = length(locW)

        n_laminaW, n_pliesW, tW, thetaW, mat_idxW = readsectorsfromfile(f, n_sector)

        return locU, n_laminaU, n_pliesU, tU, thetaU, mat_idxU, locL, n_laminaL,
        n_pliesL, tL, thetaL, mat_idxL, locW, n_laminaW, n_pliesW, tW, thetaW, mat_idxW
    end
end

"""
readsectorsfromfile(f::IOStream, n_sector::Int64)

Reads OWENSPreComp sector. Returns n_lamina,n_plies, t, theta, mat_idx
"""
function readsectorsfromfile(f::IOStream, n_sector::Int64)
    n_lamina = Int64[]
    n_plies = Int64[]
    t = Float64[]
    theta = Float64[]
    mat_idx = Int64[]

    for i = 1:n_sector
        readline(f)
        readline(f)

        line = chomp(readline(f))

        # if line == ""
        #   return []
        # end
        push!(n_lamina, Int64(Meta.parse(split(line)[2])))

        readline(f)
        readline(f)
        readline(f)
        readline(f)

        n_plies_S = zeros(Int64,n_lamina[i])
        t_S = zeros(Float64,n_lamina[i])
        theta_S = zeros(Float64,n_lamina[i])
        mat_idx_S = zeros(Int64,n_lamina[i])

        for j = 1:n_lamina[i]
            array = split(chomp(readline(f)))
            n_plies_S[j] = Int64(Meta.parse(array[2]))
            t_S[j] = Float64(Meta.parse(array[3]))
            theta_S[j] = Float64(Meta.parse(array[4]))
            mat_idx_S[j] = Int64(Meta.parse(array[5]))
        end

        append!(n_plies,n_plies_S)
        append!(t,t_S)
        append!(theta,theta_S)
        append!(mat_idx,mat_idx_S)
    end

    return n_lamina,n_plies, t, theta, mat_idx
end

"""
readprecompprofile(filename::String)
Reads precomp profile file. Returns xu, yu, xl, yl
"""
function readprecompprofile(filename::String)
    return readprofile(filename,4,true)
end

"""
readprofile(filename::String, numHeaderlines::Int64, LEtoLE::Bool)
Reads precomp profile file. Returns xu, yu, xl, yl
"""
function readprofile(filename::String, numHeaderlines::Int64, LEtoLE::Bool)
    x = Float64[]
    y = Float64[]

    # open file
    open(filename, "r") do f

        # skip through header
        for i = 1:numHeaderlines
            readline(f)
        end

        # loop through
        while !eof(f)
            line = chomp(readline(f))
            if strip(line) == ""
                break  # break if empty line
            end
            data = split(line)
            push!(x, Float64(Meta.parse(data[1])))
            push!(y, Float64(Meta.parse(data[2])))
        end
    end

    # close nose if LE to LE
    if LEtoLE
        push!(x,x[1])
        push!(y,y[1])
        return LEtoLEdata(x, y)
    else
        return TEtoTEdata(x, y)
    end
end

function TEtoTEdata(x, y)
    # separate into 2 halves
    local le_loc
    le_loc = argmin(x)

    xu = x[le_loc:-1:1]
    yu = y[le_loc:-1:1]
    xl = x[le_loc:end]
    yl = y[le_loc:end]

    # check if coordinates were input in other direction
    if y[1] < y[0]
        temp = yu
        yu = yl
        yl = temp

        temp = xu
        xu = xl
        xl = temp
    end

    return (xu, yu, xl, yl)
end

function LEtoLEdata(x, y)
    # separate into 2 halves
    local iulast
    local illast
    iulast = argmax(x)
    illast = argmax(x)
    if x[iulast] == x[iulast+1]  # blunt t.e.
        illast = iulast+1
    else
        illast = iulast
    end

    xu = x[1:iulast]
    yu = y[1:iulast]
    xl = x[end:-1:illast]
    yl = y[end:-1:illast]

    # check if coordinates were input in other direction
    if y[2] < y[1]
        temp = yu
        yu = yl
        yl = temp

        temp = xu
        xu = xl
        xl = temp
    end

    return (xu, yu, xl, yl)
end

"""
readmaterials(fname = "materials.inp")

reads material properties from OWENSPreComp materials input file `fname`
returns e1,e2,g12,nu12,rho,name
"""
function readmaterials(fname::String = "materials.inp")
    open(fname) do f

        readline(f)
        readline(f)
        readline(f)

        mat_id = Int64[]
        e1 =  Float64[]
        e2 = Float64[]
        g12 = Float64[]
        nu12 = Float64[]
        rho = Float64[]
        name = String[]
        while !eof(f)
            array = split(chomp(readline(f)))
            push!(mat_id,Meta.parse(array[1]))
            push!(e1,Meta.parse(array[2]))
            push!(e2,Meta.parse(array[3]))
            push!(g12,Meta.parse(array[4]))
            push!(nu12,Meta.parse(array[5]))
            push!(rho,Meta.parse(array[6]))
            if length(array) >= 7
                cname = array[7]
                for i = 8:length(array)
                    cname = string(cname," ",array[i])
                end
                cname = replace(cname,"("=>"")
                cname = replace(cname,")"=>"")
                push!(name,cname)
            else
                push!(name,"")
            end
        end
        order = sortperm(mat_id)
        e1 = e1[order]
        e2 = e2[order]
        g12 = g12[order]
        nu12 = nu12[order]
        rho = rho[order]
        name = name[order]
        return e1,e2,g12,nu12,rho,name
    end
end
