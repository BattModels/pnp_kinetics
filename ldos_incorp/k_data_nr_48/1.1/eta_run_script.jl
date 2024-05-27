using Pkg; Pkg.activate("/ocean/projects/cts180021p/mbabar/PhD/Gerischer/")
loc = "/ocean/projects/cts180021p/mbabar/PhD/Gerischer/ElectrochemicalKinetics.jl/src/"

using Distributed
using SharedArrays
using ClusterManagers

addprocs(61; exeflags="--project=/ocean/projects/cts180021p/mbabar/PhD/Gerischer/Project.toml")
@everywhere push!(LOAD_PATH, $loc)

# Load libs
@everywhere begin
    loc = LOAD_PATH[end]
    using Pkg; Pkg.activate("/ocean/projects/cts180021p/mbabar/PhD/Gerischer/")
    using MAT, CSV, DelimitedFiles, Glob, QuadGK, Interpolations
    import YAML
    include(loc*"ElectrochemicalKinetics.jl")
    MarcusHushChidseyDOS = ElectrochemicalKinetics.MarcusHushChidseyDOS
    calculate_Vdl_interp = ElectrochemicalKinetics.calculate_Vdl_interp
    fermi_dirac = ElectrochemicalKinetics.fermi_dirac
end

# Load functions
@everywhere function integrand(
    mhcd::MarcusHushChidseyDOS,
    ox::Bool;
    kT::Real = 0.026,
    η::Real = 0.0,
    V_q::Real = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(( E .+ mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        else
            exp_arg = -(( E .- mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        end
        exp.(exp_arg)
    end
    fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    E -> mhcd.A .* ((mhcd.dos.interp_func.(E .+ V_q)).^ 1) .* marcus_term(E .+ η) .* fd(E)
end

@everywhere function compute_k_cq(
    η,
    model::MarcusHushChidseyDOS,
    ox::Bool;
    Eo = -0.07, # E_f,red (solvent) - E_f,vac (bilayer)
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kT = 0.026,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    #Vappl_data, Vdl_data = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    #v_interp = LinearInterpolation(Vappl_data, Vdl_data)
    V_t = (Eo + η)
    V_dl = V_dl_interp(V_t)
    
    #V_dl = v_interp(V_t)
    ## Snippet to handle integration limit crossing dos energy range
    V_q = V_t - V_dl
    if V_q < 0
        E_max = E_max .- 0.05
        E_min = E_min .- V_q .+ 0.05
    elseif V_q > 0
        E_min = E_min .+ 0.05
        E_max = E_max  .- V_q .- 0.05
    end
    #print(E_min, E_max)
    k_rate = quadgk(integrand(model, ox; kT = kT, η = η, V_q = V_q), E_min, E_max)[1]
    return k_rate, V_q
end

# Modify string
@everywhere function chop_str(str::String)
         while str[length(str)] == '0'
               str = chop(str)
         end
         if str[length(str)] == '.'
            str = chop(str)
         end
         return str
end

# Load config data
config_loc = "config.yml" # ARGS[1]
constants = YAML.load_file(config_loc)

# Load matlab dos data
file = matopen(constants["dos_file"])
nsamps = floor(Int, read(file, "nsamps"));
rscx = read(file, "rscx") 
rscy = read(file, "rscy")
theta = read(file, "theta");
E_list = read(file, "E_list");
tdos = read(file, "dos");
E_f = read(file, "E_f");
print("Fermi level of system = ", E_f," eV\n")
ldos = read(file, "data");
flat_ldos = reshape(ldos, length(E_list), floor(Int, nsamps*nsamps))
if size(E_list)[1]==1 
    E_list = transpose(E_list)    
end

# Get rate k
Vq_max = 0.6
Vq_min = -0.6
C_dl = constants["C_dl"]
lambda = constants["lambda"]
A = constants["MHC_prefactor"]
Eo_rel = constants["Eo_redoxcouple"] - constants["Eo_electrode"] # eV, Relative Eo between Ruhex and tBLG
#η = constants["Vapp"] # eV, absolute overpotential, not wrt AgCl
λ = constants["lambda"] #eV
kT = constants["kbT"] #eV

η_list = LinRange(-0.5, 0.3, 41)
η_list = convert(SharedArray,collect(η_list))
location = ""

for η in η_list
    kox_data = zeros(Float64, (1,nsamps*nsamps))
    kred_data = zeros(Float64, (1,nsamps*nsamps))
    Vdl_data = zeros(Float64, (1,nsamps*nsamps))
    kox_data = convert(SharedArray,kox_data)
    kred_data = convert(SharedArray,kred_data)
    Vdl_data = convert(SharedArray,Vdl_data)
    @time begin
    print("\n","η = ",η,"\n")
    @sync @distributed for i in 1:size(flat_ldos)[2]
    	print(i," out of ",nsamps*nsamps,"\n")
        local file, dos_list, dos_data
        dos_list = flat_ldos[:,i]
        dos_data = [E_list dos_list]

        mhcd = MarcusHushChidseyDOS(A, λ, dos_data, Ef=E_f) # or assuming Ef=0 (centered)
        # Compute rates
        kred_data[i],V_q = compute_k_cq(η, mhcd, true; Eo=Eo_rel, Vq_min=Vq_min, Vq_max=Vq_max);
        kox_data[i],V_q = compute_k_cq(η, mhcd, false; Eo=Eo_rel, Vq_min=Vq_min, Vq_max=Vq_max);
        Vdl_data[i] = (Eo_rel + η) - V_q
        print("η = ",η, ", kox = ", kox_data[i],", kred = ",kred_data[i],", Vdl = ",Vdl_data[i], "\n")
    end
    end
    # Write rate data as .mat
    out_file = matopen(location*"k_data_Vapp_"*string(round(η,digits=3))*".mat", "w")
    write(out_file, "kox_data", reshape(collect(kox_data), nsamps, nsamps))
    write(out_file, "kred_data", reshape(collect(kred_data), nsamps, nsamps))
    write(out_file, "Vdl_data", reshape(collect(Vdl_data), nsamps, nsamps))
    write(out_file, "Vapp", η)
    write(out_file, "rscx", rscx)
    write(out_file, "rscy", rscy)
    write(out_file, "theta", theta)
    close(out_file)
end

##
#@time begin
#     mhcd = MarcusHushChidseyDOS(20.0, 0.82, dos)
#     k = compute_k(0.4, mhcd; calc_cq=true, Vq_min=-0.45, Vq_max=0.45)
#end


