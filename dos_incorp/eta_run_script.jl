using Pkg; Pkg.activate("/ocean/projects/cts180021p/mbabar/PhD/Gerischer/")

loc = "/ocean/projects/cts180021p/mbabar/PhD/Gerischer/ElectrochemicalKinetics.jl/src/"

using Distributed
using SharedArrays
using ClusterManagers

addprocs(128; exeflags="--project=/ocean/projects/cts180021p/mbabar/PhD/Gerischer/Project.toml")
@everywhere push!(LOAD_PATH, $loc)

# Load libs
@everywhere begin
    using Pkg; Pkg.activate("/ocean/projects/cts180021p/mbabar/PhD/Gerischer/")
    using MAT
    using CSV
    using DelimitedFiles
    using Glob
    using Interpolations
    using QuadGK
    import YAML
    using ElectrochemicalKinetics
    #calculate_Vdl_interp = ElectrochemicalKinetics.calculate_Vdl_interp
end

#include(loc*"rate_constant.jl")
#export compute_k, compute_k_cq

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

#Prefactor
@everywhere function prefactor(
    mhcd::MarcusHushChidseyDOS;
    kT::Real = 0.026,
    V_q = 0.0,
)
    #fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    Ft(E) = (1/4kT).*((sech.(E./2kT)).^2) # Thermal broadening function
    E -> ((mhcd.dos.interp_func.(E .+ V_q)).^ 1) .* Ft(E)
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

# Load angles
file = matopen("../../co_twist/q_dict.mat")
q12_list = read(file, "q12");
q23_list = read(file, "q23");
q12_list = convert(SharedArray,collect(q12_list))
q23_list = convert(SharedArray,collect(q23_list))

# Load matlab dos data
dir = "../../co_twist/";
#end_str = "_kcut_6.9282_qtype_1_nq_484_zip.mat";
#file = matopen(dir*"dos_q12_"*string(theta12)*"_q23_"*string(theta23)*end_str)
#dos_tot = read(file, "dos_tot");
#Elist = read(file, "E_list");

# Get rate k
#dos = [transpose(Elist) transpose(dos_tot)]
## Define model
Vq_max = 0.6
Vq_min = -0.6
C_dl = 10.0
Eo = 0.0
λ = 0.82
A = 1.0
kT = 0.026 #eV
#η_list = [-0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
η_list = LinRange(-0.3, 0.3, 61)
η_list = convert(SharedArray,collect(η_list))

kox_data = zeros(Float64, size(q12_list));
kred_data = zeros(Float64, size(q12_list));
kox_data = convert(SharedArray,kox_data)
kred_data = convert(SharedArray,kred_data)

for η in η_list
    @time begin
    print("\n","η = ",η,"\n")
    @sync @distributed for i in 1:size(q12_list)[2]
    	print(i,",")
        local file, theta12, theta23 
        theta12 = chop_str(string(q12_list[i]))
        theta23 = chop_str(string(q23_list[i]))
        mat_file = glob(dir*"dos_q12_"*theta12*"_q23_"*theta23*"*")[1]
        file = matopen(mat_file)
        dos_tot = read(file, "dos_tot");
        Elist = read(file, "E_list");
        dos = [transpose(Elist) transpose(dos_tot)]
        mhcd = MarcusHushChidseyDOS(A, λ, dos, Ef=0.0) # Assuming Ef=0 (centered)
        # Compute rates
        kox_data[i],V_q = compute_k_cq(η, mhcd, true; Eo=Eo, Vq_min=Vq_min, Vq_max=Vq_max);
        kred_data[i],V_q = compute_k_cq(η, mhcd, false; Eo=Eo, Vq_min=Vq_min, Vq_max=Vq_max);
        print(kox_data[i],",",kred_data[i],"\n")
    end
    end
    # Write rate data as .mat
    file = matopen("k_data_"*string(A)*"_λ_"*string(λ)*"_η_"*string(round(η,digits=2))*".mat", "w")
    write(file, "kox_list", collect(kox_data))
    write(file, "kred_list", collect(kred_data))
    write(file, "q12_list", collect(q12_list))
    write(file, "q23_list", collect(q23_list))
    #write(file, "prefactor", factor)
    close(file)
end

##
#@time begin
#     mhcd = MarcusHushChidseyDOS(20.0, 0.82, dos)
#     k = compute_k(0.4, mhcd; calc_cq=true, Vq_min=-0.45, Vq_max=0.45)
#end


