using Base
using ColorSchemes
using Plots
using LinearAlgebra
using Statistics

# config
Boltzman_constant = 1.0
N = 20
REPLICA_NUMBERS = 11
REPLICA_STEPS = 10000
INIT_RELAX_STEPS = 10000

MIN_TEMP = 0.001
MAX_TEMP = 6.001

METHOD = "REM"

# initialize
Nx = N
Ny = N
NTOTAL = Nx * Ny
REFERENCE_RANGE = 2 * Int(STEPS/5)
REPLICA_TEMPERTURE_LIST = (MIN_TEMP:(MAX_TEMP-MIN_TEMP)/(REPLICA_NUMBERS-1):MAX_TEMP)

Energy(s) = LinearAlgebra.dot(s, circshift(s, -1)) + LinearAlgebra.dot(s, circshift(s, (0, -1)))
Magnet(s) = sum(s)
Capacity(E_list, thermal) = (Statistics.mean(E_list.*E_list) - Statistics.mean(E_list)^2) / (Boltzman_constant * thermal)^2
deltaE(s_trial, i, j, s_t, s_b, s_r, s_l) = -2 * s_trial[i,j] * (s_t + s_b + s_r + s_l)

function flip(s, thermal)
    i = rand(1:Nx)
    j = rand(1:Ny)

    spin_trial = copy(s)
    spin_trial[i,j] = -1 * s[i,j]

    i == Nx ? s_t = s[1,j]  : s_t = s[i,j]
    i == 1  ? s_b = s[Nx,j] : s_b = s[i,j]
    j == Ny ? s_r = s[i,1]  : s_r = s[i,j]
    j == 1  ? s_l = s[i,Nx] : s_l = s[i,j]

    δE = deltaE(spin_trial, i, j, s_t, s_b, s_r, s_l)
    energy = Energy(s)
    energy_trial = energy + δE

    if energy_trial < energy
        return spin_trial, energy_trial
    elseif rand(Float64) < exp(-δE/thermal)
        return spin_trial, energy_trial
    else
        return s, energy
    end 
end

function create_random_state(Nx, Ny)
    spin = ones(Int8, (Nx, Ny))
    while count(x -> x == 1, spin) != Int(Nx*Ny/2)
        if count(x -> x == 1, spin) > Int(Nx*Ny/2)
            spin[rand(1:Nx),rand(1:Ny)] = -1
        else
            spin[rand(1:Nx),rand(1:Ny)] = 1
        end
    end
    return spin
end

function get_dist_name()
    name = "method_$(METHOD)-N_$(N)-STEPS_$(STEPS)-TEMPSTEPS_$(STEP_TEMP)"
    return name
end

function plot_quanity(quanity, name)
    p = plot(
        THERMAL_LIST,
        quanity,
        xlabel = "Thermal[K]",
        ylabel = name
    )
    savefig("$(name).png")
    current()
end

function main()
    Energy_list = []
    Magnet_list = []
    Capacity_list = []
    Spin_snapshot = []

    spin = create_random_state(Nx, Ny)
    
    for thermal in THERMAL_LIST
        energy = Energy(spin)
        magnet = Magnet(spin)

        energy_list = []
        magnet_list = []
        
        for k in (1:STEPS)
            spin, energy = flip(spin, thermal)
            magnet = Magnet(spin)
            push!(energy_list, energy/NTOTAL)
            push!(magnet_list, magnet/NTOTAL)
        end

        push!(Spin_snapshot, spin)
        push!(Energy_list, Statistics.mean(energy_list[REFERENCE_RANGE:end]))
        push!(Magnet_list, Statistics.mean(magnet_list[REFERENCE_RANGE:end]))
        push!(Capacity_list, Capacity(energy_list[REFERENCE_RANGE:end], thermal))
    end

    # save result systems
    dist_name = get_dist_name()
    Base.Filesystem.mkdir(dist_name)
    Base.Filesystem.cd(dist_name)

    energy_label = "Energy"
    magnet_label = "Magnet"
    capacity_label = "HeatCapacity"

    plot_quanity(Energy_list, energy_label)
    plot_quanity(Magnet_list, magnet_label)
    plot_quanity(Capacity_list, capacity_label)

    anim = Animation()
    for (thermal, snapshot) in zip(THERMAL_LIST, Spin_snapshot)
        plt = heatmap(snapshot, aspect_ratio = 1, c=:grays, title = "T = $(thermal)")
        frame(anim, plt)
    end
    gif(anim, "10fps_with_thermal.gif", fps = 10)

    # back to project root directory
    Base.Filesystem.cd("..")
end

main()
