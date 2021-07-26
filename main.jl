using ColorSchemes: thermal
using Plots: push!
# using Core: Vector
using LinearAlgebra
using Statistics
using Plots
using ColorSchemes

gr()

Boltzman_constant = 1.380649e-23
Nx = 10
Ny = 10
NTOTAL = Nx * Ny
STEPS = 40000
REFERENCE_RANGE = 2 * Int(STEPS/5)

Thermal_list = Vector(0.001:0.03:0.901)

function cal_energy(spin)
    _energy = 0
    _energy += LinearAlgebra.dot(spin, circshift(spin, -1))
    _energy += LinearAlgebra.dot(spin, circshift(spin, (0, -1)))
    return _energy
end

function create_ones_state(Nx, Ny)
    return ones(Int8, (Nx, Ny))
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

function plot_quanity(thermal, quanity, name)
    p = plot(
        thermal,
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

    for thermal in Thermal_list
        energy = - cal_energy(spin)
        magnet = sum(spin)
    
        energy_list = Vector()
        magnet_list = Vector()
        push!(energy_list, energy/NTOTAL)
        push!(magnet_list, magnet/NTOTAL)
    
        for k in (1:STEPS)
            i = rand(1:Nx)
            j = rand(1:Ny)
    
            spin_trial = copy(spin)
            spin_trial[i,j] = -1 * spin[i,j]
    
            i == Nx ? s_top = spin[1,j] : s_top = spin[i,j]
            i == 1  ? s_bottom = spin[Nx,j] : s_bottom = spin[i,j]
            j == Ny ? s_right = spin[i,1] : s_right = spin[i,j]
            j == 1 ? s_left = spin[i,Nx] : s_left = spin[i,j]
    
            δE = -2 * spin_trial[i,j] * (s_top + s_bottom + s_right + s_left)
    
            energy_trial = energy + δE
            if energy_trial < energy
                spin = spin_trial
                energy = energy_trial
                magnet = sum(spin)
            #elseif rand(Float64) < exp(-δE/Boltzman_constant*thermal)
            elseif rand(Float64) < exp(-δE/thermal)
                spin = spin_trial
                energy = energy_trial
                magnet = sum(spin)
            end
    
            push!(energy_list, energy/NTOTAL)
            push!(magnet_list, magnet/NTOTAL)
        end
    
        datarange = energy_list[REFERENCE_RANGE:end]
        push!(Energy_list, Statistics.mean(datarange))
        push!(Magnet_list, Statistics.mean(magnet_list[REFERENCE_RANGE:end]))
        push!(Capacity_list, 
            (
            Statistics.mean(datarange.*datarange)
            - Statistics.mean(datarange)^2
            #)/ (thermal * Boltzman_constant)^2
            )/thermal^2
        )
        push!(Spin_snapshot, spin)
    end

    energy_label = "Energy"
    magnet_label = "Magnet"
    capacity_label = "HeatCapacity"

    plot_quanity(Thermal_list, Energy_list, energy_label)
    plot_quanity(Thermal_list, Magnet_list, magnet_label)
    plot_quanity(Thermal_list, Capacity_list, capacity_label)

    anim = Animation()
    for (thermal, snapshot) in zip(Thermal_list, Spin_snapshot)
        plt = heatmap(snapshot, aspect_ratio = 1, c=:grays, title = "T = $(thermal)")
        frame(anim, plt)
    end

    gif(anim, "10fps_with_thermal.gif", fps = 10)
end

main()
