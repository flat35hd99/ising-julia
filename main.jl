# using Core: Vector
using LinearAlgebra
using Statistics
using Plots

gr()

Boltzman_constant = 1.380649e-23
Nx = 20
Ny = 20
NTOTAL = Nx * Ny
STEPS = 40000
FOUR_FIFTHS = 4 * Int(STEPS/5)

Thermal_list = Vector(0.001:0.03:6.001)
Energy_list = []
Magnet_list = []
Capacity_list = []

function cal_energy(spin)
    _energy = 0
    Nx = size(spin, 1)
    Ny = size(spin, 2)
    _energy += LinearAlgebra.dot(spin, circshift(spin, -1))
    _energy += LinearAlgebra.dot(spin, circshift(spin, (0, -1)))
    return _energy
end

function create_initial_state(Nx, Ny)
    spin = reshape(rand(0:1, Nx*Ny),(Nx, Ny))
    while count(x -> x == 1, spin) != Int(Nx*Ny/2)
        if count(x -> x == 1, spin) > Int(Nx*Ny/2)
            spin[rand(1:Nx),rand(1:Ny)] = 0
        else
            spin[rand(1:Nx),rand(1:Ny)] = 1
        end
    end
    return spin
end

for thermal in Thermal_list
    spin = create_initial_state(Nx, Ny)
    energy = - cal_energy(spin)
    magnet = sum(spin)

    energy_list = []
    magnet_list = []
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

    datarange = energy_list[FOUR_FIFTHS:end]
    push!(Energy_list, Statistics.mean(datarange))
    push!(Magnet_list, Statistics.mean(magnet_list[FOUR_FIFTHS:end]))
    push!(Capacity_list, 
        (
        Statistics.mean(datarange.*datarange)
        - Statistics.mean(datarange)^2
        #)/ (thermal * Boltzman_constant)^2
        )/thermal^2
    )
end

thermal_label = "Thermal[K]"
energy_label = "Energy[J]"
magnet_label = "Magnetic momoment per particle"
capacity_label = "Heat capacity"

energy_plot = plot(Thermal_list, Energy_list,
    xlabel = thermal_label,
    ylabel = energy_label
    )
savefig("energy2.png")
current()

magnet_plot = plot(Thermal_list, Magnet_list,
    xlabel = thermal_label,
    ylabel = magnet_label,
    )
savefig("magnet2.png")
current()

capacity_plot = plot(Thermal_list, Capacity_list,
    xlabel = thermal_label,
    ylabel = capacity_label
    )
savefig("capacity2.png")
current()