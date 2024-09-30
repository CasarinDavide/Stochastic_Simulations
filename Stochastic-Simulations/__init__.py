import Sampling as sampling
import SimulationsAlgorithms as simalgo
import numpy as np


#--------- MAIN -------------

# Lists to store sampled points
x_samples = []
y_samples = []


uniformDist = sampling.ProbilityDistribuition(lambda x:6,(0,100))
toSampleFunc = lambda x: 6 * np.exp(-x)

for i in range(0,100):
    x_point, y_point = sampling.rejectionSampling(uniformDist,toSampleFunc)
    x_samples.append(x_point)
    y_samples.append(y_point)

# Plot the results
sampling.plot_function_and_samples(toSampleFunc, uniformDist, x_samples, y_samples,(0,30))



# ------------- MAIN -------------

# Running the Gillespie algorithm with a max time limit
# Example instantiation of MolecularSystem
system = simalgo.MolecularSystem(
    states=["A", "B", "C", "D"],
    init_quantity=[10,10 ,10,10],
    reactions={
        "R1": simalgo.MolecularSystem.Reaction([-2, -2, 1], 0.2),
        "R2": simalgo.MolecularSystem.Reaction([3, 0, 1], 0.45),
        "R3": simalgo.MolecularSystem.Reaction([0, -4, 0], 0.20),
        "R4": simalgo.MolecularSystem.Reaction([1, 2, 0], 0.55),
        "R4": simalgo.MolecularSystem.Reaction([1, 2, 0,2], 0.40),
        "R5": simalgo.MolecularSystem.Reaction([-1, 2, 1,-2], 0.44),
    }
)

# Printing the system and reactions for debugging
print("Molecular System:")
print(f"States: {system.states}")
print(f"Initial Quantities: {system.state_quantity}")
print("Reactions:")
for reaction in system.reactions.values():
    print(reaction)
max_time = 100

x_points, y_points = simalgo.Gillespie_algorithmFRM(system, max_time)

simalgo.plotGillipieSystemVariation(x_points,y_points)



