import sampling as sampling
import simulationsalgorithms as simalgo
import numpy as np
import matplotlib.pyplot as plt
import math as math
import plotsimulations as plot_simalgo
import plotsamples as plot_samples
import ode_solver_plot as ode_plot
import ode_solver
import diffiusion_simulation_algo as diffucsion_algo

def rejection_sampling():
    x_samples = []
    y_samples = []

    uniform_dist = sampling.ProbabilityDistribution(lambda x: 6, (0, 100))

    def to_sample_func(x):
        return 6 * np.exp(-x)

    for i in range(0, 100):
        x_point, y_point = sampling.rejection_sampling(uniform_dist, to_sample_func)
        x_samples.append(x_point)
        y_samples.append(y_point)

    # Plot the results
    plot_samples.plot_function_and_samples(to_sample_func, uniform_dist, x_samples, y_samples, (0, 30))


def rejection_sampling_poisson():
    x_samples = []
    y_samples = []

    max = (10 ** 10 * np.exp(-10)) / math.factorial(int(10))

    uniform_dist = sampling.ProbabilityDistributionDiscrete(lambda x: max, (0, 3 * max))
    for i in range(0, 100):
        x_point, y_point = sampling.get_random_from_poisson(10)
        x_samples.append(x_point)
        y_samples.append(y_point)

    # Plot the results
    plot_samples.plot_function_and_samples(lambda x: (10 ** x * np.exp(-10)) / math.factorial(int(x)), uniform_dist,
                                       x_samples, y_samples, (0, 30))


def rejection_sampling_gaussian():
    x_samples = []
    y_samples = []

    max = 2

    uniform_dist = sampling.ProbabilityDistributionDiscrete(lambda x: max, (0, 3 * max))
    for i in range(0, 100):
        x_point, y_point = sampling.get_random_from_gaussian(10,5)
        x_samples.append(x_point)
        y_samples.append(y_point)

    # Plot the results

    plot_samples.plot_function_and_samples(lambda x: 1 / np.sqrt(2 * np.pi) * 5 * np.exp(-1 / 2 * ((x - 10) / 5) ** 2), uniform_dist,
                                           x_samples, y_samples, (0, 30))
def michaelis_menten_model_frm(init_quantity, max_time):
    system = simalgo.MolecularSystem(
        states=["S", "E", "ES", "P"],
        init_quantity=init_quantity,
        reactions={
            "R1": simalgo.MolecularSystem.Reaction([-1, -1, 1, 0], 0.025, lambda x: x[0] * x[1]),
            "R2": simalgo.MolecularSystem.Reaction([1, 1, -1, 0], 0.1, lambda x: x[2]),
            "R3": simalgo.MolecularSystem.Reaction([0, 1, -1, 1], 5, lambda x: x[2])
        }
    )

    # Printing the system and reactions for debugging
    print("Molecular System:")
    print(f"States: {system.states}")
    print(f"Initial Quantities: {system.state_quantity}")
    print("Reactions:")
    for reaction in system.reactions.values():
        print(reaction)

    x_points, y_points = simalgo.euler_maragama_molecular_system(system, 0.05, max_time)
    plot_simalgo.plot_system_variation(x_points, y_points)


def michaelis_menten_model_ssa_dm(init_quantity, max_time):
    system = simalgo.MolecularSystem(
        states=["S", "E", "ES", "P"],
        init_quantity=init_quantity,
        reactions={
            "R1": simalgo.MolecularSystem.Reaction([-1, -1, 1, 0], 0.025, lambda x: x[0] * x[1]),
            "R2": simalgo.MolecularSystem.Reaction([1, 1, -1, 0], 0.1, lambda x: x[2]),
            "R3": simalgo.MolecularSystem.Reaction([0, 1, -1, 1], 5, lambda x: x[2])
        }
    )

    # Printing the system and reactions for debugging
    print("Molecular System:")
    print(f"States: {system.states}")
    print(f"Initial Quantities: {system.state_quantity}")
    print("Reactions:")
    for reaction in system.reactions.values():
        print(reaction)

    x_points, y_points = simalgo.gillespie_algorithm_dm(system, max_time)
    plot_simalgo.plot_system_variation(x_points, y_points)


def lotka_volterra_reaction_plot(y):
    # Extract x-axis (values at position [1] in each sublist)
    x_values = [sublist[1] for sublist in y]

    # Extract y-axis (values at position [2] in each sublist)
    y_values = [sublist[2] for sublist in y]

    # Plot the data
    plt.plot(x_values, y_values, marker='x', markersize=1)

    # Add labels and title
    plt.xlabel("Values at position [1]")
    plt.ylabel("Values at position [2]")
    plt.title("Plot of values from position [1] vs position [2]")

    # Show the plot
    plt.show()


def lotka_volterra_reaction(init_quantity, max_time):
    system = simalgo.MolecularSystem(
        states=["x", "y1", "y2"],
        init_quantity=init_quantity,
        reactions={
            "R1": simalgo.MolecularSystem.Reaction([0, 1, 1], 10, lambda x: x[0] * x[1]),
            "R2": simalgo.MolecularSystem.Reaction([0, -1, 1], 0.01, lambda x: x[2] * x[1]),
            "R3": simalgo.MolecularSystem.Reaction([0, 0, -1], 10, lambda x: x[2])
        }
    )

    # Printing the system and reactions for debugging
    print("Molecular System:")
    print(f"States: {system.states}")
    print(f"Initial Quantities: {system.state_quantity}")
    print("Reactions:")
    for reaction in system.reactions.values():
        print(reaction)

    x_points, y_points = simalgo.gillespie_algorithm_frm(system, max_time)
    lotka_volterra_reaction_plot(y_points)
    # plot_simalgo.plotSystemVariation(x_points,y_points)


# ------------- MAIN -------------

#lotka_volterra_reaction([1,1000,0,0],10)
def f_y(t, y):
    return np.sin(y)


y0 = 1  # Condizione iniziale
dt = 0.1  # Passo temporale
T = 5  # Tempo totale di simulazione
steps = int(T / dt)

# Chiamata al metodo di Eulero
#x_values, y_values = ode_solver.runge_kutta_rk4(f_y,0 ,y0, dt, steps)
#ode_plot.plot_eulero_method(x_values,y_values)

#plot_simalgo.plot_eulero_method(x_values, y_values)
#rejection_sampling_gaussian()
# lotka_volterra_reaction([1,1000,1000],200)
#rejection_sampling_poisson()

state1 = diffucsion_algo.State(1,0.2)
state2 = diffucsion_algo.State(2,0.1)
state3 = diffucsion_algo.State(3,0.01)
state4 = diffucsion_algo.State(4,0.2)

states = diffucsion_algo.States([state1,state2,state3,state4],[100,20,8,40])
states2 = diffucsion_algo.States([state1,state2,state3,state4],[41,70,8,50])
states3 = diffucsion_algo.States([state1,state2,state3,state4],[12,12,1,1])
states4 = diffucsion_algo.States([state1,state2,state3,state4],[0,0,0,0])

subsystem1 = diffucsion_algo.Subsystem(1,10,states)
subsystem2 = diffucsion_algo.Subsystem(2,10,states2)
subsystem3 = diffucsion_algo.Subsystem(3,10,states3)
subsystem4 = diffucsion_algo.Subsystem(4,10,states4)

subsystem1.add_neighbors([subsystem2,subsystem3])
subsystem2.add_neighbors([subsystem1,subsystem4])
subsystem3.add_neighbors([subsystem1,subsystem4])
subsystem4.add_neighbors([subsystem2,subsystem3])


r1 = diffucsion_algo.Reaction(1,[-10,1,1,0],0.2,[state1],lambda x: 10 * x[0])
r2 = diffucsion_algo.Reaction(2,[10,0,0,-5],0.11,[state4],lambda x:x[3])
r3 = diffucsion_algo.Reaction(3,[0,-1,1,-1],0.3,[state2,state4],lambda x:x[2] * x[3])
r4 = diffucsion_algo.Reaction(4,[-3,-4,1,1],0.34,[state1,state2],lambda x: 3 * x[0] * 4 * x[1])

connectivity_matrix = diffucsion_algo.ConnectivityMatrix([subsystem1,subsystem2,subsystem3,subsystem4,subsystem4],[r1,r2,r3,r4])

system = diffucsion_algo.DynamicSystem(connectivity_matrix=connectivity_matrix)

diffucsion_algo.plot_quantities(diffucsion_algo.gillespie_nsm(1,system),system)