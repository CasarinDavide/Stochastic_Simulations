import sampling as sampling
import simulationsalgorithms as simalgo
import numpy as np
import matplotlib.pyplot as plt
import math as math
import plotsimulations as plot_simalgo
import plotsamples as plot_samples


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

michaelis_menten_model_frm([10000,1000,0,0],1000)
    #def f_y(t, y):
#return -2 * y * t

    # Parametri del metodo



#y0 = 1  # Condizione iniziale
#dt = 0.1  # Passo temporale
#T = 5  # Tempo totale di simulazione
#steps = int(T / dt)

# Chiamata al metodo di Eulero
#x_values, y_values = simalgo.eulero_method(f_y, y0, dt, steps)
#plot_simalgo.plot_eulero_method(x_values, y_values)
#rejection_sampling_gaussian()
# lotka_volterra_reaction([1,1000,1000],200)

# rejection_sampling_poisson()
