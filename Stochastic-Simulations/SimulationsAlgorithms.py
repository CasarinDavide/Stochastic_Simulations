import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import Sampling as sampling
class MolecularSystem:
    class Reaction:
        def __init__(self, vector_change, stochastic_constant, h_x=None):
            self.vector_change = vector_change
            self.stochastic_constant = stochastic_constant
            self.h_x = h_x if h_x is not None else self.create_h_x(vector_change)

        def create_h_x(self, vector_change):
            # Example implementation: a simple scaling function based on the vector change
            return lambda x: sum(v * q for v, q in zip(vector_change, x))

        def apply_h_x(self, quantities):
            return self.h_x(quantities)

        def a_x(self,):
            return self.h_x(self.vector_change) * self.stochastic_constant

        def get_change_vector(self):
            return self.vector_change

        def __str__(self):
            return f"Reaction(vector_change={self.vector_change}, stochastic_constant={self.stochastic_constant})"

    def __init__(self, states, init_quantity, reactions):
        self.states = states
        self.state_quantity = init_quantity
        self.reactions = reactions

    def update_init_state(self, vector_change_state):
        if len(vector_change_state) <= len(self.state_quantity):
            for i, x_i in enumerate(vector_change_state):
                if self.state_quantity[i] + x_i < 0:
                    return False  # Not enough quantity to apply change

            for i, x_i in enumerate(vector_change_state):
                self.state_quantity[i] += x_i  # Update quantities

            return True  # Successfully updated state
        else:
            print("Error: Length of vector change does not match state quantities.")
            return False

# SSA direct method
def Gillespie_algorithmDM(mol_system: MolecularSystem, max_time_elapsed):

    time_elapsed = 0

    x_points = []
    y_points_for_states = []

    while time_elapsed < max_time_elapsed:
        a_0 = 0
        a_0_weighted = []

        # Calculate total propensity and weighted propensities
        for reaction in mol_system.reactions.values():
            a_0 += reaction.a_x()  # Total propensity
            pre_sum = a_0_weighted[-1] if a_0_weighted else 0
            a_0_weighted.append(pre_sum + reaction.a_x())


        if a_0 == 0:  # No reactions possible
            print("No possible reactions.")
            break

        # Time to next reaction
        dt = (1 / a_0) * np.log(1 / rnd.random())

        # Select the next reaction based on random choice
        random_2 = rnd.random() * a_0
        next_picked_reaction = None

        for i, reaction_weighted in enumerate(a_0_weighted):
            if random_2 < reaction_weighted:
                next_picked_reaction = mol_system.reactions[list(mol_system.reactions.keys())[i]]
                break

        # Update the state of the system based on the selected reaction
        possible_state = mol_system.update_init_state(next_picked_reaction.get_change_vector())

        time_elapsed += dt  # Update elapsed time

        if possible_state:
            print("State changed successfully.")
            x_points.append(time_elapsed)
            y_points_for_states.append(list.copy(system.state_quantity))
        else:
            print("Time elapsed without possible state change.")

    return x_points,y_points_for_states
def Gillespie_algorithmFRM(mol_system: MolecularSystem, max_time_elapsed):

    time_elapsed = 0
    x_points = []
    y_points_for_states = []

    while time_elapsed < max_time_elapsed:

        next_picked_reaction = None
        min_dt = np.inf

        for reaction in mol_system.reactions.values():
            a_xi = reaction.a_x()
            dt = np.inf
            if a_xi > 0:
                dt = (1 / a_xi) * np.log(1/ rnd.random())

            if min_dt > dt:
                min_dt = dt
                next_picked_reaction = reaction

        if min_dt == np.inf:
            return

        # Update the state of the system based on the selected reaction
        possible_state = mol_system.update_init_state(next_picked_reaction.get_change_vector())

        time_elapsed += dt  # Update elapsed time

        if possible_state:
            print("State changed successfully.")
            x_points.append(time_elapsed)
            y_points_for_states.append(list.copy(mol_system.state_quantity))
        else:
            print("Time elapsed without possible state change.")

    return x_points,y_points_for_states
def plotGillipieSystemVariation(x_points,y_points):
    # Plotting the results
    plt.figure(figsize=(10, 6))

    for i, y in enumerate(y_points):
        plt.plot(x_points, y_points, marker='o', linestyle='-', label=f'Series {i + 1}')


    plt.title('Molecular System Simulation')
    plt.xlabel('Time')
    plt.ylabel('Quantity')
    plt.grid()
    plt.show()

def tau_leaping(mol_system:MolecularSystem,dt,max_time_elapsed,max_percentual_changes = 0.3):

    time_elapsed = 0
    x_points = []
    y_points_for_states = []

    while time_elapsed < max_time_elapsed:
        a_0 = 0
        a_0_weighted = []

        # Calculate total propensity and weighted propensities
        for reaction in mol_system.reactions.values():
            a_0 += reaction.a_x()  # Total propensity
            pre_sum = a_0_weighted[-1] if a_0_weighted else 0
            a_0_weighted.append(pre_sum + reaction.a_x())


        if a_0 == 0:  # No reactions possible
            print("No possible reactions.")
            break

        # Time to next reaction
        dt = (1 / a_0) * np.log(1 / rnd.random())

        # Select the next reaction based on random choice
        random_2 = rnd.random() * a_0
        next_picked_reaction = None

        for i, reaction_weighted in enumerate(a_0_weighted):
            if random_2 < reaction_weighted:
                next_picked_reaction = mol_system.reactions[list(mol_system.reactions.keys())[i]]
                break

        # Update the state of the system based on the selected reaction
        possible_state = mol_system.update_init_state(next_picked_reaction.get_change_vector())

        time_elapsed += dt  # Update elapsed time

        if possible_state:
            print("State changed successfully.")
            x_points.append(time_elapsed)
            y_points_for_states.append(list.copy(system.state_quantity))
        else:
            print("Time elapsed without possible state change.")

    return x_points,y_points_for_states



