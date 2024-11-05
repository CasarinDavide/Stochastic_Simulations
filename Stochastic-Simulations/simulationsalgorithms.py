import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import sampling as sampling


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

        def a_x(self, quantities):
            return self.apply_h_x(quantities) * self.stochastic_constant

        def get_change_vector(self):
            return self.vector_change

        def try_to_fire(self, reaction_fired_number, state_quantity, epsilon: float = 0.03):

            tau_quantities = list.copy(state_quantity)

            tau_quantities = tau_quantities + [x * reaction_fired_number for x in self.vector_change]

            a0_init = self.a_x(state_quantity)
            a0_t_tau = self.a_x(tau_quantities)

            if a0_init == 0:
                return True

            if np.abs(a0_init - a0_t_tau) / a0_init > epsilon:
                return False
            else:
                return True

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
def gillespie_algorithm_dm(mol_system: MolecularSystem, max_time_elapsed):
    time_elapsed = 0

    x_points = []
    y_points_for_states = []

    while time_elapsed < max_time_elapsed:
        a_0 = 0
        a_0_weighted = []

        # Calculate total propensity and weighted propensities
        for reaction in mol_system.reactions.values():
            a_xi = reaction.a_x(mol_system.state_quantity)
            a_0 += a_xi  # Total propensity
            pre_sum = a_0_weighted[-1] if a_0_weighted else 0
            a_0_weighted.append(pre_sum + a_xi)

        if a_0 == 0:  # No reactions possible
            print("No possible reactions.")
            break

        # Time to next reaction
        dt = (1 / a_0) * np.log(1 / np.random.uniform())

        # Select the next reaction based on random choice
        random_2 = np.random.uniform() * a_0
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
            y_points_for_states.append(list.copy(mol_system.state_quantity))
        else:
            print("Time elapsed without possible state change.")

    return x_points, y_points_for_states


def gillespie_algorithm_frm(mol_system: MolecularSystem, max_time_elapsed):
    time_elapsed = 0
    x_points = []
    y_points_for_states = []

    while time_elapsed < max_time_elapsed:

        next_picked_reaction = None
        min_dt = np.inf

        for reaction in mol_system.reactions.values():
            a_xi = reaction.a_x(mol_system.state_quantity)
            dt = np.inf
            if a_xi > 0:
                dt = (1 / a_xi) * np.log(1 / np.random.uniform())

            if min_dt > dt:
                min_dt = dt
                next_picked_reaction = reaction

        if min_dt == np.inf:
            break

        # Update the state of the system based on the selected reaction
        possible_state = mol_system.update_init_state(next_picked_reaction.get_change_vector())

        time_elapsed += min_dt  # Update elapsed time

        if possible_state:
            print("State changed successfully.")
            x_points.append(time_elapsed)
            y_points_for_states.append(list.copy(mol_system.state_quantity))
        else:
            print("Time elapsed without possible state change.")

    return x_points, y_points_for_states


def tau_leaping(mol_system: MolecularSystem, tau, max_time_elapsed, epsilon=0.03):
    possible_react = True
    fired_up_reactions = []

    for reaction in mol_system.reactions.values():
        a0_i = reaction.a_x(mol_system.state_quantity)
        if a0_i == 0:
            continue

        fired_up_reaction, y = sampling.get_random_from_poisson(tau * reaction.a_x(mol_system.state_quantity))

        possible_react = possible_react and reaction.try_to_fire(fired_up_reaction, mol_system.state_quantity, epsilon)

        if not possible_react:
            break

        fired_up_reactions.append((reaction, fired_up_reaction))

    if not possible_react:
        # retry with tau/2
        tau_leaping(mol_system, tau / 2, epsilon)
    else:

        if len(fired_up_reactions) > 0:
            for reaction in fired_up_reactions:
                react, number = reaction
                mol_system.update_init_state([x * number for x in react.vector_change])
        else:
            tau = max_time_elapsed
    return tau


def simulate_tau_leaping(mol_system: MolecularSystem, tau, max_time_elapsed, epsilon=0.03):
    time_elapsed = 0
    x_points = []
    y_points_for_states = []
    while time_elapsed < max_time_elapsed:
        dt = tau_leaping(mol_system, tau, max_time_elapsed, epsilon)
        time_elapsed += dt
        x_points.append(time_elapsed)
        y_points_for_states.append(list.copy(mol_system.state_quantity))
    return x_points, y_points_for_states


# TO FIX
def euler_maragama_molecular_system(molecular_system: MolecularSystem, dt, steps):
    x_values = []
    y_values = []

    y_values.append(molecular_system.state_quantity)
    x_values.append(0)

    for i in range(1, steps):

        quantity_variance_dt = [0 for _ in molecular_system.states]

        for reaction in molecular_system.reactions.values():
            a_xi = reaction.a_x(molecular_system.state_quantity)

            for j in range(0, len(reaction.vector_change)):
                quantity_variance_dt[j] += a_xi * reaction.vector_change[j]

        for j in range(0, len(molecular_system.states)):
            epsilon, y = sampling.get_random_from_gaussian(0,1)
            molecular_system.state_quantity[j] += quantity_variance_dt[j] * np.sqrt(dt) * epsilon

        y_values.append(list.copy(molecular_system.state_quantity))
        x_values.append(i)

    return x_values, y_values


def euler_maragama(f_y: callable, epsilon: callable, y0, dt, steps):
    # Inizializzazione della lista per memorizzare i risultati
    y_i = []
    x_i = []
    # Condizione iniziale
    y_i.append(y0)
    x_i.append(0)

    # Ciclo sul numero di passi temporali
    for i in range(1, steps):
        # Calcola il prossimo valore usando il metodo di Eulero
        t = i * dt  # Tempo attuale
        y_prev = y_i[-1]  # Ultimo valore di y
        dy = f_y(t, y_prev) * dt + epsilon(t, y_prev) * np.sqrt(
            dt) * sampling.get_random_from_gaussian(0, 1,
                                                    (0, 100))  # Derivata moltiplicata per dt piu la parte stocastica
        y_new = y_prev + dy  # Aggiorna il valore di y
        y_i.append(y_new)
        x_i.append(i)

    return x_i, y_i

