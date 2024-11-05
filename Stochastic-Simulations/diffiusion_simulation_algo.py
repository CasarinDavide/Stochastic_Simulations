import math
import random
import string
from dataclasses import dataclass, field
from tokenize import String
from typing import List, Callable

from heapq import heapify, heappush, heappop


import matplotlib.pyplot as plt
import numpy as np


class State:
    def __init__(self, label: str, diffusion_rate: float):
        self.label = label
        self.diffusion_rate = diffusion_rate

    def get_label(self):
        return self.label

    def get_diffusion_rate(self):
        return self.diffusion_rate


class Reaction:
    def __init__(self, id: int, vector_change: List[float], stochastic_constant: float, dependency: List[State],
                 h_x: Callable):
        self.id = id
        self.vector_change = vector_change
        self.stochastic_constant = stochastic_constant
        self.h_x = h_x
        self.dependency = dependency

    def a_x(self, quantities):
        return self.h_x(quantities) * self.stochastic_constant

    def get_change_vector(self):
        return self.vector_change


class States:
    def __init__(self, states: List[State], init_quantity: List[float]):
        self.states = states
        self.quantity = init_quantity
        self.map_position_to_label = {}

        for pos, state in enumerate(states):
            self.map_position_to_label[state.get_label()] = pos

    def get_states(self):
        return self.states

    def get_quantity(self):
        return self.quantity

    def update_quantity(self, vector_changes: List[float]):
        self.quantity = [q + v for q, v in zip(self.quantity, vector_changes)]

    def get_quantity_for_state(self, state_label: str):
        return self.quantity[self.map_position_to_label[state_label]]

    def update_quantity_for_state(self, state_label: str, quantity: float):
        self.quantity[self.map_position_to_label[state_label]] = quantity


class Subsystem:
    def __init__(self, id: int, area: float, reaction_system: States, neighbors: List['Subsystem'] = list):
        self.id = id
        self.area = area
        self.reaction_system = reaction_system
        self.neighbors = neighbors
        self.propensity_list = {}
        self.diffusion_list = {}
        self.a0 = 0
        self.s0 = 0

    def get_total_directions(self):
        return len(self.neighbors)

    def get_area(self):
        return self.area

    def get_id(self):
        return self.id

    def subtract_quantity(self, state: State, lost_quantities: float):
        actual_quantity = self.reaction_system.get_quantity_for_state(state.label)
        new_quantity = actual_quantity - lost_quantities
        self.reaction_system.update_quantity_for_state(state_label=state.label, quantity=new_quantity)

    def add_quantity(self, state: State, add_quantities: float):
        actual_quantity = self.reaction_system.get_quantity_for_state(state.label)
        new_quantity = actual_quantity + add_quantities
        self.reaction_system.update_quantity_for_state(state_label=state.label, quantity=new_quantity)

    def update_propensity_list(self, reactions: List[Reaction]):
        self.a0 = 0
        for reaction in reactions:
            rate = reaction.a_x(self.reaction_system.quantity)
            self.propensity_list[reaction.id] = (reaction, rate)
            self.a0 += rate

    def update_propensity_when_state_change(self, state: State, reactions: List[Reaction]):
        for reaction in reactions:
            if state in reaction.dependency:
                self.a0 -= self.propensity_list.get(reaction.id, 0)[1]
                rate = reaction.a_x(self.reaction_system.quantity)
                self.propensity_list[reaction.id] = (reaction, rate)
                self.a0 += rate

    def diffusion_constant(self, state: State):
        return state.diffusion_rate / self.area * self.reaction_system.get_quantity_for_state(
            state_label=state.label) * self.get_total_directions()

    def update_diffusion_list(self):
        self.s0 = 0
        for state in self.reaction_system.states:
            rate = self.diffusion_constant(state)
            self.diffusion_list[state.label] = (state, rate)
            self.s0 += rate

    def update_diffusion_when_state_change(self, state: State):

        _, last_rate = self.diffusion_list.get(state.label)
        self.s0 -= last_rate
        rate = self.diffusion_constant(state)
        self.diffusion_list[state.label] = (state, rate)
        self.s0 += rate

    def get_total(self):
        return self.a0 + self.s0

    def add_neighbors(self, neighbors: List['Subsystem']):
        self.neighbors = neighbors

    def occur_reaction_or_diffuse(self):

        total = self.get_total()
        rnd = np.random.uniform() * total

        for key, value in self.propensity_list.items():
            reaction, ai = value
            if rnd < ai:
                return (reaction, None)

        for key, value in self.diffusion_list.items():
            diffusion, si = value
            if rnd < self.a0 + si:
                pos = math.floor(rnd) % self.get_total_directions()
                return (diffusion,self.neighbors[pos] )

        return None

    def trigger_reaction(self, reaction: Reaction):
        self.reaction_system.update_quantity(reaction.get_change_vector())
        for state in reaction.dependency:
            self.update_diffusion_when_state_change(state)

    def trigger_diffuse(self, state: State, reactions: list[Reaction]):
        # random to neighbor
        self.update_diffusion_when_state_change(state)
        self.update_propensity_when_state_change(state, reactions)


@dataclass(order=True)
class Item:
    priority: int
    subsystem: Subsystem = field(compare=False)


class ConnectivityMatrix:

    def __init__(self, subsystems: List[Subsystem], possible_reactions: List[Reaction]):
        self.subsystems = subsystems
        self.possible_reactions = possible_reactions
        self.heap = []
        heapify(self.heap)

    def move_from_subsystem_to_subsystem(self, state: State, subsystem_from: Subsystem, subsystem_to: Subsystem):
        subsystem_from.subtract_quantity(state, 1)
        subsystem_to.add_quantity(state, 1)

        subsystem_from.update_propensity_when_state_change(state, self.possible_reactions)
        subsystem_from.update_diffusion_when_state_change(state)

        subsystem_to.update_propensity_when_state_change(state, self.possible_reactions)
        subsystem_to.update_diffusion_when_state_change(state)

    def update_propensities_for_subsystem(self, subsystem: Subsystem):
        subsystem.update_propensity_list(self.possible_reactions)

    def update_propensities(self):
        for subsystem in self.subsystems:
            subsystem.update_propensity_list(self.possible_reactions)

    def update_diffusion(self):
        for subsystem in self.subsystems:
            subsystem.update_diffusion_list()

    def update_diffusion_for_subsystem(self, subsystem: Subsystem):
        subsystem.update_diffusion_list()

    def generate_all_tau(self):
        for subsystem in self.subsystems:
            total = subsystem.get_total()
            tau = np.inf
            if total > 0:
                tau = 1 / subsystem.get_total() * np.log(1 / np.random.uniform())
            heappush(self.heap, Item(tau, subsystem))


class DynamicSystem:
    def __init__(self, connectivity_matrix: ConnectivityMatrix):
        self.connectivity_matrix = connectivity_matrix
        self.connectivity_matrix.update_propensities()
        self.connectivity_matrix.update_diffusion()


# next subsystem method
def gillespie_nsm(max_time, dynamic_system: DynamicSystem):
    elapsed_time = 0
    quantities_history = []  # To store the history of quantities

    while elapsed_time < max_time:
        dynamic_system.connectivity_matrix.generate_all_tau()
        item = heappop(dynamic_system.connectivity_matrix.heap)

        elapsed_time += item.priority
        subsystem = item.subsystem
        # Guardo che evento avviene
        event, parameters = subsystem.occur_reaction_or_diffuse()

        if isinstance(event, Reaction):
            subsystem.trigger_reaction(event)
        elif isinstance(event, State):
            dynamic_system.connectivity_matrix.move_from_subsystem_to_subsystem(event, subsystem, parameters)

        # Record the current quantities
        quantities_snapshot = [subsystem.reaction_system.get_quantity() for subsystem in dynamic_system.connectivity_matrix.subsystems]
        quantities_history.append(quantities_snapshot)

    return quantities_history


# Assume quantities_history is obtained from gillespie_nsm
def plot_quantities(quantities_history,dynamic_system):
    time_steps = range(len(quantities_history))
    state_labels = [state.get_label() for state in dynamic_system.connectivity_matrix.subsystems[0].reaction_system.get_states()]

    for i, state_label in enumerate(state_labels):
        plt.figure(figsize=(10, 6))
        for j, subsystem in enumerate(dynamic_system.connectivity_matrix.subsystems):
            quantities = [snapshot[j][i] for snapshot in quantities_history]
            plt.plot(time_steps, quantities, label=f'Subsystem {subsystem.get_id()} - {state_label}')

        plt.title(f'Quantity of {state_label} over time')
        plt.xlabel('Time steps')
        plt.ylabel('Quantity')
        plt.legend()
        plt.show()