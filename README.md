# Stochastic Algorithm Library

Welcome to the **Stochastic Algorithm Library**, a comprehensive library that provides implementations of various stochastic algorithms. Stochastic algorithms are widely used in optimization, machine learning, and probabilistic modeling due to their ability to handle complex problems where deterministic methods might fail or be inefficient.

This library offers both classical and modern stochastic algorithms, providing users with tools to solve optimization problems, perform simulations, and analyze probabilistic systems.

## Features

- **Random Sampling Methods**:
  - Methods implemented to extract samples from widely-used probability distributions using Rejection Sampling techniques.    

- **Stochastic Simulation Algorithms**:
  - A popular approach for simulating molecular systems based on chemical reactions. This library offers two implementations of the **Gillespie Algorithm**:
    - First Reaction Method: A highly optimized variant that calculates reaction times one at a time.
    - Direct Method: An alternative that delivers equivalent results to the First Reaction Method, albeit slightly less optimized.

  - **Tau-Leaping Algortihm** for simulating molecular system based on Reactions. This algorithm simulates stochastic events in molecular systems by assuming that reaction propensities remain relatively unchanged over small time intervals. It allows users to adjust the  
 tolerance level (epsilon) for changes in propensities.
  - A numerical method for solving stochastic differential equations (SDEs), useful for modeling systems influenced by randomness.

  

## Installation

To use the **Stochastic Algorithm Library**, you can install it cloning this repository.

