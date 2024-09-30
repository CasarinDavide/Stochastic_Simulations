import random as rnd
import numpy as np
import matplotlib.pyplot as plt

class ProbilityDistribuition:

    def __init__(self,func:callable,interval: tuple):
        self.func = func
        self.interval = interval

    def getRandomX(self):
        a,b = self.interval
        return rnd.uniform(a,b)

    def getRandomY(self,x,start_codomain: int = 0 ):
        return rnd.uniform(start_codomain,self.f(x))

    def getRandomYInInterval(self,start_codomain: tuple ):
        a,b = start_codomain
        return rnd.uniform(a,b)


    def f(self,x):
        return self.func(x)
def rejectionSampling(g_proxy:ProbilityDistribuition,f:callable) -> tuple:

    while(True):

        x_sampling = g_proxy.getRandomX()

        f_xo = f(x_sampling)
        g_xo = g_proxy.f(x_sampling)

        y_sampling = g_proxy.getRandomYInInterval((0,g_xo))

        if y_sampling <= f_xo:
            return (x_sampling,y_sampling)
def plot_function_and_samples(target_function, g_proxy, x_points, y_points,interval : tuple):
    """
    Plots the target function, proxy distribution, and the sampled points from rejection sampling.

    target_function: The function f(x) to plot.
    g_proxy: The proxy distribution g(x).
    x_points: The x points sampled from rejection sampling.
    y_points: The y points corresponding to the x points.
    """
    # Create a range of x values for plotting the functions
    start , end = interval
    x_vals = np.linspace(start, end, 1000)

    # Calculate the y values for the target function and the proxy distribution
    f_vals = [target_function(x) for x in x_vals]
    g_vals = [g_proxy.f(x) for x in x_vals]

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Plot the target function f(x)
    plt.plot(x_vals, f_vals, label="Target Function f(x)", color='blue', linewidth=2)

    # Plot the proxy function g(x)
    plt.plot(x_vals, g_vals, label="Proxy Function g(x)", color='green', linestyle='--', linewidth=2)

    # Plot the sampled points
    plt.scatter(x_points, y_points, color='red', label="Sampled Points", zorder=5)

    # Add labels and title
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Rejection Sampling: Target Function and Sampled Points")
    plt.legend()
    plt.grid(True)
    plt.show()








