from matplotlib import pyplot as plt


def plot_system_variation(x_points, y_points):
    # Plotting the results
    plt.figure(figsize=(10, 6))

    for i, y in enumerate(y_points):
        plt.plot(x_points, y_points, marker='o', linestyle='-', label=f'Series {i + 1}')

    plt.title('Molecular System Simulation')
    plt.xlabel('Time')
    plt.ylabel('Quantity')
    plt.grid()
    plt.show()

