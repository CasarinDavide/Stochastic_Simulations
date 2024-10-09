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


def plot_eulero_method(x_values, y_values):
    # Grafico dei risultati
    plt.plot(x_values, y_values, label="Metodo di Eulero")
    plt.xlabel("Tempo")
    plt.ylabel("y(t)")
    plt.title("Soluzione EDO con Metodo di Eulero")
    plt.legend()
    plt.show()
