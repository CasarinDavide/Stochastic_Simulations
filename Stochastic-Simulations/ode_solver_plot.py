from matplotlib import pyplot as plt


def plot_eulero_method(x_values, y_values):
    # Grafico dei risultati
    plt.plot(x_values, y_values, label="Metodo di Eulero")
    plt.xlabel("Tempo")
    plt.ylabel("y(t)")
    plt.title("Soluzione EDO con Metodo di Eulero")
    plt.legend()
    plt.show()
