# problem with eulero method: it's the wrost for approximating function, if dt -> 0 then we can have a good results otherwise
# we are overstimating that functio
def eulero_method(f_y: callable, y0, dt, simulation_time):
    # Inizializzazione della lista per memorizzare i risultati
    y_i = []
    x_i = []
    # Condizione iniziale
    y_i.append(y0)
    x_i.append(0)

    # Ciclo sul numero di passi temporali
    for i in range(1, simulation_time):
        # Calcola il prossimo valore usando il metodo di Eulero
        t = i * dt  # Tempo attuale
        y_prev = y_i[-1]  # Ultimo valore di y
        dy = f_y(t, y_prev) * dt  # Derivata moltiplicata per dt
        y_new = y_prev + dy  # Aggiorna il valore di y
        y_i.append(y_new)
        x_i.append(i)

    return x_i, y_i


# better than eulero's method, but we can improve it
# hit, what would happen if we can estimate the error of our approximation?
# if our Err(solution) >> epsilon then we can try to estimate using dt/2
# usually we can use runge kutta febherg
def runge_kutta_rk4(f_y: callable, t_0, y_0, dt, simulation_time):
    y_i = []
    x_i = []

    dt_mean = dt / 2

    while t_0 < simulation_time:
        k1 = f_y(t_0, y_0)
        k2 = f_y(t_0 + dt_mean, y_0) + k1 * dt_mean
        k3 = f_y(t_0 + dt_mean, y_0) + k2 * dt_mean
        k4 = f_y(t_0 + dt, y_0 + dt) + k3 * dt
        y_new = y_0 + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4) * dt
        t_0 = t_0 + dt
        y_0 = y_new
        y_i.append(y_new)
        x_i.append(t_0)

    return x_i, y_i


# real world problem necessity way more accuracy so dt can be that small that can not even fit on 64bit
# we use implicit method for simulating ode


