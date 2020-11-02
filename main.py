import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def initialize(x_max, y_max, numx, numy):
    dx = x_max / (numx - 1)
    dy = y_max / (numy - 1)
    x = np.linspace(0, x_max, num=numx)
    y = np.linspace(0, y_max, num=numy)
    return dx, dy, x, y

def S_update(x, y, mu_s, dt, S, I, beta, N, dx, dy, numx, numy):
    S_old = S[x][y]
    s_sink = beta * I[x][y] * S[x][y] / N
    if y == 0: #Bottom No Flux BC
        S[x][y] = S[x][y+1]
        return S
    if y == numy - 1: #Top No Flux BC
        S[x][y] = S[x][y-1]
        return S
    if x == 0: #Left No FLux BC
        S[x][y] = S[x+1][y]
        return S
    if x == numx - 1: #Right No Flux BC
        S[x][y] = S[x-1][y]
        return S
    sd_up = mu_s / dy ** 2 * S[x][y + 1]
    sd_low = mu_s / dy ** 2 * S[x][y - 1]
    sd_right = mu_s / dx ** 2 * S[x + 1][y]
    sd_left = mu_s / dx ** 2 * S[x - 1][y]
    S[x][y] = S_old + dt*(sd_up + sd_low + sd_right + sd_left - s_sink)
    S[x][y] = S[x][y] - mu_s*dt*2*S_old*(1/dx**2 + 1/dy**2)
    return S

def R_update(x, y, mu_r, dt, R, I, gamma, dx, dy, numx, numy):
    R_old = R[x][y]
    r_source = gamma*I[x][y]
    if y == 0: #Bottom No Flux BC
        R[x][y] = R[x][y+1]
        return R
    if y == numy - 1: #Top No Flux BC
        R[x][y] = R[x][y-1]
        return R
    if x == 0: #Left No FLux BC
        R[x][y] = R[x+1][y]
        return R
    if x == numx - 1: #Right No Flux BC
        R[x][y] = R[x-1][y]
        return R
    rd_up = mu_r / dy ** 2 * R[x][y + 1]
    rd_low = mu_r / dy ** 2 * R[x][y - 1]
    rd_right = mu_r / dx ** 2 * R[x + 1][y]
    rd_left = mu_r / dx ** 2 * R[x - 1][y]
    R[x][y] = R_old + dt*(rd_up + rd_low + rd_right + rd_left + r_source)
    R[x][y] = R[x][y] - mu_r*dt*2*R_old*(1/(dx**2) + 1/(dy**2))
    return R

def I_update(x, y, mu_i, dt, I, S, beta, N, gamma, dx, dy, numx, numy):
    I_old = I[x][y]
    i_source = beta * I[x][y] * S[x][y] / N
    i_sink = gamma*I[x][y]
    if y == 0:  # Bottom No Flux BC
        I[x][y] = I[x][y + 1]
        return I
    if y == numy - 1:  # Top No Flux BC
        I[x][y] = I[x][y - 1]
        return I
    if x == 0:  # Left No FLux BC
        I[x][y] = I[x + 1][y]
        return I
    if x == numx - 1:  # Right No Flux BC
        I[x][y] = I[x - 1][y]
        return I
    id_up = mu_i / dy ** 2 * I[x][y + 1]
    id_low = mu_i / dy ** 2 * I[x][y - 1]
    id_right = mu_i / dx ** 2 * I[x + 1][y]
    id_left = mu_i / dx ** 2 * I[x - 1][y]
    I[x][y] = I[x][y] + dt*(id_up + id_low + id_right + id_left + i_source - i_sink)
    I[x][y] = I[x][y] - mu_i*dt*2*I_old*(1/(dx**2) + 1/(dy**2))
    return I

def solve():
    x_max = 20 #1
    y_max = 20 #1
    num_x = 11
    num_y = 11
    num_t = 100
    dt = 0.5
    dx, dy, x, y = initialize(x_max, y_max, num_x, num_y)
    mu_s = 5e-2 #0.0002
    mu_r = 5e-2 #0.0
    mu_i = 5e-2 #0.0002
    beta = 3.2/6.5
    gamma = 1/6.5
    N = 1 #5
    S = np.ones((num_x, num_y))
    S[5][5] = 0.7 #3
    I = N - S
    R = np.zeros((num_x, num_y))
    t_vec = np.linspace(0, num_t, int(num_t/dt + 1))
    S_t = np.ones((num_x, num_y,len(t_vec)))
    I_t = np.zeros((num_x, num_y,len(t_vec)))
    R_t = np.zeros((num_x, num_y, len(t_vec)))
    S_t[:, :, 0] = S
    I_t[:, :, 0] = I
    R_t[:, :, 0] = R
    for t_ind in range(0, len(t_vec)):
        t = t_vec[t_ind]
        print(t)
        for i in range(0, num_x):
            for j in range(0, num_y):
                S = S_update(i, j, mu_s, dt, S, I, beta, N, dx, dy, num_x, num_y)
                R = R_update(i, j, mu_r, dt, R, I, gamma, dx, dy, num_x, num_y)
                I = I_update(i, j, mu_i, dt, I, S, beta, N, gamma, dx, dy, num_x, num_y)
        S_t[:, :, t_ind] = S
        I_t[:, :, t_ind] = I
        R_t[:, :, t_ind] = R
    return S_t, I_t, R_t, t_vec


if __name__ == "__main__":
   S_t, I_t, R_t, t_vec  =  solve()
   print(S_t[:, :, -1].sum() + I_t[:, :,-1].sum() + R_t[:, :, -1].sum())
   plt.plot(t_vec, S_t.sum(0).sum(0), 'yo:')
   plt.plot(t_vec, I_t.sum(0).sum(0), 'ro:')
   plt.plot(t_vec, R_t.sum(0).sum(0), 'bo:')
   plt.title("Trajectories")
   plt.show()


