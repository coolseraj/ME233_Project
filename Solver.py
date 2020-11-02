import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


class Solver:
    def __init__(self, x_dom, y_dom, t_dom, num_x, num_y, num_t, gamma, N, LD):
        self.x_min = x_dom[0]
        self.x_max = x_dom[1]
        self.y_min = y_dom[0]
        self.y_max = y_dom[1]
        self.t_min = t_dom[0]
        self.t_max = t_dom[1]
        self.num_x = num_x
        self.num_y = num_y
        self.num_t = num_t
        self.gamma = gamma
        self.N = N

        self.x, self.y, self.t, self.dx, self.dy, self.dt, self.X, self.Y = self.gen_grid()
        self.mu_s_t, self.mu_i_t, self.mu_r_t, self.beta_t = LD.mu_s_t, LD.mu_i_t, LD.mu_r_t, LD.beta_t

        self.s_t = np.ones((self.num_x, self.num_y, self.num_t)) * N
        self.i_t = np.zeros((self.num_x, self.num_y, self.num_t))
        self.r_t = np.zeros((self.num_x, self.num_y, self.num_t))

    def initialize_state(self, rand_row, rand_col, init_sus):
        s = np.ones((num_y, num_x)) * self.N
        s[rand_row][rand_col] = init_sus*self.N
        i = self.N - s
        r = 0 * i
        self.s_t[:, :, 0] = s
        self.i_t[:, :, 0] = i
        self.r_t[:, :, 0] = r
        return


    def gen_grid(self):
        x = np.linspace(self.x_min, self.x_max, self.num_x)
        y = np.linspace(self.y_min, self.y_max, self.num_y)
        t = np.linspace(self.t_min, self.t_max, self.num_t)
        dx = (self.x_max - self.x_min)/(self.num_x-1)
        dy = (self.y_max - self.y_min) / (self.num_y - 1)
        dt = (self.t_max - self.t_min) / (self.num_t - 1)
        X, Y = np.meshgrid(x,y)
        return x, y, t, dx, dy, dt, X, Y

    def update_s(self, i, j, k):
        s_t, i_t, r_t = self.s_t, self.i_t, self.r_t
        dx, dy, dt = self.dx, self.dy, self.dt
        N = self.N
        beta = self.beta_t[k]
        mu_s = self.mu_s_t[k]

        s_old = s_t[i][j][k-1]
        s_sink = beta * i_t[i][j][k-1] * s_t[i][j][k-1] / N
        if i == 0:  # Top No Flux BC
            self.s_t[i][j][k] = s_t[i + 1][j][k - 1]
            return
        if i == self.num_y - 1:  # Bot No Flux BC
            self.s_t[i][j][k] = s_t[i - 1][j][k - 1]
            return
        if j == 0:  # Left No Flux BC
            self.s_t[i][j][k] = s_t[i][j + 1][k - 1]
            return
        if j == self.num_x - 1:  # Right No Flux BC
            self.s_t[i][j][k] = s_t[i][j - 1][k - 1]
            return
        sd_up = mu_s / dy ** 2 * s_t[i - 1][j][k - 1]
        sd_low = mu_s / dy ** 2 * s_t[i + 1][j][k - 1]
        sd_right = mu_s / dx ** 2 * s_t[i][j + 1][k - 1]
        sd_left = mu_s / dx ** 2 * s_t[i][j - 1][k - 1]
        sd_center = 2 * mu_s * (1 / (dx ** 2) + 1 / (dy ** 2)) * s_t[i][j][k - 1]
        s_t[i][j][k] = s_old + dt * (sd_up + sd_low + sd_right + sd_left - sd_center - s_sink)
        self.s_t = s_t
        return

    def update_i(self, i, j, k):
        s_t, i_t, r_t = self.s_t, self.i_t, self.r_t
        dx, dy, dt = self.dx, self.dy, self.dt

        N = self.N
        beta = self.beta_t[k]
        mu_i = self.mu_i_t[k]

        i_old = i_t[i][j][k - 1]
        i_source = beta * i_t[i][j][k - 1] * s_t[i][j][k - 1] / N
        i_sink = self.gamma * i_t[i][j][k - 1]
        if i == 0:  # Top No Flux BC
            self.i_t[i][j][k] = i_t[i + 1][j][k - 1]
            return
        if i == self.num_y - 1:  # Bot No Flux BC
            self.i_t[i][j][k] = i_t[i - 1][j][k - 1]
            return
        if j == 0:  # Left No Flux BC
            self.i_t[i][j][k] = i_t[i][j + 1][k - 1]
            return
        if j == self.num_x - 1:  # Right No Flux BC
            self.i_t[i][j][k] = i_t[i][j - 1][k - 1]
            return
        id_up = mu_i / dy ** 2 * i_t[i - 1][j][k - 1]
        id_low = mu_i / dy ** 2 * i_t[i + 1][j][k - 1]
        id_right = mu_i / dx ** 2 * i_t[i][j + 1][k - 1]
        id_left = mu_i / dx ** 2 * i_t[i][j - 1][k - 1]
        id_center = 2 * mu_i * (1 / (dx ** 2) + 1 / (dy ** 2)) * i_t[i][j][k - 1]
        i_t[i][j][k] = i_old + dt * (id_up + id_low + id_right + id_left - id_center + i_source - i_sink)
        self.i_t = i_t
        return

    def update_r(self, i, j, k):
        s_t, i_t, r_t = self.s_t, self.i_t, self.r_t
        dx, dy, dt = self.dx, self.dy, self.dt

        N = self.N
        mu_r = self.mu_r_t[k]

        r_old = r_t[i][j][k - 1]
        r_source = self.gamma * i_t[i][j][k - 1]
        if i == 0:  # Top No Flux BC
            self.r_t[i][j][k] = r_t[i + 1][j][k - 1]
            return
        if i == self.num_y - 1:  # Bot No Flux BC
            self.r_t[i][j][k] = r_t[i - 1][j][k - 1]
            return
        if j == 0:  # Left No Flux BC
            self.r_t[i][j][k] = r_t[i][j + 1][k - 1]
            return
        if j == self.num_x - 1:  # Right No Flux BC
            self.r_t[i][j][k] = r_t[i][j - 1][k - 1]
            return
        rd_up = mu_r / dy ** 2 * r_t[i - 1][j][k - 1]
        rd_low = mu_r / dy ** 2 * r_t[i + 1][j][k - 1]
        rd_right = mu_r / dx ** 2 * r_t[i][j + 1][k - 1]
        rd_left = mu_r / dx ** 2 * r_t[i][j - 1][k - 1]
        rd_center = 2 * mu_r * (1 / (dx ** 2) + 1 / (dy ** 2)) * r_t[i][j][k - 1]
        r_t[i][j][k] = r_old + dt * (rd_up + rd_low + rd_right + rd_left - rd_center + r_source)
        self.r_t = r_t
        return

    def solve(self):
        rand_row = np.random.randint(0, self.num_y - 1)
        rand_col = np.random.randint(0, self.num_x-1)
        rand_row = 2
        rand_col = 2
        init_sus = 0.99

        self.initialize_state(rand_row, rand_col, init_sus)

        for t_ind in range(1, len(self.t)):
            for i in range(0, num_y):
                for j in range(0, num_x):
                    self.update_s(i, j, t_ind)
                    self.update_i(i, j, t_ind)
                    self.update_r(i, j, t_ind)
        return

    def plot(self):
        s_t, i_t, r_t, t = self.s_t, self.i_t, self.r_t, self.t
        print(s_t[:, :, -1].sum() + i_t[:, :, -1].sum() + r_t[:, :, -1].sum())
        plt.figure()
        plt.plot(t, s_t.sum(0).sum(0), 'y')
        plt.plot(t, i_t.sum(0).sum(0), 'r')
        plt.plot(t, r_t.sum(0).sum(0), 'b')
        plt.title("Trajectories")
        plt.legend(("S", "I", "R"))
        plt.show()


class Lockdown:
    def __init__(self, b0, mu_s0, mu_i0, mu_r0, t, timings, b0_scaling, mu_scaling):
        self.timings = timings
        self.b0_scaling = b0_scaling
        self.mu0_scaling = mu_scaling
        self.b0 = b0
        self.mu_s0 = mu_s0
        self.mu_i0 = mu_i0
        self.mu_r0 = mu_r0
        self.t = t
        self.beta_t = self.beta_dynamic()
        self.mu_s_t, self.mu_i_t, self.mu_r_t = self.mu_dynamic()

    def beta_dynamic(self):
        beta_t = np.ones(len(self.t)) * self.b0
        b0 = self.b0
        for i in range(0, self.timings.shape[0]):
            start_time, end_time = self.timings[i][0], self.timings[i][1]
            start_scale, end_scale = self.b0_scaling[i][0], self.b0_scaling[i][1]
            bool = (self.t >= start_time).astype(int) * (self.t < end_time).astype(int)
            beta_t[np.where(bool)] = np.linspace(start_scale*  self.b0, end_scale * self.b0, bool.sum())
        return beta_t

    def mu_dynamic(self):
        mu_s_t, mu_i_t, mu_r_t = np.ones(len(self.t)) * self.mu_s0, np.ones(len(self.t)) * self.mu_i0, np.ones(len(self.t)) * self.mu_r0
        for i in range(0, self.timings.shape[0]):
            start_time, end_time = self.timings[i][0], self.timings[i][1]
            start_scale, end_scale = self.mu0_scaling[i][0], self.mu0_scaling[i][1]
            bool = (self.t >= start_time).astype(int) * (self.t < end_time).astype(int)
            mu_s_t[np.where(bool)] = np.linspace(start_scale * self.mu_s0, end_scale * self.mu_s0, bool.sum())
            mu_i_t[np.where(bool)] = np.linspace(start_scale * self.mu_i0, end_scale * self.mu_i0, bool.sum())
            mu_r_t[np.where(bool)] = np.linspace(start_scale * self.mu_r0, end_scale * self.mu_r0, bool.sum())
        return mu_s_t, mu_i_t, mu_r_t

    def plot(self):
        t = self.t
        plt.figure()
        plt.subplot(211)
        plt.plot(t, self.beta_t, 'g')
        plt.title('Beta(t)')

        plt.subplot(212)
        plt.plot(t, self.mu_s_t, 'y')
        plt.plot(t, self.mu_i_t, 'r')
        plt.plot(t, self.mu_r_t, 'b')
        plt.title("\mu(t)")
        plt.legend(("mu_S", "mu_I", "mu_R"))
        plt.show()









if __name__ == "__main__":
    x_dom = [0, 1]
    y_dom = [0, 1]
    t_dom = [0, 200]
    num_x = 11
    num_y = 11
    num_t = 10001
    mu_s0, mu_i0, mu_r0 = 0.5e-1, 1e-2, 1e-1
    gamma = 1/6.5
    b0 = 7.5/6.5
    N = 1

    LD_timings = np.array([[0, 10], [10, 40], [40, 100]])  # Needs to connect dates. Esp ending for lift lockdown
    LD_b0_scaling = np.array([[1, 0.2], [0.2, 0.2], [0.2, 1]])  # Needs to be continuous
    LD_mu_scaling = np.array([[1, 0.05], [0.05, 0.05], [0.05, 1]])
    t = np.linspace(t_dom[0], t_dom[1], num_t)
    LD = Lockdown(b0, mu_s0, mu_i0, mu_r0, t, LD_timings, LD_b0_scaling, LD_mu_scaling)
    LD.plot()
    sol = Solver(x_dom, y_dom, t_dom, num_x, num_y, num_t, gamma, N, LD)


    sol.solve()
    sol.plot()
