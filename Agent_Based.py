import numpy as np
import matplotlib.pyplot as plt


# This class represents each person
class Person:
    def __init__(self, year, housing_loc, current_state, P_type, params_list):
        self.year = year
        self.housing_loc = housing_loc #Location class
        self.current_state = current_state #int of 0 (susceptible), 1 (infected), or 2 (recovered)
        self.P_type = P_type
        self.P = self.define_probability(params_list) # Assigns a probability distribution (all located at end of this file)
        self.past_locations = [] # Stores all the prior locations visited at each time step
        self.past_states = []  # Stores the status of the person (0 = susceptible, 1 = infected, 2 = recovered) at each prior time step

    def define_probability(self, params): #This function will calculate the person's probability
        if self.P_type == "Gaussian":
            return Gaussian(params[0], params[1])

        if self.P_type == "Poisson":
            return Poisson(params[0])

        if self.P_type == "Beta":
            return Beta(params[0], params[1])

        if self.P_type == "Lognormal":
            return Lognormal(params[0], params[1])

        raise Exception("P_type was not in selected list")

# This class represents each location (an associated ID and name)
class Location:
    def __init__(self, ind, name, protocol):
        self.ind = ind #int of (0, 100) for their location
        self.name = name #string that represents the name
        self.protocol = protocol #Restrictions/policies put in place

# This class represents protocols for each location (i.e. max people, spacing, masks. etc.)
# Still need to implement a function that takes into account these protocols and defined probability parameters
    # For example: Masks might decrease the mean of the gaussian that models the infection probability
class Protocol:
    def __init__(self, max_people, spacing, masks):
        self.max_people = max_people
        self.spacing = spacing
        self.masks = masks


# The below classes are probability distributions
class Gaussian:
    def __init__(self, mu, st_dev):
        self.mu = mu
        self.st_dev = st_dev
        self.type = "Gaussian"

    def sample(self, N):
        samples = np.random.normal(self.mu, self.st_dev, N)
        return samples

    def plot(self, num_samples, bins, title, num_rows, num_cols, ind):
        samples = self.sample(num_samples)
        plt.subplot(num_rows, num_cols, ind)
        plt.hist(samples, bins, density=True)
        plt.title(title)

class Poisson:
    def __init__(self, lam):
        self.lam = lam
        self.type = "Poisson"

    def sample(self, N):
        samples = np.random.poisson(self.lam, N)
        return samples

    def plot(self, num_samples, bins, title, num_rows, num_cols, ind):
        samples = self.sample(num_samples)
        plt.subplot(num_rows, num_cols, ind)
        plt.hist(samples, bins, density=True)
        plt.title(title)

class Beta:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def sample(self, N):
        samples = np.random.beta(self.a, self.b, N)
        return samples

    def plot(self, num_samples, bins, title, num_rows, num_cols, ind):
        samples = self.sample(num_samples)
        plt.subplot(num_rows, num_cols, ind)
        plt.hist(samples, bins, density=True)
        plt.title(title)

class Lognormal:
    def __init__(self, mu, st_dev):
        self.mu = mu
        self.st_dev = st_dev

    def sample(self, N):
        samples = np.random.lognormal(self.mu, self.st_dev, N)
        return samples

    def plot(self, num_samples, bins, title, num_rows, num_cols, ind):
        samples = self.sample(num_samples)
        plt.subplot(num_rows, num_cols, ind)
        plt.hist(samples, bins, density=True)
        plt.title(title)


# This is just to test the above classes
loc_1 = Location(ind=0, name="Wilbur Dorm", protocol=Protocol(max_people=100, spacing=6, masks=True))
person_1  = Person(year="Freshman", housing_loc=loc_1, current_state=0, P_type="Gaussian", params_list=[0.5, 0.1])



g = Gaussian(0, 0.1)
p = Poisson(2)
b = Beta(8, 2)
l = Lognormal(0, 0.1)
g.plot(10000, 100, 'Gauss', 2, 2, 1)
p.plot(10000, 100, 'Poisson', 2, 2, 2);
b.plot(10000, 100, 'Beta', 2, 2, 3);
l.plot(10000, 100, 'Lognormal', 2, 2, 4);
plt.show()