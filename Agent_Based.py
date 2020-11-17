import numpy as np
import matplotlib.pyplot as plt
import random


#Assumes that the infector determines the rate of infection. Does not depend on location. Only number of people and year
# This class runs the simulation

class Simulator:
    def __init__(self, t_domain, dt, locations, people, beta, gamma):
        self.t = np.linspace(t_domain[0], t_domain[1], int((t_domain[1] - t_domain[0])/dt + 1))
        self.locations = locations
        self.dt = dt
        self.beta = beta
        self.gamma = gamma
        self.people = people
        self.states = np.zeros((len(self.t), len(self.people)))
        self.S, self.I, self.R, self.Loc_Inf = self.simulate()


    def simulate(self):
        recovery =  abs(Gaussian(1/self.gamma, 2).sample(len(self.people)))
        loc_inf = np.zeros((len(self.locations)))
        S, I, R = np.zeros((len(self.t))), np.zeros((len(self.t))), np.zeros((len(self.t)))
        num_infected = 0
        for person_ind in range(0,len(self.people)):
            if self.people[person_ind].current_state == 1:
                num_infected = num_infected + 1
                self.states[0][person_ind] = 1

        I[0] = num_infected
        S[0] = len(self.people)-num_infected

        for t_ind in range(0, len(self.t)-1):
            for loc_ind in range(0, len(self.locations)):
                location = self.locations[loc_ind]
                ids = [person.id for person in location.people]
                potential_infs, potential_inf_ids = self.precompute_infections(location)
                for i in range(0, len(ids)):
                    potential_infector_id = ids[i]
                    #Update the location just to store
                    self.people[potential_infector_id].past_locations.append(location)
                    self.update_state(t_ind, potential_infector_id, potential_inf_ids[i], recovery)
                #Update person_1's state itself?

            S[t_ind + 1] = len(np.where(self.states[t_ind + 1 ,:] == 0)[0])
            I[t_ind + 1] = len(np.where(self.states[t_ind + 1, :] == 1)[0])
            R[t_ind + 1] = len(np.where(self.states[t_ind + 1, :] == 2)[0])

            total_pop_size = len(self.people)  # Store total size of people array
            lin_list = np.linspace(0, total_pop_size - 1, total_pop_size)
            #draw_sample = np.random.choice(lin_list, total_pop_size)  # Make linear list random for sampling
            np.random.shuffle(lin_list)
            draw_sample = lin_list
            for loc_ind in range(0, len(self.locations)):
                self.locations[loc_ind].people = []
            # Go through randomized list of people
            for person_id in draw_sample:  # pretend there is no max for each room, invert location and persons so can choose 100# person for a room
                person_id = int(person_id)
                new_loc = np.random.choice(self.locations, 1, p=self.people[person_id].P_transition)[0]
                if len(new_loc.people) > new_loc.protocol.max_people:
                    raise TypeError("Maximum Capacity exceeded")
                self.locations[new_loc.ind].people.append(self.people[person_id])


        return S, I, R, loc_inf

    def precompute_infections(self, location):
        ids = [person.id for person in location.people]
        infections = location.protocol.P.sample(len(ids))
        inf_mat = []
        # Assign who you will infect
        for i in range(0, len(ids)):
            infector_id = ids[i]
            inf_ids = random.choices(ids, k=infections[i])
            inf_mat.append(list(inf_ids))
            # Might pick the person themselves, or the same person twice. But this shouldn't be a problem
        return infections, inf_mat

    def update_state(self, t_ind, person_1_id, person_2_ids, recovery):
        #Update the state by checking for the collision

        person_1_state = self.states[t_ind, person_1_id]
        if person_1_state == 0:
            return
        if person_1_state == 2:
            self.states[t_ind+1][person_1_id] = 2
            return

        for j in range(0, len(person_2_ids)):
            person_2_id = person_2_ids[j]
            person_2_state = self.states[t_ind, person_2_id]

            if person_1_id == person_2_id:
                continue
            if person_2_state == 1:
                continue
            if person_2_state == 2:
                self.states[t_ind + 1][person_2_id] = 2
                continue
            if self.states[t_ind+1, person_2_id] == 1:
                continue
            #None of the other if statements were triggered, so this is a successful infection
            self.states[t_ind + 1,person_2_id] = 1
        #Check recovery

        days_sick = np.sum(self.states[0:t_ind, person_1_id])
        if person_1_state == 1:
            if days_sick >= recovery[person_1_id]:
                self.states[t_ind + 1, person_1_id] = 2
            else:
                self.states[t_ind + 1, person_1_id] = 1




# This class represents each person
class Person:
    def __init__(self, year, current_state, P_transition, id):
        self.year = year
        self.current_state = current_state #int of 0 (susceptible), 1 (infected), or 2 (recovered)
        self.past_locations = [] # Stores all the prior locations visited at each time step
        self.P_transition = P_transition # vector of length locations
        self.id = id


# This class represents each location (an associated ID and name)
class Location:
    def __init__(self, ind, name, protocol, people):
        self.ind = ind #int of (0, 100) for their location
        self.name = name #string that represents the name
        self.protocol = protocol #Restrictions/policies put in place
        self.people = people

# This class represents protocols for each location (i.e. max people, spacing, masks. etc.)
# Still need to implement a function that takes into account these protocols and defined probability parameters
    # For example: Masks might decrease the mean of the gaussian that models the infection probability
class Protocol:
    def __init__(self, max_people, masks, outdoors, P_type):
        self.max_people = max_people
        self.masks = masks
        self.outdoors = outdoors
        self.P_type = P_type
        if masks and outdoors:
            mu, sigma, lam = 0.6, 0.1, dt*(beta*0.8)
        if masks and not outdoors:
            mu, sigma, lam = 0.8, 0.1, beta *dt
        if not masks and outdoors:
            mu, sigma, lam = 0.5, 0.1, beta * dt
        if not masks and not outdoors:
            mu, sigma, lam = 1, 0.1, dt*(beta *1.5)
        if P_type == "Gaussian":
            parameters = [mu, sigma]#Factor in how indoors and outdoors affects parameters
        if P_type == "Poisson":
            parameters = [lam] #Poisson
        self.P = self.define_probability(parameters)

    def define_probability(self, params):
        if self.P_type == "Gaussian":
            return Gaussian(params[0], params[1])

        if self.P_type == "Poisson":
            return Poisson(params[0])

        if self.P_type == "Beta":
            return Beta(params[0], params[1])

        if self.P_type == "Lognormal":
            return Lognormal(params[0], params[1])

        raise Exception("P_type was not in selected list")


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
        plt.hist(samples, bins, density=False)
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


beta = 2/6.5
gamma = 1/6.5
dt = 1


# This is just to test the above classes
people = []
for p in range(0, 1000):
    people.append(Person(year="Grad", current_state=0, P_transition=[1, 0, 0, 0], id=p))
people[5].current_state = 1
people[55].current_state = 1
people[65].current_state = 1
people[85].current_state = 1


type = "Poisson"
loc_1 = Location(ind=0, name="EVGR", protocol=Protocol(max_people=10e5, masks=True, outdoors=False, P_type=type), people=people[0:250])
loc_2 = Location(ind=1, name="Gym", protocol=Protocol(max_people=10e5, masks=False, outdoors=False, P_type=type), people=people[250:500])
loc_3 = Location(ind=2, name="Dining Hall", protocol=Protocol(max_people=10e5, masks=True, outdoors=False, P_type=type), people=people[500:750])
loc_4 = Location(ind=3, name="Oval", protocol=Protocol(max_people=10e5, masks=True, outdoors=True, P_type=type), people=people[750:1000])
locations = [loc_1 , loc_2, loc_3, loc_4] ##### MUST ADD LOCATIONS IN ORDER OF IND




sim = Simulator(t_domain=[0, 100], dt=dt, locations=locations, people=people, beta=3/6.5, gamma=1/6.5)

plt.figure()
plt.plot(sim.t, sim.S, 'y')
plt.plot(sim.t, sim.I, 'r')
plt.plot(sim.t, sim.R, 'b')
plt.legend(("S", "I", "R"))
plt.show()

g = Gaussian(0, 0.1)
p = Poisson(0.5)
b = Beta(8, 2)
l = Lognormal(0, 0.1)
g.plot(num_samples=10000, bins=100, title='Gaussian', num_rows=2, num_cols=2, ind=1);
p.plot(num_samples=10000, bins=100, title='Poisson', num_rows=2, num_cols=2, ind=2);
b.plot(num_samples=10000, bins=100, title='Beta', num_rows=2, num_cols=2, ind=3);
l.plot(num_samples=10000, bins=100, title='Lognormal', num_rows=2, num_cols=2, ind=4);
plt.show()