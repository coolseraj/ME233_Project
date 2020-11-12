import numpy as np
import matplotlib.pyplot as plt


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
                for i in range(0, len(self.locations[loc_ind].people)):
                    person_1_id = location.people[i].id
                    person_1 = self.states[t_ind][person_1_id] #0 1 or 2
                    self.locations[loc_ind].people[i].past_locations.append(location)
                    self.people[person_1_id].past_locations.append(location)
                    if person_1 == 0:
                        continue #This might fuck up
                    if person_1 == 2:
                        self.states[t_ind+1][person_1_id] = 2
                        continue
                    for j in range(0, len(self.locations[loc_ind].people)):
                        person_2_id = location.people[j].id
                        if person_1_id == person_2_id:
                            continue
                        person_2 = self.states[t_ind][person_2_id]
                        if person_2 == 1:
                            continue
                        if person_2 == 2:
                            self.states[t_ind+1][person_2_id] = 2
                            continue
                        #Collision
                        val_infect_location = location.protocol.P.sample(1)
                        thresh = 1.3
                        """ if self.t[t_ind] > 2 and self.t[t_ind] < 8:
                            thresh = 2
                        if self.t[t_ind] > 8 and self.t[t_ind] < 13:
                            thresh = 1.4
                        if self.t[t_ind] > 13 and self.t[t_ind] < 19:
                            thresh = 2
                        if self.t[t_ind] > 19:
                            thresh = 1.4
                        """

                        if val_infect_location> thresh and self.states[t_ind+1][person_2_id] != 1: #Collision
                            self.states[t_ind+1][person_2_id] = 1
                            print('Person  ' + str(location.people[j].id) + ' infected by person ' + str(location.people[i].id)  + 'at ' + location.name)
                            loc_inf[loc_ind] = loc_inf[loc_ind] + 1

                    #Check if person i recovers
                    days_sick = np.sum(self.states[0:t_ind, person_1_id])
                    if days_sick >= recovery[person_1_id]:
                        self.states[t_ind+1, person_1_id] = 2
                    else:
                        self.states[t_ind+1, person_1_id] = 1

            S[t_ind + 1] = len(np.where(self.states[t_ind + 1 ,:] == 0)[0])
            I[t_ind + 1] = len(np.where(self.states[t_ind + 1, :] == 1)[0])
            R[t_ind + 1] = len(np.where(self.states[t_ind + 1, :] == 2)[0])
            #Reassign people
            #Zero out
            for loc_ind in range(0, len(self.locations)):
                self.locations[loc_ind].people = []
            for person in self.people:
                new_loc = np.random.choice(self.locations, 1, p=person.P_transition)[0]
                self.locations[new_loc.ind].people.append(person)
        return S, I, R, loc_inf















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
            mu, sigma = 0.6, 0.1
        if masks and not outdoors:
            mu, sigma = 0.8, 0.1
        if not masks and outdoors:
            mu, sigma = 0.5, 0.1
        if not masks and not outdoors:
            mu, sigma = 1, 0.1
        parameters = [mu, sigma]#Factor in how indoors and outdoors affects parameters
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
people = []
for p in range(0, 1000):
    people.append(Person(year="Grad", current_state=0, P_transition=[0.1, 0.5, 0.1, 0.3], id=p))
people[5].current_state = 1
people[55].current_state = 1
people[65].current_state = 1
people[85].current_state = 1



loc_1 = Location(ind=0, name="EVGR", protocol=Protocol(max_people=100, masks=True, outdoors=False, P_type="Gaussian"), people=people[0:100])
loc_2 = Location(ind=1, name="Gym", protocol=Protocol(max_people=15, masks=False, outdoors=False, P_type="Gaussian"), people=people[100:200])
loc_3 = Location(ind=2, name="Dining Hall", protocol=Protocol(max_people=10, masks=True, outdoors=False, P_type="Gaussian"), people=people[200:300])
loc_4 = Location(ind=3, name="Oval", protocol=Protocol(max_people=50, masks=True, outdoors=True, P_type="Gaussian"), people=people[300:1000])
locations = [loc_1 , loc_2, loc_3, loc_4] ##### MUST ADD LOCATIONS IN ORDER OF IND




sim = Simulator(t_domain=[0, 100], dt=1, locations=locations, people=people, beta=3/6.5, gamma=1/6.5)
plt.figure()
plt.plot(sim.t, sim.S, 'y')
plt.plot(sim.t, sim.I, 'r')
plt.plot(sim.t, sim.R, 'b')
plt.legend(("S", "I", "R"))
plt.show()

g = Gaussian(0, 0.1)
p = Poisson(2)
b = Beta(8, 2)
l = Lognormal(0, 0.1)
g.plot(10000, 100, 'Gauss', 2, 2, 1)
p.plot(10000, 100, 'Poisson', 2, 2, 2);
b.plot(10000, 100, 'Beta', 2, 2, 3);
l.plot(10000, 100, 'Lognormal', 2, 2, 4);
plt.show()