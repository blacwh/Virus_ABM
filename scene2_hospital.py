import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

STATUS = ['susceptible', 'infected', 'recovered', 'death', 'hospital']
COLORS = ['lightskyblue', 'orangered', 'green', 'gray','yellow']
contact_distance = 0.02
sns.set_theme(context="paper", style="whitegrid")

def random_vel():
    vr = 0.1 * np.sqrt(np.random.random()) + 0.05
    vphi = 2*np.pi * np.random.random()
    vx, vy = vr * np.cos(vphi), vr * np.sin(vphi) * 0.5  
    return vx, vy

class Particle:
    """A class representing a two-dimensional particle."""

    def __init__(self, x, y, vx, vy, radius, status, color):
        """Initialize the particle's position, velocity, and radius.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor.

        """

        self.r = np.array((x, y))
        self.v = np.array((vx, vy))
        self.radius = radius
        self.status = status
        self.color = color
        self.infect_time = 0
        self.hospital_time = 0
        self.time = 0
        self.wander_step_duration = int(np.random.random() * 10)
        self.last_change_time = 0

    # For convenience, map the components of the particle's position and
    # velocity vector onto the attributes x, y, vx and vy.
    @property
    def x(self):
        return self.r[0]
    @x.setter
    def x(self, value):
        self.r[0] = value
    @property
    def y(self):
        return self.r[1]
    @y.setter
    def y(self, value):
        self.r[1] = value
    @property
    def vx(self):
        return self.v[0]
    @vx.setter
    def vx(self, value):
        self.v[0] = value
    @property
    def vy(self):
        return self.v[1]
    @vy.setter
    def vy(self, value):
        self.v[1] = value

    def contact(self, other):
        """Does the circle of this agent contact that of other one?"""

        return np.hypot(*(self.r - other.r)) < self.radius + other.radius + contact_distance

    def draw(self, ax):
        """Add this Particle's Circle patch to the Matplotlib Axes ax."""

        self.circle = Circle(xy=self.r, radius=self.radius, color=self.color)
        ax.add_patch(self.circle)
        return self.circle

    def remove(self):
        self.circle.remove()

    def advance(self, dt):#move function
        """Advance the Particle's position forward in time by dt."""
        if self.status == STATUS[3] or self.status == STATUS[4]:
            return

        self.time += 1
        # it = int(np.random.random() * 10)
        if (self.time - self.last_change_time) > self.wander_step_duration:
            # self.vx, self.vy = random_vel() 
            self.vx = int(np.random.randn() * 7)/100
            self.vy = int(np.random.randn() * 7)/100
            self.last_change_time = self.time

        self.r += self.v * dt    


class Simulation:
    """A class for a simple hard-circle molecular dynamics simulation.

    The simulation is carried out on a square domain: 0 <= x < 1, 0 <= y < 1.

    """

    ParticleClass = Particle

    def __init__(self, n, infectn, radius, iterations, hospital_capacity, quarantine_probability,
                        quarantine_efficiency, quarantine_time, hospital=False):
        """Initialize the simulation with n Particles with radii radius.

        radius can be a single value or a sequence with n values.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor when drawing
        the Particles.

        """

        self.init_particles(n, infectn, radius)

        self.dt = 0.15
        self.death_ratio = 0.001
        self.infect_ratio = 0.6
        # self.contact_distance = 0.03
        self.iterations = iterations
        self.total_num = n
        self.current_susceptible = n - infectn
        self.current_infections = infectn
        self.current_recovered = 0
        self.current_death = 0
        self.hospital = hospital
        if self.hospital:
            self.hospital_capacity = hospital_capacity * n
            self.quarantine_probability = quarantine_probability
            self.quarantine_efficiency = quarantine_efficiency
            self.quarantine_time = quarantine_time
            self.hospital_num = 0


    def place_particle(self, rad, status, color):
        # # Choose x, y so that the Particle is entirely inside the
        # # domain of the simulation.
        x, y = rad + (1 - 2*rad) * np.random.random(2)

        # Choose a random velocity (within some reasonable range of
        # values) for the Particle.
        # vr = 0.1 * np.sqrt(np.random.random()) + 0.05
        # vphi = 2*np.pi * np.random.random()
        # vx, vy = vr * np.cos(vphi), vr * np.sin(vphi) * 0.5

        vx, vy = random_vel()


        particle = self.ParticleClass(x, y, vx, vy, rad, status, color)

        self.particles.append(particle)



    def init_particles(self, n, infectn, radius):
        """Initialize the n Particles of the simulation."""




        self.n = n
        self.infectn = infectn
        self.particles = []


        for i in np.arange(self.infectn):
            self.place_particle(radius, STATUS[1], COLORS[1])
        for i in np.arange(self.n - len(self.particles)):
            self.place_particle(radius, STATUS[0], COLORS[0])
            
               
    def infectionDetect(self, a1, a2):
        
        
        if a1.status == STATUS[0] and a2.status == STATUS[1]:
            infect_test = np.random.random()
            if infect_test < self.infect_ratio:
                a1.status = STATUS[1]
                a1.color = COLORS[1]
                a1.draw(self.ax[0])
                self.current_susceptible -= 1
                self.current_infections += 1

        
    def handle_collisions(self):
        """Detect and handle any collisions between the Particles.

        When two Particles collide, they do so elastically: their velocities
        change such that both energy and momentum are conserved.

        """ 

        # We're going to need a sequence of all of the pairs of particles when
        # we are detecting collisions. combinations generates pairs of indexes
        # into the self.particles list of Particles on the fly.
        pairs = combinations(range(self.n), 2)
        for i,j in pairs:
            if self.particles[i].contact(self.particles[j]):
                self.infectionDetect(self.particles[i], self.particles[j])
                self.infectionDetect(self.particles[j], self.particles[i])


    def handle_boundary_collisions(self, p):
        """Bounce the particles off the walls elastically."""

        if p.x - p.radius < 0:
            p.x = p.radius
            p.vx = -p.vx
        if p.x + p.radius > 1:
            p.x = 1-p.radius
            p.vx = -p.vx
        if p.y - p.radius < 0:
            p.y = p.radius
            p.vy = -p.vy
        if p.y + p.radius > 1:
            p.y = 1-p.radius
            p.vy = -p.vy
        


    def update(self, p):

        if p.status == STATUS[3]:
            return

        if self.hospital:
            if p.status == STATUS[1]:
                p.infect_time += 1
                death_test = np.random.random()
                if death_test < self.death_ratio:
                    p.status = STATUS[3]
                    p.color = COLORS[3]
                    p.draw(self.ax[0])
                    self.current_death += 1
                    self.current_infections -= 1
                    return
                
                if p.infect_time > 24:
                    p.status = STATUS[2]
                    p.color = COLORS[2]
                    p.draw(self.ax[0])
                    p.infect_time = 0
                    self.current_infections -= 1
                    self.current_recovered += 1

                if p.infect_time > self.quarantine_efficiency and self.hospital_num < self.hospital_capacity:
                    q_test = np.random.random()
                    if q_test < self.quarantine_probability:
                        p.status = STATUS[4]
                        p.color = COLORS[4]
                        # add into the quarantine compartment
                        p.draw(self.ax[0])
                        # p.remove()
                        p.draw(self.ax[1])
                        self.hospital_num += 1
                        

            if p.status == STATUS[4]:
                p.hospital_time += 1
                p.infect_time += 1
                deathq_test = np.random.random()
                if deathq_test < self.death_ratio / 4:
                    p.status = STATUS[3]
                    p.color = COLORS[3]
                    p.draw(self.ax[0])
                    p.infect_time = 0
                    p.hospital_time = 0
                    self.current_death += 1
                    self.current_infections -= 1
                    self.hospital_num -= 1
                    return

                if p.infect_time > 24 or p.hospital_time > self.quarantine_time:
                    p.status = STATUS[2]
                    p.color = COLORS[2]
                    p.infect_time = 0
                    p.hospital_time = 0
                    # go back to the main compartment
                    p.remove()
                    p.draw(self.ax[0])
                    self.current_recovered += 1
                    self.current_infections -= 1
                    self.hospital_num -= 1

        # else:
        #     if p.status == STATUS[1]:
        #         p.infect_time += 1
        #         death_test = np.random.random()
        #         if death_test < self.death_ratio:
        #             p.status = STATUS[3]
        #             p.color = COLORS[3]
        #             p.draw(self.ax[0])
        #             self.current_death += 1
        #             self.current_infections -= 1
        #             return
                
        #         if p.infect_time > 24:
        #             p.status = STATUS[2]
        #             p.color = COLORS[2]
        #             p.draw(self.ax[0])
        #             p.infect_time = 0
        #             self.current_infections -= 1
        #             self.current_recovered += 1



    def advance_animation(self):
        """Advance the animation by dt, returning the updated Circles list."""

        for p in self.particles:
            p.advance(self.dt)
            self.handle_boundary_collisions(p)
            self.update(p)
            

        self.handle_collisions()
        return self.circles


    def init(self):
        """Initialize the Matplotlib animation."""

        self.circles = []
        for particle in self.particles:
            self.circles.append(particle.draw(self.ax[0]))
        return self.circles, self.line_sus, self.line_inf, self.line_rec, self.line_dea,

    def create_splot(self, iterations, tnumber, line1, line2, line3, line4, line5, title):
        x = np.arange(0, iterations)
        fig, ax = plt.subplots()
        ax.plot(x, line1, label='S', color='royalblue')
        ax.plot(x, line2, label='I', color='tomato')
        ax.plot(x, line3, label='R', color='springgreen')
        ax.plot(x, line4, label='D', color='black')
        ax.plot(x, line5, label='H', color='yellow')
        ax.axhline(y=10, ls="--", c="navy",label='Cap')
        # ax.axhline(y=10, ls='-.', c="purple",label='release')
        ax.set_xlim(0, iterations)
        ax.set_ylim(0, tnumber)
        ax.set_xlabel('No. of Days', fontdict={'family': 'Calibri', 'size': 13, 'weight': "bold"})
        ax.set_ylabel('Total Number', fontdict={'family': 'Calibri', 'size': 13, 'weight': "bold"})
        ax.set_title(title, fontdict={'family': 'Times New Roman', 'size': 16, 'weight': "bold"})
        ax.legend()
        plt.show()

    def animate(self, i):
        """The function passed to Matplotlib's FuncAnimation routine."""

        self.advance_animation()
        self.xdata.append(i)
        self.ydata1.append(self.current_susceptible)
        self.ydata2.append(self.current_infections)
        self.ydata3.append(self.current_recovered)
        self.ydata4.append(self.current_death)
        self.ydata5.append(self.hospital_num)

        if len(self.xdata) == self.iterations:
            self.create_splot(100, 200, self.ydata1, self.ydata2, self.ydata3, self.ydata4, self.ydata5, 'Scenario-6')
            # with open('s2.txt','w') as f:
            #     f.write('test: ')
            #     f.write('\nx: ')
            #     f.write(str(self.xdata))
            #     f.write('\nc_sus: ')
            #     f.write(str(self.ydata1))
            #     f.write('\nc_inf: ')
            #     f.write(str(self.ydata2))
            #     f.write('\nc_rec: ')
            #     f.write(str(self.ydata3))
            #     f.write('\nc_dea: ')
            #     f.write(str(self.ydata4))
            #     f.write('\nc_h: ')
            #     f.write(str(self.ydata5))

        self.line_sus.set_data(self.xdata, self.ydata1)
        self.line_inf.set_data(self.xdata, self.ydata2)
        self.line_rec.set_data(self.xdata, self.ydata3)
        self.line_dea.set_data(self.xdata, self.ydata4)
        self.line_hos.set_data(self.xdata, self.ydata5)

        # self.line.set_data(self.xdata,self.ydata)
        return self.circles, self.line_sus, self.line_inf, self.line_rec, self.line_dea, self.line_hos,

    def setup_animation(self):
        self.fig, self.ax = plt.subplots(nrows=1, ncols=3, figsize=[18, 5])

        for s in ['top','bottom','left','right']:
            self.ax[0].spines[s].set_linewidth(2)
        self.ax[0].set_aspect('equal', 'box')
        self.ax[0].set_xlim(0, 1)
        self.ax[0].set_ylim(0, 1)
        self.ax[0].set_title('ABM Simulation')
        # create new axes on the right of the current axes
        # divider = make_axes_locatable(self.ax[0])
        # below height and pad are in inches
        # self.ax_h = divider.append_axes('right',2.0, pad=0.2, sharey=self.ax[0])
        # make some labels invisible
        # self.ax_h.xaxis.set_tick_params(labelbottom=False)
        # self.ax_h.yaxis.set_tick_params(labelleft=False)

        self.ax[1].set_aspect('equal', 'box')
        self.ax[1].set_title('hopistal compartment')

        self.ax[2].set_xlim(0,self.iterations)
        self.ax[2].set_ylim(0,self.total_num)
        self.ax[2].set_title('Statistic Data')
        self.ax[2].set_xlabel("Nº of Days")
        self.ax[2].set_ylabel("Total Number")
        self.ax[2].axhline(y=int(self.hospital_capacity), ls="--", c="black",label='capacity')
        self.xdata = []
        self.ydata1 = []
        self.ydata2 = []
        self.ydata3 = []
        self.ydata4 = []
        self.ydata5 = []

        self.lines_d = {
            'susceptible':self.current_susceptible,
            'infected':self.current_infections,
            'recovered':self.current_recovered,
            'death':self.current_death
        }

        # self.lines = []
        self.line_sus, = self.ax[2].plot([], [], c=COLORS[0], label=STATUS[0])
        self.line_inf, = self.ax[2].plot([], [], c=COLORS[1], label=STATUS[1])
        self.line_rec, = self.ax[2].plot([], [], c=COLORS[2], label=STATUS[2])
        self.line_dea, = self.ax[2].plot([], [], c=COLORS[3], label=STATUS[3])
        self.line_hos, = self.ax[2].plot([], [], c=COLORS[4], label=STATUS[4])

        # self.line, = self.ax[1].plot([], [], 'r-', animated=False, label='infections')  
        handles, labels = self.ax[2].get_legend_handles_labels()
        lgd = self.ax[2].legend(handles, labels, loc='upper right')



    def save_or_show_animation(self, anim, save, filename='hospital.mp4'):
        if save:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=4, bitrate=1800)
            anim.save(filename, writer=writer)
        else:
            plt.show()

    def do_animation(self, save=False, filename='hospital.mp4'):
        """Set up and carry out the animation of the molecular dynamics.

        To save the animation as a MP4 movie, set save=True.
        """

        self.setup_animation()
        anim = animation.FuncAnimation(self.fig, self.animate,
                init_func=self.init, frames=self.iterations, interval=100, blit=False)
        self.save_or_show_animation(anim, save, filename)


if __name__ == '__main__':
    '''
    Scenario 2 (hospital):
    para1: hospital capacity (5%--20% of total population) (10% default)
    para2: quarantine probability (50%--100%) (80% default)
    para3: quarantine efficiency (system time 7--21) (14 default)
    para4: quarantine time (system time 10--20) (14 default)
    '''
    hospital_capacity = float(input("Please input the hospital capacity (5%--20% is recommended): "))
    quarantine_probability = float(input("Please input the quarantine probability (50%--100% is recommended): "))
    quarantine_efficiency = int(input("Please input the time for the infected to be quarantined (7--21 is recommended): "))
    quarantine_time = int(input("Please input the quarantine time (10--20 is recommended): "))
    nparticles = 200
    radius = 0.012
    infectn = 2
    iterations = 100
    sim = Simulation(nparticles, infectn, radius, iterations, hospital_capacity, quarantine_probability,
                        quarantine_efficiency, quarantine_time, hospital=True)
    sim.do_animation(save=False)
