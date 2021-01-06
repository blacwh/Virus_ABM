import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations
from mpl_toolkits.axes_grid1 import make_axes_locatable


STATUS = ['susceptible', 'infected', 'recovered', 'death', 'hospital']
COLORS = ['lightskyblue', 'orangered', 'green', 'gray','yellow']
contact_distance = 0.03


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
        # move_test = np.random.random()
        # if move_test > 0.8:
        #     return
        self.vx, self.vy = random_vel() 
        self.r += self.v * dt


class Simulation:
    """A class for a simple hard-circle molecular dynamics simulation.

    The simulation is carried out on a square domain: 0 <= x < 1, 0 <= y < 1.

    """

    ParticleClass = Particle

    def __init__(self, n, infectn, radius, iterations, hospital=False):
        """Initialize the simulation with n Particles with radii radius.

        radius can be a single value or a sequence with n values.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor when drawing
        the Particles.

        """

        self.init_particles(n, infectn, radius)

        self.dt = 0.5
        self.death_ratio = 0.002
        self.infect_ratio = 0.8
        # self.contact_distance = 0.03
        self.iterations = iterations
        self.total_num = n
        self.current_susceptible = n - infectn
        self.current_infections = infectn
        self.current_recovered = 0
        self.current_death = 0
        self.hospital = hospital
        if self.hospital:
            self.hospital_capacity = 15
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

                if self.current_infections > 15 and self.hospital_num < self.hospital_capacity:
                    q_test = np.random.random()
                    if q_test < 0.8:
                        p.status = STATUS[4]
                        p.color = COLORS[4]
                        # add into the quarantine compartment
                        p.remove()
                        p.draw(self.ax_h)
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

                if p.infect_time > 24 or p.hospital_time > 14:
                    p.status = STATUS[2]
                    p.color = COLORS[2]
                    p.infect_time = 0
                    p.hospital_time = 0
                    # go back to the main compartment
                    # p.x = p.radius + (1-2*p.radius)*np.random.random()
                    # p.y = p.radius + (1-2*p.radius)*np.random.random()
                    p.draw(self.ax_h)
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

    def animate(self, i):
        """The function passed to Matplotlib's FuncAnimation routine."""

        self.advance_animation()
        self.xdata.append(i)
        self.ydata1.append(self.current_susceptible)
        self.ydata2.append(self.current_infections)
        self.ydata3.append(self.current_recovered)
        self.ydata4.append(self.current_death)
        self.ydata5.append(self.hospital_num)

        self.line_sus.set_data(self.xdata, self.ydata1)
        self.line_inf.set_data(self.xdata, self.ydata2)
        self.line_rec.set_data(self.xdata, self.ydata3)
        self.line_dea.set_data(self.xdata, self.ydata4)
        self.line_hos.set_data(self.xdata, self.ydata5)

        # self.line.set_data(self.xdata,self.ydata)
        return self.circles, self.line_sus, self.line_inf, self.line_rec, self.line_dea, self.line_hos,

    def setup_animation(self):
        self.fig, self.ax = plt.subplots(nrows=1, ncols=2, figsize=[14, 5])

        for s in ['top','bottom','left','right']:
            self.ax[0].spines[s].set_linewidth(2)
        self.ax[0].set_aspect('equal', 'box')
        self.ax[0].set_xlim(0, 1)
        self.ax[0].set_ylim(0, 1)
        self.ax[0].set_title('ABM Simulation')
        # create new axes on the right of the current axes
        divider = make_axes_locatable(self.ax[0])
        # below height and pad are in inches
        self.ax_h = divider.append_axes('right',1, pad=0.2, sharey=self.ax[0], sharex=self.ax[0])
        # make some labels invisible
        self.ax_h.xaxis.set_tick_params(labelbottom=False)
        self.ax_h.yaxis.set_tick_params(labelleft=False)


        self.ax[1].set_xlim(0,self.iterations)
        self.ax[1].set_ylim(0,self.total_num)
        self.ax[1].set_title('Statistic Data')
        self.ax[1].set_xlabel("Nº of Days")
        self.ax[1].set_ylabel("Total Number")
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
        self.line_sus, = self.ax[1].plot([], [], c=COLORS[0], label=STATUS[0])
        self.line_inf, = self.ax[1].plot([], [], c=COLORS[1], label=STATUS[1])
        self.line_rec, = self.ax[1].plot([], [], c=COLORS[2], label=STATUS[2])
        self.line_dea, = self.ax[1].plot([], [], c=COLORS[3], label=STATUS[3])
        self.line_hos, = self.ax[1].plot([], [], c=COLORS[4], label=STATUS[4])

        # self.line, = self.ax[1].plot([], [], 'r-', animated=False, label='infections')  
        handles, labels = self.ax[1].get_legend_handles_labels()
        lgd = self.ax[1].legend(handles, labels, loc='upper right')



    def save_or_show_animation(self, anim, save, filename='abm.mp4'):
        if save:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=4, bitrate=1800)
            anim.save(filename, writer=writer)
        else:
            plt.show()

    def do_animation(self, save=False, filename='abm.mp4'):
        """Set up and carry out the animation of the molecular dynamics.

        To save the animation as a MP4 movie, set save=True.
        """

        self.setup_animation()
        anim = animation.FuncAnimation(self.fig, self.animate,
                init_func=self.init, frames=self.iterations, interval=250, blit=False)
        self.save_or_show_animation(anim, save, filename)


if __name__ == '__main__':
    nparticles = 100
    radius = 0.012
    infectn = 2
    iterations = 80
    sim = Simulation(nparticles, infectn, radius, iterations, hospital=True)
    sim.do_animation(save=False)
