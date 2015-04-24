from __future__ import division
import pygame   #Imports methods from a file called 'pygame'
from vector import *   #Rename your solution to homework 2 as vector.py 
import random      #Needed to generate random numbers
from itertools import combinations
from random import randint
import sys
import math

author_name = "Erik Polsean"

class particle:
    def __init__(self,pos, vel, mass = 1.0, size = 3, shape = "circle", color = (0,0,255)):
        self.pos = pos   #Position vector
        self.vel = vel   #Velocity vectorB
        self.m = mass    #Particle mass
        self.m_inv = 1/mass
        self.size = size #Radius for circle and width and height for rectangle
        self.color = color
        self.shape = shape  #String with the name of the shape
        self.f_net = vect2d(0,0)   #Property to keep track of net force on particle
        self.pos_old = vect2d(0,0)
        self.vel_old = vect2d(0,0)
        self.f_net_old = vect2d(0,0)
        self.k1 = vect2d(0,0)
        self.k2 = vect2d(0,0)
        self.k3 = vect2d(0,0)
        self.k4 = vect2d(0,0)
        self.L1 = vect2d(0,0)
        self.L2 = vect2d(0,0)
        self.L3 = vect2d(0,0)
        self.L4 = vect2d(0,0)
        self.drag_coeff = 0.000001   #Scaling number for drag force
        self.cor = 0.99     #Coefficient of restitution
        self.rod_connection={}  #NEW: Dictionary of all rod constraint connections
        self.x_body = 0
        self.y_body = 0
        self.x_c = 0
        self.y_c = 0


    def add_force(self,force):
        """Add force to the net force vector for each particle"""
        self.f_net += force


    def reset_force(self):
        """Set net force to zero
        
        Call this at the start of each time step to reset net forces for each particle
        """
        self.f_net = vect2d(0,0) 


    def get_force(self):
        """Return f_net property"""
        return self.f_net

    def set_pos(self,new_pos):
        """Set the position of the particle
        
        Make sure that new_pos is a vect2d object
        """
        self.pos = new_pos


class world:
    def __init__(self):
        self.particle_list = []
        self.all_particles = []
        self.background_color = (255,255,255)  #World background set to white
        self.dt = 0.1     #time step size - Fixed for now
        self.width = 300  
        self.height = 300
        self.force_que = []    #Keeps track of forces, particles, and things like spring constants, etc
        self.g = vect2d(0,1)  #Constant gravitational field constant
        self.numerical = 'rk4'
        self.running = True   #Determines if while look should continue
        self.selected = None   #Particle selected with the mouse
        self.mouse_force = 1    #Coefficient for force exerted by mouse
        self.damping = 0.5    #Damping coefficient
        self.vel_max = 50     #Velocity at which the drag force kicks in
        self.rod_color=(120,120,0)  #NEW: Used for color of connecting rods
        self.setup_world()   #Create the pygame window


    def setup_world(self):
        """Create a pygame window"""
        pygame.init()
        self.screen = pygame.display.set_mode((self.width,self.height))
        pygame.display.set_caption('PHYS-360 Homework 7')
        self.screen.fill(self.background_color)
        self.clock = pygame.time.Clock()
        self.collision = collision_engine(self,self.particle_list)   #This creates an instance of the collision engine


    def set_numerical(self,method):
        """Change the self.numerical variable to use a different update method
        
        Note that 'method' must be passed as a string
        """
        self.numerical = method

    def add_particle(self,particle):
        """Add particle to particle_list"""
        self.particle_list.append(particle)
        self.all_particles.append(particle)


    def add_rod(self,part1,part2):
        """NEW: Adds a rod constraing between two particles with length equal to separation at initial time
        
        The length of the rod is set to the original distance between the two particles.  Particle class
        has new property called rod_connection.  This is a dictionary that uses the name of the particle as a key and the length of the rod 
        as the associated value.  You may want to look up python dictionary. 
        """
        L = (part1.pos - part2.pos).mag()  #Length of rod
        part1.rod_connection[part2] = L
        part2.rod_connection[part1] = L
        

    def new_force(self,force,particles):
        """Add a force to the list of forces calculated
        
        'particles' should be a list of all particles experiencing the force.  Some
        forces require other parameters and those should be included in the 'particles' list
        (e.g. see spring force defintion below)
        """
        self.force_que.append([force,particles])


    def update(self):
        """Update the positions and velocities of all particles
        
        Should also include boundary(), self.game_controls(), self.v_max_check(),
        and should use getattr() to call the numerical method used to solve
        the equations of motion.
        """
        self.game_controls()
        getattr(self,self.numerical)()  #Solve equations of motion
        self.collision.find_collision()  #Check to see if object hits the walls
        self.display()


    def display(self):
        """Draw all particles onto the surface and then flip to computer screen"""
        self.screen.fill(self.background_color)  #Wipe clean the surface in memory
        
        for particle in self.particle_list: #Loops through all of the particles in the world
            x = int(particle.pos.x)
            y = int(particle.pos.y)
            if(particle.shape == "circle"):
                pygame.draw.circle(self.screen, particle.color, (x, y) , int(particle.size*6), 0)  #Draws the particle to the screen
        
        for p in combinations(self.particle_list,2):
            if p[0] in p[1].rod_connection:
                pygame.draw.line(self.screen, self.rod_color, (p[0].pos.x,p[0].pos.y), (p[1].pos.x,p[1].pos.y), 1)
                
        pygame.display.flip()   #Displays the world
        

    def move(self):
        """Updates the position and velocity of particles under influence of a changing force"""
        self.net_force()
        for i in self.particle_list:
            i.vel += (1/i.m)*i.f_net*self.dt
            i.pos += i.vel*self.dt
            if i.shape == "rectangle":
                i.theta += i.omega*self.dt
                i.rotationUpdate()


    def euler(self):
        """Update the position and velocity using the Euler method"""
        self.net_force()
        for i in self.particle_list:
            i.vel = i.vel + (i.f_net/i.m)* self.dt
            i.pos = i.pos + i.vel * self.dt
            if i.shape == "rectangle":
                i.omega *= i.particle.cor
                i.theta += i.omega*self.dt
                i.rotationUpdate()


    def midpoint(self):
        """Update the position and velocity using the Midpoint method"""
        self.net_force()
        for i in self.particle_list:
            i.vel_old = i.vel
            i.vel = i.vel + (i.f_net/i.m)*self.dt
            i.pos = i.pos + 1/2 * (i.vel_old+i.vel)*self.dt
            if i.shape == "rectangle":
                i.theta += i.omega*self.dt
                i.rotationUpdate()


    def verlet(self):
        """Use the Velocity Verlet method to update position and velocity"""
        self.net_force()
        for i in self.particle_list:
            i.pos_old = i.pos
            i.vel_old = i.vel
            i.f_net_old = i.f_net
            
            i.pos = i.pos_old+i.vel_old*self.dt+1/2*i.f_net_old/i.m*self.dt**2
        
        self.net_force()
        for i in self.particle_list:
            i.vel = i.vel_old+self.dt/(2*i.m)*(i.f_net+i.f_net_old)
        
            if i.shape == "rectangle":
                i.theta += i.omega*self.dt
                i.rotationUpdate()


    def projectile(self):
        """Updates the position of a particle moving with a constant net force
        
        Make sure that other position update methods are not called at the same time
        as this one (e.g. either projectile() or move(), not both).
        """
        for i in self.particle_list:
            i.pos += i.vel*self.dt + 1/2*self.g*self.dt**2
            i.vel += self.g*self.dt
            
            if i.shape == "rectangle":
                i.theta += i.omega*self.dt
                i.rotationUpdate()


    def rk4(self):
        """Update the position and velocity using rk4
        
        Will need to use particle.get_force() method and particle.set_pos() at different points to 
        return the net force on each particle and set the positions before calling self.net_force().
        Store k and L results (use capital L so it isn't confused with number one)"""
        self.net_force()
        a1 = 1/6
        a2 = 1/3
        a3 = 1/3
        a4 = 1/6
        for i in self.particle_list:
            i.pos_old = i.pos
            i.vel_old = i.vel
            i.f_net_old = i.f_net
            
            i.K1 = i.vel * self.dt
            i.l1 = (1/i.m) * i.f_net * self.dt
            i.pos = i.pos_old + i.K1
        
        self.net_force()
        for i in self.particle_list:
            
            i.K2 = (i.vel + 1/2 * i.l1) * self.dt
            i.l2 = (1/i.m) * i.f_net * self.dt
            i.pos = i.pos_old + i.K2
        
        self.net_force()
        for i in self.particle_list:
            
            i.K3 = (i.vel + 1/2 * i.l2) * self.dt
            i.l3 = (1/i.m) * i.f_net * self.dt
            i.pos = i.pos_old + i.K3
        
        self.net_force()
        for i in self.particle_list:
            
            i.K4 = (i.vel + i.l3) * self.dt
            i.l4 = (1/i.m) * i.f_net * self.dt
            
            i.pos = i.pos_old + a1*i.K1 + a2*i.K2 + a3*i.K3 + a4*i.K4
            i.vel = i.vel_old + a1*i.l1 + a2*i.l2 + a3*i.l3 + a4*i.l4
            
            if i.shape == "rectangle":
                i.theta += i.omega*self.dt
                i.rotationUpdate()


    def net_force(self):
        """Find net force for all particles
        
        Set net force for each particle to zero at start of the time setup_world
        and then run through the force_que to find net force on each particle
        """
        for i in self.particle_list:
            i.reset_force()   #Reset all net forces to zero
        for f in self.force_que:
            getattr(self,f[0])(f[1])  #Calls all force methods listed in force_que


    def const_grav(self,particles):
        """Constant gravitational field
        
        Don't forget that the gravitational FORCE depends on the mass
        """
        for i in particles:
            i.add_force(i.m*3*self.g)


    def drag(self,particles):
        """Apply drag force to particles
        
        Drag force only kicks in when particle speed exceeds self.vel_max.
        Drag force is opposite direction of velocity and should scale as
        the square of the speed
        """
        for i in particles:
            if i.vel.mag() > self.vel_max:
                i.add_force(-i.drag_coeff*i.vel.mag2()*i.vel)

    def spring(self,particles):
        """Spring force between two particles
        
        Assumes only two particles are passed.  If you have more than
        two particles you will need to midify the code.  The argument 'particles' should be a
        list containing the two particles, the spring constant, and the relaxed separation distance
        """
        L = particles[0].pos - particles[1].pos
        
        sForce = -particles[2] * (L.norm()) * (L.mag() - particles[3])
        
        particles[0].add_force(sForce)
        particles[1].add_force(-sForce)


    def mouse_pull(self,particle):
        """Force applied by selecting particle with mouse"""
        (pick_x,pick_y) = pygame.mouse.get_pos()
        mouse_pos = vect2d(pick_x,pick_y)
        dx = self.selected.pos - mouse_pos
        F = -self.mouse_force*dx - self.selected.vel*self.damping
        self.selected.add_force(F)

    def create_water(self):
        x, y = pygame.mouse.get_pos()
        part = particle(vect2d(x,y),vect2d(0,0))
        self.new_force('const_grav',[part])
        self.particle_list.append(part)


    def game_controls(self):
        """Handles all mouse and keyboard events
        
        Call this in the update() method
        """
        for event in pygame.event.get(): 
            if event.type == pygame.QUIT: #If red 'x' is clicked
                self.running = False
            if event.type == pygame.MOUSEBUTTONDOWN:
                (pick_x,pick_y) = pygame.mouse.get_pos()   #Get mouse position on screen
                picked = vect2d(pick_x,pick_y)  #Turn mouse position into a vector
                for i in self.particle_list:
                    dist = i.pos - picked   #How far the mouse is from the center of a particle
                    if dist.mag() < i.size:  #If mouse click is inside circle
                        self.selected = i
                        self.new_force('mouse_pull',[self.selected])  #Add mouse force to force_que
                        self.selected.color = (255,0,0)   #Change color
            if event.type == pygame.MOUSEBUTTONUP and self.selected != None:
                self.selected.color = (0,0,255)  #change color back to blue
                self.force_que.remove(['mouse_pull',[self.selected]])   #Remove force from force_que
                self.selected = None    #No particles are selected
        if pygame.mouse.get_pressed()[0]:
            self.create_water()


class collision_engine():
    def __init__(self,world, particle_list=[]):
        self.world = world  #World instance used by this collision engine
        self.particle_list = particle_list   #List of all particles in the world
        self.has_collided = []    #List of particles that have collided with each other
        self.cor = 0.45   #Coefficient of restitution for hitting walls


    def find_collision(self):
        """EDIT: Runs through all particles in particle_list to look for overlapping particles
        
        Needs to check and see if particle has a rod connection and if it does, call check_rod_constraint
        """
        self.has_collided = []   #Reset collisions
        self.boundary()    #See if any particles collide with the wall
        for p in combinations(self.particle_list,2):
            self.bounding_sphere(p[0],p[1])
            if p[0] in p[1].rod_connection:
                self.check_rod_constraint(p[0],p[1],p[0].rod_connection[p[1]])
        if len(self.has_collided)> 0:  #If any particles have collided
            self.resolve_collision()

    def bounding_sphere(self,particle1,particle2):
        """Check for overlap using bounding spheres"""
        colDis = particle1.pos-particle2.pos
        if (colDis.mag()) < (particle1.size+particle2.size):
            newList = [particle1,particle2]
            if particle1.shape == "circle" and particle2.shape == "circle":
                self.has_collided.append(newList)
                


    def resolve_collision(self):
        """Changes velocities of point particles that have collided"""
        for i in self.has_collided:
            if(i[0].shape == "circle" and i[1].shape == "circle"):
                colNorm = self.collision_normal(i[0],i[1])
                dx = (i[0].size+i[1].size) - (i[1].pos - i[0].pos).mag()
                self.resolve_interpenetration(i[0],i[1],dx,colNorm)
                vel1_old = i[0].vel
                vel2_old = i[1].vel
                dV1 = ((i[0].m_inv)*(1+self.cor)*((vel2_old-vel1_old)*colNorm))/(i[1].m_inv+i[0].m_inv)
                dV2 = ((i[1].m_inv)*(1+self.cor)*((vel2_old-vel1_old)*colNorm))/(i[1].m_inv+i[0].m_inv)
                i[0].vel += dV1*colNorm
                i[1].vel -= dV2*colNorm
            elif(i[0].shape == "rectangle"):
                i[0].r_c = i[1].x_c*i[0].xPrimeHat+i[1].y_c*i[0].yPrimeHat
                hat = i[0].n_hat.x*i[0].xPrimeHat+i[0].n_hat.y*i[0].yPrimeHat
                dx = (i[0].r_c.mag()+i[1].size) - (i[1].pos - i[0].pos).mag()
                self.resolve_interpenetration(i[1],i[0],dx,hat)
                
                i[0].r_c_perp = vect2d(-i[0].r_c.y,i[0].r_c.x)
                v_c = i[0].vel+vect2d(i[0].omega*i[0].r_c_perp.x,-i[0].omega*i[0].r_c_perp.y)
                v_rel = v_c - i[1].vel
                v_close = v_rel.dot(i[0].n_hat)
                I_rect = (1/12*i[0].m*(i[0].width**2+i[0].height**2))
                J = ((1+self.cor)*v_close)/(i[0].m_inv+i[1].m_inv+(i[0].r_c_perp.dot(i[0].n_hat))**2/I_rect)
                i[0].omega += -(i[0].r_c_perp.dot(i[0].n_hat)*J)/I_rect
                dV1 = i[0].m_inv*J
                dV2 = i[1].m_inv*J
                i[0].vel -= dV1*i[0].n_hat
                i[1].vel += dV2*i[0].n_hat

    def collision_normal(self,particle1,particle2):
        """Returns the collision normal vector"""
        disVector = particle1.pos-particle2.pos
        unitDisVector = disVector.norm()
        return unitDisVector


    def resolve_interpenetration(self,particle1,particle2,dx,hat):
        """Move objects that overlap far enough apart that they don't continue to collide"""
        """
        particle1.pos -= particle1.m_inv / (particle2.m_inv + particle1.m_inv) * (hat)
        particle2.pos += particle2.m_inv / (particle2.m_inv + particle1.m_inv) * (hat)
        """
        particle1.pos += (dx/2*hat)
        particle2.pos -= (dx/2*hat)


    def boundary(self):
        """Checks to see if part of a particle leaves the screen.  Should shift particle back on to screen and reverse component of velocity
        
        This method has been moved from the world class.  It makes more sense to have it be part of the collision engine
        """
        for i in self.particle_list:
            if i.shape == "circle":
                if i.pos.x <= i.size:
                    i.pos.x = i.size + 0.01
                    i.vel.x = -self.cor*i.vel.x
                if i.pos.x >= self.world.width-i.size:
                    i.pos.x = self.world.width-i.size-0.01
                    i.vel.x = -self.cor*i.vel.x
                if i.pos.y <= i.size:
                    i.pos.y = i.size + 0.01
                    i.vel.y = -self.cor*i.vel.y
                if i.pos.y >= self.world.height-i.size:
                    i.pos.y = self.world.height-i.size-0.01
                    i.vel.y = -self.cor*i.vel.y

if __name__ == '__main__': 
    #The following is just to test your code out and make sure it works
    earth = world()
    
    earth.set_numerical('euler')  #Change to the Velocity Verlet method
    
    ###
    #While loop
    ###
    
    while earth.running:   
        time_passed = earth.clock.tick(120)   #Add pygame.time.Clock() to world setup.  This will allow you to set fps
        earth.update()