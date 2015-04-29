
"""
Created on Mon, Feb 6 2015

@author: Todd Zimmerman

Created by Todd Zimmerman for PHYS-360 Spring 2015 at UW-Stout

Template for introducing numerical methods.  Implements euler, midpoint, and the velocity verlet methods
"""

from __future__ import division
import ode
import pygame   #Imports methods from a file called 'pygame'
from vector import *   #Rename your solution to homework 2 as vector.py 
import random      #Needed to generate random numbers
from itertools import combinations
from random import randint
import math

author_name = "Jacob Gregor"

class particle:
	def __init__(self,pos, vel, mass = 1.0, size = 10, shape = "circle", color = (0,0,255)):
		self.pos = pos   #Position vector
		self.vel = vel   #Velocity vectorB
		self.m = mass    #Particle mass
		self.m_inv = 1/mass
		self.d = 0
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
		self.inLiquid = False
		self.submerged = False
		self.area = 1.3333*3.14*size**3#1.3333333 * self.pi * math.pow(size/2, 3) #m^3
		self.drag_coeff = 0.000001   #Scaling number for drag force
		self.cor = .99     #Coefficient of restitution
		self.rod_connection={}  #NEW: Dictionary of all rod constraint connections
		#self.screen,(0,0,255),True,[((particle.pos.x-particle.width/2)*math.cos(particle.theta),particle.pos.y-particle.height/2),(particle.pos.x-particle.width/2,particle.pos.y+particle.height/2),(particle.pos.x+particle.width/2,particle.pos.y+particle.height/2),(particle.pos.x+particle.width/2,particle.pos.y-particle.height/2)


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

class liquid(particle):
    def __init__(self,world,pos= vect2d(0,0) ,volume = 0.0001,height=100,width=100,shape = "liquid", color = (0,0,255)):
        self.width = world.width-10
	self.height = world.height/2
        self.pos = vect2d(world.width/2,world.height*7/8)
	
	
	
	#self.size = size #Radius for circle and width and height for rectangle
	self.color = color
	self.shape = shape  #String with the name of the shape
	self.f_net = vect2d(0,0)   #Property to keep track of net force on particle
	
	
	self.f_net = vect2d(0,0)   #Property to keep track of net force on particle
	self.pos_old = vect2d(0,0)
	self.vel_old = vect2d(0,0)
	self.k1 = vect2d(0,0)
	self.k2 = vect2d(0,0)
	self.k3 = vect2d(0,0)
	self.k4 = vect2d(0,0)
	self.L1 = vect2d(0,0)
	self.L2 = vect2d(0,0)
	self.L3 = vect2d(0,0)
	self.L4 = vect2d(0,0)
	#(particle.pos.x-particle.width/2,particle.pos.y-particle.height/2),(particle.pos.x-particle.width/2,particle.pos.y+particle.height/2),(particle.pos.x+particle.width/2,particle.pos.y+particle.height/2),(particle.pos.x+particle.width/2,particle.pos.y-particle.height/2)
	self.UL = vect2d(int(self.pos.x - self.width/2),int(self.pos.y - self.height/2))
	self.UR = vect2d(int(self.pos.x + self.width/2),int(self.pos.y - self.height/2))
	self.LL = vect2d(int(self.pos.x - self.width/2),int(self.pos.y + self.height/2))
	self.LR = vect2d(int(self.pos.x + self.width/2),int(self.pos.y + self.height/2))
	self.r_c = vect2d(0,0)
	self.xPrime = 0
        self.yPrime = 0
        self.xPrimeHat = 0
        self.yPrimeHat = 0
        self.n_hat = 0
	self.size = (self.UL-self.LR).mag()/2
	self.drag_coeff = 0.000001   #Scaling number for drag force
	self.cor = .99     #Coefficient of restitution
	self.rod_connection={}  #NEW: Dictionary of all rod constraint connections
	
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

class rectangle(particle):
    def __init__(self,pos, vel, theta=0,omega=0, mass = 1.0,height=100,width=100,shape = "rectangle", color = (0,0,255)):
        self.pos = pos
	self.vel = vel   #Velocity vectorB
	self.m = mass    #Particle mass
	self.m_inv = 1/mass
	#self.size = size #Radius for circle and width and height for rectangle
	self.color = color
	self.shape = shape  #String with the name of the shape
	self.f_net = vect2d(0,0)   #Property to keep track of net force on particle
	self.theta = theta
	self.omega = omega
	self.width = width
	self.height = height
	

	self.f_net = vect2d(0,0)   #Property to keep track of net force on particle
	self.pos_old = vect2d(0,0)
	self.vel_old = vect2d(0,0)
	self.k1 = vect2d(0,0)
	self.k2 = vect2d(0,0)
	self.k3 = vect2d(0,0)
	self.k4 = vect2d(0,0)
	self.L1 = vect2d(0,0)
	self.L2 = vect2d(0,0)
	self.L3 = vect2d(0,0)
	self.L4 = vect2d(0,0)
	#(particle.pos.x-particle.width/2,particle.pos.y-particle.height/2),(particle.pos.x-particle.width/2,particle.pos.y+particle.height/2),(particle.pos.x+particle.width/2,particle.pos.y+particle.height/2),(particle.pos.x+particle.width/2,particle.pos.y-particle.height/2)
	self.UL = vect2d(int(self.pos.x - self.width/2),int(self.pos.y - self.height/2))
	self.UR = vect2d(int(self.pos.x + self.width/2),int(self.pos.y - self.height/2))
	self.LL = vect2d(int(self.pos.x - self.width/2),int(self.pos.y + self.height/2))
	self.LR = vect2d(int(self.pos.x + self.width/2),int(self.pos.y + self.height/2))
	self.r_c = vect2d(0,0)
	self.xPrime = 0
        self.yPrime = 0
        self.xPrimeHat = 0
        self.yPrimeHat = 0
        self.n_hat = 0
	self.size = (self.UL-self.LR).mag()/2
	self.drag_coeff = 0.000001   #Scaling number for drag force
	self.cor = .99     #Coefficient of restitution
	self.rod_connection={}  #NEW: Dictionary of all rod constraint connections
	
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
		self.liquid_list = []
		self.all_particles = []
		self.background_color = (255,255,255)  #World background set to white
		self.dt = 0.1     #time step size - Fixed for now
		self.width = 800  
		self.height = 600
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
		self.collision = collision_engine(self,self.particle_list,self.all_particles)   #This creates an instance of the collision engine



	def set_numerical(self,method):
		"""Change the self.numerical variable to use a different update method

		Note that 'method' must be passed as a string
		"""
		self.numerical = method

	def add_particle(self,particle):
		"""Add particle to particle_list"""
		if particle.shape == "liquid":
		    self.liquid_list.append(particle)
		    self.all_particles.append(particle)
		else:    
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
		self.screen.fill(earth.background_color)  


	def display(self):
		"""Draw all particles onto the surface and then flip to computer screen"""
                for particle in self.particle_list:
                    if particle.shape == "circle":
		      pygame.draw.circle(self.screen, (0,0,255), (int(particle.pos.x), int(particle.pos.y)), particle.size, 1)
		    if particle.shape =="rectangle":
		      particle.theta = particle.theta + particle.omega
		      #particle.theta = particle.theta+particle.omega*self.dt  
		      print particle.pos.x+particle.width/2, particle.UL.x
		      
		      #pygame.draw.aalines(self.screen,(0,0,255),True,[(particle.UL.x,particle.UL.y),(particle.UR.x,particle.UR.y),(particle.LR.x,particle.LR.y),(particle.LL.x,particle.LL.y)])
		
		for particle in self.liquid_list:
  		    if particle.shape =="liquid":
  		      #pygame.draw.rect(self.screen,(0,0,255),[(int(particle.UL.x),int(particle.UL.y)),(int(particle.UR.x),int(particle.UR.y)),(int(particle.LR.x),int(particle.LR.y)),(int(particle.LL.x),int(particle.LL.y)),0])  
    		      pygame.draw.aalines(self.screen,(0,0,255),True,[(particle.UL.x,particle.UL.y),(particle.UR.x,particle.UR.y),(particle.LR.x,particle.LR.y),(particle.LL.x,particle.LL.y)],False)
		    #cool beans
		   # width = 100
		    #cool beans
		   # width = 100
	            #height = 100
	            
	            #end of coolbeans
	        
		for p in combinations(self.particle_list,2):
		    
		    if p[0] in p[1].rod_connection:
		     pygame.draw.line(self.screen,self.rod_color,(p[0].pos.x,p[0].pos.y),(p[1].pos.x,p[1].pos.y),1)
		     #pygame.draw.aalines(self.screen,self.rod_color,)
		pygame.display.flip()



	def move(self):
		"""Updates the position and velocity of particles under influence of a changing force"""
            

	def euler(self):
		"""Update the position and velocity using the Euler method"""
                self.net_force()
		for particle in self.particle_list:
		    particle.pos = particle.pos+particle.vel*self.dt
		    particle.vel =particle.vel+ particle.f_net/particle.m *self.dt
		    if particle.shape == "rectangle":
		        #print particle.UL.x,particle.UL.y
		        particle.UL = vect2d(particle.pos.x+particle.width/2*math.cos(particle.theta)-(particle.height/2)*math.sin(particle.theta),particle.pos.y+particle.height/2*math.cos(particle.theta)+(particle.width/2)*math.sin(particle.theta))
	                particle.UR = vect2d(particle.pos.x-particle.width/2*math.cos(particle.theta)-(particle.height/2)*math.sin(particle.theta),particle.pos.y+particle.height/2*math.cos(particle.theta)-(particle.width/2)*math.sin(particle.theta))
	                particle.LL = vect2d(particle.pos.x+particle.width/2*math.cos(particle.theta)+(particle.height/2)*math.sin(particle.theta),particle.pos.y-particle.height/2*math.cos(particle.theta)+(particle.width/2)*math.sin(particle.theta))
	                particle.LR = vect2d(particle.pos.x-particle.width/2*math.cos(particle.theta)+(particle.height/2)*math.sin(particle.theta),particle.pos.y-particle.height/2*math.cos(particle.theta)-(particle.width/2)*math.sin(particle.theta))
	                
		        #particle.pos.x=particle.pos.x+particle.pos.x*math.cos(particle.theta)+(-particle.pos.y*math.sin(particle.theta))
		        #particle.pos.y=particle.pos.y+particle.pos.x*math.sin(particle.theta)+(particle.pos.y*math.cos(particle.theta))
		        #particle.theta = particle.theta+particle.omega*self.dt
		        print 1
                

	def midpoint(self):
		"""Update the position and velocity using the Midpoint method"""
                self.net_force()
		for particle in self.particle_list:
		    particle.vel_Old = particle.vel
		    particle.vel = particle.vel+ particle.f_net/particle.m *self.dt
		    particle.pos = particle.pos+(1/2)*(particle.vel+particle.vel_Old)*self.dt

                
	def verlet(self):
		"""Use the Velocity Verlet method to update position and velocity"""
                self.net_force()
		for particle in self.particle_list:
		    f_Old = particle.f_net
		    particle.pos = particle.pos + particle.vel*self.dt+ 1/2*particle.f_net/particle.m*self.dt**2
		    self.net_force()
		    f_New = particle.f_net
		    particle.vel = particle.vel + self.dt/(2*particle.m)*(f_New+f_Old)
		    if particle.shape == "rectangle":
		        #print particle.UL.x,particle.UL.y
		        particle.UL = vect2d(particle.pos.x+particle.width/2*math.cos(particle.theta)-(particle.height/2)*math.sin(particle.theta),particle.pos.y+particle.height/2*math.cos(particle.theta)+(particle.width/2)*math.sin(particle.theta))
	                particle.UR = vect2d(particle.pos.x-particle.width/2*math.cos(particle.theta)-(particle.height/2)*math.sin(particle.theta),particle.pos.y+particle.height/2*math.cos(particle.theta)-(particle.width/2)*math.sin(particle.theta))
	                particle.LL = vect2d(particle.pos.x+particle.width/2*math.cos(particle.theta)+(particle.height/2)*math.sin(particle.theta),particle.pos.y-particle.height/2*math.cos(particle.theta)+(particle.width/2)*math.sin(particle.theta))
	                particle.LR = vect2d(particle.pos.x-particle.width/2*math.cos(particle.theta)+(particle.height/2)*math.sin(particle.theta),particle.pos.y-particle.height/2*math.cos(particle.theta)-(particle.width/2)*math.sin(particle.theta))

	def projectile(self):
		"""Updates the position of a particle moving with a constant net force

		Make sure that other position update methods are not called at the same time
		as this one (e.g. either projectile() or move(), not both).
		"""
                for i in self.particle_list:
			i.pos += i.vel*self.dt + 1/2*self.g*self.dt**2
			i.vel += self.g*self.dt

	def rk4(self):
		"""Update the position and velocity using rk4

		Will need to use particle.get_force() method and particle.set_pos() at different points to 
		return the net force on each particle and set the positions before calling self.net_force().
		Store k and L results (use capital L so it isn't confused with number one)"""
		
                for particle in self.particle_list:
		   particle.pos_Old = particle.pos
		   self.net_force()
		   
		   particle.L1 = self.dt/particle.m*particle.f_net
		   particle.k1 = particle.vel*self.dt
		   
		   particle.pos = particle.pos_Old+particle.L1/2
		   self.net_force()
		   
		   particle.L2 = self.dt/particle.m*particle.f_net
		   particle.k2 = (particle.vel+(1/2)*particle.L1)*self.dt
		   
		   particle.pos = particle.pos_Old+particle.L2/2
		   self.net_force()
		   
		   particle.L3 = self.dt/particle.m*particle.f_net
		   particle.k3 = (particle.vel+(1/2)*particle.L2)*self.dt
		   
		   particle.pos = particle.pos_Old+particle.L3/2
		   self.net_force()
		   
		   particle.L4 = self.dt/particle.m*particle.f_net
		   particle.k4 = (particle.vel+(1/2)*particle.L3)*self.dt
		   
		   particle.pos = particle.pos + (1/6)*(particle.k1+2*particle.k2+2*particle.k3+particle.k4)
		   particle.vel = particle.vel + (1/6)*(particle.L1+2*particle.L2+2*particle.L3+particle.L4)
		   if particle.shape == "rectangle":
		        #print particle.UL.x,particle.UL.y
		        particle.UL = vect2d(particle.pos.x+particle.width/2*math.cos(particle.theta)-(particle.height/2)*math.sin(particle.theta),particle.pos.y+particle.height/2*math.cos(particle.theta)+(particle.width/2)*math.sin(particle.theta))
	                particle.UR = vect2d(particle.pos.x-particle.width/2*math.cos(particle.theta)-(particle.height/2)*math.sin(particle.theta),particle.pos.y+particle.height/2*math.cos(particle.theta)-(particle.width/2)*math.sin(particle.theta))
	                particle.LL = vect2d(particle.pos.x+particle.width/2*math.cos(particle.theta)+(particle.height/2)*math.sin(particle.theta),particle.pos.y-particle.height/2*math.cos(particle.theta)+(particle.width/2)*math.sin(particle.theta))
	                particle.LR = vect2d(particle.pos.x-particle.width/2*math.cos(particle.theta)+(particle.height/2)*math.sin(particle.theta),particle.pos.y-particle.height/2*math.cos(particle.theta)-(particle.width/2)*math.sin(particle.theta))
                   #print particle.pos.x,particle.shape



	def net_force(self):
		"""Find net force for all particles

		Set net force for each particle to zero at start of the time setup_world
		and then run through the force_que to find net force on each particle
		"""
		for i in self.particle_list:
			i.reset_force()   #Reset all net forces to zero
		for f in self.force_que:
		        #print f[0]
			getattr(self,f[0])(f[1])  #Calls all force methods listed in force_que


	def const_grav(self,particles):
		"""Constant gravitational field

		Don't forget that the gravitational FORCE depends on the mass
		"""
                for particle in particles:
                    #print"grav"
		    force = vect2d(0,0)
		    force = particle.m * self.g
		    particle.f_net= particle.f_net + force
        def object_buoyancy(self,particles):
         
            for particle in particles:
               
                #print particle.vel
                vol = 1.3333333 * 3.14 * math.pow(particle.size, 3)
                #print particle.d
                buoyancy = .0005 * particle.area * self.g 
                if particle.inLiquid == True:    
                    particle.f_net= particle.f_net - buoyancy
                    #print particle.f_net,"after" 
                elif particle.inLiquid == False:
                    particle.f_net = particle.f_net
                #return buoyancy
            
            
        def object_gravity(self, mass):
            gravity_force = mass * self.gravity #N
            return gravity_force
    
        
        def surface_area(self, radius):
            area = 4 * self.pi * math.pow(radius, 2)
            return area
	def drag(self,particles):
		"""Apply drag force to particles

		Drag force only kicks in when particle speed exceeds self.vel_max.
		Drag force is opposite direction of velocity and should scale as
		the square of the speed
		"""
		for i in particles:
		  if i.inLiquid == True:
		      p = .5
		      i.drag_coeff = .05
		      dragF = p*i.vel.mag2()*-i.drag_coeff*i.vel.norm()
		      print dragF,"drag"
		      i.f_net = i.f_net + dragF
		  elif i.vel.mag() > self.vel_max:
		      i.drag_coeff = .000001
		      i.add_force(-i.drag_coeff*i.vel.mag2()*i.vel)

	def spring(self,particles):
		"""Spring force between two particles

		Assumes only two particles are passed.  If you have more than
		two particles you will need to midify the code.  The argument 'particles' should be a
		list containing the two particles, the spring constant, and the relaxed separation distance
		"""
                k = particles[2]
		l0 = particles[3]
		
		length = vect2d(0,0)
		length = particles[0].pos - particles[1].pos
		
		lMag = length.mag()
		lHat = vect2d(0,0)
		lHat.x = length.x/lMag
		lHat.y = length.y/lMag
		s = length.mag() - l0
		
		particles[0].f_net.x += -k * s*lHat.x
		particles[0].f_net.y += -k * s*lHat.y
		
		particles[1].f_net.x += k * s*lHat.x
		particles[1].f_net.y += k * s*lHat.y

	def mouse_pull(self,particle):
		"""Force applied by selecting particle with mouse"""
		(pick_x,pick_y) = pygame.mouse.get_pos()
		mouse_pos = vect2d(pick_x,pick_y)
		dx = self.selected.pos - mouse_pos
		F = -self.mouse_force*dx - self.selected.vel*self.damping
		self.selected.add_force(F)


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

class collision_engine():
	def __init__(self,world, particle_list=[],all_particles=[]):
		self.world = world  #World instance used by this collision engine
		self.particle_list = particle_list   #List of all particles in the world
		self.all_particles = all_particles
		self.has_collided = []    #List of particles that have collided with each other
		self.has_collidedR = []
		self.cor = 0.99   #Coefficient of restitution for hitting walls


	def find_collision(self):
		"""EDIT: Runs through all particles in particle_list to look for overlapping particles

		Needs to check and see if particle has a rod connection and if it does, call check_rod_constraint
		"""
		self.has_collided = []   #Reset collisions
		self.has_collidedR = []
		self.boundary()    #See if any particles collide with the wall
		for p in combinations(self.all_particles,2):
		    if p[0].shape =="liquid":
		      print p[1].pos.y - p[0].UR.y
		      if p[0].UR.y<=p[1].pos.y:
		          #print p[0].UR.y,"<- liquid  Particle ->",p[1].pos.y
		          p[1].inLiquid = True
		          #print p[1].pos.y - p[0].UR.y
		          if p[1].pos.y - p[0].UR.y>= p[1].size:
		           
		              p[1].submerged = True
		          p[1].d = ((p[1].pos.y - p[0].UR.y)/(self.world.height-p[0].UR.y))
		          #print p[1].d, "sswag"
		          #print p[1].d, "depth"
		          #print p[1].inLiquid
		      else:
		          p[1].inLiquid = False
		          p[1].submerged = True
		          #print p[1].inLiquid
		    '''elif p[1].shape == "liquid":
		        print"lose"
		    else:
		        print"FUCK"'''
		        
		for p in combinations(self.particle_list,2):
			self.bounding_sphere(p[0],p[1])
			if p[0] in p[1].rod_connection:
				self.check_rod_constraint(p[0],p[1],p[0].rod_connection[p[1]])
		
		if len(self.has_collided)> 0:  #If any particles have collided
			self.resolve_collision()

	def bounding_sphere(self,particle1,particle2):
		"""Check for overlap using bounding spheres"""
		if particle1.shape == "circle" and particle2.shape == "circle":
                    combinedSize = particle1.size + particle2.size
  		    if(abs(particle1.pos.x-particle2.pos.x) <= combinedSize):
  		        if(abs(particle1.pos.y-particle2.pos.y) <= combinedSize):
  		            self.has_collided.append([particle1,particle2])
                        self.resolve_collision()
                elif particle1.shape == "rectangle":
                    particle1.xPrimeHat = vect2d(math.cos(particle1.theta)+0*-math.sin(particle1.theta),math.sin(particle1.theta)+0*math.cos(particle1.theta))
                    particle1.yPrimeHat = vect2d(0*math.cos(particle1.theta)+-math.sin(particle1.theta),0*math.sin(particle1.theta)+math.cos(particle1.theta))
                    rRel = particle2.pos - particle1.pos
                    particle2.x_body = rRel.dot(particle1.xPrimeHat)
                    particle2.y_body = rRel.dot(particle1.yPrimeHat)
                    values = [abs(particle2.x_body-particle1.width/2),abs(particle2.x_body+particle1.width/2),abs(particle2.y_body-particle1.height/2),abs(particle2.y_body+particle1.height/2)]
                    if min(values)<=particle2.size:
                        
                        if(-particle1.width/2 < particle2.x_body and particle2.x_body < particle1.width/2):
                            particle2.x_c = particle2.x_body
                        elif (particle2.x_body > 0):
                            particle2.x_c = particle1.width/2
                        elif (particle2.x_body < 0):
                            particle2.x_c = -particle1.width/2
                        
                        if(-particle1.height/2 < particle2.y_body and particle2.y_body < particle1.height/2):
                            particle2.y_c = particle2.y_body
                        elif (particle2.y_body > 0):
                            particle2.y_c = particle1.height/2
                        elif (particle2.y_body < 0):
                            particle2.y_c = -particle1.height/2
                        
                        if values.index(min(values)) == 0:
                            particle1.n_hat = vect2d(1,0)
                        elif values.index(min(values)) == 1:
                            particle1.n_hat = vect2d(-1,0)
                        elif values.index(min(values)) == 2:
                            particle1.n_hat = vect2d(0,1)
                        elif values.index(min(values)) == 3:
                            particle1.n_hat = vect2d(0,-1)
                        
                        self.has_collided.append([particle1,particle2])
                    
	def resolve_collision(self):
		"""Changes velocities of point particles that have collided"""
                for particles in self.has_collided:
                    if particles[0].shape == "rectangle" or particles[1].shape == "rectangle":
                        self.rectangle_collision()
                    else:
      		    #print particles[0],particles[1]
          		    rHat = self.collision_normal(particles[0],particles[1])
          		    
          		    difference = particles[1].vel-particles[0].vel
          		    dv1 = (particles[1].m*(particles[0].cor)*(difference.dot(rHat)))/(particles[0].m+particles[1].m)
            		  
          		    dv2 = (-particles[0].m*(particles[1].cor)*(difference.dot(rHat)))/(particles[0].m+particles[1].m)
          		    
          		    particles[0].vel.x = particles[0].vel.x +dv1*rHat.x
          		    particles[0].vel.y = particles[0].vel.y +dv1*rHat.y
          		    
          		    particles[1].vel.x = particles[1].vel.x +dv2*rHat.x
          		    particles[1].vel.y = particles[1].vel.y +dv2*rHat.y
          		    
          		    self.resolve_interpenetration(particles[1],particles[0],rHat)
        def rectangle_collision(self):
		"""Changes velocities of point particles that have collided"""
                for i in self.has_collided:
                    i[0].r_c = i[1].x_c*i[0].xPrimeHat+i[1].y_c*i[0].yPrimeHat
                    hat = i[0].n_hat.x*i[0].xPrimeHat+i[0].n_hat.y*i[0].yPrimeHat
                    dx = (i[0].r_c.mag()+i[1].size) - (i[1].pos - i[0].pos).mag()
                    #self.resolve_interpenetration(i[1],i[0],dx,hat)
                    
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
                    '''
                    r_C = vect2d(r_c[0],r_c[1])
                    n_Hat = vect2d(n_hat[0],n_hat[1])
                    
                    if p1.shape == "rectangle":
                        contact1 = r_C-p1.pos
                    else:
                        contact1 = r_C-p2.pos
                    r_perp = vect2d(contact1.y, contact1.x)
                    v_rel = p1.vel - p2.vel
                    v_close=v_rel.dot(n_Hat)
                    C=(p1.cor*p2.cor)
                    
                    if p1.shape == "rectangle": 
                        I = (1/12) * p1.m*(p1.width**2+p1.height**2)
                    else:
                        I = (1/12) * p2.m*(p1.width**2+p1.height**2)
                    #print v_rel.dot(n_Hat),"swaggg"
                    J=((1+C)*v_close)/((1/p1.m)+(1/p2.m)+(1/I)*(r_perp.dot(n_Hat))**2)
                
                    delta_v1 = J/p1.m
                    delta_v2 = -J/p2.m
                    p1.vel = p1.vel +delta_v1
                    p2.vel = -p2.vel +delta_v2
                    
                    
                    delta_v1 = (-J/p1.m)
                    delta_v2 = (J/p2.m)
                    
                    print "test", J
                    
                    p1.vel.x = p1.vel.x +delta_v1*n_hat[0]
                    p1.vel.y = p1.vel.y +delta_v1*n_hat[1]
                    p2.vel.x = p2.vel.x +delta_v2*n_hat[0]
                    p2.vel.y = p2.vel.y +delta_v2*n_hat[1]
                    
                    
		    #print particles[0],particles[1]
		    r_rel = vect2d(0,0)
		    r_rel.x = particles[1].pos.x-particles[0].pos.x
		    r_rel.y = particles[1].pos.y-particles[0].pos.y
		    nHat = self.collision_normal(particles[0],particles[1])
		    xhat=vect2d(math.cos(particles[0].theta),-math.sin(particles[0].theta))
		    yhat=vect2d(math.sin(particles[0].theta),math.cos(particles[0].theta))
                    
                                        
                    
		    I = (1/12)*particles[0].m*(particles[0].width+particles[0].height)**2
		    
		    rHat = self.collision_normal(particles[0],particles[1])
		    
		    difference = particles[1].vel-particles[0].vel
		    dv1 = (particles[1].m*(particles[0].cor)*(difference.dot(rHat)))/(particles[0].m+particles[1].m)
		  
		    dv2 = (-particles[0].m*(particles[1].cor)*(difference.dot(rHat)))/(particles[0].m+particles[1].m)
		    
		    #print "Dvx =",dv1.x,dv2.x
		    
		    #print "dot =",dv1.dot(rHat)
		    particles[0].vel.x = particles[0].vel.x +dv1*rHat.x
		    particles[0].vel.y = particles[0].vel.y +dv1*rHat.y
		    #particles[0].vel = particles[0].vel+dv1.dot(rHat)
		    #particles[1].vel = particles[1].vel+dv2.dot(rHat)
		    particles[1].vel.x = particles[1].vel.x +dv2*rHat.x
		    particles[1].vel.y = particles[1].vel.y +dv2*rHat.y
		    
		    self.resolve_interpenetration(particles[1],particles[0],rHat)
		    '''
	def collision_normal(self,particle1,particle2):
		"""Returns the collision normal vector"""
                r = particle1.pos -particle2.pos
		
		#if r.x == 0 and r.y == 0:
		 #   return vect2d(1,0)
		rHat= r.norm()
		
		return r.norm()

	def resolve_interpenetration(self,particle1,particle2,hat):
		"""Move objects that overlap far enough apart that they don't continue to collide"""
                mag = (particle1.pos -particle2.pos).mag()
		dx = particle1.size + particle2.size - mag
		#print "normhat and norm norm",hat.x, norm.x
		dP1 = (-particle2.m*dx)/(particle1.m+particle2.m)
		dP2 = (particle1.m*dx/(particle1.m+particle2.m))
		
		particle1.pos.x = particle1.pos.x + dP1*hat.x
		particle1.pos.y = particle1.pos.y + dP1*hat.y
		
		particle2.pos.x = particle2.pos.x + dP2*hat.x
		particle2.pos.y = particle2.pos.y + dP2*hat.y

	def boundary(self):
		"""Checks to see if part of a particle leaves the screen.  Should shift particle back on to screen and reverse component of velocity

		This method has been moved from the world class.  It makes more sense to have it be part of the collision engine
		"""
		for particle in self.particle_list:
		    if particle.shape == "circle" :
              		if particle.pos.x > self.world.width - particle.size:
              		    particle.vel.x *= -1
              		    particle.pos.x = self.world.width - particle.size
              		elif particle.pos.x < 0+particle.size:
             		     particle.vel.x *= -1
             		     particle.pos.x = particle.size
             		     
          		if particle.pos.y > self.world.height - particle.size:
          		    particle.pos.y = self.world.height-particle.size
              		    particle.vel.y *= -1
              		elif particle.pos.y < 0+particle.size:
             		     particle.pos.y = particle.size
             		     particle.vel.y *= -1
             		     
                    elif particle.shape == "rectangle":
                        
                        maxW = max(particle.UL.x,particle.UR.x,particle.LL.x,particle.LR.x)
                        minW = min(particle.UL.x,particle.UR.x,particle.LL.x,particle.LR.x)
                        maxY = max(particle.UL.y,particle.UR.y,particle.LL.y,particle.LR.y)
                        minY = min(particle.UL.y,particle.UR.y,particle.LL.y,particle.LR.y)
                        bWidth = (maxW-minW)/2
                        bHeight = (maxY -minY)/2
                        
                        #print "Max", particle.width/2,bWidth/2
                        if particle.pos.x > self.world.width - bWidth:
          		    particle.vel.x *= -1*particle.cor
          		    particle.pos.x = self.world.width - bWidth
          		elif particle.pos.x < 0+bWidth:
         		     particle.vel.x *= -1*particle.cor
         		     particle.pos.x = bWidth
         		     
         		if particle.pos.y > self.world.height - bHeight:
         		    particle.pos.y = self.world.height-bHeight
          		    particle.vel.y *= -1*particle.cor
          		elif particle.pos.y < 0+bHeight:
          		     particle.pos.y = bHeight
         		     particle.vel.y *= -1*particle.cor
         	       


	def check_rod_constraint(self,part0, part1,L):
		"""NEW: Check to see if particles are too close or too far apart.  Call fix_rod_constraint if Needed

		L is the length of the rod.  To get the length L you will need to use part0.rod_connection[part1].  This uses the rod_connection dictionary to look
		up the associated value, which is the rod length L.
		You will need to calculate the vector along the rod (hat) and
		the current distance between the particles (dist) to pass to fix_rod_constraint.
		"""
	        if  L != (part0.pos-part1.pos).mag():
           	        hat = (part0.pos-part1.pos).norm()
           	        dist = (part0.pos-part1.pos).mag()
           	        self.fix_rod_constraint(part0,part1,L, hat, dist)
           	        #print dist, L
		    

	def fix_rod_constraint(self,part0,part1,L,hat,dist):
		"""NEW: Move particle and change velocity so rod constraints are maintained

		L is length of rod, hat is the vector pointing along the rod, and dist is the
		current distance between the particles.  You can use the same impulse equations
		used during collisions to change the velocities.  To move the particles use the 
		interpentration resolution formula.
		"""
		
        	difference = part1.vel-part0.vel
		dv1 = (part1.m*(part0.cor)*(difference.dot(hat)))/(part0.m+part1.m)
		  
		dv2 = (-part0.m*(part1.cor)*(difference.dot(hat)))/(part0.m+part1.m)
		    
		#print "Dvx =",dv1.x,dv2.x
		    
		#print "dot =",dv1.dot(rHat)
		part0.vel.x = part0.vel.x +dv1*hat.x
		part0.vel.y = part0.vel.y +dv1*hat.y
		
		part1.vel.x = part1.vel.x +dv2*hat.x
		part1.vel.y = part1.vel.y +dv2*hat.y
		
		#Find difference between original length and current length
		dx = L - dist
		
		#find out the change in position
		dP1 = (part1.m*dx)/(part0.m+part1.m)
		dP2 = (-part0.m*dx/(part0.m+part1.m))
		
		#update change in position by multiplying it by direction.
		part0.pos.x = part0.pos.x + dP1*hat.x
		part0.pos.y = part0.pos.y + dP1*hat.y
		
		part1.pos.x = part1.pos.x + dP2*hat.x
		part1.pos.y = part1.pos.y + dP2*hat.y


if __name__ == '__main__': 
	#The following is just to test your code out and make sure it works
	earth = world()

	earth.set_numerical('rk4')  #Change to the Velocity Verlet method
       
	part1 = liquid(earth)
	
	r4 = vect2d(300,100)
	v4 = vect2d(0,0)
	part4 = particle(r4,v4)   #Create another particle
	
	r5 = vect2d(310,500)
	v5 = vect2d(0,0)
	part5 = particle(r5,v5)   #Create another particle
	
	earth.add_particle(part1)
        earth.add_particle(part4)
        earth.add_particle(part5)
	num_part = 4  #Set number of particles to add to world
	part_list = []
	'''for i in range(num_part):
		r = vect2d(randint(0,700),randint(0,600))
		v = vect2d(randint(-5,5),randint(-5,5))
		part = particle(r,v)   #Create new particle
		part_list.append(part)    #Add to list of particles that have been randomly placed
		earth.add_particle(part)  #Add particle to the world
	earth.new_force('const_grav',part_list)   #Particles experience gravity
	earth.new_force('drag',part_list)   #particles also experience drag'''
        earth.new_force('const_grav',[part4,part5])  #All particles experience gravity
	#earth.new_force('spring',[part1,part2,.1,30])  #Connct both particles by a spring
	
	earth.new_force('drag',[part4,part5])
	earth.new_force("object_buoyancy",[part4,part5])
	
	
	    
	'''earth.add_rod(part1,part2)
	earth.add_rod(part2,part3)
	earth.add_rod(part3,part4)
	earth.add_rod(part4,part1)
	
	earth.add_rod(part3,part1)
	earth.add_rod(part4,part2)
	
	earth.add_rod(part5,part6)'''
	
	
	


	###
	#While loop
	###

	while earth.running:   
	  
	    
	   time_passed = earth.clock.tick(120)   #Add pygame.time.Clock() to world setup.  This will allow you to set fps
	   earth.update()

