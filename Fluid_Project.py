from __future__ import division
import os
import pygame   #Imports methods from a file called 'pygame'
from pygame.locals import *
from vector import *   #Rename your solution to homework 2 as vector.py 
import random      #Needed to generate random numbers
from itertools import combinations
from random import randint
import sys
import math


class LiquidTest(object):
    '''
    Class used to create water particles, the grid they are based on, and determine how they interact with eachother
    and most of their surroundings.
    '''
    def __init__(self, gsizeX, gsizeY, particlesX, particlesY):
        self.particles = []
        self.gsizeX = gsizeX
        self.gsizeY = gsizeY
        self.particlesX = particlesX
        self.particlesY = particlesY
        self.active = []
        self.water = Material(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.pressed = False
        self.pressedprev = False
        self.mx = 0
        self.my = 0
        self.mxprev = 0
        self.myprev = 0
        i = j = 0
        self.grid = []

        for i in xrange(gsizeX):
            self.grid.append([])

            for j in xrange(gsizeY):
                self.grid[i].append(Node())

        for i in xrange(particlesX):
            for j in xrange(particlesY):
                p = Water_Particle(self.water, i + 4, j + 64, 0.0, 0.0,"water_particle")
                self.particles.append(p)

    def draw(self, screen):
        """Method to draw blue water particles on screen"""
        for p in self.particles:
            start = [4.0 * p.x, 4.0 * p.y]
            end = [4.0 * (p.x - p.u), 4.0 * (p.y - p.v)]
            mid = [(start[0]+end[0])/2,(start[1]+end[1])/2]
            start = map(int, start)
            end = map(int, end)
            mid = map(int, mid)
            pygame.draw.circle(screen, Color('blue'), mid, 5, 0)
    
    def simulate(self):
        '''
        Simulates the movement of water in the nodes that are in grid[].
        Checks to see how each water particle should move, where it is on the grid and how it should move to its surrounding nodes.
        '''
        drag = False
        mdx = mdy = 0.0

        if self.pressed and self.pressedprev:
            drag = True
            mdx = 0.25 * (self.mx - self.mxprev)
            mdy = 0.25 * (self.my - self.myprev)

        self.pressedprev = self.pressed
        self.mxprev = self.mx
        self.myprev = self.my

        for a in self.active:
            a.clear()
        self.active = []
        
        fx = fy = 0.0
        
        for p in self.particles:
            p.cx = int(p.x - 0.5)
            p.cy = int(p.y - 0.5)

            x = p.cx - p.x
            p.px[0] = (0.5 * x * x + 1.5 * x + 1.125)
            p.gx[0] = (x + 1.5)
            x += 1.0
            p.px[1] = (-x * x + 0.75)
            p.gx[1] = (-2.0 * x)
            x += 1.0
            p.px[2] = (0.5 * x * x - 1.5 * x + 1.125)
            p.gx[2] = (x - 1.5)

            y = p.cy - p.y
            p.py[0] = (0.5 * y * y + 1.5 * y + 1.125)
            p.gy[0] = (y + 1.5)
            y += 1.0
            p.py[1] = (-y * y + 0.75)
            p.gy[1] = (-2.0 * y)
            y += 1.0
            p.py[2] = (0.5 * y * y - 1.5 * y + 1.125)
            p.gy[2] = (y - 1.5)


            for i in xrange(3):
                for j in xrange(3):
                    n = self.grid[p.cx + i][p.cy + j]
                    if not n.active:
                        self.active.append(n)
                        n.active = True

                    phi = p.px[i] * p.py[j]
                    n.m += phi * p.mat.m
                    n.d += phi
                    n.gx += p.gx[i] * p.py[j]
                    n.gy += p.px[i] * p.gy[j]

        
        for p in self.particles:

            cx = int(p.x)
            cy = int(p.y)
            cxi = cx + 1
            cyi = cy + 1

            n01 = self.grid[cx][cy]
            n02 = self.grid[cx][cyi]
            n11 = self.grid[cxi][cy]
            n12 = self.grid[cxi][cyi]

            pdx = n11.d - n01.d
            pdy = n02.d - n01.d
            C20 = 3.0 * pdx - n11.gx - 2.0 * n01.gx
            C02 = 3.0 * pdy - n02.gy - 2.0 * n01.gy
            C30 = -2.0 * pdx + n11.gx + n01.gx
            C03 = -2.0 * pdy + n02.gy + n01.gy
            csum1 = n01.d + n01.gy + C02 + C03
            csum2 = n01.d + n01.gx + C20 + C30
            C21 = 3.0 * n12.d - 2.0 * n02.gx - n12.gx - 3.0 * csum1 - C20
            C31 = -2.0 * n12.d + n02.gx + n12.gx + 2.0 * csum1 - C30
            C12 = 3.0 * n12.d - 2.0 * n11.gy - n12.gy - 3.0 * csum2 - C02
            C13 = -2.0 * n12.d + n11.gy + n12.gy + 2.0 * csum2 - C03
            C11 = n02.gx - C13 - C12 - n01.gx
 
            u = p.x - cx
            u2 = u * u
            u3 = u * u2
            v = p.y - cy
            v2 = v * v
            v3 = v * v2
            density = n01.d + n01.gx * u + \
                      n01.gy * v + C20 * u2 + \
                      C02 * v2 + C30 * u3 + \
                      C03 * v3 + C21 * u2 * v + \
                      C31 * u3 * v + C12 * u * v2 + \
                      C13 * u * v3 + C11 * u * v
 
            pressure = density - 1.0
            if pressure > 2.0:
                pressure = 2.0
 
            fx = 0.0
            fy = 0.0
 
            if p.x < 1.5:
                fx += p.mat.m * (1.5 - p.x);
            elif p.x > self.gsizeX-5:
                fx += p.mat.m * (self.gsizeX - 5 - p.x)
 
            if p.y < 4.0:
                fy += p.mat.m * (4.0 - p.y)
            elif p.y > self.gsizeY - 5:
                fy += p.mat.m * (self.gsizeY - 5 - p.y)
 
            if drag:
                vx = abs(p.x - 0.25 * self.mx)
                vy = abs(p.y - 0.25 * self.my)
                if vx < 10.0 and vy < 10.0:
                    weight = p.mat.m * (1.0 - vx * 0.10) * (1.0 - vy * 0.10)
                    fx += weight * (mdx - p.u)
                    fy += weight * (mdy - p.v)
 
            for i in xrange(3):
                for j in xrange(3):
                    n = self.grid[(p.cx + i)][(p.cy + j)]
                    phi = p.px[i] * p.py[j]
                    n.ax += -((p.gx[i] * p.py[j]) * pressure) + fx * phi
                    n.ay += -((p.px[i] * p.gy[j]) * pressure) + fy * phi
 
        for n in self.active:
            if n.m > 0.0:
                n.ax /= n.m
                n.ay /= n.m
                n.ay += 0.03
 
        for p in self.particles:
            for i in xrange(3):
                for j in xrange(3):
                    n = self.grid[(p.cx + i)][(p.cy + j)]
                    phi = p.px[i] * p.py[j]
                    p.u += phi * n.ax
                    p.v += phi * n.ay

            mu = p.mat.m * p.u
            mv = p.mat.m * p.v
            for i in xrange(3):
                for j in xrange(3):
                    n = self.grid[(p.cx + i)][(p.cy + j)]
                    phi = p.px[i] * p.py[j]
                    n.u += phi * mu
                    n.v += phi * mv
 
        for n in self.active:
            if n.m > 0.0:
                n.u /= n.m
                n.v /= n.m
 
        for p in self.particles:
            gu = 0.0
            gv = 0.0
            for i in xrange(3):
                for j in xrange(3):
                    n = self.grid[(p.cx + i)][(p.cy + j)]
                    phi = p.px[i] * p.py[j]
                    gu += phi * n.u
                    gv += phi * n.v
            p.x += gu
            p.y += gv
            p.u += 1.0 * (gu - p.u)
            p.v += 1.0 * (gv - p.v)
            if p.x < 1.0:
                p.x = (1.0 + random.randint(0, 100) / 10000.0)
                p.u = 0.0
            elif p.x > self.gsizeX - 2:
                p.x = (self.gsizeX - 2 - random.randint(0, 100) / 10000.0)
                p.u = 0.0
            if p.y < 1.0:
                p.y = (1.0 + random.randint(0, 100) / 10000.0)
                p.v = 0.0
            elif p.y > self.gsizeY - 2:
                p.y = (self.gsizeY - 2 - random.randint(0, 100) / 10000.0)
                p.v = 0.0
                
class Node(object):
    '''
    Used to define individual cell in grid.
    '''
    __slots__ = ['m', 'd', 'gx', 'gy', 'u', 'v', 'ax', 'ay', 'active']
    def __init__(self):
        self.m = 0
        self.d = 0
        self.gx = 0
        self.gy = 0
        self.u = 0
        self.v = 0
        self.ax = 0
        self.ay = 0
        self.active = False
    

    def clear(self):
        self.m = self.d = self.gx = self.gy = self.u = \
        self.v = self.ax = self.ay = 0.0
        self.active = False

   
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
        self.area = 1.3333*3.14*size**3#1.3333333 * self.pi * math.pow(size/2, 3) #m^3
        self.drag_coeff = 0.000001   #Scaling number for drag force
        self.cor = .99     #Coefficient of restitution


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
    '''
    Create instance of body of water. Used to replicate a body of water and used to determine if game objects need to have the buoyancy
    force applied and the Drag Force changed accordingly.
    '''
    def __init__(self,world,pos= vect2d(0,0) ,volume = 0.0001,height=100,width=100,shape = "liquid", color = (0,0,255)):
        self.width = world.width
        self.height = world.height/4
        self.pos = vect2d(self.width/2,self.height*3.5)
        self.color = color
        self.shape = shape  #String with the name of the shape
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
        self.UL = vect2d(int(self.pos.x - self.width/2),int(self.pos.y - self.height/2))
        self.UR = vect2d(int(self.pos.x + self.width/2),int(self.pos.y - self.height/2))
        self.LL = vect2d(int(self.pos.x - self.width/2),int(self.pos.y + self.height/2))
        self.LR = vect2d(int(self.pos.x + self.width/2),int(self.pos.y + self.height/2))
        self.secretHeight = self.UR.y
        self.r_c = vect2d(0,0)
        self.xPrime = 0
        self.yPrime = 0
        self.xPrimeHat = 0
        self.yPrimeHat = 0
        self.n_hat = 0
        self.size = (self.UL-self.LR).mag()/2
        self.drag_coeff = 0.000001   #Scaling number for drag force
        self.cor = .99     #Coefficient of restitution

	
    
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
		
class Water_Particle(object):
    __slots__ = ['mat', 'x', 'y', 'u', 'v', 'dudx', 'dudy', 'dvdx',
                 'dvdy', 'cx', 'cy', 'px', 'py', 'gx', 'gy','shape','inLiquid','m']
    
    def __init__(self, mat, x, y, u, v,shape):
        self.m = 1
        self.shape = shape
        self.mat = mat
        self.x = x #Position x
        self.y = y #Position y
        self.u = u #Velocity x
        self.v = v #Velocity y
        self.dudx = 0 #Change in x, and so on with the other 3
        self.dudy = 0 
        self.dvdx = 0
        self.dvdy = 0
        self.cx = 0
        self.cy = 0
        self.px = [0, 0, 0]
        self.py = [0, 0, 0]
        self.gx = [0, 0, 0]
        self.gy = [0, 0, 0]

#Determines the property of water particles
class Material(object):
    __slots__ = ['m', 'rd', 'k', 'v', 'd', 'g']
    def __init__(self, m, rd, k, v, d, g):
        self.m = m;
        self.rd = rd;
        self.k = k;
        self.v = v;
        self.d = d;
        self.g = g;

class Visual(object):
    def __init__(self, screenSize, liquidTest):
        self.liquidTest = liquidTest
        self.screenSize = screenSize
        self.screen = pygame.display.set_mode(screenSize)
        self.clock = pygame.time.Clock()
        self.pop = False
        self.spawn = False
        
        self.particle_list = []
        self.liquid_list = []
        self.all_particles = []
        self.background_color = (255,255,255)  #World background set to white
        self.dt = 0.1     #time step size - Fixed for now
        self.width = 200  
        self.height = 400
        self.force_que = []    #Keeps track of forces, particles, and things like spring constants, etc
        self.g = vect2d(0,1)  #Constant gravitational field constant
        self.numerical = 'rk4'
        self.running = True   #Determines if while look should continue
        self.selected = None   #Particle selected with the mouse
        self.mouse_force = 1    #Coefficient for force exerted by mouse
        self.damping = 0.5    #Damping coefficient
        self.vel_max = 50     #Velocity at which the drag force kicks in
        self.collision = collision_engine(self,self.particle_list,self.all_particles) 
        self.solid = False
        
    
    def run(self):
        self.collision = collision_engine(self,self.particle_list,self.all_particles)   #This creates an instance of the collision engine
        while True:
    
            for event in pygame.event.get():
                
                if event.type == QUIT:
                    pygame.quit()
                    raise SystemExit
                #Move game_controls() into run, ran into less issues.
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
               
                    
                        
                    
                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_c: #Press and hold c to set pop to true
                        self.pop =True
                    if event.key == pygame.K_z: #Press and hold hold z to set spawn to true
                        self.spawn = True
                    elif event.key == pygame.K_RETURN: #Press enter to toggle liquid from being see-through to solid.
                        if self.solid == True:
                            self.solid = False
                        else:
                            self.solid = True
                if event.type == pygame.KEYUP:
                    self.pop =False
                    self.spawn = False;
                
            if self.pop == True: #If pop is true, delete oldest particles
                if(len(self.liquidTest.particles)>0):
                    self.liquidTest.particles.pop(0)
                if(len(self.all_particles)>2):
                    self.all_particles.pop(2)
    	       
            if self.spawn == True: #If spawn is true, spawn particles at pipe's location
                p = Water_Particle(self.liquidTest.water, (16+random.randrange(-4,4))/4, (16+random.randrange(-4,4))/4, 0.0, 0.0,"water_particle")
                self.liquidTest.particles.append(p)
                self.all_particles.append(p)  
            
                
            getattr(self,self.numerical)()
            self.collision.find_collision()  #Check to see if object hits the walls
            if pygame.mouse.get_pressed()[2]: #If right mouse is pressed, spawn particles at mouse's location
                x, y = pygame.mouse.get_pos()
                p = Water_Particle(self.liquidTest.water, (x+random.randrange(-4,4))/4, (y+random.randrange(-4,4))/4, 0.0, 0.0,"water_particle")
                self.liquidTest.particles.append(p)
                self.all_particles.append(p)

            self.screen.fill(Color('white'))
            
            #Draw static lines two static lines at different positions
            pygame.draw.line(self.screen, (0,0,0), (0,50), (125,100), 1)
            pygame.draw.line(self.screen, (0,0,0), (75,200), (200,150), 1)
            
            #Draw Liquid/particles
            for particle in self.liquid_list:
                if particle.shape =="liquid":
                    if self.solid == False:
                        pygame.draw.rect(self.screen, (0,0,255), [particle.UL.x, particle.secretHeight+5, particle.LR.x, particle.LR.y],2)
                    else:
                        pygame.draw.rect(self.screen, (0,0,255), [particle.UL.x, particle.secretHeight+5, particle.LR.x, particle.LR.y],0)
                    		  	  
            
            self.liquidTest.simulate()
            print len(self.liquidTest.particles) #Prints out number of particles in world
            
            #Set boundary for static lines
            self.line_boundary(vect2d(0,50/4),vect2d(125/4,100/4))
            self.line_boundary(vect2d(75/4,200/4),vect2d(200/4,150/4))
            
            self.liquidTest.draw(self.screen)
            
            #Creates red buoyant circle 
            for particle in self.particle_list:
                if particle.shape == "circle":
                    pygame.draw.circle(self.screen, (255,0,0), (int(particle.pos.x), int(particle.pos.y)), particle.size, 0)
                    
            #Draws two rectangles to create a pipe at the top left corner of the world
            pygame.draw.rect(self.screen,(200,200,235),(7,0,17,18),0)
            pygame.draw.rect(self.screen,(200,200,235),(5,12,21,6),0)
            
            pygame.display.update()
            
            #Displays FPS as the caption
            pygame.display.set_caption('fps: %d' % self.clock.get_fps()," len %d" % len(self.liquidTest.particles))

            self.clock.tick(60)
            
    #Determine if collision happens against lines
    def line_boundary(self,start,end):
        #Distance between start and end in grid coordinates
        r = start-end 
        m = r.y/r.x
        
        #Distance between start and end in world coordinates
        r2 = (start*4)-(end*4)
        m2 = r2.y/r2.x
        
        #Test if water particles are colliding with boundary lines.
        for i in self.liquidTest.particles:
            if i.x >= start.x and i.x < end.x:
                if start.x == 0:
                    y = m*(i.x)+start.y
                elif start.x == 75/4:
                    y = m*(i.x)+start.y+30/4
                if i.y <= y+1 and i.y > y-1:
                    i.y = y-1.001
                    i.v = -i.v/3
                    i.u += 1/(1000*m)
                    
        #Check if particles are colliding with boundary lines.
        for i in self.particle_list:
            if i.pos.x >=start.x*4 and i.pos.x<end.x*4:
                if start.x == 0:
                    y = m2*(i.pos.x)+start.y*4
                elif start.x*4 == 75:
                    y = m2*(i.pos.x)+start.y*4+30
                if i.pos.y <= y+1 and i.pos.y > y-1:
            
                    i.pos.y = y-1.001
                    #Si.vel.x *= -1
                    i.vel.y *= -1


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
            
        for particle in self.liquid_list:
            if particle.shape =="liquid":
                pygame.draw.aalines(self.screen,(0,0,255),True,[(particle.UL.x,particle.UL.y),(particle.UR.x,particle.UR.y),(particle.LR.x,particle.LR.y),(particle.LL.x,particle.LL.y)],False)
                
        pygame.display.flip()
        
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
        for particle in particles:
            force = vect2d(0,0)
            force = particle.m * self.g
            particle.f_net= particle.f_net + force
    
    #Caluclate buoyancy of object using object's size, density, area, and gravity of world
    def object_buoyancy(self,particles):
        for particle in particles:
            vol = 1.3333333 * 3.14 * math.pow(particle.size, 3)
            vol2 = particle.size* (particle.size*particle.d)*particle.size
            buoyancy = .0005 * particle.area*particle.d * self.g
            if particle.shape == "circle": 
                if particle.inLiquid == True:    
                    particle.f_net= particle.f_net - buoyancy
                elif particle.inLiquid == False:
                    particle.f_net = particle.f_net
                elif particle.shape == "water particle":
                    print "Test"
    
    #Calculate little g.
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
                i.f_net = i.f_net + dragF
            elif i.vel.mag() > self.vel_max:
                i.drag_coeff = .000001
                i.add_force(-i.drag_coeff*i.vel.mag2()*i.vel)

    def mouse_pull(self,particle):
        """Force applied by selecting particle with mouse"""
        (pick_x,pick_y) = pygame.mouse.get_pos()
        mouse_pos = vect2d(pick_x,pick_y)
        dx = self.selected.pos - mouse_pos
        F = -self.mouse_force*dx - self.selected.vel*self.damping
        self.selected.add_force(F)

class collision_engine():
    def __init__(self,world, particle_list=[],all_particles=[]):
        self.world = world  #World instance used by this collision engine
        self.particle_list = particle_list   #List of all particles in the world
        self.all_particles = all_particles
        self.has_collided = []    #List of particles that have collided with each other
        self.has_collidedR = []
        self.cor = 0.99   #Coefficient of restitution for hitting walls
		

    #Find what type of collisions are happening
    def find_collision(self):
        """EDIT: Runs through all particles in particle_list to look for overlapping particles

        Needs to check and see if particle has a rod connection and if it does, call check_rod_constraint
        """
        self.has_collided = []   #Reset collisions
        self.has_collidedR = []
        self.boundary()    #See if any particles collide with the wall
        
        #Check if water particle is colliding with liquid body. 
        for p in combinations(self.all_particles,2):
            #If so, raise water level.
            if p[0].shape =="liquid" and p[1].shape == "water_particle":	       
                if p[1].y/.2475>= p[0].secretHeight:
                    if(len(self.world.liquidTest.particles)>200):	               
                        self.world.liquidTest.particles.pop(0)
                        self.world.all_particles.pop(2)  
                    p[0].secretHeight = p[1].y/.255
            # If a particle collides with liquid, apply buoyancy force and determine how much of the particle is submerged
            # for buoyancy formula.
            elif p[0].shape =="liquid" and p[1].shape == "circle":
                if p[0].secretHeight<=p[1].pos.y:
                    p[1].inLiquid = True
                    if(abs(p[1].pos.y-p[0].secretHeight)>p[1].size):
                        p[1].d = 1
                    else:
                        p[1].d=abs(p[1].pos.y-p[0].secretHeight)/p[1].size
                else:
                    p[1].inLiquid = False


    def bounding_sphere(self,particle1,particle2):
        """Check for overlap using bounding spheres"""
        if particle1.shape == "circle" and particle2.shape == "circle":
            combinedSize = particle1.size + particle2.size
            if(abs(particle1.pos.x-particle2.pos.x) <= combinedSize):
                if(abs(particle1.pos.y-particle2.pos.y) <= combinedSize):
                    self.has_collided.append([particle1,particle2])
                self.resolve_collision()  
	      

    def resolve_collision(self):
        """Changes velocities of point particles that have collided"""
        for particles in self.has_collided:
            rHat = self.collision_normal(particles[0],particles[1])
            
            difference = particles[1].vel-particles[0].vel
            
            dv1 = (particles[1].m*(particles[0].cor)*(difference.dot(rHat)))/(particles[0].m+particles[1].m)
            dv2 = (-particles[0].m*(particles[1].cor)*(difference.dot(rHat)))/(particles[0].m+particles[1].m)
            
            particles[0].vel.x = particles[0].vel.x +dv1*rHat.x
            prticles[0].vel.y = particles[0].vel.y +dv1*rHat.y
            
            particles[1].vel.x = particles[1].vel.x +dv2*rHat.x
            particles[1].vel.y = particles[1].vel.y +dv2*rHat.y	
	

    def boundary(self):
        """
        World boundary for buoyant object
        """
        for particle in self.particle_list:
            if particle.shape == "circle" :
                if particle.pos.x > self.world.width - particle.size:
                    particle.vel.x *= -1*self.cor
                    particle.pos.x = self.world.width - particle.size
                elif particle.pos.x < 0+particle.size:
                    particle.vel.x *= -1*self.cor
                    particle.pos.x = particle.size
             		     
            if particle.pos.y > self.world.height - particle.size:
                particle.pos.y = self.world.height-particle.size
                particle.vel.y *= -1*self.cor
            elif particle.pos.y < 0+particle.size:
                particle.pos.y = particle.size
                particle.vel.y *= -1*self.cor



def main():

    os.environ['SDL_VIDEO_CENTERED'] = '1'

    liquidTest = LiquidTest(54, 81, 4, 15)
    visual = Visual((200, 400), liquidTest) #Sets window of world
    visual.set_numerical('rk4')  

    water = liquid(visual,vect2d(100,0))
    r1 = vect2d(35,200)
    v1 = vect2d(0,10)
    part1 = particle(r1,v1)
    visual.add_particle(water)
    visual.add_particle(part1)
    visual.new_force('const_grav',[part1])
    visual.new_force("object_buoyancy",[part1])
    visual.new_force('drag',[part1])
    
    visual.run()


if __name__ == '__main__':
    main()