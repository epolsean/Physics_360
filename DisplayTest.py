# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 20:19:16 2015

@author: student
"""

from __future__ import division
import pygame
import random
from vector import *   #Rename your solution to homework 2 as vector.py
import math


class canvas:
    def __init__(self):
        self.running = True
        self.width = 600
        self.height = 600
        self.backgroundColor = (0, 0, 0)
        self.setupScreen()
        self.x = 0
        
    def setupScreen(self):
        pygame.init()
        pygame.display.set_caption('Fluid Dynamics')
        self.screen = pygame.display.set_mode((self.width, self.height))
        self.clock = pygame.time.Clock()
    
    def display(self):
        self.screen.fill(self.backgroundColor)  #Wipe clean the surface in memory
        rand1 = random.randrange(1,4)
        rand2 = random.randrange(-1,1)
        size = 15
        length = 4
        self.x += rand2/10
        for i in range(-20,20):
            for j in range(-20, 20):
                U = math.sin(j/8 + self.x)
                V = math.cos(j/2 + self.x)
                #print "("+str(U)+","+str(V)+")"
                line = vect2d(V,U)
                norm = line.norm()
                pygame.draw.line(self.screen, (255,255,255), ((i*size+300-norm.x*length),(j*size+300-norm.y*length)), ((i*size+300+norm.x*length),(j*size+300+norm.y*length)), 1)
                
        pygame.display.flip()   #Displays the world
    
    def update(self):
        self.display()
        self.controls()
        
    def controls(self):
        for event in pygame.event.get(): 
            if event.type == pygame.QUIT:
                self.running = False

if __name__ == '__main__' :
    world = canvas()
    
    while world.running:
        time_passed = world.clock.tick(120)
        world.update()
        