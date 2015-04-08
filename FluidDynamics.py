# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 20:19:16 2015

@author: student
"""

import pygame
import random


class fluid:
    def __init__(self, pos, vel, density, vicosity, energy):
        self.pos = pos
        self.vel = vel
        self.density = density
        self.vicosity = vicosity
        self.energy = energy

class canvas:
    def __init__(self):
        self.running = True
        self.width = 600
        self.height = 600
        self.backgroundColor = (0, 0, 0)
        self.setupScreen()
        
    def setupScreen(self):
        pygame.init()
        pygame.display.set_caption('Fluid Dynamics')
        self.screen = pygame.display.set_mode((self.width, self.height))
        self.screen.fill(self.backgroundColor)
        
    def update(self):
        self.controls()
        
    def controls(self):
        for event in pygame.event.get(): 
            if event.type == pygame.QUIT:
                self.running = False

if __name__ == '__main__' :
    world = canvas()
    
    while world.running:
        world.update()
        pygame.display.flip()
        