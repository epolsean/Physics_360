import sys
import random
import math

import pygame
from pygame.locals import *
from pygame.color import *

import pymunk as pm

def add_ball(space, position):
    radius = random.randint(5, 10)
    mass = radius
    inertia = pm.moment_for_circle(mass, 0, radius, (0, 0))

    body = pm.Body(mass, inertia)

    x, y = position

    x = random.randint(x-2, x+2)

    y = random.randint(y-2, y+2)

    body.position = x, y
    shape = pm.Circle(body, radius, (0, 0))
    space.add(body, shape)
    shape.color = random.randrange(50, 200), 0, 0
    return shape

def draw_ball(screen, ball):
    p = int(ball.body.position.x), 600-int(ball.body.position.y)

#    rad = ball.radius
#
#    rect = pygame.rect.Rect(p[0]-rad, p[1]-rad, rad*2, rad*2)
#
#    pygame.draw.rect(screen, ball.color, rect)

    pygame.draw.circle(screen, ball.color, p, int(ball.radius), 0)

def add_l(space):
    body = pm.Body(pm.inf, pm.inf)
    body.position = (0, 0)
    l1 = pm.Segment(body, (400, 50), (590, 150.0), 3.0)
    l1b = pm.Segment(body, (590, 150), (590, 590.0), 3.0)

    l2 = pm.Segment(body, (10, 150), (200, 50), 3.0)
    l2b = pm.Segment(body, (10, 150), (10, 590.0), 3.0)

    lc = pm.Segment(body, (10, 590), (590, 590), 3.0)

    lt1 = pm.Segment(body, (250, 0), (300, 65), 3.0)
    lt2 = pm.Segment(body, (300, 65), (350, 0), 3.0)

    space.add(l1, l2, l1b, l2b, lc, lt1, lt2)
    return l1, l2, l1b, l2b, lc, lt1, lt2

def draw_lines(screen, lines):
    for line in lines:
        body = line.body
        pv1 = body.position + line.a.rotated(math.degrees(body.angle))
        pv2 = body.position + line.b.rotated(math.degrees(body.angle))
        p1 = to_pygame(pv1)
        p2 = to_pygame(pv2)
        pygame.draw.lines(screen, THECOLORS["lightgray"], False, [p1, p2])

def to_pygame(p):
    return int(p.x), int(-p.y+600)

def to_pygame_tup(p):
    return int(p[0]), int(-p[1]+600)

def main():
    pygame.init()

    screen = pygame.display.set_mode((600, 600))
    pygame.display.set_caption("PyMunk")
    clock = pygame.time.Clock()

    running = True

    pm.init_pymunk()
    space = pm.Space()
    gravity = 0.0
    space.gravity = (0.0, gravity)

    balls = []

    lines = add_l(space)
    
    font = pygame.font.Font(None, 50)

    while running:
        for event in pygame.event.get():
            if event.type == QUIT:
                running = False
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                running = False
            if event.type == MOUSEBUTTONDOWN:
                if event.button == 4:
                    gravity += 50.0
                    space.gravity = (0.0, gravity)
                elif event.button == 5:
                    gravity -= 50.0
                    space.gravity = (0.0, gravity)
                
        
        if pygame.mouse.get_pressed()[0]:
            x, y = pygame.mouse.get_pos()
            ball_shape = add_ball(space, to_pygame_tup((x, y)))
            balls.append(ball_shape)
    
            
        if pygame.mouse.get_pressed()[2]:
            x, y = pygame.mouse.get_pos()
            
            for ball in balls:
                ball.body.apply_impulse(ball.body.position,
                                        to_pygame_tup((x, y)))
            
            
            
        for ball in balls:

            #ball.body.apply_impulse((50, 50), (50,50))
            if ball.body.position.y < 0-ball.radius:
                balls.remove(ball)
            elif ball.body.position.y > 600+ball.radius:
                balls.remove(ball)

        gravnum = font.render("""Gravity = %d"""%gravity, True,
                              (159, 182, 205), (255, 255, 255))

        gravrect = gravnum.get_rect()
        gravrect.centerx = screen.get_rect().centerx
        gravrect.centery = screen.get_rect().centery-20

        ballnum = font.render("""Balls = %d"""%len(balls), True,
                              (159, 182, 205), (255, 255, 255))

        ballrect = ballnum.get_rect()
        ballrect.centerx = screen.get_rect().centerx
        ballrect.centery = screen.get_rect().centery+20

        #PAINTING

        screen.fill(THECOLORS["white"])

        screen.blit(gravnum, gravrect)
        screen.blit(ballnum, ballrect)

        for ball in balls:
            draw_ball(screen, ball)

        draw_lines(screen, lines)

        space.step(1/50.0)

        pygame.display.flip()
        clock.tick(50)

if __name__ == '__main__':
    sys.exit(main())