# This is a simulation of a double pendulum using the Runge-Kutta 4th order method
# The Runge-Kutta 4th order method is used to solve ordinary differential equations
# that represent the motion of the double pendulum. The simulation is displayed using
# the pygame library.

from math import sin, cos
import pygame
from pygame.locals import *
from scipy.constants import pi
import sys

# All of the colour RGB values
BLUE = (0, 0, 255)
GREEN = (0, 255, 0)
RED = (255, 0, 0)
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)

# Sets the width and height of the window
WIDTH, HEIGHT = 1280, 960

# Sets the FPS for the pygame window
FPS = 60.0

# Sets the pivot point for the pendulum
PIVOT = (WIDTH // 2, (HEIGHT // 2) - 50)

# Sets the accelerating force of gravity (is minus because pygame is flipped weirdly)
g = -9.8

# Sets the time step for the simulation
dt = 0.5


def accel1(m1, m2, a1, a2, v1, v2, l1, l2):
    """Calculates the acceleration of the first pendulum mass"""
    num = -g * (2 * m1 + m2) * sin(a1) - m2 * g * sin(a1 - 2 * a2) - 2 * sin(a1 - a2) * m2 * (v2 * v2 * l2 + v1 * v1 * l1 * cos(a1 - a2))
    den = l1 * (2 * m1 + m2 - m2 * cos(2 * a1 - 2 * a2))
    return num / den

def accel2(m1, m2, a1, a2, v1, v2, l1, l2):
    """Calculates the acceleration of the second pendulum mass"""
    num = 2 * sin(a1 - a2) * (v1 * v1 * l1 * (m1 + m2) + g * (m1 + m2) * cos(a1) + v2 * v2 * l2 * m2 * cos(a1 - a2))
    den = l2 * (2 * m1 + m2 - m2 * cos(2 * a1 - 2 * a2))
    return num / den


def draw(screen, coords1, coords2, angle1, angle2, mass1, mass2):
    """Draws the pendulum on the screen"""
    screen.fill(BLACK)

    # Converts the angles to degrees
    deg1 = round((angle1 * 180 / pi) % 360)
    deg2 = round((angle2 * 180 / pi) % 360)

    # Draws the lines and circles that make up the pendulum
    pygame.draw.line(screen, WHITE, PIVOT, coords1, 5)
    pygame.draw.line(screen, WHITE, coords1, coords2, 5)
    pygame.draw.circle(screen, RED, PIVOT, 5)
    pygame.draw.circle(screen, GREEN, coords1, mass1 / 3)
    pygame.draw.circle(screen, BLUE, coords2, mass2 / 3)
    
    # Renders the angles on the screen
    font = pygame.font.Font(None, 36)
    a1 = font.render(str(deg1) + "\u00B0", True, WHITE)
    a2 = font.render(str(deg2) + "\u00B0", True, WHITE)
    screen.blit(a1, (10, 10))
    screen.blit(a2, (10, 40))

    pygame.display.flip()


def run():
    pygame.init()
    screen = pygame.display.set_mode((WIDTH, HEIGHT))

    # Initial conditions
    m1, m2 = 60, 60
    l1, l2 = 200, 200
    a1, a2 = pi / 2, pi / 2
    v1, v2 = 0, 0

    while True:
        # Runge-Kutta 4th order method
        
        # dt * accel = change in velocity over time (dv/dt)
        # dt * v = change in angle over time (dθ/dt)
        
        # The method starts of with the initial known values of the system, which in this case are the angles and velocities,
        # and calculates the first slope (acceleration/velocity) using the initial values.
        # It then uses the next values of the system to calculate the next slope, and so on.
        # The average slope is then calculated and used to update the values of the system.
        
        # k1_v = dt * accel(t, v)
        # k1_a = dt * v
        k1_v1 = dt * accel1(m1, m2, a1, a2, v1, v2, l1, l2)
        k1_v2 = dt * accel2(m1, m2, a1, a2, v1, v2, l1, l2)
        k1_a1 = dt * v1
        k1_a2 = dt * v2
        
        # k2_v = dt * accel(t + 0.5 * dt, v + 0.5 * k1)
        # k2_a = dt * (v + 0.5 * k1)
        k2_v1 = dt * accel1(m1, m2, a1 + 0.5 * k1_a1, a2 + 0.5 * k1_a2, v1 + 0.5 * k1_v1, v2 + 0.5 * k1_v2, l1, l2)
        k2_v2 = dt * accel2(m1, m2, a1 + 0.5 * k1_a1, a2 + 0.5 * k1_a2, v1 + 0.5 * k1_v1, v2 + 0.5 * k1_v2, l1, l2)
        k2_a1 = dt * (v1 + 0.5 * k1_v1)
        k2_a2 = dt * (v2 + 0.5 * k1_v2)
        
        # k3 = dt * accel(t + 0.5 * dt, v + 0.5 * k2)
        # k3_a = dt * (v + 0.5 * k2)
        k3_v1 = dt * accel1(m1, m2, a1 + 0.5 * k2_a1, a2 + 0.5 * k2_a2, v1 + 0.5 * k2_v1, v2 + 0.5 * k2_v2, l1, l2)
        k3_v2 = dt * accel2(m1, m2, a1 + 0.5 * k2_a1, a2 + 0.5 * k2_a2, v1 + 0.5 * k2_v1, v2 + 0.5 * k2_v2, l1, l2)
        k3_a1 = dt * (v1 + 0.5 * k2_v1)
        k3_a2 = dt * (v2 + 0.5 * k2_v2)
        
        # k4 = dt * accel(t + dt, v + k3)
        # k4_a = dt * (v + k3)
        k4_v1 = dt * accel1(m1, m2, a1 + k3_a1, a2 + k3_a2, v1 + k3_v1, v2 + k3_v2, l1, l2)
        k4_v2 = dt * accel2(m1, m2, a1 + k3_a1, a2 + k3_a2, v1 + k3_v1, v2 + k3_v2, l1, l2)
        k4_a1 = dt * (v1 + k3_v1)
        k4_a2 = dt * (v2 + k3_v2)
        
        # Average slope (acceleration/velocity) = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        a1 += (k1_a1 + 2 * k2_a1 + 2 * k3_a1 + k4_a1) / 6
        a2 += (k1_a2 + 2 * k2_a2 + 2 * k3_a2 + k4_a2) / 6
        v1 += (k1_v1 + 2 * k2_v1 + 2 * k3_v1 + k4_v1) / 6
        v2 += (k1_v2 + 2 * k2_v2 + 2 * k3_v2 + k4_v2) / 6
        
        # Calculate the coordinates of the pendulum masses
        c1 = (PIVOT[0] + (l1 * sin(a1)), PIVOT[1] - (l1 * cos(a1)))
        c2 = (c1[0] + (l2 * sin(a2)), c1[1] - (l2 * cos(a2)))

        # Event handling to prevent the window from crashing/freezing
        for event in pygame.event.get():
            if event.type == QUIT:
                pygame.quit()
                sys.exit()
        draw(screen, c1, c2, a1, a2, m1, m2)

        pygame.time.Clock().tick(60.0)


if __name__ == "__main__":
    run()
