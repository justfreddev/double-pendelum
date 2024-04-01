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
WIDTH, HEIGHT = 1920, 1080

# Sets the FPS for the pygame window
FPS = 60.0

# Sets all the minimum and maximum values for the pendulum parameters
MIN_MASS = 10
MAX_MASS = 100
MIN_G = 1
MAX_G = 20
MIN_L = 100
MAX_L = 300

# Sets the pivot point for the pendulum
PIVOT = (WIDTH // 2, (HEIGHT // 2) - 50)

# Sets the time step for the simulation
dt = 0.5


def accel1(m1, m2, a1, a2, v1, v2, l1, l2, g):
    """Calculates the acceleration of the first pendulum mass"""
    num = -g * (2 * m1 + m2) * sin(a1) - m2 * g * sin(a1 - 2 * a2) - 2 * sin(a1 - a2) * m2 * (v2 * v2 * l2 + v1 * v1 * l1 * cos(a1 - a2))
    den = l1 * (2 * m1 + m2 - m2 * cos(2 * a1 - 2 * a2))
    return num / den

def accel2(m1, m2, a1, a2, v1, v2, l1, l2, g):
    """Calculates the acceleration of the second pendulum mass"""
    num = 2 * sin(a1 - a2) * (v1 * v1 * l1 * (m1 + m2) + g * (m1 + m2) * cos(a1) + v2 * v2 * l2 * m2 * cos(a1 - a2))
    den = l2 * (2 * m1 + m2 - m2 * cos(2 * a1 - 2 * a2))
    return num / den


class Slider:
    """A class that creates a slider for the pygame window"""
    def __init__(self, label, min_val, max_val, pos, width):
        self.label = label
        self.min = min_val
        self.max = max_val
        self.value = min_val + (max_val - min_val) / 2
        self.pos = pos
        self.slider_width = width
        self.dragging = False

    def draw(self, screen):
        """Some really confusing pygame stuff"""
        font = pygame.font.Font(None, 24)
        text = font.render(self.label, True, WHITE)
        screen.blit(text, (self.pos[0] - 100, self.pos[1] - 8))
        
        pygame.draw.line(screen, WHITE, (self.pos[0] - 5, self.pos[1]), (self.pos[0] + self.slider_width - 10, self.pos[1]), 5)
        knob_x = self.pos[0] + (self.value - self.min) / (self.max - self.min) * (self.slider_width - 10)
        pygame.draw.circle(screen, WHITE, (knob_x, self.pos[1]), 10)
        
        text = font.render(str(int(self.value)), True, WHITE)
        screen.blit(text, (knob_x - 10, self.pos[1] - 30))

    def handle_event(self, event):
        """Handles the event of the slider being clicked and dragged"""
        if event.type == pygame.MOUSEBUTTONDOWN:
            if self.pos[0] <= event.pos[0] <= self.pos[0] + self.slider_width and self.pos[1] <= event.pos[1] <= self.pos[1] + 20:
                self.dragging = True
        elif event.type == pygame.MOUSEBUTTONUP:
            self.dragging = False
        elif event.type == pygame.MOUSEMOTION and self.dragging:
            mouse_x = event.pos[0]
            self.value = (mouse_x - self.pos[0]) / (self.slider_width - 10) * (self.max - self.min) + self.min
            self.value = max(min(self.value, self.max), self.min)
        
        return self.value


# Creates the sliders for the pendulum parameters
mass1_slider = Slider("Mass 1:", MIN_MASS, MAX_MASS, (WIDTH - 250, 50), 200)
mass2_slider = Slider("Mass 2:", MIN_MASS, MAX_MASS, (WIDTH - 250, 100), 200)

g_slider = Slider("Gravity:", MIN_G, MAX_G, (WIDTH - 250, 150), 200)

l1_slider = Slider("Length 1:", MIN_L, MAX_L, (WIDTH - 250, 200), 200)
l2_slider = Slider("Length 2:", MIN_L, MAX_L, (WIDTH - 250, 250), 200)


def draw(screen, coords1, coords2, mass1, mass2):
    """Draws the pendulum on the screen"""
    screen.fill(BLACK)

    # Draws the lines and circles that make up the pendulum
    pygame.draw.line(screen, WHITE, PIVOT, coords1, 5)
    pygame.draw.line(screen, WHITE, coords1, coords2, 5)
    pygame.draw.circle(screen, RED, PIVOT, 5)
    pygame.draw.circle(screen, GREEN, coords1, mass1 / 3)
    pygame.draw.circle(screen, BLUE, coords2, mass2 / 3)
    
    # Draws the sliders on the screen
    mass1_slider.draw(screen)
    mass2_slider.draw(screen)
    g_slider.draw(screen)
    l1_slider.draw(screen)
    l2_slider.draw(screen)

    pygame.display.flip()


def run():
    pygame.init()
    screen = pygame.display.set_mode((WIDTH, HEIGHT))
    pygame.display.set_caption("Double Pendulum Simulation")

    # Initial conditions
    m1, m2 = 60, 60
    l1, l2 = 200, 200
    a1, a2 = pi / 2, pi / 2
    v1, v2 = 0, 0
    g = (MAX_G + MIN_G) / 2
    

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
        k1_v1 = dt * accel1(m1, m2, a1, a2, v1, v2, l1, l2, g)
        k1_v2 = dt * accel2(m1, m2, a1, a2, v1, v2, l1, l2, g)
        k1_a1 = dt * v1
        k1_a2 = dt * v2
        
        # k2_v = dt * accel(t + 0.5 * dt, v + 0.5 * k1)
        # k2_a = dt * (v + 0.5 * k1)
        k2_v1 = dt * accel1(m1, m2, a1 + 0.5 * k1_a1, a2 + 0.5 * k1_a2, v1 + 0.5 * k1_v1, v2 + 0.5 * k1_v2, l1, l2, g)
        k2_v2 = dt * accel2(m1, m2, a1 + 0.5 * k1_a1, a2 + 0.5 * k1_a2, v1 + 0.5 * k1_v1, v2 + 0.5 * k1_v2, l1, l2, g)
        k2_a1 = dt * (v1 + 0.5 * k1_v1)
        k2_a2 = dt * (v2 + 0.5 * k1_v2)
        
        # k3 = dt * accel(t + 0.5 * dt, v + 0.5 * k2)
        # k3_a = dt * (v + 0.5 * k2)
        k3_v1 = dt * accel1(m1, m2, a1 + 0.5 * k2_a1, a2 + 0.5 * k2_a2, v1 + 0.5 * k2_v1, v2 + 0.5 * k2_v2, l1, l2, g)
        k3_v2 = dt * accel2(m1, m2, a1 + 0.5 * k2_a1, a2 + 0.5 * k2_a2, v1 + 0.5 * k2_v1, v2 + 0.5 * k2_v2, l1, l2, g)
        k3_a1 = dt * (v1 + 0.5 * k2_v1)
        k3_a2 = dt * (v2 + 0.5 * k2_v2)
        
        # k4 = dt * accel(t + dt, v + k3)
        # k4_a = dt * (v + k3)
        k4_v1 = dt * accel1(m1, m2, a1 + k3_a1, a2 + k3_a2, v1 + k3_v1, v2 + k3_v2, l1, l2, g)
        k4_v2 = dt * accel2(m1, m2, a1 + k3_a1, a2 + k3_a2, v1 + k3_v1, v2 + k3_v2, l1, l2, g)
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
        
        # Event handling to prevent the window from crashing/freezing (pygame things)
        for event in pygame.event.get():
            if event.type == QUIT:
                pygame.quit()
                sys.exit()
            # Updates the values of the pendulum parameters
            m1 = mass1_slider.handle_event(event)
            m2 = mass2_slider.handle_event(event)
            g = -g_slider.handle_event(event)
            l1 = l1_slider.handle_event(event)
            l2 = l2_slider.handle_event(event)
        draw(screen, c1, c2, m1, m2)

        pygame.time.Clock().tick(60.0)


if __name__ == "__main__":
    run()
