import numpy as np
from Particle import Particle
from Boundary import Boundary
from Timestep import TimeFuncs
from Gfx import Display
import time

disp = Display()
p = Particle(disp = disp, radius = 25, speed = [50,0])
b1 = Boundary(disp = disp, angle = 0, offset=[0,-60])
b2 = Boundary(disp = disp, angle = 0, offset=[0,60], direction=False)
b3 = Boundary(disp = disp, angle = np.pi/2, offset=[60,0])
b4 = Boundary(disp = disp, angle = np.pi/2, offset=[-60,0], direction=False)
stepper = TimeFuncs()

while(True):
    tStart = time.time()
    p.update(stepper)
    time.sleep(time.time() - tStart + 1/60)