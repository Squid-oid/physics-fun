import sys
import sdl2
import sdl2.ext
import W_object as W
import numpy as np
import Timestep as ts
import time as time
from W_object import Collision as Col

#################################################################ch
# Mini Classes ##################################################

class SoftwareRenderer(sdl2.ext.SoftwareSpriteRenderSystem):
    def __init__(self, window):
        super(SoftwareRenderer, self).__init__(window)

    def render(self, components):
        sdl2.ext.fill(self.surface, sdl2.ext.Color(255, 255, 255))
        super(SoftwareRenderer, self).render(components)

class r_object(sdl2.ext.Entity):
    def __init__(self, world, object:W.W_object):
        self.sprite = object.sprite

    def cloneTo(self, object:W.W_object):
        self.sprite = object.sprite

class InitFuncs():
    def create_objects(factory):
     
        objects = []

        bll = W.Ball(fact = factory, coord=[101,151], vel=[50,0], radius=50, mass = 100)        #Velocity in pixels per second
        bll2 = W.Ball(fact = factory , coord=[700,51], vel=[-250,0], radius=10, mass = 0.5)        #Velocity in pixels per second
        bll3 = W.Ball(fact = factory, coord=[300,50], vel=[-120,0], radius=30, mass= 4)        #Velocity in pixels per second
        objects.append(bll)
        objects.append(bll2)
        objects.append(bll3)

        bndt = W.Barrier(coord = [0,0], fact = factory, angle=np.pi)
        bndl = W.Barrier(coord=[0,0], fact = factory, angle=np.pi*1/2)
        bndb = W.Barrier(coord=[600,400], fact = factory, angle= 0)
        bndr = W.Barrier(coord=[600,400], fact = factory, angle=np.pi*3/2)
        objects.append(bndt)
        objects.append(bndl)
        objects.append(bndb)
        objects.append(bndr)

        tri1_coords = 350 + np.asmatrix([[10,0],[0,20],[20,20]])
        tri1 = W.Tri(fact = factory, coords=tri1_coords, vel = [10,10])
        tri2_coords = 250 + np.asmatrix([[10,0],[0,30],[40,30]])
        tri2 = W.Tri(fact = factory, coords=tri2_coords, vel = [10,10])
        tri3_coords = np.asmatrix([[0,40],[0,0],[40,40]])
        tri3 = W.Tri(fact = factory, coords=tri3_coords, vel = [20,5])
        objects.append(tri1)
        objects.append(tri2)
        objects.append(tri3)

        return objects

#######################################################
# Code begins #########################################

sdl2.ext.init()

world = sdl2.ext.World()
window = sdl2.ext.Window("Box", size=(600, 400))
window.show()
processor = sdl2.ext.TestEventProcessor()
game_running = True

spriterenderer = SoftwareRenderer(window)
world.add_system(spriterenderer)

factory = sdl2.ext.SpriteFactory(sdl2.ext.SOFTWARE)
RESOURCES = sdl2.ext.Resources(__file__, "basesprites")
bckg = factory.from_color(color=[255,255,255], size=[5000,5000])

objects = InitFuncs.create_objects(factory)
r_objects = []
for object in objects: 
    r_objects.append(r_object(world, object))
stepper = ts.Time_Funcs(t = time.time())

refresh_args = {
    W.Ball: {'stepper' : stepper, 'objs' : objects},
    W.Barrier: {}, # The barrier doesn't need args
    W.W_object: {'stepper' : stepper, 'objs' : objects},
    W.Tri: {'stepper' : stepper, 'objs' : objects}
}

while(game_running):
    stepper.find_set_stepsize(time.time())
    events = sdl2.ext.get_events()
    for event in events:
        if event.type == sdl2.SDL_QUIT:
            game_running = False
            break
        elif event.type == sdl2.SDL_KEYDOWN:
            prevE = Col.elasticity
            if event.key.keysym.sym == sdl2.SDLK_UP:   # UP arrow increases elasticity
                Col.setElasticity(min(1.0, Col.elasticity + 0.05))
            elif event.key.keysym.sym == sdl2.SDLK_DOWN:  # DOWN arrow decreases elasticity
                Col.setElasticity(max(0.1, Col.elasticity - 0.05))           
            if event.key.keysym.sym == sdl2.SDLK_g:
                W.W_object.setGravity(not W.W_object.gravity)

    for [object, r_obj] in zip(objects, r_objects):
        object.refresh(refresh_args[type(object)])
        r_obj = r_obj.cloneTo(object)

    status_line = f"Elasticity: {Col.elasticity:.2f} | Gravity: {'ON' if W.W_object.gravity else 'OFF'}"
    print(f"\r\033[K{status_line}", end="")
    sys.stdout.flush()

    world.process()
    
