import sys
import sdl2
import sdl2.ext
import W_object as W
import numpy as np
import Timestep as ts
import time as time

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
        barriers = []
        balls = []

        bll = W.Ball(fact = factory, coord=[51,51], vel=[60,0], radius=50, mass = 100)        #Velocity in pixels per second
        bll2 = W.Ball(fact = factory,coord=[50,32], vel=[-60,40])        #Velocity in pixels per second
        bll3 = W.Ball(fact = factory , coord=[700,51], vel=[-600,0])        #Velocity in pixels per second
        bll4 = W.Ball(fact = factory, coord=[300,50], vel=[-300,0])        #Velocity in pixels per second
        bll5 = W.Ball(fact = factory, coord=[58,38], vel=[-60,400])        #Velocity in pixels per second

        bndt = W.Barrier(coord = [0,0], fact = factory, angle=np.pi)
        bndl = W.Barrier(coord=[0,0], fact = factory, angle=np.pi*1/2)
        bndb = W.Barrier(coord=[800,600], fact = factory, angle= 0)
        bndr = W.Barrier(coord=[800,600], fact = factory, angle=np.pi*3/2)

        tri1_coords = np.asmatrix([[30,30],[40,30],[40,40]])
        tri1 = W.Tri(fact = factory, coords=tri1_coords, vel = [10,10])

        objects.append(tri1)

        objects.append(bndt)
        objects.append(bndl)
        objects.append(bndb)
        objects.append(bndr)

        objects.append(bll)
        objects.append(bll2)
        objects.append(bll3)
        objects.append(bll4)
        objects.append(bll5)

        balls.append(bll)
        balls.append(bll2)
        balls.append(bll3)
        balls.append(bll4)
        balls.append(bll5)

        barriers.append(bndt)
        barriers.append(bndl)
        barriers.append(bndb)
        barriers.append(bndr)

        return [objects, balls, barriers]

#######################################################
# Code begins #########################################

sdl2.ext.init()

world = sdl2.ext.World()
window = sdl2.ext.Window("Box", size=(800, 600))
window.show()
processor = sdl2.ext.TestEventProcessor()
game_running = True

spriterenderer = SoftwareRenderer(window)
world.add_system(spriterenderer)

factory = sdl2.ext.SpriteFactory(sdl2.ext.SOFTWARE)
RESOURCES = sdl2.ext.Resources(__file__, "basesprites")
bckg = factory.from_color(color=[255,255,255], size=[5000,5000])

[objects, balls, barriers] = InitFuncs.create_objects(factory)
r_objects = []
for object in objects: 
    r_objects.append(r_object(world, object))
stepper = ts.Time_Funcs(t = time.time())

refresh_args = {
    W.Ball: {'barriers': barriers, 'stepper' : stepper, 'balls' : balls},
    W.Barrier: {},
    W.W_object: {},
    W.Tri: {}
}


while(game_running):
    stepper.find_set_stepsize(time.time())
    events = sdl2.ext.get_events()
    for event in events:
        if event.type == sdl2.SDL_QUIT:
            game_running = False
            break
    for [object, r_obj] in zip(objects, r_objects):
        object.refresh(refresh_args[type(object)])
        r_obj = r_obj.cloneTo(object)
    world.process()
    
