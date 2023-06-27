import sys
import sdl2
import sdl2.ext
import W_object as W
import numpy as np
import Timestep as ts
import time as time

game_running = True
RESOURCES = sdl2.ext.Resources(__file__, "basesprites")
sdl2.ext.init()

window = sdl2.ext.Window("Test Window", size=(640, 480))
window.show()

factory = sdl2.ext.SpriteFactory(sdl2.ext.SOFTWARE)
sprite_renderer = factory.create_sprite_render_system(window)
bckg = factory.from_color(color=[255,255,255], size=[5000,5000])


processor = sdl2.ext.TestEventProcessor()

objects = []
barriers = []

bll = W.Ball(coord=[10,30], vel=[60,40])        #Velocitty in pixels per second
bndt = W.Barrier(coord = [0,0], fact = factory, angle=np.pi)
bndl = W.Barrier(coord=[0,0], fact = factory, angle=np.pi*1/2)
bndb = W.Barrier(coord=[640,480], fact = factory, angle= 0)
bndr = W.Barrier(coord=[640,480], fact = factory, angle=np.pi*3/2)
objects.append(bndt)
objects.append(bndl)
objects.append(bndb)
objects.append(bndr)

objects.append(bll)

barriers.append(bndt)
barriers.append(bndl)
barriers.append(bndb)
barriers.append(bndr)

stepper = ts.Time_Funcs(t = time.time())

refresh_args = {
    W.Ball: {'barriers': barriers, 'stepper' : stepper},
    W.Barrier: {},
    W.W_object: {},
}

while(game_running):
    stepper.find_set_stepsize(time.time())
    for obj in objects:
        args = refresh_args[type(obj)]
        obj.refresh(args)
    sprite_renderer.render(bckg)
    for obj in objects:
        obj.draw(sprite_renderer = sprite_renderer)
        sdl2.SDL_Delay(90)
    
