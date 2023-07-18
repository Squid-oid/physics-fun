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
balls = []

bll = W.Ball(coord=[51,51], vel=[60,0], radius=50)        #Velocitty in pixels per second
bll2 = W.Ball(coord=[50,32], vel=[-60,40])        #Velocitty in pixels per second
bll3 = W.Ball(coord=[700,51], vel=[-600,0])        #Velocitty in pixels per second
bll4 = W.Ball(coord=[300,50], vel=[-300,0])        #Velocitty in pixels per second
bll5 = W.Ball(coord=[58,38], vel=[-60,400])        #Velocitty in pixels per second

bndt = W.Barrier(coord = [0,0], fact = factory, angle=np.pi)
bndl = W.Barrier(coord=[0,0], fact = factory, angle=np.pi*1/2)
bndb = W.Barrier(coord=[640,480], fact = factory, angle= 0)
bndr = W.Barrier(coord=[640,480], fact = factory, angle=np.pi*3/2)

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


stepper = ts.Time_Funcs(t = time.time())

refresh_args = {
    W.Ball: {'barriers': barriers, 'stepper' : stepper, 'balls' : balls},
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
    sdl2.SDL_Delay(1)