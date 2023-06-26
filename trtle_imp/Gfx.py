import logging
import numpy as np
import turtle
import time
import sched

class Display:
    def __init__(self):
        self.window = turtle.Screen()
        self.trt = turtle.Turtle(visible=False)
        self.drawables = []
        self.trt.speed(speed = 0)
        self.updates = []
        self.sch = sched.scheduler(time.time, time.sleep)
        self.flag_update("initialized")
        self.bnds = []

    def add_drawable(self, obj):
        self.drawables.append(obj)
        self.flag_update("add_drawable")

    def add_boundary(self, bnd):
        self.bnds.append(bnd)
        self.add_drawable(bnd)

    def redraw(self):
        stime = time.time()
        logging.warning('Began Clearing Updates \n' + str(self.updates) + '\n At time:' + str(stime))
        self.updates.clear()
        self.window.clear()
        for  obj in self.drawables :
            print('drawing: ' + str(obj))
            obj.draw(self.trt)
            self.window.update()
        logging.warning('Finished Clearing Updates at time: ' + str(time.time()) + '   time elapsed: ' + str(time.time()-stime))

    def flag_update(self, update:str):
        if not self.updates:
            self.updates.append([update, time.time()])
            self.sch.enter(delay= 1/60, priority= 1, action = self.redraw)
            self.sch.run()
        else:
            self.updates.append([update, time.time()])
class Drawable:
    def draw(self, trt:turtle.Turtle):
        trt.dot(1)
        logging.warning('default draw used')
