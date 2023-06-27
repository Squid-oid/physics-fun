#The basic Boundary, a line with a direction in which things may not exist
import numpy as np
import turtle
from Gfx import Drawable
from Gfx import Display

class Boundary(Drawable):
    def __init__(self, disp:Display, angle:float = 0, direction:bool = True, offset:list = [0,0]):
        sin = np.sin(angle)
        cos = np.cos(angle)
        self.angle = angle
        self.direction = direction
        directedNormal = [v*(-1)**direction for v in [sin, -cos]]
        directedTangent = [cos, sin]
        self.mat = np.mat([directedNormal,directedTangent])
        self.offset = np.mat(offset)                                            #The offset in absolute coordinates
        local_offsetT = np.linalg.solve(np.transpose(self.mat),np.transpose(self.offset))                      #Finds the offset expressed in the boundaries direction vectors
        self.local_offset = np.transpose(local_offsetT)
        self.disp = disp
        disp.add_boundary(self)


    def draw(self, trt:turtle.Turtle):
        trt.pen(pendown = False, pencolor= "black")
        trs = trt.getscreen()
        print(trs.screensize())
        mLength = np.linalg.norm(trs.screensize())
        trt.setposition(x = self.offset[0,0], y = self.offset[0,1])
        trt.setheading(float(self.angle)*float(180/(np.pi)))
        trt.pen(pendown = True)
        trt.forward(mLength)
        trt.setposition(x = self.offset[0,0], y = self.offset[0,1])
        trt.forward(-mLength)

        trt.pen(pendown = False, pencolor= "red")
        trt.setposition(x = self.offset[0,0], y = self.offset[0,1])
        trt.pen(pendown = True) 
        trt.setheading(float(self.angle-np.pi/2)*float(180/(np.pi)))
        trt.forward(25*(-1)**self.direction)
        trt.right(165)
        trt.forward(10*(-1)**self.direction)
        trt.forward(-10*(-1)**self.direction)
        trt.right(30)
        trt.forward(10*(-1)**self.direction)
