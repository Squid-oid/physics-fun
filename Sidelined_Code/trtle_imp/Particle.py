#The basic particle, a circle with a radius and a coordinate
import numpy as np
from Gfx import Drawable
from Gfx import Display
import logging
import Boundary
import turtle
import Timestep

class Particle(Drawable):
    def __init__(self, disp:Display = None, radius:float = 1, coordinates:list = [0,0], speed:list = [0,0], acceleration:list = [0,0]):
        if(radius < 0):
            radius = 0
        self.radius = radius
        self.coordinates = coordinates
        self.speed = speed
        self.acceleration = acceleration
        self.disp = disp
        if not disp is None:
            disp.add_drawable(self)

    def update(self, stepper:Timestep.TimeFuncs):
        pre_u = [self.coordinates + self.speed + self.acceleration]
        u = np.transpose(np.mat(pre_u))
        f = np.mat([[1,0,1,0,0,0],
                    [0,1,0,1,0,0],
                    [0,0,1,0,1,0],
                    [0,0,0,1,0,1],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])
        [u,t] = stepper.timestep(f,u)
        u = np.transpose(u).tolist()[0]
        [self.coordinates, self.speed, self.acceleration] = [[u[0],u[1]],[u[2],u[3]],[u[4],u[5]]]
        for bnd in self.disp.bnds:
            overlap = self.overlapsBoundary(bnd)
            if overlap[0]:
                overlap = overlap[1]
                rebasedT = np.linalg.solve(np.transpose(bnd.mat),np.transpose(self.coordinates))
                rebased = np.transpose(rebasedT)
                rebased = rebased - [overlap[0,0],0]
                print(rebased)
                print(bnd.mat) 
                print('////////////////////////')
                self.coordinates = np.transpose(np.matmul(bnd.mat, np.transpose(rebased)))
                print(self.coordinates)
                rebasedT = np.linalg.solve(np.transpose(bnd.mat),np.transpose(self.speed))
                rebased = np.transpose(rebasedT)
                rebased = [-rebased[0],rebased[0]]
                self.speed = np.matmul(bnd.mat, rebased)
        if not self.disp is None:
            self.disp.flag_update("updated particle")

    def overlapsParticle(self, other:'Particle'):
        if self.disp == other.disp:
            delta_coords = self.coordinates - other.coordinates
            distance = np.linalg.norm(delta_coords)
            combined_radius = self.radius + other.radius
            return distance < combined_radius
        else:
            logging.warning('mismatched display comparison attempted, (particle particle) /n    ' + self + '/n    ' + other)
            return False
        
    
    def overlapsBoundary(self, boundary:Boundary):
        if self.disp == boundary.disp:
            rebasedT = np.linalg.solve(np.transpose(boundary.mat),np.transpose(self.coordinates))        # See notes
            rebased = np.transpose(rebasedT)
            rebased_delta = rebased - boundary.local_offset
            return [rebased_delta[0,0] < self.radius, rebased_delta]
        else:
            logging.warning('mismatched display comparison attempted, (particle boundary) /n    ' + self + '/n    ' + boundary)
            return False                                      

    def updateCoords(self, coordinates:list = None, delta:list = None):
        if not coordinates is None:
            self.coordinates = np.mat(coordinates)
        if not delta is None :
            self.coordinates = self.coordinates + np.mat(delta)
        if not self.disp is None:
            self.disp.flag_update("moved particle")
         

    def draw(self, trt:turtle.Turtle):
        trt.setheading(0)
        trt.pen(pendown = False)
        trt.setposition(x = self.coordinates[0], y = self.coordinates[1])
        trt.pen(pendown = True, pencolor='blue')
        trt.forward(0)
        trt.dot(self.radius*2, "blue")

