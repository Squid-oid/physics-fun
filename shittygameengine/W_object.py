import numpy as np
import sdl2.ext
import copy

###########################################################################
# World object classes                                                    #
###########################################################################
class W_object():
    def __init__(self, coord:'Coordinate' = None, vel:'Velocity' = None, sprite:sdl2.ext.sprite = None):
        if coord is None:
            coord = Coordinate()
        if vel is None:
            vel = Velocity()
        self.coord = coord
        self.vel = vel
        self.sprite = sprite

    def draw(self, sprite_renderer:sdl2.ext.SpriteRenderSystem):
        sprite_renderer.render(self.sprite, round(self.coord.get_x()), round(self.coord.get_y()))

    def refresh(self, args):
        pass


class Barrier(W_object):
    def __init__(self, coord:'Coordinate' = None, vel:'Velocity' = None, sprite:sdl2.ext.sprite = None, angle:float = 0, fact:sdl2.ext.SpriteFactory = None):
        self.angle = 0 ##Temp flat angle only
        sin = np.sin(angle)
        cos = np.cos(angle)
        directed_normal = [sin, -cos]
        directed_tangent = [cos, sin]
        
        coord = Coordinate(coord)
        vel = Velocity(vel)

        if sprite is None:
            base_sprites = sdl2.ext.Resources(__file__, "basesprites")
            sprite = fact.from_image(base_sprites.get_path("barrier.bmp"))
            
        
        self.mat = np.mat([directed_normal,directed_tangent])
        super().__init__(coord = coord, vel = vel, sprite = sprite)
            
    def draw(self, sprite_renderer:sdl2.ext.SpriteRenderSystem):
        pass

class Ball(W_object):
    def __init__(self, coord:'Coordinate' = None, vel:'Velocity' = None):
        base_sprites = sdl2.ext.Resources(__file__, "basesprites")
        sprite = sdl2.ext.SpriteFactory(sdl2.ext.SOFTWARE).from_image(base_sprites.get_path("ball.bmp"))
        
        coord = Coordinate(coord)
        vel = Velocity(vel)

        self.radius = 10
        super().__init__(coord = coord, vel = vel, sprite = sprite)

    def bar_overlap(self:'Ball', bar:Barrier):
        '''
        Finds and returns the greatest overlapping distance between a ball and a barrier then returns it in the standard basis vectors as a vector,
        args:
            bar:Barrier
                The barrier to check overlap with
        returns:
            overlap:np.mat 
                The maximum overlap in vector form, directed such that a negative overlap implies that the two are not touching yet
        '''
        delta = CVector(self.coord.delta(bar.coord))
        proj_delta = CVector([delta.transform(bar.mat)[0,0], 0])
        proj_gap = proj_delta - CVector([self.radius, 0])
        overlap = np.matmul(proj_gap.vec, bar.mat)
        return overlap
    
    def bar_collide(self:'Ball', bar:Barrier):
        '''
        Handles a bar collision
        args:
            bar:Barrier
                The barrier to check overlap with
        returns:
            overlap:np.mat 
                The maximum overlap in vector form, directed such that a negative overlap implies that the two are not touching yet
        '''
        delta = CVector(self.coord.delta(bar.coord))
        proj_delta = CVector([delta.transform(bar.mat)[0,0], 0])
        proj_gap = proj_delta - CVector([self.radius, 0])
        if proj_gap.vec[0,0] <= 0 :
            proj_vel = CVector(self.vel.transform(bar.mat))
            proj_vel.vec[0,0] = -proj_vel.vec[0,0]
            self.vel = Velocity(np.matmul(proj_vel.vec, bar.mat))

            proj_coord = CVector(self.coord.transform(bar.mat))
            proj_coord = proj_coord - proj_gap
            self.coord = Coordinate(np.matmul(proj_coord.vec, bar.mat))
            return True
        
        else:
            return False
    
    def refresh(self, args):
        self.coord = self.coord + self.vel
        bars = args["barriers"]
        for bar in bars:
            self.bar_collide(bar)
##########################################################################
# Collider Class, used in Objects for generic collision detection        #
##########################################################################


##########################################################################
# Vector types below, eg Coordinates w/ some special functions and       #
# a generic CVector parent for all the rest, more may be added as needed #
##########################################################################
class CVector():
    def __init__(self, vec:np.matrix = [0,0]):
        if isinstance(vec, np.matrix) or isinstance(vec, list):
            self.vec = np.mat(vec)
        elif issubclass(type(vec), CVector) or isinstance(vec, CVector):
            self.vec = np.mat(vec.vec)
        else:
            self.vec = np.mat([0,0])
        
    def __add__(self, other:'CVector'):
        return type(self)(self.vec + other.vec)

    def __sub__(self, other:'CVector'):
        return type(self)(self.vec - other.vec)
    
    def delta(self, other:'CVector'):
        return self.vec - other.vec

    def dist(self, other:'CVector'):
        return np.linalg.norm(self.vec - other.vec)
    
    def transform(self, base:np.mat):
        return np.transpose(np.linalg.solve(np.transpose(base), np.transpose(self.vec)))            #Go study linalg if you don't get this
    

class Coordinate(CVector):
    def __init__(self, coords = [0,0]):
        super().__init__(coords)

    def coords(self):
        return self.vec
    
    def get_x(self):
        return self.vec[0,0]
    
    def get_y(self):
        return self.vec[0,1]


class Velocity(CVector):
    def __init__(self, vel = [0,0]):
        super().__init__(vel)

    def vel(self):
        return self.vec