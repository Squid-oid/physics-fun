import numpy as np
import sdl2.ext
import copy

###########################################################################
# World object classes                                                    #
###########################################################################
class W_object():
    """
    Superclass for all game/physics objects handles coordinate, velocity and sprite storage by default
    """
    def __init__(self, coord:'Coordinate' = None, vel:'Velocity' = None, sprite:sdl2.ext.sprite = None):
        """
        Initializes a new W_object object at the coordinates, velocity and sprite provided

        Parameters
        ----------
        coord : W_object.Coordinate
            The coordinates of the object, if none is given the object will be at the origin
        vel : W_object.Velocity
            The velocity of the object, if none is given the object will have no velocity
        sprite : sdl2.ext.sprite
            The sprite of the object, if none is given the object has no sprite and wil not draw
        """
        if coord is None:
            coord = Coordinate()
        if vel is None:
            vel = Velocity()
        self.coord = coord
        self.vel = vel
        self.sprite = sprite

    def draw(self, sprite_renderer:sdl2.ext.SpriteRenderSystem):
        """
        Draws the object using the sprite renderer provided

        Parameters
        ----------
        sprite_renderer : sdl2.ext.SpriteRenderSystem
            The Sprite renderer with which to draw the sprite
        """
        sprite_renderer.render(self.sprite, round(self.coord.get_x()), round(self.coord.get_y()))

    def refresh(self, args):
        """
        Provides a default implementation of refreshing the object (taking a new timestep)

        Parameters
        ----------
        args : List
            A list of args, should usually be empty when the default refresh is used
        """
        pass


class Barrier(W_object):
    """
    A W_object subclass that implements Barriers which collide with everything behind them
    """
    def __init__(self, coord:'Coordinate' = None, vel:'Velocity' = None, sprite:sdl2.ext.sprite = None, angle:float = 0, fact:sdl2.ext.SpriteFactory = None):
        """
        Initializes a new Barrier object at the coordinates coord, the barrier extends infinetly in each direction and collides with everything behind it 
        
        Parameters
        ----------
        coord : Coordinate
            The coordinate at which the Barrier is centered
        vel : Velocity
            The velocity at which the Barrier would move if it could move
        sprite : sdl2.ext.sprite
            The sprite which the barrier wil be represented by
        angle : float
            The angle of the Barrier
        fact : sdl2.ext.SpriteFactory
            A sprite factory with which to create a default sprite
        """
        self.angle = angle
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
        """
        Draws the Barrier using the SpriteRenderSystem provided (Does not actually do this right now)

        Parameters
        ----------
        sprite_renderer : sdl2.ext.SpriteRenderSystem
            The renderer with which to draw the Barrier
        """
        pass

class Ball(W_object):
    """
    A W_object subclass that moves around and bounces off Barriers
    """
    def __init__(self, coord:'Coordinate' = None, vel:'Velocity' = None):
        """
        Initializes a new Ball object at the provided Coordinates with the provided Velocity

        Parameters
        ----------
        coord : Coordinate
            The coordinate at which to initiliaze the Ball, if none is provided it is created at the origin
        vel : Velocity
            The initial velocity, if none is provided the Ball is stationary
        """
        base_sprites = sdl2.ext.Resources(__file__, "basesprites")
        sprite = sdl2.ext.SpriteFactory(sdl2.ext.SOFTWARE).from_image(base_sprites.get_path("ball.bmp"))
        
        coord = Coordinate(coord)
        vel = Velocity(vel)

        self.radius = 10
        super().__init__(coord = coord, vel = vel, sprite = sprite)

    def bar_overlap(self:'Ball', bar:Barrier):
        """
        Finds and returns the greatest overlapping distance between a ball and a barrier then returns it in the standard basis vectors as a vector

        Parameters
        ----------
        bar : W_object.Barrier
            The barrier to check overlap with
        Returns
        -------
        overlap : numpy.mat 
            The maximum overlap in vector form, directed such that a negative overlap implies that the two are not touching yet
        """
        delta = CVector(self.coord.delta(bar.coord))
        proj_delta = CVector([delta.transform(bar.mat)[0,0], 0])
        proj_gap = proj_delta - CVector([self.radius, 0])
        overlap = np.matmul(proj_gap.vec, bar.mat)
        return overlap
    
    def bar_collide(self:'Ball', bar:Barrier):
        """
        Handles a ball bar collision

        Parameters
        ----------
        bar : W_object.Barrier
            The barrier to check overlap with
        Returns
        -------
        - : bool 
            Wether or not the Ball collided with the Barrier
        """
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
        """
        Overrides the defualt refresh to create a new one which moves the Ball and checks if it has collided with any barriers

        Parameters
        ----------
        args : List of Dictionaries
            A list, if passed correctly it should contain a dictionary with an entry keyed "barriers" containing a list of Barriers  
        """
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