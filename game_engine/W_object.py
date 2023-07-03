import numpy as np
import sdl2.ext
import Timestep as ts
from typing import cast

###########################################################################
# World object classes                                                    #
###########################################################################
class W_object():
    """
    Superclass for all game/physics objects handles coordinate, velocity and sprite storage by default
    """
    def __init__(self, coord:np.mat = None, vel:np.mat = None, sprite:sdl2.ext.sprite = None):
        """
        Initializes a new W_object object at the coordinates, velocity and sprite provided

        Parameters
        ----------
        coord : np.mat/List[float]
            The coordinates of the object, if none is given the object will be at the origin
        vel : np.mat/List[float]
            The velocity of the object, if none is given the object will have no velocity
        sprite : sdl2.ext.sprite
            The sprite of the object, if none is given the object has no sprite and wil not draw
        """
        if coord is None:
            coord = np.mat([0,0])
        if vel is None:
            vel = np.mat([0,0])
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
        sprite_renderer.render(self.sprite, round(self.coord[0,0]), round(self.coord[0,1]))

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
    def __init__(self, coord:np.mat = None, vel:np.mat = None, sprite:sdl2.ext.sprite = None, angle:float = 0, fact:sdl2.ext.SpriteFactory = None):
        """
        Initializes a new Barrier object at the coordinates coord, the barrier extends infinetly in each direction and collides with everything behind it 
        
        Parameters
        ----------
        coord : numpy.mat/List[float]
            The coordinate at which the Barrier is centered
        vel : numpy.mat/List[float]
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
        
        coord = np.mat(coord)
        vel = np.mat(vel)

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
    def __init__(self, coord:np.mat = None, vel:np.mat = None, radius:int = 5):
        """
        Initializes a new Ball object at the provided Coordinates with the provided Velocity

        Parameters
        ----------
        coord : np.mat/List[int]
            The coordinate at which to initiliaze the Ball, if none is provided it is created at the origin
        vel : np.mat/List[int]
            The initial velocity, if none is provided the Ball is stationary
        radius : int
            The radius of the ball
        """
        base_sprites = sdl2.ext.Resources(__file__, "basesprites")
        sprite = sdl2.ext.SpriteFactory(sdl2.ext.SOFTWARE).from_image(base_sprites.get_path("ball.bmp"))
        
        coord = np.mat(coord)
        vel = np.mat(vel)

        self.radius = radius
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
        delta = self.coord - bar.coord
        proj_delta = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(delta)))
        proj_gap = proj_delta - np.mat([self.radius, 0])
        overlap = np.matmul(proj_gap, bar.mat)
        return overlap
    
    def bar_collide(self:'Ball', bar:Barrier):
        """
        Handles a ball bar collision

        Parameters
        ----------
        bar : W_object.Barrier
            The barrier to collide with
            
        Returns
        -------
        - : bool 
            Wether or not the Ball collided with the Barrier
        """
        delta = self.coord - bar.coord
        proj_delta = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(delta)))
        proj_gap = np.mat([proj_delta[0,0] - self.radius, 0])
        if proj_gap[0,0] <= 0 :
            proj_vel = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(self.vel)))
            proj_vel[0,0] = -proj_vel[0,0]
            self.vel = np.matmul(proj_vel, bar.mat)

            proj_coord = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(self.coord)))
            proj_coord = proj_coord - 2*proj_gap                    #The minus 2 factor here makes the ball bounce off the wall continuing it's post reflection path
            self.coord = np.matmul(proj_coord, bar.mat)             #This means that the balls path will progress exactly the same regardless of what frame rate the program runs at
            return True
        
        else:
            return False
    

    def ball_collide(self:'Ball', other:'Ball'):
        """
        Handles a ball on ball collision

        Parameters
        ----------
        other : W_object.Ball
            The Ball to collide  with
            
        Returns
        -------
        - : bool 
            Wether or not the Ball collided with the Barrier
        """
        if other is self:
            return False
        else:
            delta = -self.coord + other.coord
            dist = np.linalg.norm(delta)
            overlap = dist - self.radius - other.radius
            if overlap < 0 : 
                print(overlap)
                base_e_1 = delta/dist                           #The normalized basis vector in the direction of the other balls center from this balls center
                base_e_2 = np.mat([base_e_1[0,1], -base_e_1[0,0]])      #Arbitrary normalized right angle basis vector to e1
                base = np.concatenate((base_e_1, base_e_2))
                               
                proj_vel = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(self.vel)))
                proj_o_vel = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(other.vel)))
                combo_vel = proj_vel[0,0] + proj_o_vel[0,0]
                [proj_vel[0,0], proj_o_vel[0,0]] = [proj_o_vel[0,0], proj_vel[0,0]]

                elapsed_time = overlap/combo_vel

                other.vel = np.matmul(proj_o_vel, base)
                self.vel = np.matmul(proj_vel, base)

                print(self.vel)
                print(other.vel)

                proj_coord = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(self.coord)))
                proj_coord = proj_coord + np.mat([overlap,0]) + np.mat([elapsed_time*proj_vel[0,0]]) #The minus 2 factor here makes the ball bounce off the wall continuing it's post reflection path
                proj_o_coord = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(other.coord)))
                proj_o_coord = proj_o_coord - np.mat([overlap,0]) + np.mat([elapsed_time*proj_o_vel[0,0]]) #The minus 2 factor here makes the ball bounce off the wall continuing it's post reflection path

                other.coord = np.matmul(proj_o_coord, base)
                self.coord = np.matmul(proj_coord, base)             #This means that the balls path will progress exactly the same regardless of what frame rate the program runs at 
                return True
            return False  

    def refresh(self, args):
        """
        Overrides the defualt refresh to create a new one which moves the Ball and checks if it has collided with any barriers

        Parameters
        ----------
        args : List of Dictionaries
            A list, if passed correctly it should contain a dictionary with an entry keyed "barriers" containing a list of Barriers  
        """
        stepper = cast(ts.Time_Funcs, args['stepper'])
        state = np.transpose(np.concatenate((self.coord, self.vel), axis = 1))
        f_mat = np.zeros([4,4])
        f_mat[0:2,2:4] = np.identity(2)
        state = np.transpose(stepper.time_step(f  = f_mat, u = state))
        [self.coord, self.vel] = np.split(state,[2], axis = 1)
        bars = args["barriers"]
        balls = args["balls"]

        collided = True
        while collided is True :
            collided = False
            for bar in bars:
                if self.bar_collide(bar):                   #If a collision occurs a recheck happens in case another collision is caused by the first one resolving
                    collided = True
            for ball in balls:
                if self.ball_collide(ball):                   #If a collision occurs a recheck happens in case another collision is caused by the first one resolving
                    collided = True



##########################################################################
# Collider Class, used in Objects for generic collision detection        #
##########################################################################

## Not yet implemented lmaoski ##