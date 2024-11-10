import numpy as np
import sdl2.ext
import Timestep as ts
from typing import cast
from PIL import Image


###########################################################################
# World object classes                                                    #
###########################################################################
class W_object():
    """
    Superclass for all game/physics objects handles coordinate, velocity and sprite storage by default
    """
    def __init__(self, coord:np.matrix = None, vel:np.matrix = None, sprite:sdl2.ext.sprite = None):
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
            coord = np.asmatrix([0,0])
        if vel is None:
            vel = np.asmatrix([0,0])
        self.coord = coord
        self.vel = vel
        self.sprite = sprite

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
    def __init__(self, coord:np.asmatrix = None, vel:np.asmatrix = None, sprite:sdl2.ext.sprite = None, angle:float = 0, fact:sdl2.ext.SpriteFactory = None):
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
        
        coord = np.asmatrix(coord)
        vel = np.asmatrix(vel)

        if sprite is None:
            base_sprites = sdl2.ext.Resources(__file__, "basesprites")
            sprite = fact.from_image(base_sprites.get_path("barrier.png"))
            sprite.position = -10000,-10000
        self.mat = np.asmatrix([directed_normal,directed_tangent])
        super().__init__(coord = coord, vel = vel, sprite = sprite)
            
class Ball(W_object):
    """
    A W_object subclass that moves around and bounces off Barriers
    """
    def __init__(self, fact:sdl2.ext.SpriteFactory, coord:np.asmatrix = None, vel:np.asmatrix = None, radius:int = 5, mass:float = 1.00):
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
        coord = np.asmatrix(coord)
        vel = np.asmatrix(vel)

        self.mass = mass
        self.radius = radius
        
        with Image.open('game_engine\\basesprites\\ball.png') as im:
            im = im.resize(size = (2*radius,2*radius), resample= Image.BILINEAR)
            im.save('game_engine\\tempsprites\\temp.png')
        temp_sprites = sdl2.ext.Resources("game_engine\\tempsprites")
        sprite = fact.from_image(temp_sprites.get_path("temp.png"))
        super().__init__(coord = coord, vel = vel, sprite = sprite)
        sprite.coord = coord - np.asmatrix([radius, radius])

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
        proj_gap = proj_delta - np.asmatrix([self.radius, 0])
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
        proj_gap = np.asmatrix([proj_delta[0,0] - self.radius, 0])
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
                base_e_1 = delta/dist                           #The normalized basis vector in the direction of the other balls center from this balls center
                base_e_2 = np.asmatrix([base_e_1[0,1], -base_e_1[0,0]])      #Arbitrary normalized right angle basis vector to e1
                base = np.concatenate((base_e_1, base_e_2))

                m_a = self.mass #self.mass
                m_b = other.mass #other.mass

                proj_vel = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(self.vel)))
                proj_o_vel = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(other.vel)))
                u_a = proj_vel[0,0]
                u_b = proj_o_vel[0,0]
                combo_vel = u_a + u_b
                elapsed_time = overlap/combo_vel

                coeff_sq = m_a**2 + m_a*m_b
                coeff_li = -2*m_a**2*u_a - 2*m_a*m_b*u_b
                coeff_co = m_a**2*u_a**2 + 2*m_a*m_b*u_a*u_b - m_a*m_b*u_a**2
                
                roots = np.roots([coeff_sq, coeff_li, coeff_co])
                v_a = min(roots)
                v_b =  (-m_a*v_a + m_a*u_a + m_b*u_b)/(m_b)

                proj_vel[0,0] = v_a
                proj_o_vel[0,0] = v_b
                
                other.vel = np.matmul(proj_o_vel, base)
                self.vel = np.matmul(proj_vel, base)               
                
                proj_coord = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(self.coord)))
                proj_coord = proj_coord - np.asmatrix([u_a*elapsed_time, 0]) + np.asmatrix([v_a*elapsed_time, 0])
                proj_o_coord = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(other.coord)))
                proj_o_coord = proj_o_coord - np.asmatrix([u_b*elapsed_time, 0]) + np.asmatrix([v_b*elapsed_time, 0])

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
                    collided = False
        
        self.sprite.position = round(self.coord[0,0] - self.radius), round(self.coord[0,1] - self.radius)


class Tri(W_object):
    @staticmethod
    def generateSprite(coords, fact):
        """
        Generates a sprite from the coordinates of the corners of a triangle
        """
        ## Finds the size of the bounding rectangle
        bounding_box = np.max(coords,0) - np.min(coords,0)
        
        ## Creates an array of overscaled but correct size that we will then scale down(naive way to improve aliasing)
        sz1 = int(np.rint(bounding_box[0,0]*3))
        sz2 = int(np.rint(bounding_box[0,1]*3))
        arr = np.asmatrix(np.zeros((sz1,sz2)))

        # Generates Clockwise coords
        c1 = coords[1,:] - coords[0,:]
        c2 = coords[2,:] - coords[1,:]
        cross = c1[0,0]*c2[0,1] - c1[0,1]*c2[0,0]
        print(cross)
        if cross > 0:
            coords = np.flip(coords,1)

        # Generates counterclockwise edges from vertices
        side1 = coords[0,:] - coords[1,:]
        side2 = coords[1,:] - coords[2,:]
        side3 = coords[2,:] - coords[0,:]

        ## Finds orthogonal outwardsfacing vector
        orth1 = np.asmatrix([side1[0,1], -side1[0,0]])
        orth2 = np.asmatrix([side2[0,1], -side2[0,0]])
        orth3 = np.asmatrix([side3[0,1], -side3[0,0]])

        ## Now we wish to shade all pixels that are within the triangle, we do this by checking each pixel against the polytope criterion
        c_const = np.asmatrix([1/2,1/2]).T
        for i in range(0,sz1):
            for j in range(0,sz2):
                x = (np.asmatrix([i,j]).T + c_const) / 3 + np.min(coords,0).T
                cond1 = all(orth1*x <= orth1*coords[1,:].T)
                cond2 = all(orth2*x <= orth2*coords[2,:].T)
                cond3 = all(orth3*x <= orth3*coords[0,:].T)

                if all([cond1,cond2,cond3]):
                    arr[i,j] = 255

        ## Now we use the array together with pysdl to create a sprite
        arr = np.flip(arr,0)
        arr = arr.reshape(sz1*sz2)
        arr = np.pad(arr.T, ((0,0),(3,0)))
        arr.shape = (sz1,sz2,4)
        upscaled_img = Image.fromarray(arr.astype(np.uint8))
                
        im = upscaled_img.resize(size = (sz1,sz2), resample= Image.BILINEAR)

        sprite = fact.from_image(upscaled_img)
        return sprite


    def __init__(self, fact:sdl2.ext.SpriteFactory, coords:np.asmatrix = None, vel:np.asmatrix = None, mass:float = 1.00):
        """
        Initializes a new Triangle object at the provided Coordinates with the provided Velocity

        Parameters
        ----------
        coords : np.mat
            The coordinates of the corners of the triangle, each row is a corner
        vel : np.mat/List[int]
            The initial velocity, if none is provided the Ball is stationary
        """
        self.coords = np.asmatrix(coords)
        vel = np.asmatrix(vel)

        self.mass = mass

        avg_point = np.ones((1,3)) @ self.coords / 3

        print("Attempting to Generate Sprite")
        sprite = Tri.generateSprite(coords, fact)
        super().__init__(coord = avg_point, vel = vel, sprite = sprite)
        self.sprite.position = round(np.min(self.coords,1)[0,0]), round(np.min(self.coords,1)[1,0])
        
            
    def bar_overlap(self:'Tri', bar:Barrier):
        """
        Finds how far along the triangles velocity vector it has travelled into a barrier

        First finds deepest point then projects the velocity vector onto the overlap vector.

        Parameters
        ----------
        bar : W_object.Barrier
            The barrier to check overlap with

        Returns
        -------
        overlap : numpy.mat 
            The maximum overlap in vector form, directed such that a negative overlap implies that the two are not touching yet
        """
        
        ## The deepest point will always be one of the corners
        max_proj_delta = np.zeros((1,2))
        for coord in self.coords:
            delta = coord - bar.coord
            proj_delta = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(delta)))
            if max_proj_delta[0,0] < proj_delta[0,0]:
                max_proj_delta = proj_delta

        ## If the item has velocity, ie could have clipped into the barrier of it's own accord we want to push it out along it's own vector
        if not self.vel is None:
            ## Find the velocity expressed in the bar's coordinate system
            proj_vel = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(self.vel)))
            proj_vel_hat = proj_vel/np.linalg.norm(proj_vel)

            ## Find the vector giving how far into the barrier the point has travelled, in the bar coordinate system
            proj_gap = max_proj_delta[0,0]/proj_vel_hat[0,0] * proj_vel_hat

            ## Transform back into global coordinate system
            overlap = np.matmul(proj_gap, bar.mat)
        ## If instead we have no velocity, just shove it out the shortest path
        else:
            overlap = np.matmul(max_proj_delta, bar.mat)

        return overlap


    def bar_collide(self:'Tri', bar:Barrier):
        """
        Handles a triangle bar collision

        Parameters
        ----------
        bar : W_object.Barrier
            The barrier to collide with
            
        Returns
        -------
        - : bool 
            Wether or not the Tri collided with the Barrier
        """
        max_proj_delta = np.zeros((1,2))
        for coord in self.coords:
            delta = coord - bar.coord
            proj_delta = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(delta)))
            if max_proj_delta[0,0] < proj_delta[0,0]:
                max_proj_delta = proj_delta

        ## Check if there is an overlap
        if max_proj_delta[0,0] < 0:
            ## If the item has velocity, ie could have clipped into the barrier of it's own accord we want to push it out along it's own vector
            if not self.vel is None:
                ## Find the velocity expressed in the bar's coordinate system
                proj_vel = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(self.vel)))
                proj_vel[0,0] = -proj_vel[0,0]

                self.vel = np.matmul(proj_vel, bar.mat)
                
            ## If instead we have no velocity, just shove it out the shortest path
            else:
                pass
            
            projected_offset = - 2*np.asmatrix([max_proj_delta[0,0],0])
            for coord in self.coords:
                coord = coord + np.matmul(projected_offset, bar.mat)

            return True
        else:
            return False

    def refresh(self, args):
        pass

##########################################################################
# Collider Class, used in Objects for generic collision detection        #
##########################################################################
