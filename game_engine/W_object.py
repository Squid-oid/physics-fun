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

    def collide(self, other):
        return Collision.collide(self, other)

    def refresh(self, args):
        """
        Provides a default implementation of refreshing the object (taking a new timestep)

        Parameters
        ----------
        args : List
            A list of args, includes the world stepper and a list of all world objects
        """
        stepper = cast(ts.Time_Funcs, args['stepper'])
        objs = args['objs']

        state = np.transpose(np.concatenate((self.coord, self.vel), axis = 1))
        f_mat = np.zeros([4,4])
        f_mat[0:2,2:4] = np.identity(2)
        state = np.transpose(stepper.time_step(f  = f_mat, u = state))
        [self.coord, self.vel] = np.split(state,[2], axis = 1)

        if isinstance(self, Tri):
            self.coords = self.coord + self.offsets # Update Coords to match center

        collided = True
        while collided is True :
            collided = False
            for obj in objs:
                if self.collide(obj):                   #If a collision occurs a recheck happens in case another collision is caused by the first one resolving
                    collided = True


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

    
    def refresh(self, args):
        """ Since we know that barriers/half planes never move, we can leave refreshing and colliding to Tris and Balls """
        pass
            
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
       
    def refresh(self, args):
        """
        Overrides the defualt refresh to create a new one which moves the BallSprite

        Parameters
        ----------
        args : Dcit
            A dict, if passed correctly it should contain a list of all objects under the key 'objs' and a stepper under the key 'stepper'
        """
        # Run Update
        super().refresh(args)        
        # Update the sprite
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
            coords = np.flip(coords,0)

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

        # Transpose the array to get it to y,x format for image
        arr = arr.T

        ## Now we use the array together with pysdl to create a sprite
        arr = arr.reshape(sz1*sz2)
        arr = np.pad(arr.T, ((0,0),(3,0)))
        arr.shape = (sz1,sz2,4)
        upscaled_img = Image.fromarray(arr.astype(np.uint8))
                
        im = upscaled_img.resize(size = (bounding_box[0,0], bounding_box[0,1]), resample= Image.BILINEAR)

        sprite = fact.from_image(im)
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

        # Find 'center' of triangle
        avg_point = np.ones((1,3)) @ self.coords / 3
        # Find mapping from center to vertices
        self.offsets = self.coords - avg_point

        print("Attempting to Generate Sprite")
        sprite = Tri.generateSprite(coords, fact)
        super().__init__(coord = avg_point, vel = vel, sprite = sprite)
        self.sprite_offset = [-avg_point[0,0] + np.min(self.coords[0,:]) , -avg_point[0,1] + np.min(self.coords[1,:])]
        self.sprite.position = round(self.coord[0,0] + self.sprite_offset[0]) , round(self.coord[0,1] + self.sprite_offset[1])
        
    def refresh(self, args):
        """
        Overrides the default refresh to create a new one which moves the TriSprite

        Parameters
        ----------
        args : Dcit
            A dict, if passed correctly it should contain a list of all objects under the key 'objs' and a stepper under the key 'stepper'
        """
        # Update pos and collide
        super().refresh(args)
        # Move Sprite
        self.sprite.position = round(self.coord[0,0] + self.sprite_offset[0]) , round(self.coord[0,1] + self.sprite_offset[1])
        pass

##########################################################################
# Collider Class, used in Objects for generic collision detection        #
##########################################################################
class Collision:
    """
    Static class for handling collision and overlap calculations
    between world objects like Ball, Barrier, and Tri.
    """

    @ staticmethod 
    def collide(obj1:'W_object', obj2:'W_object') -> bool:
        """
        Dispatches to the correct collision function depending on types.
        Returns True if a collision occurred, False otherwise.
        """
        from W_object import Ball, Barrier, Tri  

        # Balls and Walls
        if isinstance(obj1, Ball) and isinstance(obj2, Barrier):
            return Collision.ball_bar_collide(obj1, obj2)
        elif isinstance(obj1, Barrier) and isinstance(obj2, Ball):
            return Collision.ball_bar_collide(obj2, obj1)  # reverse order

        # Balls and Balls
        elif isinstance(obj1, Ball) and isinstance(obj2, Ball):
            return Collision.ball_ball_collide(obj1, obj2)

        # Tris and Walls
        elif isinstance(obj1, Tri) and isinstance(obj2, Barrier):
            return Collision.tri_bar_collide(obj1, obj2)
        elif isinstance(obj1, Barrier) and isinstance(obj2, Tri):
            return Collision.tri_bar_collide(obj2, obj1)

        # Tris and Balls
        elif isinstance(obj1, Tri) and isinstance(obj2, Ball):
            return Collision.tri_ball_collide(obj1, obj2)
        elif isinstance(obj1, Ball) and isinstance(obj2, Tri):
            return Collision.tri_ball_collide(obj2, obj1)

        # Tris and Tris
        elif isinstance(obj1, Tri) and isinstance(obj2, Tri):
            return Collision.tri_tri_collide(obj1, obj2)  
         
        return False



    @staticmethod
    def ball_bar_overlap(ball:'Ball', bar:'Barrier') -> np.matrix:
        """
        Finds and returns the greatest overlapping distance between a ball and a barrier.
        Negative overlap means not touching yet.
        """
        delta = ball.coord - bar.coord
        proj_delta = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(delta)))
        proj_gap = proj_delta - np.asmatrix([ball.radius, 0])
        overlap = np.matmul(proj_gap, bar.mat)
        return overlap

    @staticmethod
    def ball_bar_collide(ball:'Ball', bar:'Barrier') -> bool:
        """
        Handles a ball-barrier collision.
        """
        delta = ball.coord - bar.coord
        proj_delta = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(delta)))
        proj_gap = np.asmatrix([proj_delta[0,0] - ball.radius, 0])

        if proj_gap[0,0] <= 0:
            proj_vel = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(ball.vel)))
            proj_vel[0,0] = -proj_vel[0,0]
            ball.vel = np.matmul(proj_vel, bar.mat)

            proj_coord = np.transpose(np.linalg.solve(np.transpose(bar.mat), np.transpose(ball.coord)))
            proj_coord = proj_coord - 2*proj_gap
            ball.coord = np.matmul(proj_coord, bar.mat)
            return True

        return False

    @staticmethod
    def ball_ball_collide(ball:'Ball', other:'Ball') -> bool:
        """
        Handles a ball-on-ball collision.
        """
        if other is ball:
            return False

        delta = -ball.coord + other.coord
        dist = np.linalg.norm(delta)
        overlap = dist - ball.radius - other.radius

        if overlap < 0:
            base_e_1 = delta/dist
            base_e_2 = np.asmatrix([base_e_1[0,1], -base_e_1[0,0]])
            base = np.concatenate((base_e_1, base_e_2))

            m_a, m_b = ball.mass, other.mass

            proj_vel = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(ball.vel)))
            proj_o_vel = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(other.vel)))
            u_a, u_b = proj_vel[0,0], proj_o_vel[0,0]
            combo_vel = u_a + u_b
            elapsed_time = overlap/combo_vel

            coeff_sq = m_a**2 + m_a*m_b
            coeff_li = -2*m_a**2*u_a - 2*m_a*m_b*u_b
            coeff_co = m_a**2*u_a**2 + 2*m_a*m_b*u_a*u_b - m_a*m_b*u_a**2

            roots = np.roots([coeff_sq, coeff_li, coeff_co])
            v_a = min(roots)
            v_b = (-m_a*v_a + m_a*u_a + m_b*u_b)/(m_b)

            proj_vel[0,0] = v_a
            proj_o_vel[0,0] = v_b

            other.vel = np.matmul(proj_o_vel, base)
            ball.vel = np.matmul(proj_vel, base)

            proj_coord = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(ball.coord)))
            proj_coord = proj_coord - np.asmatrix([u_a*elapsed_time, 0]) + np.asmatrix([v_a*elapsed_time, 0])
            proj_o_coord = np.transpose(np.linalg.solve(np.transpose(base), np.transpose(other.coord)))
            proj_o_coord = proj_o_coord - np.asmatrix([u_b*elapsed_time, 0]) + np.asmatrix([v_b*elapsed_time, 0])

            other.coord = np.matmul(proj_o_coord, base)
            ball.coord = np.matmul(proj_coord, base)
            return True

        return False

    @staticmethod
    def tri_bar_collide(tri:'Tri', bar:'Barrier') -> bool:
        """
        Handles a triangle-barrier collision.
        Checks all vertices and resolves penetration along the barrier normal.
        Reflects velocity along the barrier normal.
        """
        eps = 1e-6

        # Barrier normal (first row of bar.mat)
        bar_normal = np.asmatrix(bar.mat[0,:]) 

        # Project all triangle vertices onto barrier normal
        deltas = tri.coords - bar.coord  
        projections = deltas @ bar_normal.T 

        # Find the deepest penetration (most negative projection)
        max_penetration = np.min(projections)

        if max_penetration < 0:  # collision occurred
            if tri.vel is not None:
                # Reflect velocity along the barrier normal
                v_dot_n = float(tri.vel @ bar_normal.T)
                tri.vel = tri.vel - 2 * v_dot_n * bar_normal

            # Move the triangle out along the barrier normal
            correction = -(max_penetration - eps)* bar_normal  # positive vector, we add a tiny offset since this collision is prone to getting stuck
            tri.coords = tri.coords + correction
            tri.coord = tri.coord + correction
            print("Correction applied:", correction)  # DEBUG
            print("New triangle coords:\n", tri.coords)  # DEBUG

            return True

        return False
    
    @staticmethod
    def tri_ball_collide(tri:'Tri', ball:'Ball') -> bool:
        """
        Handles a Triange Ball Collision
        """
        # Find closest triangle vertex to ball center
        closest_point = None
        min_dist = float('inf')
        for pt in tri.coords:
            delta = ball.coord - pt
            dist = np.linalg.norm(delta)
            if dist < min_dist:
                min_dist = dist
                closest_point = pt

        # Check if ball overlaps with the point
        if min_dist < ball.radius:
            #  Compute normal vector from point to ball
            normal = ball.coord - closest_point
            normal_hat = normal / np.linalg.norm(normal)

            #  Reflect ball velocity along normal
            v_dot_n = np.dot(ball.vel, normal_hat.T)[0,0]  # scalar projection
            ball.vel = ball.vel - 2*v_dot_n*normal_hat

            # Push ball out of triangle along normal
            penetration = ball.radius - min_dist
            ball.coord = ball.coord + normal_hat * penetration

            return True

        return False

    # -------------------
    # Tri–Tri (placeholder)
    # -------------------
    @staticmethod
    def tri_tri_collide(tri1:'Tri', tri2:'Tri') -> bool:
        """
        Placeholder for triangle–triangle collision.
        Currently unimplemented.
        """
        #print("Tri–Tri collision not implemented yet.")
        return False