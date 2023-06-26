import numpy as np

class TimeFuncs:
        def __init__(self, t = 0,  h = 1/60):
              self.t = t
              self.h = h

        def setTime(self, t):
              self.t = t

        def setStepSize(self, h):
              self.h = h

        def timestep(self, f:np.matrix, u:np.matrix):
            """
            Takes a time step forward using the framework u_y+1 = u_y + du/dt_y+1
            f is a matrix such that fu gives du\dt so u_y+1 = u_y
            
            Parameters
            ----------
            f : np.matrix or array that can be used as matrix
                Derivative matrix
            u : np.matrix or array that can be used as matrix
                Initial state for some variable or collection of variables
            
            Returns
            -------
            u : np.matrix
                Updated state
            t : float
                New time
            """
            f = np.mat(f)
            u = np.mat(u)
            unew = np.linalg.solve(np.eye(f.shape[0]) - self.h*f, u)
            self.t = self.t + self.h
            return [unew, self.t]