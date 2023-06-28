import numpy as np
import logging

class Time_Funcs:
    def __init__(self, t:float = 0,  h:float = 0):
        """
        Instantiates a new TimeFuncs object with time t and stepsize h

        Parameters
        ----------
        t : float
            The current time in seconds
        h : float
            The default stepsize
        """
        self.t = t
        self.h = h

    def set_time(self, t:float):
        """
        Sets the current time

        Parameters
        ----------
        t : float
            Current time in seconds, reccomend using time.time or manually counting time from 0
        """
        self.t = t

    def find_set_stepsize(self, t:float):
        h = t - self.t
        if h <= 0 : 
            h = 0
            logging.warning('Negative or zero stepsize found at time: ' + str(t))
        self.h = h
        self.t = t
        return h

    def setStepSize(self, h:float):
        """
        Manually sets the current stepsize to a value h

        Parameters
        ----------
        h : float
            The new stepsize
        """
        self.h = h

    def time_step(self, f:np.matrix, u:np.matrix, h:float = None):
        """
        Takes a time step forward using the framework u_y+1 = u_y + du/dt_y+1
        f is a matrix such that fu gives du\dt so u_y+1 = u_y
        
        Parameters
        ----------
        f : numpy.matrix or array that can be used as matrix
            Derivative matrix
        u : numpy.matrix or array that can be used as matrix
            Initial state for some variable or collection of variables
            
        Returns
        -------
        u : np.matrix
            Updated state
        """
        if h is None: h = self.h
         
        f = np.mat(f)
        u = np.mat(u)
        unew = np.linalg.solve(np.eye(f.shape[0]) - self.h*f, u)
        return unew
