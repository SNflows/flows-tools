"""
This type stub file was generated by pyright.
"""

import abc

"""
Optimization algorithms used in `~astropy.modeling.fitting`.
"""
__all__ = ["Optimization", "SLSQP", "Simplex"]
DEFAULT_MAXITER = ...
DEFAULT_EPS = ...
DEFAULT_ACC = ...
DEFAULT_BOUNDS = ...
class Optimization(metaclass=abc.ABCMeta):
    """
    Base class for optimizers.

    Parameters
    ----------
    opt_method : callable
        Implements optimization method

    Notes
    -----
    The base Optimizer does not support any constraints by default; individual
    optimizers should explicitly set this list to the specific constraints
    it supports.

    """
    supported_constraints = ...
    def __init__(self, opt_method) -> None:
        ...
    
    @property
    def maxiter(self): # -> int:
        """Maximum number of iterations"""
        ...
    
    @maxiter.setter
    def maxiter(self, val): # -> None:
        """Set maxiter"""
        ...
    
    @property
    def eps(self): # -> DEFAULT_EPS:
        """Step for the forward difference approximation of the Jacobian"""
        ...
    
    @eps.setter
    def eps(self, val): # -> None:
        """Set eps value"""
        ...
    
    @property
    def acc(self): # -> float:
        """Requested accuracy"""
        ...
    
    @acc.setter
    def acc(self, val): # -> None:
        """Set accuracy"""
        ...
    
    def __repr__(self): # -> str:
        ...
    
    @property
    def opt_method(self): # -> Unknown:
        """ Return the optimization method."""
        ...
    
    @abc.abstractmethod
    def __call__(self):
        ...
    


class SLSQP(Optimization):
    """
    Sequential Least Squares Programming optimization algorithm.

    The algorithm is described in [1]_. It supports tied and fixed
    parameters, as well as bounded constraints. Uses
    `scipy.optimize.fmin_slsqp`.

    References
    ----------
    .. [1] http://www.netlib.org/toms/733
    """
    supported_constraints = ...
    def __init__(self) -> None:
        ...
    
    def __call__(self, objfunc, initval, fargs, **kwargs): # -> tuple[Unknown, dict[str, None]]:
        """
        Run the solver.

        Parameters
        ----------
        objfunc : callable
            objection function
        initval : iterable
            initial guess for the parameter values
        fargs : tuple
            other arguments to be passed to the statistic function
        kwargs : dict
            other keyword arguments to be passed to the solver

        """
        ...
    


class Simplex(Optimization):
    """
    Neald-Mead (downhill simplex) algorithm.

    This algorithm [1]_ only uses function values, not derivatives.
    Uses `scipy.optimize.fmin`.

    References
    ----------
    .. [1] Nelder, J.A. and Mead, R. (1965), "A simplex method for function
       minimization", The Computer Journal, 7, pp. 308-313
    """
    supported_constraints = ...
    def __init__(self) -> None:
        ...
    
    def __call__(self, objfunc, initval, fargs, **kwargs): # -> tuple[Unknown, dict[str, None]]:
        """
        Run the solver.

        Parameters
        ----------
        objfunc : callable
            objection function
        initval : iterable
            initial guess for the parameter values
        fargs : tuple
            other arguments to be passed to the statistic function
        kwargs : dict
            other keyword arguments to be passed to the solver

        """
        ...
    


