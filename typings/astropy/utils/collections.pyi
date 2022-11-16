"""
This type stub file was generated by pyright.
"""

"""
A module containing specialized collection classes.
"""
class HomogeneousList(list):
    """
    A subclass of list that contains only elements of a given type or
    types.  If an item that is not of the specified type is added to
    the list, a `TypeError` is raised.
    """
    def __init__(self, types, values=...) -> None:
        """
        Parameters
        ----------
        types : sequence of types
            The types to accept.

        values : sequence, optional
            An initial set of values.
        """
        ...
    
    def __iadd__(self, other): # -> Self@HomogeneousList:
        ...
    
    def __setitem__(self, idx, value): # -> None:
        ...
    
    def append(self, x): # -> None:
        ...
    
    def insert(self, i, x): # -> None:
        ...
    
    def extend(self, x): # -> None:
        ...
    


