import astropy.units as u
import numpy as np


class Corner:

    def __init__(self, x: u.Quantity, y: u.Quantity, hw: u.Quantity):
        self.x = x
        self.y = y
        self.hw = hw
        self.corner_xy = self.set_corners()
        self.corners = self.process_corner()

    def set_corners(self):
        """Get corners of a rectangle for a given ra dec and side length
        :return: list[corner 1, corner 2, ..]
        """
        return [(self.x - self.hw, self.y - self.hw), (self.x - self.hw, self.y + self.hw),
                (self.x + self.hw, self.y + self.hw), (self.x + self.hw, self.y - self.hw)]

    def process_corner(self):
        """Given corners of a rectangle defined as astropy Quantity objects, return it as an np array of floats"""
        _points = [u.quantity.Quantity(corner) for corner in self.corner_xy]
        return np.array(_points)
