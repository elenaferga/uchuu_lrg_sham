#! /usr/bin/env python
import numpy as np
from scipy.interpolate import splev, splrep


class Coordinates(object):
    
    def __init__(self, cosmology):
        """
        Class containing useful methods for converting between cartesian and equitorial coordinates

        Args:
            cosmology: instance of astropy.cosmology class
        """
        self.cosmology = cosmology
        self.__initialize_redshift_interpolator()


    def __initialize_redshift_interpolator(self):
        # cubic spline interpolation for converting comoving distance to redshift
        z = np.arange(0, 10, 0.0001)
        rcom = self.comoving_distance(z)
        self.__tck = splrep(rcom, z)
    
    
    def redshift(self, comoving_distance):
        """
        Convert comoving distance to redshift

        Args:
            distance: comoving distance in units [Mpc/h]
        Returns:
            array of redshift
        """
        return splev(comoving_distance, self.__tck)
    
    
    def comoving_distance(self, redshift):
        """
        Convert redshift to comoving distance

        Args:
            array of redshift
        Returns:
            array of comoving distance in units [Mpc/h]
        """
        return self.cosmology.comoving_distance(redshift).value*self.cosmology.h


    def equitorial_to_pos3d(self, ra, dec, z):
        """
        Convert ra, dec, z to 3d cartesian coordinates

        Args:
            ra:  array of ra [deg]
            dec: array of dec [deg]
            z:   array of redshift
        Returns:
            2d array of position vectors [Mpc/h]
        """
        # convert degrees to radians
        ra = ra.copy()*np.pi / 180
        dec = dec.copy()*np.pi / 180

        # comoving distance to redshift z
        r_com = self.comoving_distance(z)

        pos = np.zeros((len(r_com),3))
        pos[:,0] = r_com * np.cos(ra) * np.cos(dec) # x coord
        pos[:,1] = r_com * np.sin(ra) * np.cos(dec) # y coord
        pos[:,2] = r_com * np.sin(dec)              # z coord

        return pos


    def pos3d_to_equitorial(self, pos):
        """
        Convert 3d cartesian coordinates to ra, dec, z

        Args:
            pos: 2d array of position vectors [Mpc/h]
        Returns:
            ra:  array of ra [deg]
            dec: array of dec [deg]
            z:   array of redshift
        """

        # get ra
        ra = np.arctan(pos[:,1] / pos[:,0])
        ind = np.logical_and(pos[:,1] < 0, pos[:,0] >= 0)
        ra[ind] += 2*np.pi
        ind = pos[:,0] < 0
        ra[ind] += np.pi

        # get z from comoving distance
        r_com = np.sqrt(np.sum(pos**2, axis=1))
        z = self.redshift(r_com)
        
        # get dec
        dec = (np.pi/2) - np.arccos(pos[:,2] / r_com)

        # convert radians to degrees
        ra *= 180 / np.pi
        dec *= 180 / np.pi

        return ra, dec, z

    
    def vel_to_vlos(self, pos, vel):
        """
        Projects velocity vector along line of sight, with observer positioned at the origin

        Args:
            pos: 2d array of comoving position vectors [Mpc/h]
            vel: 2d array of velocity vectors [km/s]
        Returns:
            array of line of sight velocity [km/s]
        """
        # comoving distance to each object
        distance = np.sum(pos**2, axis=1)**0.5

        # normalize postion vectors
        pos_norm = pos.copy()
        for i in range(3):
            pos_norm[:,i] = pos_norm[:,i] / distance

        # project velocity along position vectors
        v_los = np.sum(pos_norm*vel,axis=1)

        return v_los
    
    
    def zcos_to_zobs(self, z_cos, v_los):
        """
        Adds the effect of the los velocity to the cosmological redshift

        Args:
            z_cos: array of cosmological redshift
            v_los: array of line of sight velocity [km/s]
        Returns:
            array of observed redshift
        """
        z_obs = ((1 + z_cos) * (1 + v_los/3e5)) - 1.
        return z_obs

    
    def zobs_to_vel(self, z_cos, z_obs):
        """
        Convert line of sight velocity to observed redshift

        Args:
            z_cos: array of cosmological redshift
            z_obs: array of observed redshift
        Returns:
            array of line of sight velocity [km/s]
        """
        v_los = 3e5 * ((1. + z_obs)/(1. + z_cos) - 1)
        return v_los