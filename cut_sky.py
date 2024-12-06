import numpy as np
from coordinates import Coordinates


def cut_sky(position, velocity, cosmology, Lbox=2000, zmin=None, zmax=None, replication=(0,0,0)):
    """
    Creates a cut sky mock by converting the cartesian coordiantes of a cubic box mock to ra, dec, z
    
    Args:
        position:  array of comoving position vectors (Mpc/h), in the range -Lbox/2 < pos < Lbox/2
        velocity:  array of proper velocity vectors (km/s)
        cosmology: instance of astropy.cosmology class
        [Lbox]:    comoving box length of simulation (Mpc/h), default is 2000 Mpc/h
        [zmin]:    If provided, will only return galaxies with z>=zmin.
                   By default will return all galaxies.
        [zmax]:    If provided, will only return galaxies with z<=zmax. 
                   By default will return all galaxies.
        [replication]: tuple indicating which periodic replication to use. 
                   Default value is (0,0,0) (no replications).
    """
    
    # array of the index of each galaxy in the original cube
    index = np.arange(position.shape[0])
    
    coords = Coordinates(cosmology)
    
    position_rep = position.copy()
    if replication==(0,0,0):
        print("No periodic replications")
    else:
        print("Applying periodic replications")
        for i in range(3):
            print("%.1f < %s < %.1f"%((-1+2*replication[i])*Lbox/2., chr(120+i), (1+2*replication[i])*Lbox/2.))
            position_rep[:,i] += Lbox*replication[i]
    
    ra, dec, zcos = coords.pos3d_to_equitorial(position_rep)
    print(len(ra))
    vlos = coords.vel_to_vlos(position_rep, velocity)
    zobs = coords.zcos_to_zobs(zcos, vlos)
    print(min(zobs), max(zobs))
    
    if not zmin is None:
        print("Applying redshift cut z > %.1f"%zmin)
        keep = zobs >= zmin
        print(sum(keep))
        ra, dec, zcos, zobs, index = ra[keep], dec[keep], zcos[keep], zobs[keep], index[keep]

    print(len(ra))  
        
    if not zmax is None:
        print("Applying redshift cut z < %.1f"%zmax)
        keep = zobs <= zmax
        ra, dec, zcos, zobs, index = ra[keep], dec[keep], zcos[keep], zobs[keep], index[keep]
                      
        
    return ra, dec, zcos, zobs, index


