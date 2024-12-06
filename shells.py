import numpy as np
import h5py
from cut_sky import cut_sky
from astropy.cosmology import FlatLambdaCDM
import pandas as pd

from coordinates import Coordinates
cosmology = FlatLambdaCDM(H0=67.74, Om0=0.3089, Ob0=0.0486)
coords = Coordinates(cosmology)


def make_shell(cosmology, Lbox=2000, zmin=0, zmax=2):
    """
    Converts a cubic box mock into a cut-sky mock, covering the full sky
    
    Args:
        input_path:  path of the cubic box SHAM mocks
        output_path: path where to store output files
        zsnap:     redshift of the snapshot
        cosmology: instance of astropy.cosmology class
        [Lbox]:    comoving box length of simulation (Mpc/h), default is 2000 Mpc/h
        [zmin]:    If provided, will only return galaxies with z>=zmin. Default is 0
        [zmax]:    If provided, will only return galaxies with z<=zmax. Default is 2
    """
    
    box_size = Lbox
    # convert zmin and zmax of the shell into comoving distances
    rmin_shell, rmax_shell = cosmology.comoving_distance(np.array([zmin, zmax])).value * cosmology.h
    column_names = ['id', 'pid', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mvir_all', 'Vpeak0', 'Vpeak', 'galmass']
    df = pd.read_csv('/gals_boxes/Uchuu_galaxies_complete_34.txt', delim_whitespace=True, names=column_names, dtype={'galmass': np.float64})
    pos = np.column_stack((df['x'].values, df['y'].values, df['z'].values))
    vel = np.column_stack((df['vx'].values, df['vy'].values, df['vz'].values))

    for j in range(3):
        # shift positions so they are in the range -Lbox/2 < x < +Lbox/2
        pos[:,j] =pos[:,j] - Lbox/2.

    nn=0
    # loop through periodic replications of the box in the x, y and z-directions
    for i in range(-5,6):
        for j in range(-5,6):
            for k in range(-5,6):
                rep = (i,j,k)
                print(rep)
                
                # check if corner of box nearest the observer is beyond the max z of the shell
                # if this is true, we can skip this replication
                X = np.clip(abs(i)-0.5, 0, None) 
                Y = np.clip(abs(j)-0.5, 0, None) 
                Z = np.clip(abs(k)-0.5, 0, None) 
                rmin_rep = (X**2 + Y**2 + Z**2)**0.5 * box_size
                
                if rmin_rep>rmax_shell: continue 
                
                # check if corner of box furthest from the observer is within the min z of the shell
                # if this is true, we can skip this replication
                X = abs(i)+0.5
                Y = abs(j)+0.5
                Z = abs(k)+0.5
                rmax_rep = (X**2 + Y**2 + Z**2)**0.5 * box_size
                
                if rmax_rep<rmin_shell: continue 

                # convert the cartesian coordinates to ra, dec and z
                # zcos is cosmological redshift
                # zobs is observed redshift (including velocities)
                ra, dec, zcos, zobs, index = \
                    cut_sky(pos, vel, cosmology, Lbox=Lbox, zmin=zmin, zmax=zmax, replication=rep)
                # index gives the index in the original cubic box
                # can use this if there are other properties in the mock you want to keep
                
                print("NGAL:", np.count_nonzero(ra))


                id = df['id'].values[index]
                pid = df['pid'].values[index]
                x = df['x'].values[index]
                y = df['y'].values[index]
                z = df['z'].values[index]
                vx = df['vx'].values[index]
                vy = df['vy'].values[index]
                vz = df['vz'].values[index]
                mass = df['mvir_all'].values[index]
                vpeak = df['Vpeak0'].values[index]
                vpeak_s = df['Vpeak'].values[index]
                galmass = df['galmass'].values[index]
                out = np.column_stack((id, pid, x, y, z, vx, vy, vz, mass, vpeak, vpeak_s, galmass, ra, dec, zcos, zobs))

                np.savetxt("shells/z0p94/cut_sky_%i_%i_%i.txt"%(i,j,k), out)
                
                # these files will then all need to be joined together afterwards
        
        


# Uchuu cosmology
cosmology = FlatLambdaCDM(H0=67.74, Om0=0.3089, Ob0=0.0486)

Lbox = 2000. #box size in Mpc/h

zmin = 0.9 # min redshift of shell
zmax = 0.95 # max redshift of shell


make_shell(cosmology, Lbox=2000, zmin=zmin, zmax=zmax)
