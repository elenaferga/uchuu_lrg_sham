import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.interpolate import interp2d
from coordinates import Coordinates
import math
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
# Inicializamos el cosmology
cosmology = FlatLambdaCDM(H0=67.74, Om0=0.3089, Ob0=0.0486)
coords = Coordinates(cosmology)
cte = 4 * math.pi / 3

# Cargar el archivo de datos
column_names = ['id', 'pid', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mvir_all', 'Vpeak0', 'Vpeak', 'galmass', 'ra', 'dec', 'zcos', 'zobs']
df = pd.read_csv('/home/users/dae/efernandez/SHAM/LRGs/shells/cut_sky.txt', delim_whitespace=True, names=column_names)

mass = df['galmass']

#redshift bins we use to do random downsampling
zrange = np.arange(0.4, 1.1+0.025, 0.025)
zrange = np.round(zrange, 3)

logMstll_range = np.linspace(9.5, max(mass), 200)

# We read the SMF from files and append it in an array
smf_values = []
for i in range(len(zrange) - 1):
    logMstll, smf, smf_err = np.loadtxt(f'smf_lrg_y1_z{zrange[i]}_{zrange[i+1]}.txt', delimiter=' ', unpack=True)
    smf_interp = interp1d(logMstll, smf, kind='linear', fill_value='extrapolate')
    smf_values.append(smf_interp(logMstll_range))

smf_values = np.array(smf_values)

# we interpolate smf in each redshift bin (and stellar mass bin=
smf_interp_2d = interp2d(logMstll_range, (zrange[:-1] + zrange[1:]) / 2, smf_values, kind='linear')


df2 = pd.DataFrame()
#we use narrower redshift bins once we have the interpolator!
zrange = np.arange(0.4, 1.1+0.01, 0.01)

for i in range(len(zrange) - 1):
    print(zrange[i])
    keep = (df['zobs'] >= zrange[i]) & (df['zobs'] < zrange[i + 1])
    dfaux = df[keep]

    rmin = coords.comoving_distance(zrange[i])
    rmax = coords.comoving_distance(zrange[i + 1])
    vol = cte * (rmax**3 - rmin**3) / (0.6774**3)

    binwidth = logMstll_range[1] - logMstll_range[0]

    for j in range(len(logMstll_range)):
        mask = (dfaux['galmass'] >= logMstll_range[j] - 0.5 * binwidth) & (dfaux['galmass'] < logMstll_range[j] + 0.5 * binwidth)

        # Interpolación 2D para obtener el valor de SMF interpolado
        zmean = (zrange[i] + zrange[i + 1]) / 2
        smf_value = smf_interp_2d(logMstll_range[j], zmean)

        smf_complete = len(dfaux[mask])
        excess = smf_complete - round(float(smf_value) * binwidth * vol)

        if excess > 0:
            indices_to_keep = dfaux[mask].sample(n=len(dfaux[mask]) - excess).index
            df2 = df2.append(dfaux.loc[indices_to_keep], ignore_index=True)
        elif excess < 0:
            df2 = df2.append(dfaux[mask], ignore_index=True)


dfaux = df[df['galmass']>12.4]
df2 = df2.append(dfaux, ignore_index=True)

df2.columns = column_names
out = df2.to_numpy()
np.savetxt('uchuu-desi-y3.txt', out, fmt="%i %i %.5f %.5f %.5f %.5f %.5f %.5f %.15f %.15f %.15f %.15f %.5f %.5f %.5f %.5f")

