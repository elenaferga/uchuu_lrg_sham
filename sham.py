import matplotlib.pyplot as plt
import numpy as np
import h5py
import pandas as pd
from numba import njit, prange, config
from scipy import interpolate
from scipy.optimize import root_scalar

def process_file(file_path):
    with h5py.File(file_path, 'r') as f:
        return {
            'pid': f['pid'][:],
            'id': f['id'][:],
            'Mvir_all': f['Mvir_all'][:],
            'Vpeak': f['Vpeak'][:],
            'x': f['x'][:],
            'y': f['y'][:],
            'z': f['z'][:],
            'vx': f['vx'][:],
            'vy': f['vy'][:],
            'vz': f['vz'][:],
        }
    
def find_a_for_nvpeak(f, nvpeak_value):
    def objective(logMstll):
        return f(logMstll) - nvpeak_value
    result = root_scalar(objective, bracket=[minmass, maxmass], method='brentq')
    return result.root


file_pattern = "/home/users/dae/ishiyama/Uchuu/RockstarExtendedAll/RockstarExtended/halodir_034/halolist_z0p94_{}.h5"
#file_pattern = "/home/users/dae/ishiyama/Uchuu/RockstarExtendedAll/RockstarExtended/halodir_036/halolist_z0p78_{}.h5"
chunk_size = 10  # Adjust the chunk size based on your memory constraints
dfs = []

for i in range(0, 100, chunk_size):
    print(f"Processing files {i} to {i + chunk_size - 1}")
    files = [file_pattern.format(j) for j in range(i, i + chunk_size)]
    data_chunks = [process_file(file) for file in files]
    df_chunk = pd.DataFrame({
        'pid': np.concatenate([chunk['pid'] for chunk in data_chunks]),
        'id': np.concatenate([chunk['id'] for chunk in data_chunks]),
        'Vpeak': np.concatenate([chunk['Vpeak'] for chunk in data_chunks]),
        'Mvir_all': np.concatenate([chunk['Mvir_all'] for chunk in data_chunks]),
	'x': np.concatenate([chunk['x'] for chunk in data_chunks]),
	'y': np.concatenate([chunk['y'] for chunk in data_chunks]),
	'z': np.concatenate([chunk['z'] for chunk in data_chunks]),
	'vx': np.concatenate([chunk['vx'] for chunk in data_chunks]),
	'vy': np.concatenate([chunk['vy'] for chunk in data_chunks]),
	'vz': np.concatenate([chunk['vz'] for chunk in data_chunks]),
    })
    dfs.append(df_chunk)

df = pd.concat(dfs, ignore_index=True)
print(len(df))

df = df[df.Mvir_all>pow(10, 10.4)].reset_index()

vol = pow(2000, 3)/pow(0.6774,3)

minmass = 8.5
maxmass = 12.3

df['Vpeak_scatter'] = df['Vpeak']*(1+np.random.normal(loc=0, scale=0.25, size=len(df)))
print(len(df[df['Vpeak_scatter']>150]), df[df['Vpeak_scatter']>120])

df = df[df['Vpeak_scatter']>150]

print(len(df))
print('Scatter added to Vpeak')

Vpeak = df['Vpeak_scatter'].values
x = df['x'].values
y = df['y'].values
z = df['z'].values
vx = df['vx'].values
vy = df['vy'].values
vz = df['vz'].values
mvir_all= df['Mvir_all'].values
pid= df['pid'].values
id = df['id'].values
Vpeak0 = df['Vpeak'].values


bins = np.linspace(min(Vpeak), max(Vpeak), 100000)
nvpeak, bins_vpeak, patches = plt.hist(Vpeak, bins=bins, cumulative=-1)
print('Histogram of Vpeak calculated')
nvpeak = nvpeak/vol


#I will interpolate nvpeak just in case it is necessary
nvpeak_interpolation = interpolate.interp1d(bins_vpeak[:-1], nvpeak, fill_value='extrapolate')

df['nvpeak'] = nvpeak_interpolation(Vpeak)    

print('nVpeak interpolated')

logMstll, ngals = np.loadtxt('SMF_cumulative_high_tres_trozos.csv', skiprows=1, delimiter=',', unpack=True)

#I will	interpolate the complete smf just in case it is necessary
f = interpolate.interp1d(logMstll, ngals, kind='cubic', fill_value='extrapolate')

print('ngals read')

#Now we have to host haloes with galaxies!
nvpeak = df['nvpeak'].values
galmass = np.zeros(len(Vpeak))
for i in range(len(nvpeak)):
    # Interpolar el valor acumulado de nvpeak para Vpeak
    nvpeak_value = nvpeak[i]
    galmass_value = find_a_for_nvpeak(f, nvpeak_value)
    #print(i, 'of', len(nvpeak)-1, Vpeak[i], nvpeak_value, galmass_value)
    galmass[i] = galmass_value

print('Haloes hosted with galaxies')


keep = galmass>9

out = np.column_stack((id[keep], pid[keep], x[keep], y[keep], z[keep], vx[keep], vy[keep], vz[keep], mvir_all[keep], Vpeak0[keep], Vpeak[keep], galmass[keep]))
np.savetxt('/home/users/dae/efernandez/SHAM/LRGs/Uchuu-DESI_complete_Y1/tres_trozos/scatter_variable/Uchuu_galaxies_complete_34.txt', out, fmt="%i %i %.5f %.5f %.5f %.5f %.5f %.5f %.15f %.15f %.15f %.15f")


