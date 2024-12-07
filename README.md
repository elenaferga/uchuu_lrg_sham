# SHAM method for Uchuu-LRG 🌌
This repository contains all the needed codes in order to construct a Uchuu-LRG mock which reproduces DESI Y1-Y3 2pcf. 

## Files of the respository 📁

### **smf_Y1/** 
This folder contains files with the observed SMF from DESI LRG Y1 using CIGALE stellar masses estimations. They are calculated in redshift bins of width 0.01 from z=0.4 to z=1.1.

### **SMF_cumulative_{}_redshift.csv** 
{} can be *high* or *low*. It is the complete SMF we use in the SHAM method.

### **sham.py** 
This code applies the SHAM method to each box. One has to change the snapshot (file_pattern, line 31) and run it once per snapshot

### **shells.py** 
Each of the codes generated from sham.py is an input from shells.py. This code generates shells using the method from [Smith, A. et al. 2022](https://academic.oup.com/mnras/article/516/3/4529/6694100).  From each box, you generate a shell in a redshift range. 

### **concatenate.sh** 
This is a shell code that concatenates all the shells generated from shells.py in order to have a cut sky lightcone. 

### **lightcone_downsampling.py** 
This code reads the cut sky lightcone and does a random downsampling to an observed SMF. We use Y1 SMF in redshift bins of width 0.01 starting in z=0.4 until z=1.1.

Until this point, we have a full sky lightcone. We can apply DESI footprint to it and have a good version of the code to be compared with data!

