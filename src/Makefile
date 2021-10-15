all: mts
mts: src/driver.o src/mtserr.o src/preprocessing.o src/vecfunc.o src/em.o src/mtsio.o
g++ -lm src/driver.o src/mtserr.o src/preprocessing.o src/vecfunc.o src/em.o src/mtsio.o -o
mts
driver.o: src/driver.c src/mts.h src/preprocessing.h src/mtserr.h src/em.h
g++ -c drive.c
mtserr.o: src/mtserr.c src/mts.h src/mtserr.h
g++ -c src/mtserr.c
preprocessing.o: src/preprocessing.c src/mts.h src/preprocessing.h src/mtsio.h src/mtserr.h src/em.h
src/vecfunc.h
g++ -lm -c src/preprocessing.c
vecfunc.o: src/vecfunc.c src/vecfunc.h src/mts.h
g++ -lm -c src/vecfunc.c
em.o: src/em.c src/mts.h src/em.h src/mtsio.h
g++ -lm -c src/em.c
mtsio.o: src/mtsio.c src/mts.h src/mtsio.h src/mtserr.h src/vecfunc.c
g++ -c src/mtsio.c
clean:
rm src/*.o
### Configuration file for MTS ###
### The '#' symbol declare the proceeding line as a comment ###
### To run MTS, first enter 'make' in the same directory of makefile ###
### Then enter command - './mts config.tst' ###
### Input File - Required (output has same name but with .em extension) ###
input input/test.txt
### Voxel Sizes - Millimeters - Optional (default is 100) ###
x_voxel_size 50
y_voxel_size 50
z_voxel_size 50
#### Em - 3 Options (average, median, 3D) - All may be run concurrently ###
em average
em median
#em 3D
#### Bin Info for Median EM (default number of bins are 100 and default size per bin is 10000) ###
bins 200
bin_size 100000
### Iterations for EM - Optional (default is 100) ###
iterations 200
### Iterations for EM - Optional (default is 100) ###
#online 100000
### Initial lambda value for EM in units of mrad^2/cm - Optional (default is 0.1) ###
### If STD option is chosen, intial lambda is based off of poca results ###
lambda 0.1
### Milliradians option uses milliradians instead of radians in EM calculations ###
#milliradians
### Nominal option lets user set the nomimal momentum used in EM (default is 3) ###
#nominal 3
### Weighted EM - Optional selection will turn on weighting voxels in EM based on POCA (default
is turned off) ###
### Option 1 is pure POCA weight; i.e. voxel where POCA is found is weighted 1, all others 0###
### Option 2 is linear POCA weight; weight is determined by the formula (n - |c - p|) /n, n = voxels
before or after poca voxel, c = current voxel, p = poca voxel ###
#weight 2
### Output file - Optional ###
### Contains all information from EM preprocessing ###
#output pleasework.txt
### Standard deviation - Optional (default is off) - Can provide filename to output to separate file
(Used as intial value for EM lambda regardless) ###
#std test.std
### Angle distribution - Optional (default is off) - Can provide filename to output to separate file
###
#dist EMDIST.txt
### POCA - Optional (default is off) - Can provide filename to output to separate file ###
#poca VertSmallNew.txt
### Precise L & T - Use these options to estimate the L and T for every voxel rather than for every
event ###
precise_l
precise_t
#### Units - the metric units of length to run EM in relative to mm (i.e. 1m is 1000mm, so for meter
units use 1000 as the parameter) - default value is 10 (i.e. cm)
#units 10