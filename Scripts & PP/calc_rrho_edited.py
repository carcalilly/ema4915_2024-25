#import modules and functions 
import numpy  as np
import sys
import matplotlib.pyplot as plt
from   scipy.interpolate import interp1d, make_interp_spline   
from   scipy.integrate   import simpson, trapezoid
from   scipy.signal import savgol_filter
from   scipy.ndimage import gaussian_filter1d

"""
 @author Liping Yu
 repurposed for ema4915 2024/2025
"""

data = open(sys.argv[1],'r').readlines() #opens first argument in command line and stores it in data
fout = open(sys.argv[2],'w') #opens second argument in command line and stores it in fout



"""
read 3-5 data lines (in data[2-4])
    represents lattice vectors in the CHG file 
splits into components (.split()) (for the O example, data[2].split() = ['8.0', '0.0', '0.0'])
designates as float (float(x))

squares (**2)
takes the sum of squares (np.sum)
takes the square root of the sum (**0.5), this finds the norm/magnitude of the vector

multiplies by a (a*) (scaling factor, line 2 of the CHG file)

designates as a1-3
    will equal the lattice parameters a, b, c
calculates volume (vol)
    of the unit cell
"""
a = float(data[1].split()[0])
a1=a*(np.sum([float(x)**2 for x in data[2].split()]))**0.5 
a2=a*(np.sum([float(x)**2 for x in data[3].split()]))**0.5
a3=a*(np.sum([float(x)**2 for x in data[4].split()]))**0.5
vol = a1*a2*a3
print(f"Lattice Parameters: a1={a1:.2f}, a2={a2:.2f}, a3={a3:.2f}")
print(f"Unit Cell Volume: vol={vol:.2f}")

"""
 read string from data[6] and split into components (data[6].split())
    represents line 7 of the CHG file (number of atoms)
 elements designated as integers (int(x))
    sum of elements (np.sum)
    designated as nat (for our purposes, nat = 1 always)

data[nat+9]= data[10] (line 11 of the CHG file, which in O example: data[10].split() → ["80", "80", "80"]
    writes into nx, ny, nz
    represents the number of grid points in each direction (x, y, z)

calculate the spacing in each direction (dx, dy, dz) 
    via division of lattice parameter (a1-3) by number of grid points in each direction (nx-z)
    volume of a grid cell (dvol) = product of the three spacings
"""

nat      = np.sum([int(x) for x in data[6].split()])
nx,ny,nz = [ int(x) for x in data[nat+9].split() ]
dx = a1/nx
dy = a2/ny
dz = a3/nz
dvol = dx*dy*dz
print(f"Number of atoms: nat={nat}")
print(f"Grid Points: nx={nx}, ny={ny}, nz={nz}")
print(f"Spacing: dx={dx}, dy={dy}, dz={dz}")
print(f"Grid Cell Volume: dvol={dvol}")

"""
 read charge data
"""
"""
read charge data from data[nat+11] (line 12 in CHG) 
    calculate number of columns (ncolums) using len and split (split to components, len to count)
        O example has 10 columns of data in CHG starting from line 12
grid_lines = number of grid points (nx*ny*nz) divided by number of columns (ncolums)
    if the division is not an integer, grid_lines is rounded up
    should output approximate number of lines of data in the CHG file
"""
ncolums    = len(data[nat+11].split())
grid_lines = int(nx*ny*nz/ncolums)
if float(nx*ny*nz)/ncolums > grid_lines :
   grid_lines += 1

"""
 calc lines to store atoms, ends up being 0 all the time since nat=1?
"""
nat_lines = int(nat/ncolums)
if float(nat)/ncolums > nat_lines : 
   nat_lines  += 1

"""
declare data_1D as an empty list
for loop to read data from data[nat+10] (line 11) to data[nat+10+grid_lines] (to about the end of the file)
    split and designate as x 
append adds x to data_1D list as a float 
reshape data_1D into a numpy 3D array (data_3D) with dimensions nx, ny, nz
    order='F' is for Fortran-style indexing
    divide by volume (vol) to normalize
    this is the charge density in the unit cell
"""
data_1D = []
for l in range(nat+10, nat+10+grid_lines) : 
    for x in data[l].split() :              
        data_1D.append(float(x))
data_3D = np.reshape(np.array(data_1D),(nx,ny,nz),order='F')/vol
  
"""
total charge density (sum of all r in data_3d) multiplied by volume (dvol)
    returns total charge
minimum value (np.min) 
"""
print("# total charge (all r's) and rho_min: ", np.sum(data_3D)*dvol, np.min(data_3D))


"""
calculate radial charge density rrho, and Rs
    multiply #of grid points in each direction (nx*ny*nz) to get total grid points = N
    create two arrays of zeros with N elements, Rs_all and rrho_all
"""
N = nx*ny*nz 
Rs_all   = np.zeros(N) 
rrho_all = np.zeros(N) 
print("Total # of grid points: ", N)

"""
integer division of array (0 to N) by (ny*nz) to get xinds (x grid indices)
int. division of array by nz modulo (remainder) ny to get yinds (y grid indices)
remainder of division by nz to get zinds (z grid indices)
"""

xinds = np.arange(N)//(ny*nz) 
yinds = np.arange(N)//nz%ny
zinds = np.arange(N)%nz

"""
symmetry shift of grid indices
    if xinds > nx-xinds, xinds = nx-xinds
    if yinds > ny-yinds, yinds = ny-yinds
    if zinds > nz-zinds, zinds = nz-zinds
square
"""

xinds = np.where(xinds > nx-xinds, nx-xinds, xinds)**2
yinds = np.where(yinds > ny-yinds, ny-yinds, yinds)**2
zinds = np.where(zinds > nz-zinds, nz-zinds, zinds)**2

"""
calculate real space distance Rs_all
    transpose of array of xinds, yinds, zinds (stacked vertically) to get Nx3 array
        rows are grid points, columns are x, y, z indices
    array of dx, dy, dz (spacing in each direction) squared
    dot product of the two arrays (element-wise multiplication and sum) and square root
results in real space distance from the origin to each grid point
reshape data_3D into 1D array (rrho_all) with N elements
    charge density at each grid point
"""

Rs_all   = np.dot( np.transpose([xinds,yinds,zinds]), np.array([dx,dy,dz])**2 )**0.5
rrho_all = data_3D.reshape(N)

"""
sort Rs_all in ascending order
pairs Rs_all and rrho_all values after sorting
"""

ind      = np.argsort(Rs_all)
Rs_all   = Rs_all[ind]
rrho_all = rrho_all[ind]

"""
removes charge density values that exceed eps
    reverse rrho_all (rrho_all[::-1]) (now in descending order)
    cumulative sum of rrho_all (np.cumsum) 
    multiply by dvol (volume of a grid cell) to get charge
    true if charge is greater than eps, false if less than eps
    count number of true values (values above eps) (len)
    i_cut is cutoff index
"""

eps      = 1e-5    
i_cut    = len(Rs_all[np.cumsum(rrho_all[::-1])*dvol > eps ])

"""
truncate Rs_all and rrho_all to the index specified by i_cut
    Rs_red and rrho_red are the truncated arrays
    r_cut is the last value of Rs_red
"""

Rs_red   = Rs_all[:i_cut]
rrho_red = rrho_all[:i_cut]
r_cut    = Rs_red[-1]

"""
find total charge in all grid points and truncated grid points
    sum of charge density values (rrho_all and rrho_red) multiplied by dvol (volume of a grid cell)
"""

Qtot_all = np.sum(rrho_all)*dvol
Qtot_red = np.sum(rrho_red)*dvol

print("Residual charge/points skipped/r_cut:", Qtot_all-Qtot_red, N-len(Rs_red), r_cut)


"""
 grouping the grids at the same R's
"""
counter       = np.arange(1, len(Rs_red))
Rs_splitted   = np.split(Rs_red,   counter[Rs_red[1:]!=Rs_red[:-1]])
rrho_splitted = np.split(rrho_red, counter[Rs_red[1:]!=Rs_red[:-1]])

Rs_red    = np.array([ Rs_splitted[i][0]         for i in range(len(Rs_splitted))])
rrho_red  = np.array([ np.mean(rrho_splitted[i]) for i in range(len(Rs_splitted))])
Qtot_red  = trapezoid(rrho_red*4*np.pi*Rs_red**2, Rs_red)


# rrho interpolation over a regular mesh
f = interp1d(Rs_red, rrho_red, kind='slinear')

step  = 1e-2
Rs_reg   = np.arange(0,r_cut,step)
rrho_reg = f(Rs_reg)
Drho_reg = [rrho_reg[i]*np.pi*4.*Rs_reg[i]**2 for i in range(len(Rs_reg))]
Qtot_reg = simpson(Drho_reg, Rs_reg)

"""
Gaussian smoothing 
"""
sigma = (dvol*3/(4*np.pi))**(1/3.)*0.6

rrho_reg_sm = np.zeros(len(Rs_reg))

for i in range(len(Rs_reg)) : 
    delta_x     = Rs_reg - Rs_reg[i] 
    weights     = np.exp(-0.5*((delta_x /sigma)**2))
    rrho_reg_sm[i] = np.sum(weights*rrho_reg)/np.sum(weights)

    if i%100 == 0 : print("...Calculating rrho at regular R point: ", i, Rs_reg[i], rrho_reg[i])


Drho_reg_sm = 4*np.pi*Rs_reg**2 *rrho_reg_sm
Qtot_reg_sm = simpson(Drho_reg_sm, Rs_reg)

"""
normalized to the right total charge
"""
rrho = rrho_reg_sm/Qtot_reg_sm*round(Qtot_all)
Drho = rrho*4.*np.pi*Rs_reg**2 
Qtot_normalized = simpson(Drho, Rs_reg)

print("# total charge (all r's | regular r's | regular smoothed | normalized): ", Qtot_all, Qtot_reg, Qtot_reg_sm, Qtot_normalized)

"""
write data
"""

# write data for ML
zsigma = 0.1
trho =  round(Qtot_all) * np.exp(-0.5*(Rs_reg/zsigma)**2)/zsigma/(0.5*np.pi)**0.5 - Drho

fout.write("# total charge: %f \n" %(simpson(Drho, Rs_reg)))
fout.write("# r(A)     rho (e/A^3)     4*pi*r^2*rho(e/A)  dip  \n")
for i in range(len(Rs_reg)) :
    fout.write("%7.4f  %15.7E  %15.7E  %15.7E\n" %(Rs_reg[i], rrho[i], Drho[i], trho[i]))

fout.close()

# write source data for debugging...
sout = open('ir_rrho.dat','w') 

sout.write("# total charge (all r's | reg r's | normalized): %f %f %f" %(Qtot_all, Qtot_reg, Qtot_normalized))
sout.write("# r(A)     rho (e/A^3)\n")
sout.write("# ..........................@ all r's...................... \n")
for i in range(len(Rs_all)) :
    sout.write("%7.4f  %15.7E\n" %(Rs_all[i], rrho_all[i]))

sout.write("\n")
sout.write("# ..........................intepolatd @ even-spaced r's.............. \n")
for i in range(len(Rs_reg)) :
    sout.write("%7.4f  %15.7E\n" %(Rs_reg[i], rrho_reg[i]))

sout.write("\n")
sout.write("# ..........................smoothed @  even-spaced r's.............. \n")
for i in range(len(Rs_reg)) :
    sout.write("%7.4f  %15.7E\n" %(Rs_reg[i], rrho_reg_sm[i]))

sout.write("\n")
sout.write("# ..........................normalized @  even-spaced r's.............. \n")
for i in range(len(Rs_reg)) :
    sout.write("%7.4f  %15.7E\n" %(Rs_reg[i], rrho[i]))

sout.close()

# plot
plt.scatter(Rs_all, rrho_all)
plt.plot(Rs_reg, rrho,'r')
plt.plot(Rs_reg, Drho,'b')

plt.xlim(0,4)
plt.xlabel('r (angstrom)')
plt.ylabel('radial charge density/distribution')
plt.show() 