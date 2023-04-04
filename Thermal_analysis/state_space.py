import numpy as np

# Constants
k1 = 1.6   # Conductivity of exterior wall
k2 = 1.6   # Conductivity of middle wall
k3 = 1.6   # Conductivity of interior wall
l1 = 0.001 # Thickness of exterior wall
l2 = 0.5   # Thickness of middle wall
l3 = 0.001 # Thickness of interior wall

# Density
rho1 = 2000 # Density of exterior wall
rho2 = 2000 # Density of middle wall
rho3 = 2000 # Density of interior wall
rho4 = 1.2  # Density of air

# Specific heat, J / kg deg C
cp1 = 880 # Specific heat of exterior wall
cp2 = 880 # Specific heat of middle wall
cp3 = 880 # Specific heat of interior wall
cp4 = 1000 # Specific heat of air

# Area size
a = 10 # Area of wall

# Heat exchange coefficients
ho = 64 * 3.6 # Radiative heat transfer coefficient between exterior wall and space
hi = 11 * 3.6 # Convective heat transfer coefficient between interior wall and air

# Volume
V1 = a * l1 # Volume of exterior wall
V2 = a * l2 # Volume of middle wall
V3 = a * l3 # Volume of interior wall
V4 = 27     # Volume of air

# Emissivity of exterior wall
emissivity = 0.94

# Stefan-Boltzmann constant
stef_coeff = 5.67 * 10**(-8)

Tspace = 3
# Solar radiation
qsol = 1360 * a

# Internal heat source
qhvac = 0

# Time variables
dt = 3600   # Time step, in seconds
timestep = 1
hours = 48

# State space matrices
C = np.array([[rho1 * cp1 * V1, 0, 0, 0],
              [0, rho2 * cp2 * V2, 0, 0],
              [0, 0, rho3 * cp3 * V3, 0],
              [0, 0, 0, rho4 * cp4 * V4]])

A = np.array([[-2*k1/l1-k2/l2, k2/l2,           0,       0],
              [k2/l2,         -2*k2/l2-k3/l3, k3/l3,   0],
              [0,             k3/l3,          -ho/hi-k3/l3, ho/hi],
              [0,             0,              ho/hi,   -ho/hi]])

B = np.array([[a*emissivity*stef_coeff, 1, 0],
              [0, 0, 0],
              [0, 0, 0],
              [0, 0, 1]])

U = np.array([[Tspace**4], [qsol], [qhvac]])

T0 = np.array([[273], [273], [273], [273]])

# Simulation loop
for i in range(int(hours*dt/timestep)):
    Tn = np.linalg.inv(np.eye(4, dtype=int) - dt*np.linalg.inv(C) @ A) @ (T0 + dt*np.linalg.inv(C) @ B @ U)
    T0 = Tn
    print(Tn)

#
# figure1 = figure;
#
# axes1 = axes('Parent', figure1);
# hold(axes1, 'on');
#
# plot1 = plot(Temp, 'Parent', axes1);
# set(plot1(1), 'DisplayName', 'Texternal surface', 'Color', 'blue');
# set(plot1(2), 'DisplayName', 'T2', 'Visible', 'on');
# set(plot1(3), 'DisplayName', 'Tinternal surface', 'Color', 'green');
# set(plot1(4), 'DisplayName', 'Tindoor', 'Color', 'r', 'Linewidth', 1);
#
# xlim(axes1, [0 n - 1]);
# ylim(axes1, [-20 40]);
#
# box(axes1, 'on');
# set(axes1, 'XGrid', 'on', 'YGrid', 'on');
#
# legend(axes1, 'show');
