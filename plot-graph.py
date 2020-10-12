import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table


# The puropose of this program is to plot data from fortran scripts
# This is a prototype which will be perfectioned over time

data_inicial = Table.read("CI-advection.dat", format="ascii")
data_final = Table.read("advection-final.dat", format="ascii")
data_middle = Table.read("advection-middle.dat", format="ascii")
#Arrays with initial conditions
x_CI = data_inicial["x"]
U_CI = data_inicial["U(x,0)"]
 
#Arrays with middle state
x_mid = data_middle["x"]
U_mid = data_middle["U(x,t)"]

#Arrays with the final state
x_f = data_final["x"]
U_f = data_final["U(x,tlim)"]

#Plot everything
plt.plot(x_CI, U_CI, label="Initial Conditions (t=0)")
plt.plot(x_mid, U_mid, label=r"Intermediate state")
plt.plot(x_f, U_f, label="Final state (t=1.5)")

plt.legend()
plt.xlabel("Length x")
plt.ylabel("Conserved Quantity U")
plt.savefig("advection.pdf")
