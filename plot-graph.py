
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import argparse

# The puropose of this program is to plot data from fortran scripts
# This is a prototype which will be perfectioned over time

parser = argparse.ArgumentParser(description="""choose which simulation you want to plot""")

parser.add_argument("--data", type=str, default="advection", choices=["advection", "burguers"],
                    help="choose the simulation you want to plot")

cmd_args = parser.parse_args()
filename_start = cmd_args.data

data_inicial = Table.read("CI-"+filename_start+".dat", format="ascii")
data_final = Table.read(filename_start+"-final.dat", format="ascii")
data_middle = Table.read(filename_start+"-middle.dat", format="ascii")
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
plt.savefig(filename_start+".pdf")
