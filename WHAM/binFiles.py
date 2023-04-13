import numpy as np 

inNames = [f'omega{i}.dat' for i in range(12)]
outNames = [f'histogram{i}.dat' for i in range(12)]

start = -180
end = 180
binNumber = 120+1
Bins = np.linspace(start,end,binNumber)
binSize = Bins[1] - Bins[0]

for inp, outp in zip(inNames, outNames):
        omegaVals = np.loadtxt(inp)# raw Data
        hist, bin_edges = np.histogram(omegaVals, bins= Bins, density =True)
        np.savetxt(outp, np.column_stack((bin_edges[:-1]+binSize/2, hist)), delimiter = ' ', fmt=('%0.1f', '%0.4f'))