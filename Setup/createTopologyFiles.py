import  os

def writeParmedInputFile(parmedInputFileName, outputTopologyName, startTopology, ResidueID_Proline, newHeightOneFold, newHeightTwoFold):
    '''
    Write input files for ParmEd:
    ParmEd creates the topology files used for Hamiltonian Replica Exchange simulations
    - startTopology: AMBER topology file with original forcefield parameters
    
    '''
    with open(parmedInputFileName, 'w') as fil:


        term1 = f':{ResidueID_Proline -1}@O :{ResidueID_Proline -1}@C :{ResidueID_Proline}@N :{ResidueID_Proline}@CD'
        term2 = f':{ResidueID_Proline -1}@O :{ResidueID_Proline -1}@C :{ResidueID_Proline}@N :{ResidueID_Proline}@CA'
        term3 = f':{ResidueID_Proline -1}@CA :{ResidueID_Proline -1}@C :{ResidueID_Proline}@N :{ResidueID_Proline}@CD'
        term4 = f':{ResidueID_Proline -1}@CA :{ResidueID_Proline -1}@C :{ResidueID_Proline}@N :{ResidueID_Proline}@CA'
        parmedCommands = (f'parm {startTopology} \n'
                          f'deleteDihedral {term1} \n'
                          f'deleteDihedral {term2} \n'
                          f'deleteDihedral {term3} \n'
                          f'deleteDihedral {term4} \n'
                          f'\n'
                          f'addDihedral {term1} {newHeightTwoFold} 2.0000 180.0001 1.2000 2.0000\n'
                          f'addDihedral {term2} {newHeightTwoFold} 2.0000 180.0001 1.2000 2.0000\n'
                          f'addDihedral {term3} {newHeightTwoFold} 2.0000 180.0001 1.2000 2.0000\n'
                          f'addDihedral {term4} {newHeightTwoFold} 2.0000 180.0001 1.2000 2.0000\n'
                          f'\n'
                          f'addDihedral {term4} {newHeightOneFold} 1.0000 180.001 1.2000 2.0000\n' # destabilize trans state. for cis destabilization use: 1.0000 0.001 1.2000 2.0000
                          f'\n'
                          f'HMassRepartition \n' #optional
                          f'parmout {outputTopologyName} \n'
                          f'go \n'
                          f'quit \n'
                          f'\n')

        fil.write(parmedCommands)


######
#  MAIN
######

# exemplary setup as used in publication, customise to your needs:

replicaNum = 12

effectiveBarrierHeightTwoFold = [0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 15.0, 18.0, 20.0] #in kcal/mol
barrierParmedTwoFold = [effHeight / 8.0 for effHeight in effectiveBarrierHeightTwoFold] # explanation factor 8.0: 4 contributing terms and V(x) = k[1+cos(x+ pi/2)], details in publication

barrierParmedOneFold = [1 for i in range(replicaNum)] #amounts to 2kcal/mol penalty on trans state for all replicas

#execute barrierheight adjustment for dummy system
i=0
for barrier1F, barrier2F in zip(barrierParmedOneFold, reversed(barrierParmedTwoFold)):
    writeParmedInputFile(f'parm{i}.in', f'rep{i}.top', 'start.prmtop', 3, barrier1F, barrier2F)
    os.system(f'parmed -i parm{i}.in')
    i += 1