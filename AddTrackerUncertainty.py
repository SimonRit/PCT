#!/usr/bin/env python

import itk
import hepunits
import numpy as np
import click

@click.command()
@click.option("xOverX0", "--xOverX0", default=5e-3, help="Material budget (unitless)")
@click.option("--sp", default=0.15, help="Standard deviation of the tracker uncertainty (mm)")
@click.option("--dt", default=10, help="Distance between trackers (cm)")
@click.option("-i", "--input", required=True, help="Input file name")
@click.option("-o", "--output", required=True, help="Output file name")
@click.option("entryTranslation", "--entryTranslation", default=0., help="Amount of translation of the entry detector position in the direction perpendicular to the detector (mm)")
@click.option("exitTranslation", "--exitTranslation", default=0., help="Amount of translation of the exit detector position in the direction perpendicular to the detector (mm)")

def AddTrackerUncertainty(xOverX0, sp, dt, input, output, entryTranslation, exitTranslation):
    dt = dt * hepunits.cm
    sp = sp * hepunits.mm

    # Equation (25) and (26) from [Krah et al, PMB, 2018]
    T = np.zeros((2,2))
    T[0,1] = 1
    T[1,0] = -1/dt
    T[1,1] = 1/dt
    
    # Equation (27) and (28) from [Krah et al, PMB, 2018]
    def GetSigmaSc(energy, xOverX0):
        proton_mass_c2 = 938.272013 * hepunits.MeV
        betap = (energy + 2*proton_mass_c2)*energy / (energy + proton_mass_c2)
        sigmaSc = 13.6 * hepunits.MeV / betap * np.sqrt(xOverX0) * (1+0.038*np.log(xOverX0))
        SigmaSc = np.zeros((energy.size,2,2))
        SigmaSc[:,1,1] = sigmaSc**2
        return np.tile(sp**2*T@T.T,(energy.size,1,1)) + SigmaSc

    pairs = itk.imread(input)
    pairs = itk.GetArrayFromImage(pairs)

    # Move to entrance and exit detector (new) positions
    pairs[:,0,:] += ((entryTranslation / pairs[:,2,2]) * pairs[:,2,:].T).T
    pairs[:,1,:] += ((exitTranslation / pairs[:,3,2]) * pairs[:,3,:].T).T

    eEntry = pairs[:,4,0]
    eExit = pairs[:,4,1]
    SigmaEntry = GetSigmaSc(eEntry, xOverX0)
    SigmaExit = GetSigmaSc(eExit, xOverX0)
    wEntry, Qentry = np.linalg.eig(np.linalg.inv(SigmaEntry))
    wExit, Qexit = np.linalg.eig(np.linalg.inv(SigmaExit))
    xrEntry = np.random.randn(eEntry.size, 2, 2)
    xrExit = np.random.randn(eExit.size, 2, 2)
    #Wentry = np.diag(1./np.sqrt(wEntry[i,:]))
    Wentry = np.zeros((eEntry.size, 2, 2))
    Wentry[:,0,0] = 1./np.sqrt(wEntry[:,0])
    Wentry[:,1,1] = 1./np.sqrt(wEntry[:,1])
    #Wexit = np.diag(1./np.sqrt(wExit[i,:]))
    Wexit = np.zeros((eExit.size, 2, 2))
    Wexit[:,0,0] = 1./np.sqrt(wExit[:,0])
    Wexit[:,1,1] = 1./np.sqrt(wExit[:,1])
    #dYuncertEntry = Qentry[i,:,:].dot(Wentry).dot(xrEntry[i,:,:]).T
    dYuncertEntry = np.matmul(np.matmul(Qentry,Wentry), xrEntry)
    #dYuncertExit = Qexit[i,:,:].dot(Wexit).dot(xrExit[i,:,:]).T
    dYuncertExit = np.matmul(np.matmul(Qexit,Wexit), xrExit)
    pairs[:,0,0:2] += dYuncertEntry[:,0,:] # Entrance position X/Y
    pairs[:,2,0:2] += dYuncertEntry[:,1,:] # Entrance direction X/Y
    pairs[:,1,0:2] += dYuncertExit[:,0,:] # Exit position X/Y
    pairs[:,3,0:2] += dYuncertExit[:,1,:] # Exit direction X/Y
    
    itk.imwrite(itk.GetImageFromArray(pairs, is_vector=True), output)

if __name__ == '__main__':
    AddTrackerUncertainty()

