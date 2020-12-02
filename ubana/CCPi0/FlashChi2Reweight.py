import ROOT, math, sys, os
import uproot
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
import re
import scipy as sci
from mpl_toolkits.mplot3d import Axes3D
from ROOT import TH1, TAxis, gROOT, TCanvas
from scipy import stats, optimize

####################################################################################################
channel_or_particle = 'channel'

def chanToHistogram(channel):
    if channel == "QE":
       return 0
    elif channel == "RES":
       return 1
    elif channel == "DIS":
       return 2
    elif channel == "MEC":
       return 3
    elif channel == "Other":
       return 4
    elif channel == "OFFV":
       return 5
    else:
       return -1

def splitAndSort(inputArray, nDivisions):
  preSplit = np.sort(inputArray)
  dropEntries = preSplit.size % nDivisions
  print "Dropping %d entries prior to split" % dropEntries
  #preSplit = preSplit[dropEntries:]

  return np.array_split(preSplit, nDivisions)

def tagDuplicateEvents(inputFrame):
  inputFrame.insert(inputFrame.shape[1], "DuplicatedEvent", inputFrame.duplicated() )

def loadMCEventInfo(inputFrame):
  tagDuplicateEvents(inputFrame)
  inputFrame.insert(inputFrame.shape[1], "mc_channel", [getChan(x, y) for x, y in zip(inputFrame['mc_nu_interaction_type'], inputFrame['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
  inputFrame.eval('mc_Ehad = mc_nu_energy - mc_nu_lepton_energy', inplace=True) #Insert the true energy transfer (nu)
  inputFrame.insert(inputFrame.shape[1], "mc_expQ2", [getQ2(x, y, z) for x, y, z in zip(inputFrame['mc_nu_energy'], inputFrame['mc_nu_lepton_energy'], inputFrame['mc_nu_lepton_theta'])] )
  inputFrame.insert(inputFrame.shape[1], "mc_expW", [getW(x, y, z) for x, y, z in zip(inputFrame['mc_Ehad'], inputFrame['mc_expQ2'], inputFrame['mc_channel'] ) ] )
  inputFrame.insert(inputFrame.shape[1], "mc_expXbj", [getXbj(x, y) for x, y in zip(inputFrame['mc_Ehad'], inputFrame['mc_expQ2'] ) ] )
  inputFrame.insert(inputFrame.shape[1], "mc_expY", [getInel(x, y) for x, y in zip(inputFrame['mc_Ehad'], inputFrame['mc_nu_energy'] ) ] )
  #inputFrame.insert(inputFrame.shape[1], "template_wgt", [getChanWeight(x, y) for x, y in zip(inputFrame['mc_nu_interaction_type'], inputFrame['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
  inputFrame.insert(inputFrame.shape[1], "pot", mcPOT)
  #inputFrame.insert(inputFrame.shape[1], "isTrueCC", [isTrueCC(x, y) for x, y in zip(inputFrame['mc_pdg'], inputFrame['mc_nu_ccnc'])])

def loadMCTrackInfo(inputFrame):
  inputFrame.insert(inputFrame.shape[1], 'particle', [getParticle(x) for x in inputFrame['mc_pdg'] ])

def loadTrackInfo(inputFrame, isMC=False):
  inputFrame.insert(inputFrame.shape[1], "DuplicatedEvent", inputFrame.duplicated())
  inputFrame.insert(inputFrame.shape[1], "phi", [getPhi(x, y) for x, y in zip(inputFrame['track_diry'], inputFrame['track_dirx'] ) ] )
  inputFrame.eval('track_chi2_ratio = track_chi2_proton / track_chi2_muon', inplace=True)
  inputFrame.insert(inputFrame.shape[1], "isContained", [isContained(x, y, z, a, b, c) for x, y, z, a, b, c in zip(inputFrame['vx'], inputFrame['vy'], inputFrame['vz'], inputFrame['track_endx'], inputFrame['track_endy'], inputFrame['track_endz'])])
  inputFrame.insert(inputFrame.shape[1], 'track_mom_best', [getBestP(x, y, z) for x, y, z in zip(inputFrame['isContained'], inputFrame['track_range_mom_mu'], inputFrame['track_mcs_mom'])])
  if(isMC):
    inputFrame.insert(inputFrame.shape[1], 'particle', [getParticle(x) for x in inputFrame['mc_pdg'] ])

def AggregateFrame(inputFrame, var, stat):
  if stat not in ["count", "max", "mean"]:
    print "Cannot aggregate based on stat %s" % stat
    return
  
  statFrame = inputFrame.groupby(level=["run", "subrun", "event"]).agg({var: [stat]})

  statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]

  inputFrame = inputFrame.join(statFrame['%s_%s' % (var, stat)], on=["run", "subrun", "event"])
  
  if(stat == "max"): 
    inputFrame.eval('isMax_%s = (%s == %s_max)' % (var, var, var), inplace=True)


def makeMCHistogram(mc, channel, binRange, nBins, filename, Titles):
  dir_name = "PlotDir"
  colors = {"QE":'b', "RES":'g', "DIS":'y', "2p2h":'r', "NC / Other":'grey', "Ext":'magenta'}

  plt.hist(mc, bins=nBins, stacked=False, range=binRange, color = colors[channel])
  plt.legend([channel])
  try:
    plotTitle, xAxisTitle, yAxisTitle = Titles
  except(ValueError):
    print "Pleast provide three titles for plot name, x and y axis names" 
    plotTitle  = ""
    xAxisTitle = ""
    yAxisTitle = "" 
  plt.title(plotTitle)
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  plt.savefig("%s/%s%s.png" % ( dir_name, filename, channel.replace(" / Other", "")) )
  plt.close()

def make2DMCHistogram(mc, channel, binRange, nBins, filename, Titles):
  dir_name = "PlotDir"
  colors = {"QE":'b', "RES":'g', "DIS":'y', "2p2h":'r', "NC / Other":'grey', "Ext":'magenta'}
  zMin = 0.01

  try:
    xBins, yBins = binRange

  except(ValueError):
    print "Please provide a range of bins for each axis"
    return

  plt.hist2d(mc[0], mc[1], bins=nBins, range=binRange, cmin=zMin, normed=True)
  plt.legend([channel])

  try:
    plotTitle, xAxisTitle, yAxisTitle = Titles
  except(ValueError):
    print "Pleast provide three titles for plot name, x and y axis names" 
    plotTitle  = ""
    xAxisTitle = ""
    yAxisTitle = "" 
  plt.title(plotTitle)
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  plt.savefig("%s/%s%s.png" % ( dir_name, filename, channel.replace(" / Other", "")) )
  plt.close()

def makeDataMCHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles, xRange = [], yRange = [], rRange = []):
  dir_name = "FittingPlots"
  xLimits  = []
  yLimits  = []
  for x in xRange:
    xLimits.append(x)
  for y in yRange:
    yLimits.append(y)  
  if(len(xLimits) == 1):
    xLimits.insert(0, 0.0)
  if(len(yLimits) == 1):
    yLimits.insert(0, 0.0)

  plt.hist(mcList, bins=nBins, stacked=True, range=binRange, color = ['b', 'g', 'y', 'r', 'grey', 'gold', 'magenta'], weights = mcWeights )
  plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Dirt', 'Ext'])

  #plotTitle, xAxisTitle, yAxisTitle =  Titles
  try:
    plotTitle, xAxisTitle, yAxisTitle = Titles
  except(ValueError):
    print "Pleast provide three titles for plot name, x and y axis names" 
    plotTitle  = ""
    xAxisTitle = ""
    yAxisTitle = "" 
  plt.title(plotTitle)
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  
  if(len(xLimits) == 2):
    plt.xlim(xLimits)
  if(len(yLimits) == 2):
    plt.ylim(yLimits)
  
  data_hist = dataify(dataList, nBins, binRange)
  plt.errorbar(data_hist[0], data_hist[1], yerr=data_hist[2], fmt='o', color='black')
  plt.savefig("%s/%s.png" % ( dir_name, filename) )
  plt.close()

  makeDataMCRatioHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles, xLimits, rRange)

def makeDataMCRatioHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles, xRange = [], yRange = []):
  dir_name  = "FittingPlots"
  xLimits  = []
  yLimits  = []
  for x in xRange:
    xLimits.append(x)
  for y in yRange:
    yLimits.append(y)  
  if(len(xLimits) == 1):
    xLimits.insert(0, 0.0)
  if(len(yLimits) == 1):
    yLimits.insert(0, -1.0)
  mcSum = np.full(nBins, 0.0 )
  for mc, weight in zip(mcList, mcWeights):
     mc_hist   = np.histogram(mc, bins=nBins, range=binRange, weights = weight )
     np.add(mc_hist[0], mcSum, out=mcSum)
  data_hist = dataify(dataList, nBins, binRange)
  MCScalarSum   = np.sum(mcSum)
  DataScalarSum = np.sum(data_hist[1])
  sumRatio = DataScalarSum / MCScalarSum 
  ratio = np.divide(data_hist[1], mcSum)
  err   = np.multiply(ratio, np.divide(1.0, data_hist[2]))
  np.nan_to_num(ratio, copy=False)
  np.nan_to_num(err, copy=False)

  fig, axi = plt.subplots() #create subplots so I can put a textbox in

  plt.errorbar(data_hist[0], ratio, yerr=err, fmt='o', color='black') #This ignores MC stats.
  try:
    plotTitle, xAxisTitle, yAxisTitle = Titles
  except(ValueError):
    print "Pleast provide three titles for plot name, x and y axis names" 
    plotTitle  = ""
    xAxisTitle = ""
    yAxisTitle = "" 
  plt.title(plotTitle)
  plt.xlabel(xAxisTitle)
  plt.ylabel("Data / MC")
  
  if(len(xLimits) == 2):
    plt.xlim(xLimits)
  if(len(yLimits) == 2):
    plt.ylim(yLimits)

  text = r'$\int \frac{data}{MC} = %.3f$' % sumRatio
  props = dict(boxstyle='round', facecolor='lightsteelblue', alpha=0.5)

  axi.text(0.75, 1.1, text, transform=axi.transAxes, fontsize=14, verticalalignment='top', bbox=props)
  plt.savefig("%s/%sRatio.png" % ( dir_name, filename) )
  plt.close()
  

def KEtoMomentum(T, restMass):
    TotalE = T + restMass
    return math.sqrt(TotalE**2 - restMass**2)
    '''
    0.2 + 0.1 = 0.3 
    sqrt(0.3^2 - 0.1^2) = sqrt(0.09 - 0.01)
    sqrt(0.08) = 0.28
    '''

def getMultibin(dataSequence, weightsSequence):
  
  #Munerahs Bins in KE
  T_mu_bins         = (0.0, 0.12, 0.34, 0.46, 2.0)
  T_p_bins          = (0.0, 0.009, 0.15, 0.21, 1.5)

  #Convert these to momentum
  #P_mu_bins         = tuple(KEtoMomentum(T, 0.1) for T in T_mu_bins)
  #P_p_bins          = tuple(KEtoMomentum(T, 1.0) for T in T_p_bins)
  P_mu_bins          = (0.12, 0.34, 0.52, 0.80, 15.0)
  P_p_bins          = (0.1, 0.42, 0.56, 0.74, 3.0) 
  cos_theta_mu_bins = (-1.0, 0.48, 0.75, 0.89, 1.0)
  cos_theta_p_bins  = (-1.0, 0.48, 0.75, 0.89, 1.0)
  phi_mu_p_bins    = (0.0, 0.69, 0.95, 1.04, 2.0)
  

  bins         = [P_mu_bins, P_p_bins, cos_theta_mu_bins, cos_theta_p_bins, phi_mu_p_bins]
  #bins          = [T_mu_bins, T_p_bins]
  return sci.stats.binned_statistic_dd(dataSequence, weightsSequence, statistic="sum", bins=bins)

def dataify(array, bins, limits):
   counts, bin_edges = np.histogram(array, bins, range=limits)
   bin_centers       = (bin_edges[:-1] + bin_edges[1:]) / 2
   errs = np.sqrt(counts)
   return [bin_centers, counts, errs]

def Stack(dataframe, dirtDF, extDF, variable, longest = False):
  retlist=[]
  addons = ''
  qry = channel_or_particle
  
  if qry == 'channel':
    q_attribute = 'mc_channel'
    value_list = ['QE','RES','DIS','2p2h','NC / Other']
  elif qry == 'particle':
    q_attribute = 'particle'
    value_list = ['muon', 'proton','pion','electron','muon+','other']
  else:
    print('please enter either channel or particle')

  dirt = dirtDF[variable].to_numpy(dtype=float)
  ext = extDF[variable].to_numpy(dtype=float)

  if longest:
    addons = " & isLongestTrack == True"
    dirt = dirtDF.query('isLongestTrack == True')[variable].to_numpy(dtype=float)
    ext = extDF.query('isLongestTrack == True')[variable].to_numpy(dtype=float)

  for value in value_list:
      call = '{} == "{}"'.format(q_attribute,value) + addons
      item = dataframe.query(call)[variable].to_numpy(dtype=float)
      retlist.append(item)
  retlist.append(dirt)
  retlist.append(ext)
  return retlist
   
def getChan(interaction, isNC):
    if (isNC):
      return "NC / Other"
    elif(interaction == 0):
      return "QE"
    elif(interaction == 1):
      return "RES"
    elif(interaction == 2):
      return "DIS"
    elif(interaction == 10):
      return "2p2h"

def getParticle(pdg):
    if (pdg == 13):
      return "muon"
    elif(pdg == 2212):
      return "proton"
    elif(pdg == 211):
      return "pion"
    elif(pdg == 11):
      return "electron"
    else:
      return "other"      

def getChanWeight(interaction, isNC):
    templateFactors = [6.47365e-01, 1.20327e-07, 6.02801e-08, 2.71038e-08, 3.42514e-01, 1.00000e-02]
    if (isNC):
      return (1.0 + templateFactors[4])
    elif(interaction == 0):
      return (1.0 + templateFactors[0])
    elif(interaction == 1):
      return (1.0 + templateFactors[1])
    elif(interaction == 2):
      return (1.0 + templateFactors[2])
    elif(interaction == 10):
      return (1.0 + templateFactors[3])



def dotProduct(df1, df2, dim1, dim2, dimList):
   returnSeries = pd.Series(0.0, index=df1.index.copy())
  # print df1['%s' % dim1]*df1['%s' % "track_dirx"]
  # print df2['%s' % dim2]*df2['%s' % "track_dirx"]
  # print df1['%s' % dim1]*df1['%s' % "track_dirx"] + df2['%s' % dim2]*df2['%s' % "track_dirx"]
   denomSeries = df1['%s' % dim1]*df2['%s' % dim2]
   #print denomSeries

   for dim in dimList:    
       returnSeries += (df1['%s' % dim1]*df1['%s' % dim]*df2['%s' % dim2]*df2['%s' % dim])/denomSeries
       #print returnSeries
   #print returnSeries
   return returnSeries.to_numpy()

def getPhi(pY, pX):
    return (np.arctan2(pY, pX) / np.pi)

def getTheta(pTotal, pZ):
    return (pZ / pTotal)

def getQ2(nuE, nuMu, thetaMu):
    return 4*nuE*nuMu*math.pow(math.sin(thetaMu/2), 2)

#Containment taken from mcc8 Nproton analysis, note 1065
def isContained(xStart, yStart, zStart, xEnd, yEnd, zEnd):
    return (checkContained(xStart, yStart, zStart) and checkContained(xEnd, yEnd, zEnd))

def isFiducial(x, y, z):
  if(x > 270.0 or x < 7.5):
      return False
  elif(y > 110.0 or y < -110.0):  
      return False
  elif(z > 990.0 or z <9.0):
      return False
  else:  
      return True

def checkContained(x, y, z):
    if(x > 245.0 or x < 10.0):
      return False
    elif(y > 100.0 or y < -100.0):  
      return False
    elif(z > 1030.0 or z <10.0):
      return False
    else:  
      return True

def getBestP(isContained, pOne, pTwo):
  return (pOne if isContained else pTwo)

def isTrueCC(pdg, isNC):
  if(pdg == 14 and not isNC):
    return True
  else:
    return False

def getW2(Ehad, Q2, interaction):
    targetMass = 0.989
    #if it's a 2p 2 h, you have 2ps as your target! Wakka Wakka Wakka!
    if(interaction == "2p2h"):
      targetMass = 2*targetMass
    return 2*targetMass*Ehad + math.pow(targetMass, 2) - Q2

def getW(Ehad, Q2, interaction):
    W2 = getW2(Ehad, Q2, interaction)
    if(W2 >= 0.0):
      return math.sqrt(W2)
    else:
      return -1.0

def getXbj(Ehad, Q2):
    targetMass = 0.989
    if(Ehad > 0.0):
      return Q2 / (2*targetMass*Ehad)
    else:
      return -1.0

def getInel(Ehad, Enu):
    if(Enu > 0.0):
      return (Ehad / Enu)
    else:
      return -1.0

def chi2(data, fit, err):
    diff = np.square(np.subtract(data, fit))
    errSquared = np.square(err)
    div = np.divide(diff, errSquared)
    return np.sum(div)

def cubic(x, a, b, c, d):
   return (a*pow(x, 3) + b*pow(x, 2) + c*x + d)

def quartic(x, a, b, c, d, e):
   return (a*pow(x, 4) + b*pow(x, 3) + c*pow(x,2) + d*x + e)
def degSix(x, a, b, c, d, e, f, g):
   return (a*pow(x, 6) + b*pow(x, 5) + c*pow(x,4) + d*pow(x,3) + e*pow(x,2) + f*x + g)
def expQuad(x, a, b, c, d):
   return(np.exp(-a*x)*(b*pow(x, 2) + c*x + d))   
def expQuar(x, a, b, c, d, e, f):
   return(np.exp(-a*x)*(b*pow(x, 4) + c*pow(x,3) + d*pow(x,2) + e*x + f))
def expDegSix(x, a, b, c, d, e, f, g, h, i):
   return(np.exp(-a*x)*(b*pow(x, 6) + c*pow(x,5) + d*pow(x,4) + e*pow(x,3) + f*pow(x,2) + g*x + h) + i)      
InputFiles = ["/uboone/data/users/joelam/stv-ntuples-new/numu_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/bnb_5e19_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC1_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/dirt_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC2_run1.root"]

#InputFiles = ["/uboone/data/users/ametzler/DetVarNtuples/prodgenie_bnb_nu_overlay_DetVar_LYAttenuation.root", "/uboone/data/users/joelam/stv-ntuples-new/bnb_5e19_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC1_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/dirt_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC2_run1.root"]


#InputFiles = ["/uboone/data/users/joelam/stv-ntuples-new/numu_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/nucc_run1_bnb.root", "/uboone/data/users/joelam/stv-ntuples-new/extC1_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/dirt_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC2_run1.root"]


#OverlayScale  = 1.0
ExtScale     = 0.97
numMCTemplates = 6
empty = []


minProtonChi2 = 60.0
maxMuonChi2   = 30.0
minRatioChi2  = 7.0
minTrackScore = 0.5
minNeutrinoScore = 0.1
minMuonTrackScore = 0.85
minTrackL = 20
maxVtxDist = 4
requiredGen = 2
numupdg = 14

#Python library to read in ROOT ntuples files.
overlayEvents = uproot.open(InputFiles[0])["NuCCanalyzer"]["Event"]
bnbEvents     = uproot.open(InputFiles[1])["NuCCanalyzer"]["Event"]
extEventsC1    = uproot.open(InputFiles[2])["NuCCanalyzer"]["Event"]
extEventsC2   = uproot.open(InputFiles[4])["NuCCanalyzer"]["Event"]
dirtEvents    = uproot.open(InputFiles[3])["NuCCanalyzer"]["Event"]

overlayPOT    = uproot.open(InputFiles[0])["NuCCanalyzer"]["subruns"]
dirtPOT       = uproot.open(InputFiles[3])["NuCCanalyzer"]["subruns"]

#Scale factors, because we generate more simulation than data. We also do not take an equal ammount of on and off beam data (though it is close)
mcPOT         = pd.Series(overlayPOT.array("pot"))
sumPOT        = mcPOT.sum()

sumDirtPOT    = (pd.Series(dirtPOT.array("pot"))).sum()


useC2 = True
weightSeperately = False
#dataPOT       = 4.418e+19
dataPOT       = 4.08e+19
#dataPOT       = 1.469e+20 
bnbSpills     = 9045263.0
#bnbSpills     = 10408266.0
#bnbSpills     = 33837051.0  
extTriggersC1 = 33630174.0
extTriggersC2 = 31587147.0
extTriggers   = extTriggersC1
maxEvents     = 35000

subtractDirt = True
subtractExt  = True

'''
mc_pdg = 14
cc_nc = 0
'''


ExtTemplateWeight = (1.0 + 1.00000e-02)

#Create frames of the event tree (which has information about the interaction) and the duaghters tree (which has information about the particles within the interaction).
#Do this for "overlay" simulation, beam data, and off-beam data
overlayDaughters = uproot.open(InputFiles[0])["NuCCanalyzer"]["Daughters"]
trackOverlay   = pd.DataFrame(overlayDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "generation", "mc_pdg", "run", "subrun", "event"] ) )
filteredEvents   = pd.DataFrame(overlayEvents.arrays(["run", "subrun", "event", "mc_nu_interaction_type", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "nu_pdg", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "mc_nu_ccnc", "nu_mu_cc_selected", "mc_nu_lepton_energy", "mc_nu_energy", "mc_nu_lepton_theta"]) )
filteredEvents.insert(filteredEvents.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredEvents['nu_vx'], filteredEvents['nu_vy'], filteredEvents['nu_vz'])] )
filteredEvents.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)

overlayCVWeights = pd.read_csv("/uboone/data/users/joelam/MCWeights/AllBNBWeights.csv", names=["run", "subrun", "event", "wgt_tune"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_tune" : float})
#overlayCVWeights = pd.read_csv("/uboone/data/users/joelam/MCWeights/LYAttenCVWeights.csv", names=["run", "subrun", "event", "wgt_tune"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_tune" : float})

dirtCVWeights    = pd.read_csv("/uboone/data/users/joelam/MCWeights/DirtWeights1.csv", names=["run", "subrun", "event", "wgt_tune"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_tune" : float})

#overlaySplineWeights = pd.read_csv("/uboone/data/users/joelam/MCWeights/LYAttenSplineWeights.csv", names=["run", "subrun", "event", "wgt_spline"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_spline" : float})
overlaySplineWeights = pd.read_csv("/uboone/data/users/joelam/MCWeights/BNBSplineWeights.csv", names=["run", "subrun", "event", "wgt_spline"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_spline" : float})
dirtSplineWeights    = pd.read_csv("/uboone/data/users/joelam/MCWeights/DirtSplineWeights.csv", names=["run", "subrun", "event", "wgt_spline"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_spline" : float})

#OverlayScale = dataPOT / (sumPOT*(float(maxEvents)/trackOverlay.shape[0]))

overlayCVWeights.insert(overlayCVWeights.shape[1], "DuplicatedEvent", overlayCVWeights.duplicated() ) #Tag the events which are duplicated
overlayCVWeights.replace([np.inf], 0.0, inplace=True)

dirtCVWeights.insert(dirtCVWeights.shape[1], "DuplicatedEvent", dirtCVWeights.duplicated() ) #Tag the events which are duplicated
dirtCVWeights.replace([np.inf], 0.0, inplace=True)

overlaySplineWeights.insert(overlaySplineWeights.shape[1], "DuplicatedEvent", overlaySplineWeights.duplicated() ) #Tag the events which are duplicated
overlaySplineWeights.replace([np.inf], 0.0, inplace=True)

dirtSplineWeights.insert(dirtSplineWeights.shape[1], "DuplicatedEvent", dirtSplineWeights.duplicated() ) #Tag the events which are duplicated
dirtSplineWeights.replace([np.inf], 0.0, inplace=True)

dataDaughters = uproot.open(InputFiles[1])["NuCCanalyzer"]["Daughters"]
trackData     = pd.DataFrame(dataDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz","track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "is_track_daughter", "generation", "run", "subrun", "event"] ) )
filteredData  = pd.DataFrame(bnbEvents.arrays(["run", "subrun", "event", "nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "nu_pdg"]) )
filteredData.insert(filteredData.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredData['nu_vx'], filteredData['nu_vy'], filteredData['nu_vz'])] )
filteredData.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)

extDaughtersC1 = uproot.open(InputFiles[2])["NuCCanalyzer"]["Daughters"]
trackExtC1     = pd.DataFrame(extDaughtersC1.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "is_track_daughter", "generation", "run", "subrun", "event"] ) )
filteredExtC1  = pd.DataFrame(extEventsC1.arrays(["run", "subrun", "event", "nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "nu_pdg"]) )
filteredExtC1.insert(filteredExtC1.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredExtC1['nu_vx'], filteredExtC1['nu_vy'], filteredExtC1['nu_vz'])] )
filteredExtC1.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)

if(weightSeperately):
  filteredExtC1.insert(filteredExtC1.shape[1], "wgt", (bnbSpills / (2*extTriggers) ) )

extDaughtersC2 = uproot.open(InputFiles[4])["NuCCanalyzer"]["Daughters"]
trackExtC2     = pd.DataFrame(extDaughtersC2.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "is_track_daughter", "generation", "run", "subrun", "event"] ) )
filteredExtC2  = pd.DataFrame(extEventsC2.arrays(["run", "subrun", "event", "nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "nu_pdg"]) )
filteredExtC2.insert(filteredExtC2.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredExtC2['nu_vx'], filteredExtC2['nu_vy'], filteredExtC2['nu_vz'])] )
filteredExtC2.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)

if(weightSeperately):
   filteredExtC2.insert(filteredExtC2.shape[1], "wgt", (bnbSpills / (2*extTriggersC2) ) )

dirtDaughters = uproot.open(InputFiles[3])["NuCCanalyzer"]["Daughters"]
trackDirt     = pd.DataFrame(dirtDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "is_track_daughter", "generation", "run", "subrun", "event"] ) )
filteredDirt  = pd.DataFrame(dirtEvents.arrays(["run", "subrun", "event", "nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "nu_pdg"]) )
filteredDirt.insert(filteredDirt.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredDirt['nu_vx'], filteredDirt['nu_vy'], filteredDirt['nu_vz'])] )
filteredDirt.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)


if(useC2):
  filteredExt = pd.concat([filteredExtC1, filteredExtC2])
  trackExt    = pd.concat([trackExtC1, trackExtC2])
  extTriggers = extTriggers + extTriggersC2

else:
  filteredExt = filteredExtC1
  trackExt    = trackExtC1



#Here, we calculate some additional event information that isn't part of the input ROOT ntuple
#This is because the grad. student who created the files didn't include this information
loadMCEventInfo(filteredEvents)
loadTrackInfo(trackOverlay, True)

tagDuplicateEvents(filteredDirt)
loadTrackInfo(trackDirt)

#SAVE THIS
#print trackOverlay.loc[:100,['isContained', 'track_range_mom_mu', 'track_mcs_mom', 'track_mom_best']]

tagDuplicateEvents(filteredExt)

if not weightSeperately:
  ExtScale     = bnbSpills / extTriggers
  extWeights              = np.full(filteredExt.shape[0],  ExtScale)
  extTemplateWeights      = np.full(filteredExt.shape[0],  ExtTemplateWeight)
  filteredExt.insert(filteredExt.shape[1], "wgt", extWeights )

tagDuplicateEvents(filteredData)
loadTrackInfo(trackData)

loadTrackInfo(trackExt)

OverlayScale = dataPOT / sumPOT
DirtScale    = dataPOT / sumDirtPOT
print "MC POT: %e or %e Overlay Scale: %.3f Ext Scale: %.3f Dirt Scale: %.3f" % (sumPOT, sumPOT*(float(maxEvents)/trackOverlay.shape[0]), OverlayScale, ExtScale, DirtScale)
print "Total MC POT: %e total MC events: %d" % (mcPOT.sum(), trackOverlay.shape[0])
print "Total Dirt POT: %e" % sumDirtPOT

#Index the events and daugthers by the run, subrun, event tuple
#This is IMPORTANT. The only infomration we have to connect the two frames a priori is this set of 3 ints
#A single event can have multiple tracks (and often does!)
#Multiindexing makes our life much easier, cuz we can grab the event info for ANY track from it's multiindex

trackOverlay   =  trackOverlay.set_index(['run', 'subrun', 'event'])
filteredEvents =  filteredEvents.set_index(['run', 'subrun', 'event'])

filteredData   =  filteredData.set_index(['run', 'subrun', 'event'])
trackData      =  trackData.set_index(['run', 'subrun', 'event'])

filteredExt    = filteredExt.set_index(['run', 'subrun', 'event'])
trackExt       = trackExt.set_index(['run', 'subrun', 'event'])

filteredDirt    = filteredDirt.set_index(['run', 'subrun', 'event'])
trackDirt       = trackDirt.set_index(['run', 'subrun', 'event'])

overlayCVWeights = overlayCVWeights.set_index(['run', 'subrun', 'event'])
dirtCVWeights    = dirtCVWeights.set_index(['run', 'subrun', 'event'])

overlaySplineWeights = overlaySplineWeights.set_index(['run', 'subrun', 'event'])
dirtSplineWeights    = dirtSplineWeights.set_index(['run', 'subrun', 'event'])

#Do this to make our loops and lookups a bit more efficienct

trackOverlay.sort_index()
filteredEvents.sort_index()


filteredEvents   =  filteredEvents[filteredEvents.DuplicatedEvent == False]
filteredData     =  filteredData[filteredData.DuplicatedEvent == False]
filteredExt      =  filteredExt[filteredExt.DuplicatedEvent == False]
filteredDirt     =  filteredDirt[filteredDirt.DuplicatedEvent == False]

trackOverlay     =  trackOverlay[trackOverlay.DuplicatedEvent == False]
trackData        =  trackData[trackData.DuplicatedEvent == False]
trackExt         =  trackExt[trackExt.DuplicatedEvent == False]
trackDirt        =  trackDirt[trackDirt.DuplicatedEvent == False]


overlayCVWeights =  overlayCVWeights[overlayCVWeights.DuplicatedEvent == False]
dirtCVWeights    =  dirtCVWeights[dirtCVWeights.DuplicatedEvent == False]
overlaySplineWeights =  overlaySplineWeights[overlaySplineWeights.DuplicatedEvent == False]
dirtSplineWeights    =  dirtSplineWeights[dirtSplineWeights.DuplicatedEvent == False]



numberFiltered = 0


#create a dict of event info we want to associate with each daughter.
#by doing this, we have the complete event information for each track.
#Really what we want is to look at the particles' properties as a funciton of the underlying event information
#This is extendible to any event varaible we want to associate to a particle
interactionInfo = ("mc_channel", "nu_mu_cc_selected","nu_score", "nu_pdg", "nu_flash_chi2", "obvious_cosmic_chi2", "flash_chi2_ratio", "nu_vx", "nu_vy", "nu_vz", "daughters_start_contained", "isFiducial", "mc_Ehad", "mc_expQ2", "mc_expXbj", "mc_expY", "mc_expW") 

for field in interactionInfo:
  trackOverlay   = trackOverlay.join(filteredEvents['%s' % field], on=["run", "subrun", "event"])

trackOverlay = trackOverlay.join(overlayCVWeights['wgt_tune'], on=["run", "subrun", "event"])
trackOverlay = trackOverlay.join(overlaySplineWeights['wgt_spline'], on=["run", "subrun", "event"])

extInfo = { "nu_mu_cc_selected", "nu_score", "nu_pdg", "nu_flash_chi2", "obvious_cosmic_chi2", "flash_chi2_ratio", "nu_vx", "nu_vy", "nu_vz", "daughters_start_contained", "isFiducial", "wgt" }

for field in extInfo:
  trackExt   = trackExt.join(filteredExt['%s' % field], on=["run", "subrun", "event"])

dirtInfo = ("nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "flash_chi2_ratio", "nu_vx", "nu_vy", "nu_vz", "daughters_start_contained", "isFiducial", "nu_pdg")

for field in dirtInfo:
  trackDirt = trackDirt.join(filteredDirt['%s' % field], on=["run", "subrun", "event"])

#SAVE THIS
#print dirtCVWeights.at[(6553, 129, 6458), 'wgt_tune']

for field in dirtInfo:
  trackData  = trackData.join(filteredData['%s' % field], on=["run", "subrun", "event"])

dirtCVWeightMeans     = dirtCVWeights.groupby(level=["run", "subrun", "event"]).agg({"wgt_tune" : ["mean"]})
dirtSplineWeightMeans = dirtSplineWeights.groupby(level=["run", "subrun", "event"]).agg({"wgt_spline" : ["mean"]})

dirtCVWeightMeans.columns = ["_".join(x) for x in dirtCVWeightMeans.columns.ravel()]
dirtCVWeightMeans.rename(columns={"wgt_tune_mean" : "wgt_tune"}, inplace=True)

dirtSplineWeightMeans.columns = ["_".join(x) for x in dirtSplineWeightMeans.columns.ravel()]
dirtSplineWeightMeans.rename(columns={"wgt_spline_mean" : "wgt_spline"}, inplace=True)

trackDirt = trackDirt.join(dirtCVWeightMeans['wgt_tune'], on=["run", "subrun", "event"])
trackDirt = trackDirt.join(dirtSplineWeightMeans['wgt_spline'], on=["run", "subrun", "event"])

muonMomentumRange   = (0.0, 2.0)
protonMomentumRange = (0.0, 1.5)
phiRange = (-1.0, 1.0)
isSelectedRange = (0.0, 1.0)
phiDiffRange = (0.0, 2.0)
trkScoreRange= (0.0, 1.0)
chi2Range    = (0.0, 50.0)
chi2PRange   = (0.0, 350.0)
vtxRange     = (0.0, 6.0)
lengthRange  = (0.0, 200.0)
flashChi2Range = (0.0, 300.0)
pdgRange     = (0, 30)


overlayWeights = np.full(trackOverlay.shape[0], OverlayScale )
dirtWeights    = np.full(trackDirt.shape[0], DirtScale )




trackOverlay.insert(trackOverlay.shape[1], "pot_wgt", overlayWeights )
trackOverlay.eval('wgt = pot_wgt*wgt_tune*wgt_spline', inplace=True) 

trackDirt.insert(trackDirt.shape[1], "pot_wgt", dirtWeights )
trackDirt.eval('wgt = pot_wgt*wgt_tune*wgt_spline', inplace=True)

overlaySliceScoreStack = '''[trackOverlay.query('mc_channel == "QE"')['VAR'].to_numpy(), trackOverlay.query('mc_channel == "RES"')['VAR'].to_numpy(), trackOverlay.query('mc_channel == "DIS"')['VAR'].to_numpy(), trackOverlay.query('mc_channel == "2p2h"')['VAR'].to_numpy(), trackOverlay.query('mc_channel == "NC / Other"')['VAR'].to_numpy(), trackDirt['VAR'].to_numpy(), trackExt['VAR'].to_numpy()]'''
exec( "incSliceScoreStack   = "  + re.sub(r'VAR', 'nu_score', overlaySliceScoreStack) )
exec( "incSliceScorekWeights     = "  + re.sub(r'VAR', 'wgt',    overlaySliceScoreStack) )

#makeDataMCHistogram(incSliceScoreStack, incSliceScorekWeights, trackData['nu_score'].to_numpy(), trkScoreRange, 25, "IncSliceScore", ["Slice Score", "Score", "Number of Daughters"])

exec( "incIsSelectedStack   = "  + re.sub(r'VAR', 'nu_mu_cc_selected', overlaySliceScoreStack) )
#exec( "incSliceScorekWeights     = "  + re.sub(r'VAR', 'wgt',    overlaySliceScoreStack) )

#makeDataMCHistogram(incIsSelectedStack, incSliceScorekWeights, trackData['nu_mu_cc_selected'].to_numpy(), isSelectedRange, 2, "IncIsSelected", ["Selected", "Selected", "Number of Daughters"])


#extNuScore     = trackExt.query('DuplicatedEvent == False & nu_score > @minNeutrinoScore  & nu_pdg == @numupdg')

extMuonCandidates      = trackExt.query('DuplicatedEvent == False & track_score > @minMuonTrackScore  & vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen & track_chi2_proton > @minProtonChi2 & track_chi2_muon < @maxMuonChi2 & track_chi2_ratio > @minRatioChi2 & nu_score > @minNeutrinoScore & isContained == True')
dirtMuonCandidates     = trackDirt.query('DuplicatedEvent == False & track_score > @minMuonTrackScore &  vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen & track_chi2_proton > @minProtonChi2 & track_chi2_muon < @maxMuonChi2 & track_chi2_ratio > @minRatioChi2 & nu_score > @minNeutrinoScore & isContained == True')
overlayMuonCandidates  = trackOverlay.query('DuplicatedEvent == False & track_score > @minMuonTrackScore & vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen & track_chi2_proton > @minProtonChi2 & track_chi2_muon < @maxMuonChi2 & track_chi2_ratio > @minRatioChi2 & nu_score > @minNeutrinoScore & isContained == True')
dataMuonCandidates     = trackData.query('DuplicatedEvent == False & track_score > @minMuonTrackScore &  vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen & track_chi2_proton > @minProtonChi2 & track_chi2_muon < @maxMuonChi2 & track_chi2_ratio > @minRatioChi2 & nu_score > @minNeutrinoScore & isContained == True')


statFrame = extMuonCandidates.groupby(level=["run", "subrun", "event"]).agg({"track_length": ["max"]})
statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]
extMuonCandidates = extMuonCandidates.join(statFrame['track_length_max'], on=["run", "subrun", "event"])  
extMuonCandidates.eval('isLongestTrack = (track_length == track_length_max)', inplace=True)

statFrame = dirtMuonCandidates.groupby(level=["run", "subrun", "event"]).agg({"track_length": ["max"]})
statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]
dirtMuonCandidates = dirtMuonCandidates.join(statFrame['track_length_max'], on=["run", "subrun", "event"])  
dirtMuonCandidates.eval('isLongestTrack = (track_length == track_length_max)', inplace=True)

statFrame = overlayMuonCandidates.groupby(level=["run", "subrun", "event"]).agg({"track_length": ["max"]})
statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]
overlayMuonCandidates = overlayMuonCandidates.join(statFrame['track_length_max'], on=["run", "subrun", "event"])  
overlayMuonCandidates.eval('isLongestTrack = (track_length == track_length_max)', inplace=True)

statFrame = dataMuonCandidates.groupby(level=["run", "subrun", "event"]).agg({"track_length": ["max"]})
statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]
dataMuonCandidates = dataMuonCandidates.join(statFrame['track_length_max'], on=["run", "subrun", "event"])  
dataMuonCandidates.eval('isLongestTrack = (track_length == track_length_max)', inplace=True)



'''
maxFlashChi2 = 10
minNeutrinoScoreFlashFails = 0.25
maxFlashChi2Ratio  = 5

extInclusiveEvents = extMuonCandidates.query('isLongestTrack == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & (nu_flash_chi2 < @maxFlashChi2 | nu_score > @minNeutrinoScoreFlashFails)')

dirtInclusiveEvents = dirtMuonCandidates.query('isLongestTrack == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & (nu_flash_chi2 < @maxFlashChi2 | nu_score > @minNeutrinoScoreFlashFails)')

overlayInclusiveEvents = overlayMuonCandidates.query('isLongestTrack == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & (nu_flash_chi2 < @maxFlashChi2 | nu_score > @minNeutrinoScoreFlashFails)')

dataInclusiveEvents = dataMuonCandidates.query('isLongestTrack == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & (nu_flash_chi2 < @maxFlashChi2 | nu_score > @minNeutrinoScoreFlashFails)')
'''

#SAVE THIS
#print dataInclusiveEvnets.loc[(5774, 15,762 ),('track_chi2_muon', 'track_chi2_proton', 'track_chi2_ratio', 'isFiducial', 'nu_score', 'nu_flash_chi2', 'nu_mu_cc_selected')]



overlayPrimMuonChi2FlashInclusiveStack = Stack(overlayMuonCandidates, dirtMuonCandidates, extMuonCandidates, "nu_flash_chi2", True)
overlayPrimMuonPhiInclusiveStack = Stack(overlayMuonCandidates, dirtMuonCandidates, extMuonCandidates, "phi", True)
overlayPrimMuonScoreInclusiveStack = Stack(overlayMuonCandidates, dirtMuonCandidates, extMuonCandidates, "nu_score", True)
overlayPrimMuonWeightInclusiveStack    = Stack(overlayMuonCandidates, dirtMuonCandidates, extMuonCandidates, "wgt", True)

makeDataMCHistogram(overlayPrimMuonChi2FlashInclusiveStack, overlayPrimMuonWeightInclusiveStack, dataMuonCandidates['nu_flash_chi2'].to_numpy(), (0, 50), 50, "ContainedPrimMuonFlashChi2", ["Flash Chi2", "Chi2", "Number of Events"])
makeDataMCHistogram(overlayPrimMuonPhiInclusiveStack, overlayPrimMuonWeightInclusiveStack, dataMuonCandidates['phi'].to_numpy(), phiRange, 30, "ContainedPrimMuonPhi", ["Muon Phi Angle", "Angle / pi (radians)", "Number of Primary Muons"])
makeDataMCHistogram(overlayPrimMuonScoreInclusiveStack, overlayPrimMuonWeightInclusiveStack, dataMuonCandidates['nu_score'].to_numpy(), (0, 1.0), 50, "ContainedPrimMuonNuScore", ["Neutrino Score", "Nu Score", "Number of Primary Muons"])


flash_chi2_bins = []
#flash_chi2_bins.append(0)
firstBreak = 11
secondBreak = 20
thirdBreak = 26
fourthBreak = 37
fifthBreak  = 40
sixthBreak  = 45
seventhBreak = 50
eigthBreak  = 54
ninthBreak  = 55
maxBins = firstBreak + secondBreak + thirdBreak + fourthBreak


firstBin = 0.5
for i in range(maxBins):
  if (i == 0):
    flash_chi2_bins.append(0.0)
  elif i == 1:
    flash_chi2_bins.append(0.5)  
  elif i > 1 and i <= firstBreak:
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 0.2)
  elif i > firstBreak and i <= secondBreak:
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 0.5)
  elif i > secondBreak and i<= thirdBreak:
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 1.0)
  elif  i > thirdBreak and i<= fourthBreak:
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 2.0) 
  elif i > fourthBreak and i<= fifthBreak:    
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 5.0)
  elif i > fifthBreak and i<= sixthBreak:    
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 10.0)
  elif i > sixthBreak and i<= seventhBreak:    
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 20.0)
  elif i > seventhBreak and i<= eigthBreak:    
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 50.0)
  elif i > eigthBreak and i<= ninthBreak:    
    flash_chi2_bins.append(flash_chi2_bins[i-1] + 100.0)                  



#print flash_chi2_bins
flash_chi2_bins = [0.0, 0.5, 0.7, 0.9, 1.1, 1.5,  2.0,  2.5, 3.0,  4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 31.0, 33.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0]
#flash_chi2_bins = [0.0, 0.5, 0.9, 1.5,  2.5, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 17.0, 21.0, 25.0, 29.0, 33.0, 40.0, 50.0, 70.0, 90.0, 120.0, 160.0, 200.0, 300.0, 400.0, 500.0]
#print flash_chi2_bins

numberOfBins = len(flash_chi2_bins)-1
binRange = (0, 500.0)

dir_name  = "FittingPlots"
mcSum = np.full(numberOfBins, 0.0 )
dirtHist =  np.full(numberOfBins, 0.0 )
extHist  =  np.full(numberOfBins, 0.0 )
dirtErr =  np.full(numberOfBins, 0.0 )
extErr  =  np.full(numberOfBins, 0.0 )
#print overlayPrimMuonChi2FlashInclusiveStack[:-1]
it = 0
for mc, weight in zip(overlayPrimMuonChi2FlashInclusiveStack, overlayPrimMuonWeightInclusiveStack ):
   #print mc
   mc_hist   = np.histogram(mc, bins=flash_chi2_bins, range=binRange, weights = weight )
   if(subtractDirt and it == len(overlayPrimMuonChi2FlashInclusiveStack) - 2):
      dirtHist = mc_hist[0]
      dirtErr  = np.sqrt(dirtHist)
      it = it + 1
      continue
   elif(subtractExt and it == len(overlayPrimMuonChi2FlashInclusiveStack) - 1):
      extHist = mc_hist[0]
      extErr  = np.sqrt(extHist)
      it = it +1
      continue
   else:
      np.add(mc_hist[0], mcSum, out=mcSum)
      it = it + 1       
data_hist = dataify(dataMuonCandidates['nu_flash_chi2'].to_numpy(), flash_chi2_bins, binRange)
MCScalarSum   = np.sum(mcSum)
dataNP = np.array(data_hist[1], dtype=float)
ratioNoSub = np.divide(dataNP, mcSum)
if(subtractDirt):
  np.subtract(dataNP, dirtHist, out=dataNP)
if(subtractExt):
  np.subtract(dataNP, extHist,  out=dataNP)
DataScalarSum = np.sum(dataNP)
sumRatio = DataScalarSum / MCScalarSum 
ratio = np.divide(dataNP, mcSum)
dataErr = np.array(data_hist[2])
print dataErr
np.add(np.square(dataErr), np.square(dirtErr), out=dataErr)
np.add(dataErr, np.square(extErr), out=dataErr)
np.sqrt(dataErr, out=dataErr)
print dataErr
err   = np.multiply(ratioNoSub, np.divide(1.0, dataErr) )
np.nan_to_num(ratio, copy=False)
np.nan_to_num(err, copy=False)

bin_widths  = np.subtract(flash_chi2_bins[1:], flash_chi2_bins[:-1])
bin_widths  = np.divide(bin_widths, 2.0)
bin_centers = np.add(flash_chi2_bins[:-1], flash_chi2_bins[1:]) / 2
plt.errorbar(data_hist[0], ratio, xerr=bin_widths, yerr=err, fmt='o', color='black')
filename = "DataMCPreFit"
plt.xlabel("Flash Chi2")
plt.ylabel("Data / MC")
plt.xlim(0, 200)
plt.ylim(0, 2.0)
plt.savefig("%s/%s.png" % ( dir_name, filename) )
plt.close()


plt.errorbar(data_hist[0], ratio , yerr=err, xerr=bin_widths, fmt='o', color='black')
filename = "DataMCPreFitZoom"
plt.xlabel("Flash Chi2")
plt.ylabel("Data / MC")
plt.xlim(0, 25)
plt.ylim(0, 2.0)
plt.savefig("%s/%s.png" % ( dir_name, filename) )
plt.close()

dropBins = 5
#print bin_centers
#print bin_centers[:-dropBins]
func, cov = sci.optimize.curve_fit(expDegSix, bin_centers[:-dropBins], ratio[:-dropBins], sigma=err[:-dropBins], absolute_sigma=True)


#print len(func)
#print func
#print chi2(ratio[:-dropBins], expDegSix(bin_centers, *func)[:-dropBins], err[:-dropBins])
print "Chi2 : %.3f dof: %d" % (chi2(ratio[:-dropBins], expDegSix(bin_centers, *func)[:-dropBins], err[:-dropBins]), (len(ratio) - dropBins - len(func) ) )


plt.plot(bin_centers, expDegSix(bin_centers, *func), 'r-')
plt.errorbar(data_hist[0], ratio, xerr=bin_widths, yerr=err, fmt='o', color='black')
filename = "FitFunction"
plt.xlabel("Flash Chi2")
plt.ylabel("Data / MC")
plt.xlim(0, 200)
plt.ylim(0, 2.0)
plt.savefig("%s/%s.png" % ( dir_name, filename) )
plt.close() 

plt.plot(bin_centers, expDegSix(bin_centers, *func), 'r-')
plt.errorbar(data_hist[0], ratio, xerr=bin_widths, yerr=err, fmt='o', color='black')
filename = "FitFunctionZoom"
plt.xlabel("Flash Chi2")
plt.ylabel("Data / MC")
plt.xlim(0, 25)
plt.ylim(0, 2.0)
plt.savefig("%s/%s.png" % ( dir_name, filename) )
plt.close() 

#makeDataMCHistogram(overlayPrimMuonChi2FlashInclusiveStack, overlayPrimMuonWeightInclusiveStack , dataMuonCandidates['nu_flash_chi2'].to_numpy(), (0, 500), flash_chi2_bins, "InclusiveEventsPrimMuonFlashChi2", ["Flash Chi2", "Chi2", "Number of Events"])

#makeDataMCHistogram(overlayPrimMuonPhiInclusiveStack, overlayIsSelectedInclusiveWeights, dataInclusiveEvents['phi'].to_numpy(), phiRange, 30, "InclusiveEventsPrimMuonPhi", ["Muon Phi Angle", "Angle / pi (radians)", "Number of Primary Muons"])


#print dataInclusiveEvents.query('nu_mu_cc_selected == False')

print func

sys.exit()
