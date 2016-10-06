#!/usr/bin/python

# Imports
import ROOT
import numpy as np
import sys, os
import collections

# Setup WCSim libraries
WCSimRootLibPath = "/storage/shared/cvilela/WCSim_v151/WCSim/libWCSimRoot.so"
ROOT.gSystem.Load(WCSimRootLibPath)
#ROOT.SetMemoryPolicy( ROOT.kMemoryStrict )

# Change the following to run on a subset of files
#hornCurrentModes = ["FHC", "RHC"]
hornCurrentModes = ["FHC"]

#nuTypes = {"nue"      : 12,
#           "numu"     : 14,
#           "nutau"    : 16,
#           "nuebar"   : -12,
#           "numubar"  : -14,
#           "nutaubar" : -16}

nuTypes = {"nue"      : 12,
           "numu"     : 14}

# Number of files per horn mode per nutype to process. Maximum is 100.
nFiles = 100

baseFilePath = "/storage/shared/cvilela/DuneWC/"
OutputPath = baseFilePath+"DuneWCana/Out/"

# Incoming neutrino direction is along z. This should be corrected in a future iteration...
# Hopefully not a very significant problem with such a big tank. Look at efficiency vs angle (relative to tank) and momentum
trueNuDir = [0., 0., 1.]
# Detector radius and half-height
detR = 3700.
detHH = 3000.
dWallCut = 200.

# These functions return the paths to input files
def fqFilePath(nuType, hornMode, fileNumber) :
    return baseFilePath+"gntp_"+hornMode+"_"+nuType+"_fluxosc/"+"fitqun/"+"/gntp_"+hornMode+"_"+nuType+"_fluxosc_fQv5r4_"+str.rjust(str(fileNumber), 3, '0')+".root"

def wcFilePath(nuType, hornMode, fileNumber) :
    return baseFilePath+"gntp_"+hornMode+"_"+nuType+"_fluxosc/"+"wcsim/"+"/gntp_"+hornMode+"_"+nuType+"_fluxosc_"+str.rjust(str(fileNumber), 3, '0')+".root"

# don't think this will be necessary, but just in case
def vecFilePath(nuType, hornMode, fileNumber) :
    return baseFilePath+"gntp_"+hornMode+"_"+nuType+"_fluxosc/"+"vectors/"+"/gntp_"+hornMode+"_"+nuType+"_fluxosc_"+str.rjust(str(fileNumber), 3, '0')+".nuance"

# Selection criteria
# Single-ring mu-like cuts:
def FCmu (event):
    return event.fq1rpcflg[0*7+2] == 0

def FVmu (event):
    pos_l = [event.fq1rpos[0*7*3 + 2*3 + 0],
             event.fq1rpos[0*7*3 + 2*3 + 1],
             event.fq1rpos[0*7*3 + 2*3 + 2]]
    # FV
    if (pos_l[0]**2 + pos_l[1]**2) > (detR - dWallCut)**2 :
        return False    
    if pos_l[2]**2 > (detHH - dWallCut)**2 :
        return False
    return True

def nRingMu ( event ) :
    return event.fqmrnring[0] == 1

def pidMuE (event) :
    return event.fq1rnll[0*7 + 2] - event.fq1rnll[0*7 + 1] < 2000

def nseMu (event) :
    return True
#    return event.fqnse <= 2

# Single-ring e-like cuts:
def FCe (event):
    return event.fq1rpcflg[0*7+1] == 0

def FVe (event):
    pos_l = [event.fq1rpos[0*7*3 + 1*3 + 0],
             event.fq1rpos[0*7*3 + 1*3 + 1],
             event.fq1rpos[0*7*3 + 1*3 + 2]]
    # FV
    if (pos_l[0]**2 + pos_l[1]**2) > (detR - dWallCut)**2 :
        return False
    if pos_l[2]**2 > (detHH - dWallCut)**2 :
        return False
    return True

def nRingE ( event ) :
    return event.fqmrnring[0] == 1

def pidEmu (event) :
    return event.fq1rnll[0*7 + 2] - event.fq1rnll[0*7 + 1] >= 2000

def nseE (event) :
    return True
    #return event.fqnse == 1

def pidEpi0 (event) :
    return event.fqpi0nll[0] - event.fq1rnll[0*7+1] >= 0 

# Interaction mode category criteria
def CC0pi0pcriterion(event) :

    isCC = True

    wcEv = event.wcsimrootevent
    if abs(wcEv.GetTrigger(0).GetMode()) >= 30 :
        isCC = False

    # Loop through WCSim true particles and count pions/protons
    for par in range(0, wcEv.GetTrigger(0).GetNtrack()) :
        if abs(wcEv.GetTrigger(0).GetTracks().At(par).GetIpnu()) in [111, 211, 2122] :
            if wcEv.GetTrigger(0).GetTracks().At(par).GetFlag() == 0 :
                isCC = False
                break 
    del wcEv
    return isCC

def CCothercriterion(event) :

    isCCother = True
    
    wcEv = event.wcsimrootevent
    # Check that is CC and not CC0pi0pcriterion
    if abs(wcEv.GetTrigger(0).GetMode()) >= 30 :
        isCCother = False
    if CC0pi0pcriterion(event) :
        isCCother = False

    return isCCother

def NCcriterion(event) :
    wcEv = event.wcsimrootevent
    mode = wcEv.GetTrigger(0).GetMode()
    del wcEv
    return abs(mode) >= 30 
        
def main() :

    # Initialise event selectors
    muSel = eventSelector("SingleRingMuLike")
    muSel.addCriterion("FCmu", FCmu)
    muSel.addCriterion("FVmu", FVmu)
    muSel.addCriterion("nRingMu", nRingMu)
    muSel.addCriterion("pidMuE", pidMuE)
    muSel.addCriterion("nseMu", nseMu)

    eSel = eventSelector("SingleRingELike")
    eSel.addCriterion("FCe", FCe)
    eSel.addCriterion("FVe", FVe)
    eSel.addCriterion("nRingE", nRingE)
    eSel.addCriterion("pidEmu", pidEmu)
    eSel.addCriterion("nseE", nseE)
    eSel.addCriterion("pidEpi0", pidEpi0)

    # Initialise mode selectors
    cc0pi0pSel = intModeSelector("CC0pi0p")
    cc0pi0pSel.addCriterion("AllCuts", CC0pi0pcriterion)
    ccOtherSel = intModeSelector("CCother")
    ccOtherSel.addCriterion("AllCuts", CCothercriterion)
    ncSel = intModeSelector("NC")
    ncSel.addCriterion("AllCuts", NCcriterion)

    # Dictionaries to keep histograms
    hists1D = {}
    hists2D = {}
    
    # Loop through horn modes
    for hornMode in hornCurrentModes :
        # Loop through incoming neutrino types
        for nuType in nuTypes :

            # Open output file
            outFpath = OutputPath+"/"+hornMode+"/"+nuType
            if not os.path.isdir( outFpath ) :
                os.makedirs( outFpath );

            fOut = ROOT.TFile(outFpath+"/DuneWCHists"+".root", "RECREATE")
            # Book histograms:
            #              Name       nBins  Min Max
            histList1D = {"hErec" : [100,   0., 10000]}

            # Cut flow histograms, one for each combination of interaction mode and event sample:
            for mode in intModeSelector.modes :
                for sample in eventSelector.samples :
                    histList1D["cutFlow_"+sample.name] = [len(sample.criteria), 0, len(sample.criteria)]

            # 2D histograms
            histList2D = {"hPCosTheta" : [100, 0, 5000, 100, -1, 1],
                          "hNLLRePi0" : [100, 0, 300, 100, -5000, 5000],
                          "hNLLReMu"  : [100, 0, 5000, 100, -5000, 5000],
                          "hEtrueErec"  : [100, 0, 5000, 100, 0, 5000]}

            hists1D = merge( hists1D, bookHists(histList1D, 1, hornMode, nuType) )
            hists2D = merge( hists2D, bookHists(histList2D, 2, hornMode, nuType) )

            # File loop
            for fNum in range(0, nFiles) :
                # Open ROOT files
                fFQ = ROOT.TFile(fqFilePath(nuType = nuType, hornMode = hornMode, fileNumber = fNum))
                fWC = ROOT.TFile(wcFilePath(nuType = nuType, hornMode = hornMode, fileNumber = fNum))

                # Grab TTrees
                tFQ = fFQ.Get("fiTQun")
                tWC = fWC.Get("wcsimT")

                # Slightly paranoid, but wihtout (at least one of) these there is a bad memory leak!
                tWC.GetBranch("wcsimrootevent").SetAutoDelete(True)
                tFQ.AddFriend(tWC)
                tFQ.GetBranch("wcsimrootevent").SetAutoDelete(True)

                print "Analysing ", hornMode, nuType, fNum

                # Loop through events
                for event in tFQ :

                    sampleName = ""
                    modeName = ""

                    isMode = False
                    isSelected = False

                    for mode in intModeSelector.modes :
                        isMode =  mode.passCriteria(event)
                        if isMode :
                            modeName = mode.name
                            break

                    for sample in eventSelector.samples :
                        isSelected = sample.passCriteria(event, modeName, nuType)
                        if isSelected :
                            sampleName = sample.name
                            break

                    if isMode and isSelected :

                        # Variables for getting kinematics / Erec
                        fQindex = -1
                        m_l = -1
                        m_n = 939.6
                        m_p = 938.2
                        v_nuc = 27.

                        if sampleName == "SingleRingELike" :
                            fQindex = 1
                            m_l = 0.511
                        elif sampleName == "SingleRingMuLike" :
                            fQindex = 2
                            m_l = 105.7
                        else :
                            print "Unkown sample, not filling histograms"
                            continue

                        p_l = event.fq1rmom[0*7+fQindex]

                        E_tot = (p_l**2 + m_l**2)**0.5

                        dir_nu = trueNuDir

                        dir_l = [event.fq1rdir[0*7*3 + 1*3 + 0],
                                 event.fq1rdir[0*7*3 + 1*3 + 1],
                                 event.fq1rdir[0*7*3 + 1*3 + 2]]

                        cos_th_beam = sum([dir_nu[i]*dir_l[i] for i in 0,1,2]) / (sum([x**2 for x in dir_l])**0.5 * sum([x**2 for x in dir_nu])**0.5)

                        if hornMode == "FHC" :
                            Erec = ( (m_n - v_nuc)*E_tot - m_l**2/2 + m_n*v_nuc - v_nuc**2/2 + (m_p**2 - m_n**2)/2 ) / ( m_n - v_nuc - E_tot + p_l *cos_th_beam )
                        elif hornMode == "RHC" :
                            # Double check this... is Vnuc the same?
                            Erec = ( (m_p - v_nuc)*E_tot - m_l**2/2 + m_p*v_nuc - v_nuc**2/2 + (m_n**2 - m_p**2)/2 ) / ( m_p - v_nuc - E_tot + p_l *cos_th_beam )

                        # Fill Histograms here
                        hists1D["hErec"][hornMode][sampleName][modeName][nuType].Fill(Erec)
                        hists2D["hPCosTheta"][hornMode][sampleName][modeName][nuType].Fill(p_l, cos_th_beam)
                        hists2D["hNLLRePi0"][hornMode][sampleName][modeName][nuType].Fill(event.fqpi0mass[0], event.fqpi0nll[0] - event.fq1rnll[0*7+1])
                        hists2D["hNLLReMu"][hornMode][sampleName][modeName][nuType].Fill(event.fq1rmom[0*7+1], event.fq1rnll[0*7+2] - event.fq1rnll[0*7+1])
                        wcEv = event.wcsimrootevent
                        hists2D["hEtrueErec"][hornMode][sampleName][modeName][nuType].Fill(Erec, wcEv.GetTrigger(0).GetTracks().At(0).GetP())

                        del wcEv
                    del event

                    # Make some plots (set to false for batch)
            makePlots = True
            if makePlots :
                ErecStack = ROOT.THStack("Erec", hornMode+" "+nuType+";Erec [MeV]")
                hists1D["hErec"][hornMode]["SingleRingMuLike"]["CC0pi0p"][nuType].SetFillColor(ROOT.kRed)
                hists1D["hErec"][hornMode]["SingleRingMuLike"]["CCother"][nuType].SetFillColor(ROOT.kBlue)
                hists1D["hErec"][hornMode]["SingleRingMuLike"]["NC"][nuType].SetFillColor(ROOT.kGreen)

                ErecStack.Add(hists1D["hErec"][hornMode]["SingleRingMuLike"]["CC0pi0p"][nuType])
                ErecStack.Add(hists1D["hErec"][hornMode]["SingleRingMuLike"]["CCother"][nuType])
                ErecStack.Add(hists1D["hErec"][hornMode]["SingleRingMuLike"]["NC"][nuType])

                c1 = ROOT.TCanvas()
                ErecStack.Draw()
                c1.Draw()

                hEMuPID = hists2D["hNLLReMu"][hornMode]["SingleRingMuLike"]["CC0pi0p"][nuType].Clone()
                hEMuPID.Add(hists2D["hNLLReMu"][hornMode]["SingleRingMuLike"]["CCother"][nuType])
                hEMuPID.Add(hists2D["hNLLReMu"][hornMode]["SingleRingELike"]["CC0pi0p"][nuType])
                hEMuPID.Add(hists2D["hNLLReMu"][hornMode]["SingleRingELike"]["CCother"][nuType])

                hEMuPID.SetMarkerColor(ROOT.kRed)

                c2 = ROOT.TCanvas()
                hEMuPID.Draw()
                c2.Draw()

                print "Press enter to quit"
                raw_input() 
                # End plots

            # Fill cut flow histograms and print out summary
            for sample in eventSelector.samples :
                sample.printSummary()
                for mode in intModeSelector.modes :
                    combinedName = mode.name+"_"+nuType
                    cutCount = 1
                    if combinedName in sample.countByIntMode :
                        for cut in sample.cutOrder:
                            hists1D["cutFlow_"+sample.name][hornMode][sample.name][mode.name][nuType].SetBinContent(cutCount, sample.countByIntMode[combinedName][cut])
                            cutCount += 1
            
            fOut.Write()

    return 

class baseSelector :

    def __init__ (self, name) :
        self.name = name
        self.cutOrder = []
        self.criteria = {}
        self.count = {}
        self.countByIntMode = {}
        self.count["All"] = 0

    def addCriterion (self, cutName, criterion) :
        self.cutOrder.append(cutName)
        self.criteria[cutName] = criterion
        self.count[cutName] = 0

    # If interaction mode is given *ONLY* interaction mode counters will be incremented
    def passCriteria(self, event, intModeName = None, nuType = None) :
        
        combinedName = None
        if intModeName != None and nuType != None :
            combinedName = intModeName+"_"+nuType
            if combinedName not in self.countByIntMode :

                self.countByIntMode[combinedName] = {}
                self.countByIntMode[combinedName]["All"] = 0
                for cutName in self.cutOrder :
                    self.countByIntMode[combinedName][cutName] = 0

        if combinedName :
            self.countByIntMode[combinedName]["All"] += 1
        self.count["All"] += 1
        
        for cutName in self.cutOrder :
            if self.criteria[cutName] (event) :
                if combinedName :
                    self.countByIntMode[combinedName][cutName] += 1
                self.count[cutName] += 1
            else :
                return False
        return True

    def printSummary(self) :
        nCols = len(self.countByIntMode)+2
        colWidth = 16
        width = nCols*colWidth
        print str.center("", width, '#')
        print "#"+str.center(self.name, width-2)+"#"
        print str.center("", width, '#')

        line = ""
        line += str.ljust("Cut", colWidth)
        line += str.rjust("All modes", colWidth)
        for mode in self.countByIntMode :
            line += str.rjust(mode, colWidth)
        print line
        print str.center("", width, "-")
        
        for criterion in ["All"]+self.cutOrder :
            line = ""
            line += str.ljust(criterion, colWidth)
            line += str.rjust(str(self.count[criterion]),colWidth/2)+str.rjust("%.2f" % (self.count[criterion]/float(self.count["All"])*100.) ,colWidth/2)
            for mode in self.countByIntMode :
                line += str.rjust(str(self.countByIntMode[mode][criterion]),colWidth/2)+str.rjust( "%.2f" % (self.countByIntMode[mode][criterion]/float(self.countByIntMode[mode]["All"])*100.),colWidth/2)
            print line
        print ""
    
class eventSelector (baseSelector) :

    samples = []

    def __init__ (self, name) :
        baseSelector.__init__(self, name)
        self.samples.append(self)
                
class intModeSelector (baseSelector) :

    modes = []

    def __init__ (self, name) :
        baseSelector.__init__(self, name)
        self.modes.append(self)

def bookHists (histList, nDim) :
    
    hists = {}
    for hist in histList :
        hists[hist] = {}
        for hornMode in hornCurrentModes :
            hists[hist][hornMode] = {}
            for sample in eventSelector.samples :
                hists[hist][hornMode][sample.name] = {}
                for interactionMode in intModeSelector.modes :
                    hists[hist][hornMode][sample.name][interactionMode.name] = {} 
                    for nuType in nuTypes :
                        histName = hist+"_"+hornMode+"_"+sample.name+"_"+interactionMode.name+"_"+nuType
                        if nDim == 1:
                            hists[hist][hornMode][sample.name][interactionMode.name][nuType] = ROOT.TH1D(histName, histName, histList[hist][0], histList[hist][1], histList[hist][2])
                        elif nDim == 2:
                            hists[hist][hornMode][sample.name][interactionMode.name][nuType] = ROOT.TH2D(histName, histName, histList[hist][0], histList[hist][1], histList[hist][2], histList[hist][3], histList[hist][4], histList[hist][5])
                        else :
                            print "bookHists ERROR: Invalid number of dimentions", nDim, "Quitting"
                            exit(-11)
    return hists

def bookHists (histList, nDim, hornMode, nuType) :
    hists = {}
    for hist in histList :
        hists[hist] = {}
        hists[hist][hornMode] = {}
        for sample in eventSelector.samples :
            hists[hist][hornMode][sample.name] = {}
            for interactionMode in intModeSelector.modes :
                hists[hist][hornMode][sample.name][interactionMode.name] = {} 
                histName = hist+"_"+hornMode+"_"+sample.name+"_"+interactionMode.name+"_"+nuType
                if nDim == 1:
                    hists[hist][hornMode][sample.name][interactionMode.name][nuType] = ROOT.TH1D(histName, histName, histList[hist][0], histList[hist][1], histList[hist][2])
                elif nDim == 2:
                    hists[hist][hornMode][sample.name][interactionMode.name][nuType] = ROOT.TH2D(histName, histName, histList[hist][0], histList[hist][1], histList[hist][2], histList[hist][3], histList[hist][4], histList[hist][5])
                else :
                    print "bookHists ERROR: Invalid number of dimentions", nDim, "Quitting"
                    exit(-11)
    return hists

def merge(source, destination):
    for key, value in source.items():
        if isinstance(value, dict):
            # get node or create one
            node = destination.setdefault(key, {})
            merge(value, node)
        else:
            destination[key] = value
            
    return destination

if __name__ == '__main__' :
    main()
