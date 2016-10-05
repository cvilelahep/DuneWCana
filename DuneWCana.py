#!/us/bin/python

import ROOT
ROOT.gSystem.Load("/home/cvilela/WCSim/v151/WCSim/libWCSimRoot.so")

import numpy as np

verbosity = False

# Some definitions:
# hornCurrentModes = ["FHC", "RHC"]
hornCurrentModes = ["FHC"]

#nuTypes = {"nue"      : 12,
#           "numu"     : 14,
#           "nutau"    : 16,
#           "nuebar"   : -12,
#           "numubar"  : -14,
#           "nutaubar" : -16}

nuTypes = {"numu"     : 14}

baseFilePath = "/home/cvilela/DuneWC/"

nFiles = 5

OutputPath = baseFilePath+"Analysis/Out/"

# Incoming neutrino direction is along z. This should be corrected in a future iteration...
# Hopefully not a very significant problem with such a big tank. Look at efficiency vs angle (relative to tank) and momentum
trueNuDir = [0., 0., 1.]
# Detector radius and half-height
detR = 3700.
detHH = 3000.
dWallCut = 200.

def fqFilePath(nuType, hornMode, fileNumber) :
    return baseFilePath+"Data/gntp_"+hornMode+"_"+nuType+"_fluxosc_fQv5r4_"+str.rjust(str(fileNumber), 3, '0')+".root"

def wcFilePath(nuType, hornMode, fileNumber) :
    return baseFilePath+"Data/gntp_"+hornMode+"_"+nuType+"_fluxosc_"+str.rjust(str(fileNumber), 3, '0')+".root"

# don't think this will be necessary, but just in case
def vecFilePath(nuType, hornMode, fileNumber) :
    return baseFilePath+"Data/gntp_"+hornMode+"_"+nuType+"_fluxosc_"+str.rjust(str(fileNumber), 3, '0')+".nuance"

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
    return event.fq1rnll[0*7 + 1] - event.fq1rnll[0*7 + 2] >= -100

def nseMu (event) :
    return event.fqnse <= 2

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
    return event.fq1rnll[0*7 + 2] - event.fq1rnll[0*7 + 1] >= 100

def nseE (event) :
    return event.fqnse == 1

def pidEpi0 (event) :
    return event.fqpi0nll[0] - event.fq1rnll[0*7+1] >= 0 

# Interaction mode category criteria
def CC0pi0pcriterion(event) :

    isCC = True

    wcEv = event.wcsimrootevent
    if abs(wcEv.GetTrigger(0).GetMode()) >= 30 :
        isCC = False

    for par in range(0, wcEv.GetTrigger(0).GetNtrack()) :
#        print "Ipnu", wcEv.GetTrigger(0).GetTracks().At(par).GetIpnu(), "Flag", wcEv.GetTrigger(0).GetTracks().At(par).GetFlag(), "is piOrp", abs(wcEv.GetTrigger(0).GetTracks().At(par).GetIpnu()) in [111, 211, 2122] 
        if abs(wcEv.GetTrigger(0).GetTracks().At(par).GetIpnu()) in [111, 211, 2122] :
            if wcEv.GetTrigger(0).GetTracks().At(par).GetFlag() == 0 :
                isCC = False
                break 
    
    return isCC

def CCothercriterion(event) :

    isCCother = True

    wcEv = event.wcsimrootevent
    # Check that is CC and not CC0pi0pcriterion
    if abs(wcEv.GetTrigger(0).GetMode()) >= 30 :
        isCCother = False
        if verbosity :
            print "CCothercriterion:", "Failed is CC", abs(wcEv.GetTrigger(0).GetMode())
    if CC0pi0pcriterion(event) :
        isCCother = False
        if verbosity :
            print "CCothercriterion:", "Failed is not CC0pi0pCC"

    return isCCother

def NCcriterion(event) :
    wcEv = event.wcsimrootevent
    # Check that is CC and not CC0pi0pcriterion
    if abs(wcEv.GetTrigger(0).GetMode()) >= 30 :
        return True
    else :
        return False
        
def main() :

    hTrueNuDir = ROOT.TH3D("trueNuDir", "trueNuDir", 10, -1, 1, 10, -1, 1, 10, -1, 1)
    
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

    # Book histograms:
    #              Name       nBins  Min Max
    histList1D = {"hErec" : [100,   0., 10000]}

    histList2D = {"hPCosTheta" : [100, 0, 5000, 100, -1, 1],
                  "hNLLRePi0" : [100, 0, 300, 100, -5000, 5000],
                  "hNLLReMu"  : [100, 0, 5000, 100, -5000, 5000],
                  "hEtrueErec"  : [100, 0, 5000, 100, 0, 5000]}

    hists1D = bookHists(histList1D, 1)
    hists2D = bookHists(histList2D, 2)
   
    # File loop
    for hornMode in hornCurrentModes :
        for nuType in nuTypes :
            for fNum in range(0, nFiles) :
                # Open fiTQun file and get TTree
                fFQ = ROOT.TFile(fqFilePath(nuType = nuType, hornMode = hornMode, fileNumber = fNum))
                tree = fFQ.Get("fiTQun")

                # Add WCSim TTree as Friend
                tree.AddFriend("wcsimT", wcFilePath(nuType = nuType, hornMode = hornMode, fileNumber = fNum))

                for event in tree :

                    for sample in eventSelector.samples :
                        isSelected = sample.passCriteria(event)

                        for mode in intModeSelector.modes :
                            isMode =  mode.passCriteria(event)

                            if isMode :
                                sample.passCriteria(event, mode.name, nuType)

                            if isMode and isSelected :

                                fQindex = -1
                                m_l = -1
                                m_n = 939.6
                                m_p = 938.2
                                v_nuc = 27.
                                
                                if sample.name == "SingleRingELike" :
                                    fQindex = 1
                                    m_l = 0.511
                                elif sample.name == "SingleRingMuLike" :
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
                                
                                hists1D["hErec"][hornMode][sample.name][mode.name][nuType].Fill(Erec)
                                hists2D["hPCosTheta"][hornMode][sample.name][mode.name][nuType].Fill(p_l, cos_th_beam)
                                hists2D["hNLLRePi0"][hornMode][sample.name][mode.name][nuType].Fill(event.fqpi0mass[0], event.fqpi0nll[0] - event.fq1rnll[0*7+1])
                                hists2D["hNLLReMu"][hornMode][sample.name][mode.name][nuType].Fill(event.fq1rmom[0*7+1], event.fq1rnll[0*7+2] - event.fq1rnll[0*7+1])

                                hists2D["hEtrueErec"][hornMode][sample.name][mode.name][nuType].Fill(Erec, event.wcsimrootevent.GetTrigger(0).GetTracks().At(0).GetP())
                                
#                    print modeString
                fFQ.Close()

    # Draw some histograms
    c1 = ROOT.TCanvas()
    hists1D["hErec"]["FHC"]["SingleRingMuLike"]["CC0pi0p"]["numu"].SetFillColor(ROOT.kRed)
    hists1D["hErec"]["FHC"]["SingleRingMuLike"]["CCother"]["numu"].SetFillColor(ROOT.kBlue)
    hists1D["hErec"]["FHC"]["SingleRingMuLike"]["NC"]["numu"].SetFillColor(ROOT.kGreen)

    Erec_FHC_numu_1rmu = ROOT.THStack("Erec_FHC_numu_1rmu", "FHC numu Single-ring Mu-like;E_{rec} [MeV]")
    Erec_FHC_numu_1rmu.Add(hists1D["hErec"]["FHC"]["SingleRingMuLike"]["CC0pi0p"]["numu"])
    Erec_FHC_numu_1rmu.Add(hists1D["hErec"]["FHC"]["SingleRingMuLike"]["CCother"]["numu"])
    Erec_FHC_numu_1rmu.Add(hists1D["hErec"]["FHC"]["SingleRingMuLike"]["NC"]["numu"])
    Erec_FHC_numu_1rmu.Draw()
    
    c1.Draw()


    c2 = ROOT.TCanvas()
    hNLLReMu_numu = hists2D["hNLLReMu"]["FHC"]["SingleRingMuLike"]["CC0pi0p"]["numu"].Clone()
    hNLLReMu_numu.Add(hists2D["hNLLReMu"]["FHC"]["SingleRingELike"]["CC0pi0p"]["numu"])
    hNLLReMu_numu.Add(hists2D["hNLLReMu"]["FHC"]["SingleRingMuLike"]["CCother"]["numu"])
    hNLLReMu_numu.Add(hists2D["hNLLReMu"]["FHC"]["SingleRingELike"]["CCother"]["numu"])
    hNLLReMu_numu.Add(hists2D["hNLLReMu"]["FHC"]["SingleRingMuLike"]["NC"]["numu"])
    hNLLReMu_numu.Add(hists2D["hNLLReMu"]["FHC"]["SingleRingELike"]["NC"]["numu"])

    hNLLReMu_numu.SetTitle("FHC numu;p_{e} [MeV;]ln(L_{e}/L_{#mu})")
    hNLLReMu_numu.Draw("COLZ")
    c2.Draw()

    c3 = ROOT.TCanvas()
    hists2D["hEtrueErec"]["FHC"]["SingleRingMuLike"]["CC0pi0p"]["numu"].Draw("COLZ")
    c3.Draw()

    for sample in eventSelector.samples :
        sample.printSummary()
    
    print "Press enter to quit"
    raw_input()
    

class baseSelector :

    def __init__ (self, name) :
        self.name = name
        self.cutOrder = []
        self.criteria = {}
        self.count = {}
        self.countByIntMode = {}
#        self.cutOrder.append("All")
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
        else :
            self.count["All"] += 1
        
        for cutName in self.cutOrder :
            if self.criteria[cutName] (event) :
                if combinedName :
                    self.countByIntMode[combinedName][cutName] += 1
                else :
                    self.count[cutName] += 1
            else :
                return False
        return True

    def printSummary(self) :
        nCols = len(self.countByIntMode)+2
        colWidth = 30
        width = nCols*colWidth
        print str.center("", width, '#')
        print "#"+str.center(self.name, width-2)+"#"
        print str.center("", width, '#')

        line = ""
        line += str.ljust("Cut", colWidth)
        line += str.center("All modes", colWidth)
        for mode in self.countByIntMode :
            line += str.center(mode, colWidth)
        print line
        print str.center("", width, "-")
        
        for criterion in ["All"]+self.cutOrder :
            line = ""
            line += str.ljust(criterion, colWidth-2)
            line += str.rjust(str(self.count[criterion]),colWidth/2)+str.rjust(str(self.count[criterion]/float(self.count["All"])*100.),colWidth/2)+"  "
            for mode in self.countByIntMode :
                line += str.rjust(str(self.countByIntMode[mode][criterion]),colWidth/2)+str.rjust(str(self.countByIntMode[mode][criterion]/float(self.countByIntMode[mode]["All"])*100.),colWidth/2)+"  "
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
        
if __name__ == '__main__' :
    main()
