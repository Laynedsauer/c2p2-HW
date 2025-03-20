import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

file = uproot.open("C:/Users/layne/miniconda3/c2p2HW1data/DAOD_PHYSLITE.37620644._000012.pool.root.1")

tree = file['CollectionTree;1']
Electrons = tree.arrays(["AnalysisElectronsAuxDyn.charge", "AnalysisElectronsAuxDyn.pt","AnalysisElectronsAuxDyn.eta", "AnalysisElectronsAuxDyn.phi","AnalysisElectronsAuxDyn.m"])
Muons = tree.arrays(["AnalysisMuonsAuxDyn.charge", "AnalysisMuonsAuxDyn.pt","AnalysisMuonsAuxDyn.eta", "AnalysisMuonsAuxDyn.phi"])

e = ak.zip({
    "pt" : Electrons["AnalysisElectronsAuxDyn.pt"],
    "eta" : Electrons["AnalysisElectronsAuxDyn.eta"],
    "phi" : Electrons["AnalysisElectronsAuxDyn.phi"],
    "mass" : Electrons["AnalysisElectronsAuxDyn.m"],
    "charge" : Electrons["AnalysisElectronsAuxDyn.charge"]
})
mu = ak.zip({
    "pt" : Muons["AnalysisMuonsAuxDyn.pt"],
    "eta" : Muons["AnalysisMuonsAuxDyn.eta"],
    "phi" : Muons["AnalysisMuonsAuxDyn.phi"],
    "charge" : Muons["AnalysisMuonsAuxDyn.charge"],
    "mass" : 0
})



#  Momentum components
def Px(pt, phi):
    return pt*np.cos(phi)

def Py(pt,phi):
    return pt*np.sin(phi)

def Pz(pt,eta):
    return pt*np.sinh(eta)

def E(pt, eta, m):
    return np.sqrt(pt**2 * np.cosh(eta)**2 + m**2)

e_px = Px(e.pt,e.phi)
e_py = Py(e.pt,e.phi)
e_pz = Pz(e.pt,e.eta)
e_E = E(e.pt,e.eta,e.mass)
emasses=[]
mumasses=[]
allmasses=[]
def withmass_calc_invm(p1, p2):
    sumpx = Px(p1.pt,p1.phi)+Px(p2.pt, p2.pt)
    sumpy =  Py(p1.pt,p1.phi)+Py(p2.pt, p2.pt)
    sumpz = Pz(p1.pt,p1.eta)+Pz(p2.pt, p2.eta)
    sumE = E(p1.pt,p1.eta, p1.mass)+E(p2.pt, p2.eta,p2.mass)
    
    invm = np.sqrt(sumE**2-sumpz**2-sumpx**2-sumpy**2) # fix me
    return invm/1000.0 # in GeV



def selection_calc(lepton):
    lepton_cut = (lepton.pt>25000) & (ak.num(lepton.charge) >=2 )
    cleaned_lepton = lepton[lepton_cut]
    h = ak.combinations(cleaned_lepton, 2, axis = 1)
    h_non_empty = h[ak.num(h, axis=1) > 0]
    p1temp,p2temp=ak.unzip(h)
    
    pairs_cut=(p1temp.charge!=p2temp.charge) #only keep pairs of oppositely charged particles
    pairs=h[pairs_cut]
    p1temp2,p2temp2=ak.unzip(pairs)
    #print(p1.charge[:5],p2.charge[:5])
    #print(pairs[:5])
    nonemptypairs=pairs[ak.num(p1temp2.charge)>0]
    p1,p2=ak.unzip(nonemptypairs)
    #print(nonemptypairs[:5])
    print(ak.sum(ak.num(cleaned_lepton)))
    print(ak.sum(ak.num(pairs)))
    print(p1.mass[:5])
    #print(len(pairs))
    for i in range(len(nonemptypairs)):
        allmasses.append(withmass_calc_invm(p1[i],p2[i]))



print('Calculations started')
selection_calc(e)
selection_calc(mu)
print('Calculations done')
#allmasses=emasses+mumasses
# %%


plt.hist(allmasses, bins=200, histtype="stepfilled")
plt.xlim(0,1000)
plt.title('Histogram of allmasses')
plt.xlabel('Value')
plt.ylabel('Frequency')

# Show the plot
plt.show()
print('hist plotted')
massaverage=ak.sum(allmasses)/len(allmasses)
print(massaverage , "GeV")
