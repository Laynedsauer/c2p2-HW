
import uproot
import awkward as ak
import numpy as np
import cupy as cp
from numba import cuda, jit
import matplotlib.pyplot as plt
file = uproot.open("/lstr/sahara/act/data/DAOD_PHYSLITE.37620644._000012.pool.root.1")
tree = file['CollectionTree;1']
Electrons = tree.arrays(["AnalysisElectronsAuxDyn.charge", "AnalysisElectronsAuxDyn.pt","AnalysisElectronsAuxDyn.eta", "AnalysisElectronsAuxDyn.phi","AnalysisElectronsAuxDyn.m"])
Muons = tree.arrays(["AnalysisMuonsAuxDyn.charge", "AnalysisMuonsAuxDyn.pt","AnalysisMuonsAuxDyn.eta", "AnalysisMuonsAuxDyn.phi"])
Taus = tree.arrays([ "AnalysisTauJetsAuxDyn.pt","AnalysisTauJetsAuxDyn.eta", "AnalysisTauJetsAuxDyn.phi", "AnalysisTauJetsAuxDyn.m"])
Jets = tree.arrays([ "AnalysisJetsAuxDyn.pt","AnalysisJetsAuxDyn.eta", "AnalysisJetsAuxDyn.phi", "AnalysisJetsAuxDyn.m"])
ps=tree.arrays(["BTagging_AntiKt4EMPFlowAuxDyn.DL1dv01_pb", "BTagging_AntiKt4EMPFlowAuxDyn.DL1dv01_pc", "BTagging_AntiKt4EMPFlowAuxDyn.DL1dv01_pu"])

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


tau = ak.zip({
    "pt" : Taus["AnalysisTauJetsAuxDyn.pt"],
    "eta" : Taus["AnalysisTauJetsAuxDyn.eta"],
    "phi" : Taus["AnalysisTauJetsAuxDyn.phi"],
    "mass" : Taus["AnalysisTauJetsAuxDyn.m"]
})

jets = ak.zip({
    "pt" : Jets["AnalysisJetsAuxDyn.pt"],
    "eta" : Jets["AnalysisJetsAuxDyn.eta"],
    "phi" : Jets["AnalysisJetsAuxDyn.phi"],
    "mass" : Jets["AnalysisJetsAuxDyn.m"]
})


p = ak.zip({
    "b" : ps["BTagging_AntiKt4EMPFlowAuxDyn.DL1dv01_pb"],
    "c" : ps["BTagging_AntiKt4EMPFlowAuxDyn.DL1dv01_pc"],
    "u" : ps["BTagging_AntiKt4EMPFlowAuxDyn.DL1dv01_pu"],
})

ak.to_backend(e, "cuda")
ak.to_backend(mu, "cuda")
ak.to_backend(tau, "cuda")
ak.to_backend(jets, "cuda")
ak.to_backend(p, "cuda")

jetspt=ak.ArrayBuilder()
eeta=[]
ept=[]
mupt=[]

num_events = len(mu.pt) # fix me
print("Number of events total:", num_events)
# print(ak.num(jets.eta, axis=-1)[:10])




#tried setting function to run with GPU, but couldn't get it to work

#@cuda.jit('float32')
def selection_cut(events):

    #first cut checking leptons
    lepton_cut=(ak.any(e.pt>30000, axis=-1) & ak.any(np.abs(e.eta<2.1), axis=-1)) | (ak.any(mu.pt>30000, axis=-1) & ak.any(np.abs(mu.eta)<2.1, axis=-1)) | (ak.any(tau.pt>30000, axis=-1) & ak.any(np.abs(tau.eta)<2.1, axis=-1))  

    
    #jet cut
    jets_cut=ak.num(jets.eta, axis=-1)>=4 & ak.all(jets.pt>30000) & ak.all(np.abs(jets.eta)<2.4)

    #combine cuts
    total_cut=np.logical_and(lepton_cut, jets_cut)
    cutevents=events[total_cut]

    
     #histograms after lepton and jet cut
    print("first cut number of events ", len(events.pt[total_cut]))
    plt.hist(ak.sum(cutevents.pt, axis=-1), bins=100)
    plt.title("Histogram of first cut")
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.show()

    #  btag math
    f_c=0.018
    D_DL1=np.log(p.b/(f_c*p.c+(1-f_c)*p.c))


    #  b tag cut, checks each value in D_DL1, if greater than 2.4, add 1 to "cut",
    
    #  if in the last element of the event, append either True or False to btagcut depending if cut >=2
    #print(D_DL1)
    # print(D_DL1[0][2])
    btagcut=[]
    cut=0
    l=0
    for i in range(len(D_DL1)):
        # print(len(D_DL1[i]))
        for l in range(len(D_DL1[i])):
            # print(l)
            # cut[i][l]=D_DL1[i][l]>2.4
            # print(np.sum(int(D_DL1[0][0]>2.4))>=2)
            if D_DL1[i][l]>2.4:
                cut=cut+1
            if l>=len(D_DL1[i])-1 and cut>=2:
                btagcut.append(True)
                cut=0
                # print("btagtrue")
                
                
            elif l>=len(D_DL1[i])-1:
                btagcut.append(False)
                cut=0
                # print("btagfalse")
            
                
    finalcut=np.logical_and(total_cut,btagcut)
    print('number of events after all cuts is ', len(events[finalcut]))

    finalcutevents=events[finalcut]
    #  Histograms of pts after all cuts
    plt.hist(ak.sum(finalcutevents.pt, axis=-1), bins=100)
    # print("Final cut number of events ", len(events.pt[finalcut]))
    plt.title("Histogram of finalcut events")
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.show()
print("Selection and histogram of e events")
selection_cut(e)
print("Selection and histogram of mu events")
selection_cut(mu)
print("Selection and histogram of jet events")
selection_cut(jets)    
print("Selection and histogram of tau events")
selection_cut(tau)
#print(ak.sum(mu.pt,axis=-1))
