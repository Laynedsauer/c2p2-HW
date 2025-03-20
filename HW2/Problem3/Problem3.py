import uproot
import awkward as ak
import numpy as np
import cupy as cp
from numba import cuda, jit
import matplotlib.pyplot as plt
file = uproot.open("/lstr/sahara/act/data/TTto2L2Nu_NanoAODv12-130x_mcRun3_13p6TeV_powheg-pythia.root")
print(file.keys())
tree=file['Events;1']
#print(tree.keys())
triggerevents=tree.arrays(['HLT_Mu8'])

trigger = ak.zip({
    "cut" : triggerevents['HLT_Mu8']
})

array=(ak.ravel(triggerevents))
print('number of all events',len(triggerevents))
#print('number of triggered events',np.sum(array))
print('efficiency of trigger is therefore', round(np.sum(array)/len(triggerevents)*100, 3),'%')
mu=tree.arrays(['Muon_mvaLowPt', 'HLT_HighPtTkMu100', 'Muon_highPtId', 'L1_SingleMu8er1p5'])

mu = ak.zip({
    "lowpt" : mu['Muon_mvaLowPt'],
    "HLThighpt" : mu['HLT_HighPtTkMu100'],
    "highptId" : mu['Muon_highPtId'],
    "singmu8" : mu['L1_SingleMu8er1p5'],



})
# print('sum(ak.any(mu.HLThighpt,axis=-1)',sum(ak.count(mu.HLThighpt,axis=-1)))
# print('ak.any(mu.HLThighpt,axis=-1',sum(ak.any(mu.HLThighpt,axis=-1)))
# print(ak.count(mu.HLThighpt, axis=-1))
# print(mu.lowpt[mu.HLThighpt])
# print(trigger.cut)
# print(mu.lowpt)
# print(mu.lowpt[trigger.cut])
onemucut=np.logical_and(trigger.cut ,ak.count(mu.lowpt, axis=-1)==1)
# print(ak.count(mu.lowpt, axis=-1)==1)
# print('onemucut is ',onemucut)
print('efficiency of trigger for events with one muon', round(len(triggerevents[onemucut])/len(triggerevents)*100, 3),'%')
# print(mu.lowpt[onemucut])

twomucut=np.logical_and(trigger.cut, ak.count(mu.HLThighpt, axis=-1)==2)
# print('twomucut is ',sum(ak.count(mu.HLThighpt,axis=-1)==2))

print('efficiency of trigger for events with two muons', round(len(triggerevents[twomucut])/len(triggerevents)*100, 3),'%')
# print(triggerevents[:10])

plt.hist(mu.lowpt[onemucut])
plt.title("Histogram of mu.lowpt one mu cut")
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()



plt.hist(mu.lowpt[twomucut])
plt.title("Histogram of mu.lowpt two mu cut")
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()


print("I doubt I have done this right, but the cutoff for these histograms are -1,1 , so I would say this is when the trigger turns on")
