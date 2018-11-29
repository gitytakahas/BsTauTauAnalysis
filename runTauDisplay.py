import os, math, sys
from ROOT import TFile, TH1F, gROOT, TTree
import numpy as num
from DataFormats.FWLite import Events, Handle
from DeltaR import *

class GenTau(object):
    def __init__(self, p4):
        self.p4 = p4

    def __getattr__(self,name):
        return getattr(self.p4, name)



def daughters(p):
    ds = []
    for i_d in xrange(p.numberOfDaughters()):
        ds.append(p.daughter(i_d))
    return ds

def getNChargedAndNeutral(ds, n_charged=0, n_pizero=0):
    for d in ds:
        if d.charge():
            n_charged += 1
        elif d.pdgId() == 111:
            n_pizero += 1
        else:
            drs = daughters(d)
            if drs:
                n_charged, n_pizero = getNChargedAndNeutral(drs, n_charged, n_pizero)
            else:
                if d.pdgId() not in [310, 130, 22]:
                    pass
                    # import pdb; pdb.set_trace()

    return n_charged, n_pizero



def getNChargedAndNeutral_dr(ds, n_charged=0, n_pizero=0):

    charged_eta = 0
    charged_phi = 0
    charged_p4= 0

    neutral_eta = 0
    neutral_phi = 0
    neutral_p4= 0

    for d in ds:
        if d.charge():
            n_charged += 1
            charged_eta = d.eta()
            charged_phi = d.phi()
            charged_p4 = d.p4()
            
        elif d.pdgId() == 111:
            n_pizero += 1
            neutral_eta = d.eta()
            neutral_phi = d.phi()
            neutral_p4 = d.p4()

        else:
            drs = daughters(d)
            if drs:
                n_charged, n_pizero = getNChargedAndNeutral_dr(drs, n_charged, n_pizero)
            else:
                if d.pdgId() not in [310, 130, 22]:
                    pass
                    # import pdb; pdb.set_trace()

        dR2 = deltaR2(charged_eta,
                      charged_phi,
                      neutral_eta,
                      neutral_phi)
#
        print 'mass=', (charged_p4 + neutral_p4).mass()


    return dR2







gROOT.SetBatch(True)


tauH = Handle('vector<pat::Tau>')
muonH = Handle('vector<pat::Muon>')
vertexH = Handle('std::vector<reco::Vertex>')
genParticlesH  = Handle ('std::vector<reco::GenParticle>')
jetH = Handle('vector<pat::Jet>')
metH = Handle('vector<pat::MET>')


argvs = sys.argv
argc = len(argvs)

print argc, argvs

runtype = ''

if argc == 2:
#    print 'Please specify the runtype : python runTauDisplay.py <ZTT, ZEE, ZMM, QCD>'
#    sys.exit(0)

    runtype = argvs[1]




filelist = []

if runtype == 'Bs':
    filelist = [
        'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/cgalloni/PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test/miniAOD/miniAOD_PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test_1.root',
        'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/cgalloni/PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test/miniAOD/miniAOD_PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test_2.root',
        'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/cgalloni/PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test/miniAOD/miniAOD_PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test_3.root',
        'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/cgalloni/PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test/miniAOD/miniAOD_PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test_4.root',
        'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/cgalloni/PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test/miniAOD/miniAOD_PYTHIA8_BsPhiKKtautau_13TeV_2018_11_27_test_5.root',
        ]
    
else:
    filelist = [
        'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/GluGluHToTauTau_M125_13TeV_powheg_pythia8/1C6F3F7F-96C8-E611-A0D7-0025905A4964.root'
        ]


events = Events(filelist)
print len(filelist), 'files will be analyzed'



outputname = 'root/Myroot' + runtype + '.root'
file = TFile(outputname, 'recreate')

tree = TTree('tree', 'tree')


gen_pt = num.zeros(1, dtype=float)
gen_eta = num.zeros(1, dtype=float)
gen_phi = num.zeros(1, dtype=float)
gen_dm = num.zeros(1, dtype=float)
gen_ch_nh_dr = num.zeros(1, dtype=float)
gen_ncharged = num.zeros(1, dtype=int)
gen_npizeros = num.zeros(1, dtype=int)

pt = num.zeros(1, dtype=float)
eta = num.zeros(1, dtype=float)
phi = num.zeros(1, dtype=float)
mass = num.zeros(1, dtype=float)
rawMVAoldDM = num.zeros(1, dtype=float)
min_dr = num.zeros(1, dtype=float)
chad_min_dr = num.zeros(1, dtype=float)

iso_tight = num.zeros(1, dtype=int)
iso_loose = num.zeros(1, dtype=int)
iso_medium = num.zeros(1, dtype=int)
iso_vtight = num.zeros(1, dtype=int)
iso_vloose = num.zeros(1, dtype=int)


q = num.zeros(1, dtype=int)
decayMode = num.zeros(1, dtype=int)


againstElectronLooseMVA6 = num.zeros(1, dtype=int)
againstElectronMediumMVA6 = num.zeros(1, dtype=int)
againstElectronTightMVA6 = num.zeros(1, dtype=int)
againstElectronVLooseMVA6 = num.zeros(1, dtype=int)
againstElectronVTightMVA6 = num.zeros(1, dtype=int)
againstMuonLoose3 = num.zeros(1, dtype=int)
againstMuonTight3 = num.zeros(1, dtype=int)
chargedIsoPtSum = num.zeros(1, dtype=float)
neutralIsoPtSum = num.zeros(1, dtype=float)
puCorrPtSum = num.zeros(1, dtype=float)


tree.Branch('gen_pt', gen_pt, 'gen_pt/D')
tree.Branch('gen_eta', gen_eta, 'gen_eta/D')
tree.Branch('gen_phi', gen_phi, 'gen_phi/D')
tree.Branch('gen_dm', gen_dm, 'gen_dm/D')
tree.Branch('gen_ch_nh_dr', gen_ch_nh_dr, 'gen_ch_nh_dr/D')
tree.Branch('gen_ncharged', gen_ncharged, 'gen_ncharged/I')
tree.Branch('gen_npizeros', gen_npizeros, 'gen_npizeros/I')

tree.Branch('pt', pt, 'pt/D')
tree.Branch('eta', eta, 'eta/D')
tree.Branch('phi', phi, 'phi/D')
tree.Branch('mass', mass, 'mass/D')
tree.Branch('rawMVAoldDM', rawMVAoldDM, 'rawMVAoldDM/D')
tree.Branch('min_dr', min_dr, 'min_dr/D')
tree.Branch('chad_min_dr', chad_min_dr, 'chad_min_dr/D')

tree.Branch('iso_tight', iso_tight, 'iso_tight/I')
tree.Branch('iso_loose', iso_loose, 'iso_loose/I')
tree.Branch('iso_medium', iso_medium, 'iso_medium/I')
tree.Branch('iso_vtight', iso_vtight, 'iso_vtight/I')
tree.Branch('iso_vloose', iso_vloose, 'iso_vloose/I')


tree.Branch('q', q, 'q/I')
tree.Branch('decayMode', decayMode, 'decayMode/I')

tree.Branch('againstElectronLooseMVA6',
            againstElectronLooseMVA6, 'againstElectronLooseMVA6/I')
tree.Branch('againstElectronMediumMVA6',
            againstElectronMediumMVA6, 'againstElectronMediumMVA6/I')
tree.Branch('againstElectronTightMVA6',
            againstElectronTightMVA6, 'againstElectronTightMVA6/I')
tree.Branch('againstElectronVLooseMVA6',
            againstElectronVLooseMVA6, 'againstElectronVLooseMVA6/I')
tree.Branch('againstElectronVTightMVA6',
            againstElectronVTightMVA6, 'againstElectronVTightMVA6/I')
tree.Branch('againstMuonLoose3', againstMuonLoose3,
            'againstMuonLoose3/I')
tree.Branch('againstMuonTight3', againstMuonTight3,
            'againstMuonTight3/I')
tree.Branch('chargedIsoPtSum', chargedIsoPtSum, 'chargedIsoPtSum/D')
tree.Branch('neutralIsoPtSum', neutralIsoPtSum, 'neutralIsoPtSum/D')
tree.Branch('puCorrPtSum', puCorrPtSum, 'puCorrPtSum/D')



evtid = 0

for event in events:
    
    evtid += 1  
    
    if evtid%1000 == 0:
        print 'Event ', evtid, 'processed'

    if evtid > 100:
        pass
#        break


    event.getByLabel("slimmedMuons", muonH)

    if runtype=='Bs':
        event.getByLabel("slimmedTausLowPt", tauH)
    else:
        event.getByLabel("slimmedTaus", tauH)

    event.getByLabel("offlineSlimmedPrimaryVertices", vertexH)
    event.getByLabel('prunedGenParticles',genParticlesH)
    event.getByLabel("slimmedJets", jetH)
    event.getByLabel("slimmedMETs", metH)

    taus = tauH.product()
    muons = muonH.product()
    vertices = vertexH.product()
    gens = genParticlesH.product()
    jets = [jet for jet in jetH.product() if jet.pt() > 20 and abs(jet.eta()) < 2.3]
    met = metH.product()

    
    vis_had_taus = []
    for p in gens:
        if abs(p.pdgId()) != 15: continue

        ds = daughters(p)

#        print 'status=', p.status()


        # check mothers
        if p.numberOfMothers()==1:

            parent = p.mother(0)
            
            while abs(parent.pdgId()) == 15:
                parent = parent.mother(0)

#            if abs(parent.pdgId()) in [25,531]: 
            if abs(parent.pdgId()) in [25, 531]: 
                pass
            else:
                print 'mother is not from H nore Bs', parent.pdgId()
                continue
        else:
            print '# of mother is not 1', parent.numberOfMothers()
            continue



        # remove leptonic decays
        if any(abs(d.pdgId()) in [11, 13] for d in ds): 
            continue

        ds = [d for d in ds if abs(d.pdgId()) not in [12, 14, 16]] # keep visible products


        if not ds: continue

        vis_had_taus.append(GenTau(ds[-1].p4()))
        vis_tau = vis_had_taus[-1]
        for i_d in xrange(len(ds) - 1):
            vis_tau.p4 += ds[i_d].p4()
                
        n_charged, n_pizero = getNChargedAndNeutral(ds)

        if n_charged == 1:
            vis_tau.dm = n_pizero
        elif n_charged == 3:
            vis_tau.dm = 10 + n_pizero
        else:
            print 'Weird gen DM:', n_charged, n_pizero
            vis_tau.dm = 100
                    # import pdb; pdb.set_trace()

        vis_tau.charged = [d for d in ds if d.charge()]
#        vis_tau.pizeros = [daughters(d) for d in ds if d.pdgId()==111]
        vis_tau.pizeros = [d for d in ds if abs(d.pdgId())==111]

#        print vis_tau.pizeros, vis_tau.charged
#        import pdb; pdb.set_trace()        
#        if vis_tau.pizeros:
#            pass



    if len(vis_had_taus)==0: continue

    for i_gen_tau, gen_tau in enumerate(vis_had_taus):
#        if gen_tau.pt() > 25. and abs(gen_tau.eta()) < 2.3:
#        if abs(gen_tau.eta()) < 2.3:
#        if abs(gen_tau.eta()) < 99:
        gen_pt[0] = gen_tau.pt()
        gen_eta[0] = gen_tau.eta()
        gen_phi[0] = gen_tau.phi()
        gen_dm[0] = gen_tau.dm
        gen_ncharged[0] = len(gen_tau.charged)
        gen_npizeros[0] = len(gen_tau.pizeros)
#        import pdb; pdb.set_trace()

        chad = gen_tau.charged[-1].p4()
        for i_d in xrange(len(gen_tau.charged) - 1):
            chad += gen_tau.charged[i_d].p4()

        nhad = None


        if gen_tau.pizeros:
            nhad = gen_tau.pizeros[-1].p4()
            for i_d in xrange(len(gen_tau.pizeros) - 1):
                nhad += gen_tau.pizeros[i_d].p4()
                

        if nhad!= None:
            gen_ch_nh_dr[0] = returndR(chad.eta(), 
                                       chad.phi(),
                                       nhad.eta(),
                                       nhad.phi())
        else:
            gen_ch_nh_dr[0] = -1. 

        _min_dr = 999

#        print '-'*80

        for tau in taus:

            print 'reco. pt =', tau.pt(), 'dR =', deltaR(gen_tau, tau)

            dR = deltaR(gen_tau, tau)

            if  _min_dr > dR :
                
                pt[0] = tau.pt()
                eta[0] = tau.eta()
                phi[0] = tau.phi()
                mass[0] = tau.mass()
                rawMVAoldDM[0] = tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw")

                iso_tight[0] = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT")
                iso_loose[0] = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT")
                iso_medium[0] = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT")
                iso_vloose[0] = tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT")
                iso_vtight[0] = tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT")

                q[0] = tau.charge()
                decayMode[0] = tau.decayMode()

                againstElectronLooseMVA6[0] = tau.tauID('againstElectronLooseMVA6')
                againstElectronMediumMVA6[0] = tau.tauID('againstElectronMediumMVA6')
                againstElectronTightMVA6[0] = tau.tauID('againstElectronTightMVA6')
                againstElectronVLooseMVA6[0] = tau.tauID('againstElectronVLooseMVA6')
                againstElectronVTightMVA6[0] = tau.tauID('againstElectronVTightMVA6')
                againstMuonLoose3[0] = tau.tauID('againstMuonLoose3')
                againstMuonTight3[0] = tau.tauID('againstMuonTight3')
                chargedIsoPtSum[0] = tau.tauID('chargedIsoPtSum')
                neutralIsoPtSum[0] = tau.tauID('neutralIsoPtSum')
                puCorrPtSum[0] = tau.tauID('puCorrPtSum')
                min_dr[0] = dR



                chad_dR = returndR(chad.eta(),
                                   chad.phi(),
                                   tau.eta(),
                                   tau.phi())
                
                

                chad_min_dr[0] = chad_dR


                _min_dr = dR


        if _min_dr != 999:
            tree.Fill()


print evtid, 'events are processed !'

file.Write()
file.Close()

