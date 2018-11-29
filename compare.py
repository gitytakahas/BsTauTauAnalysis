import copy, math, os
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TTree, gROOT, gStyle
from DisplayManager import DisplayManager
from officialStyle import officialStyle

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

#from optparse import OptionParser, OptionValueError 
#usage = "usage: python BrazilianPlots.py [signal] [var]" 
#parser = OptionParser(usage) 
#
#parser.add_option(
#    "-c", "--channel",
#    default="Single",
#    type="string",
#    dest="channel"
#)
#
#(options, args) = parser.parse_args() 


colours = [1, 2, 4, 6, 8, 13, 15]
styles = [1, 2, 4, 3, 5, 1, 1]

def applyHistStyle(h, i):
    print h, i
    h.SetLineColor(colours[i])
    h.SetMarkerColor(colours[i])
    h.SetMarkerSize(0)
    h.SetLineStyle(styles[i])
    h.SetLineWidth(3)
    h.SetStats(False)

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def comparisonPlots(hists, titles, isLog=False, pname='sync.pdf', isRatio=False, isLegend=True):

    display = DisplayManager(pname, isLog, isRatio, 0.6, 0.7)
    display.draw_legend = isLegend

    display.Draw(hists, titles)


def sproducer(key, rootfile, name, ivar):

    hist = TH1F('h_' + key + '_' + name, 
                'h_' + key + '_' + name, 
                ivar['nbin'], ivar['xmin'], ivar['xmax'])

    hist.Sumw2()
    exp = '(' + ivar['sel'] + ')*' + ivar['weight']
        
    tree = rootfile.Get(ivar['tree'])

    print ivar['var'] + ' >> ' + hist.GetName(), exp
    
    tree.Draw(ivar['var'] + ' >> ' + hist.GetName(), exp)
    hist.GetXaxis().SetTitle(ivar['title'])
    hist.GetYaxis().SetTitle('a.u.')
        
    return copy.deepcopy(hist)

vardict = {

    'pt':{'tree':'tree', 'var':'pt', 'nbin':50, 'xmin':0, 'xmax':100, 'title':"tau pT (GeV)", 'sel':'1', 'weight':'1'},
    'eta':{'tree':'tree', 'var':'eta', 'nbin':30, 'xmin':-4., 'xmax':4., 'title':"tau eta", 'sel':'1', 'weight':'1'},
    'gen_pt':{'tree':'tree', 'var':'gen_pt', 'nbin':50, 'xmin':0, 'xmax':100, 'title':"gen. vis. tau pT (GeV)", 'sel':'1', 'weight':'1'},
    'gen_eta':{'tree':'tree', 'var':'gen_eta', 'nbin':30, 'xmin':-4., 'xmax':4., 'title':"gen. vis. tau eta", 'sel':'1', 'weight':'1'},
    'decayMode':{'tree':'tree', 'var':'decayMode', 'nbin':11, 'xmin':0., 'xmax':11., 'title':"reco. decay mode", 'sel':'1', 'weight':'1'},
    'gen_dm':{'tree':'tree', 'var':'gen_dm', 'nbin':11, 'xmin':0., 'xmax':11., 'title':"gen. decay mode", 'sel':'1', 'weight':'1'},
    'q':{'tree':'tree', 'var':'q', 'nbin':3, 'xmin':-1., 'xmax':2., 'title':"tau charge", 'sel':'1', 'weight':'1'},
    'againstElectronLooseMVA6':{'tree':'tree', 'var':'againstElectronLooseMVA6', 'nbin':2, 'xmin':0, 'xmax':2., 'title':"againstElectronLooseMVA6", 'sel':'1', 'weight':'1'},
    'againstMuonLoose3':{'tree':'tree', 'var':'againstMuonLoose3', 'nbin':2, 'xmin':0, 'xmax':2., 'title':"againstMuonLoose3", 'sel':'1', 'weight':'1'},
    'iso_tight':{'tree':'tree', 'var':'iso_tight', 'nbin':2, 'xmin':0, 'xmax':2., 'title':"MVA Iso. tight", 'sel':'1', 'weight':'1'},
    'chargedIsoPtSum':{'tree':'tree', 'var':'chargedIsoPtSum', 'nbin':30, 'xmin':0, 'xmax':50, 'title':"chargedIsoPtSum", 'sel':'1', 'weight':'1'},
    'neutralIsoPtSum':{'tree':'tree', 'var':'neutralIsoPtSum', 'nbin':30, 'xmin':0, 'xmax':50, 'title':"neutralIsoPtSum", 'sel':'1', 'weight':'1'},
    'puCorrPtSum':{'tree':'tree', 'var':'puCorrPtSum', 'nbin':30, 'xmin':0, 'xmax':200, 'title':"puCorrPtSum", 'sel':'1', 'weight':'1'},
    'min_dR':{'tree':'tree', 'var':'min_dr', 'nbin':30, 'xmin':0, 'xmax':2*math.pi, 'title':"min. dR (gen, reco)", 'sel':'1', 'weight':'1'},
    'min_dR_zoom':{'tree':'tree', 'var':'min_dr', 'nbin':30, 'xmin':0, 'xmax':0.2, 'title':"min. dR (gen, reco): zoom", 'sel':'1', 'weight':'1'},
    'chad_dR':{'tree':'tree', 'var':'chad_min_dr', 'nbin':30, 'xmin':0, 'xmax':2*math.pi, 'title':"min. dR (gen. ch, reco)", 'sel':'1', 'weight':'1'},
    'gen_ncharged':{'tree':'tree', 'var':'gen_ncharged', 'nbin':10, 'xmin':0, 'xmax':10, 'title':"# of chad", 'sel':'1', 'weight':'1'},
    'gen_npizeros':{'tree':'tree', 'var':'gen_npizeros', 'nbin':10, 'xmin':0, 'xmax':10, 'title':"# of nhad", 'sel':'1', 'weight':'1'},
    'resolution':{'tree':'tree', 'var':'(pt - gen_pt)/gen_pt', 'nbin':30, 'xmin':-10, 'xmax':10, 'title':"(reco. pT - gen.)/gen. (dR < 0.3)", 'sel':'min_dr < 0.3', 'weight':'1'},
    'gen_ch_nh_dr':{'tree':'tree', 'var':'gen_ch_nh_dr', 'nbin':30, 'xmin':0, 'xmax':2*math.pi, 'title':"gen. #DeltaR(ch, nh) (#nh!=0)", 'sel':'gen_npizeros!=0', 'weight':'1'},
    }




ensureDir('Plots')

for vkey, ivar in vardict.iteritems():

    hists = []
    titles = []

    print vkey, ivar

    for runtype in ['HTT', 'Bs']:
    
        sfile = TFile('root/Myroot' + runtype + '.root')

        hist_madgraph = sproducer('hist_' + runtype, sfile, vkey, ivar)
        hists.append(copy.deepcopy(hist_madgraph))
        titles.append(str(runtype))

   
    for ii, ihist in enumerate(hists):
        applyHistStyle(ihist, ii)

        ihist.Scale(1./ihist.GetSumOfWeights())
        ihist.SetMaximum(ihist.GetBinContent(ihist.GetMaximumBin())*1.2)

#        comparisonPlots(hists, titles, False, 'Pair_Plots_lo_nlo/' + str(runtype) + 'GeV/compare_' + str(runtype) + '_' + vkey + '.pdf')

    comparisonPlots(hists, titles, False, 'Plots/' + vkey + '.pdf')


