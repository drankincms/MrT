import ROOT as rt
import random
rt.gStyle.SetOptStat(0)

com = '250'
gLQ = 1.
ilumi = 250.
jetres = 0.03
metres = 0.05
elres = 0.01
mures = 0.01
muloss = 0.01
elloss = 0.01
#https://www.desy.de/f/students/2018/reports/NiallMcHugh.pdf

samples = {
    'qqvv':[0,r'e\mu \rightarrow qq\nu\nu',rt.kMagenta+2,'events_SM_qqvv_%sgev.root'%com,{'250':0.01904,'550':0.1537}],
    'LQ_m650':[10,r'LQ (M = 650 GeV)',rt.kBlue-7,'events_LQ_S1_M650_qq_%sgev.root'%com,{'250':0.0008647*gLQ*gLQ*gLQ*gLQ,'550':0.00228*gLQ*gLQ*gLQ*gLQ}],
    'LQ_m2000':[50,r'LQ (M = 2 TeV)',rt.kGreen+1,'events_LQ_S1_M2000_qq_%sgev.root'%com,{'250':1.16e-05*gLQ*gLQ*gLQ*gLQ,'550':5.227e-05*gLQ*gLQ*gLQ*gLQ}],
    'qqmv':[0,r'e\mu \rightarrow qql\nu',rt.kRed+2,'events_SM_qqmv_%sgev.root'%com,{'250':0.04965*2.,'550':0.205*2.}], #gen is for munu only, x2 for enu as well
}

sample_order = [
    'qqvv',
    'qqmv',
    'LQ_m650',
    'LQ_m2000',
]

files = {s:rt.TFile('%s'%samples[s][3]) for s in samples}
trees = {s:files[s].Get('t') for s in samples}

hists = {
    'pt0':{s:rt.TH1F('pt0_%s'%s,'pt0_%s'%s,50,0.,(float(com)/2.)+10.) for s in sample_order},
    'pt1':{s:rt.TH1F('pt1_%s'%s,'pt1_%s'%s,50,0.,(float(com)/2.)+10.) for s in sample_order},
    'eta0':{s:rt.TH1F('eta0_%s'%s,'eta0_%s'%s,60,-3.,3.) for s in sample_order},
    'eta1':{s:rt.TH1F('eta1_%s'%s,'eta1_%s'%s,60,-3.,3.) for s in sample_order},
    'theta0':{s:rt.TH1F('theta0_%s'%s,'theta0_%s'%s,40,-1.,1.) for s in sample_order},
    'theta1':{s:rt.TH1F('theta1_%s'%s,'theta1_%s'%s,40,-1.,1.) for s in sample_order},
    'mqq':{s:rt.TH1F('mqq_%s'%s,'mqq_%s'%s,50,0.,float(com)+10.) for s in sample_order},
    'met':{s:rt.TH1F('met_%s'%s,'met_%s'%s,50,0.,(float(com)/2.)+10.) for s in sample_order},
}

histnames = {
    'pt0':r'p_{T}(leading q)',
    'pt1':r'p_{T}(subleading q)',
    'eta0':r'\eta (leading j)',
    'eta1':'r\eta (subleading j)',
    'theta0':r'cos(\theta) (leading j)',
    'theta1':r'cos(\theta) (subleading j)',
    'mqq':r'm_{jj} [GeV]',
    'met':r'MET [GeV]',
}

for s in sample_order:
    print('starting %s'%s)
    iev = 0
    q1 = None
    q2 = None
    v1 = None
    v2 = None
    e1 = None
    e2 = None
    e3 = None
    m1 = None
    m2 = None
    m3 = None
    met = 0.
    for ev in trees[s]:
        if iev!=ev.event:
            iev = ev.event
            if v1 is not None:
                met = v1.Pt()
                if v2 is not None:
                    met = (v1+v2).Pt()*random.gauss(1., metres)
            hists['met'][s].Fill(met)
    
            if q2 is not None:
                lepton_veto = False
                if e1 is not None:
                    if random.uniform(0., 1.)<elloss:
                        lepton_veto = True
                if e2 is not None:
                    if random.uniform(0., 1.)<elloss:
                        lepton_veto = True
                if e3 is not None:
                    if random.uniform(0., 1.)<elloss:
                        lepton_veto = True
                if m1 is not None:
                    if random.uniform(0., 1.)<elloss:
                        lepton_veto = True
                if m2 is not None:
                    if random.uniform(0., 1.)<elloss:
                        lepton_veto = True
                if m3 is not None:
                    if random.uniform(0., 1.)<elloss:
                        lepton_veto = True
                if not lepton_veto:
                    q1.SetE(q1.E()*random.gauss(1., jetres))
                    q2.SetE(q2.E()*random.gauss(1., jetres))
                    hists['pt0'][s].Fill(q1.Pt())
                    hists['pt1'][s].Fill(q2.Pt())
                    hists['eta0'][s].Fill(q1.Eta())
                    hists['eta1'][s].Fill(q2.Eta())
                    hists['theta0'][s].Fill(rt.TMath.Cos(q1.Theta()))
                    hists['theta1'][s].Fill(rt.TMath.Cos(q2.Theta()))
                    hists['mqq'][s].Fill((q1+q2).M())
            
            q1 = None
            q2 = None
            v1 = None
            v2 = None
            e1 = None
            e2 = None
            e3 = None
            m1 = None
            m2 = None
            m3 = None
            met = 0.

        pid = ev.id
        pp4 = ev.p4

        if abs(pid)==12 or abs(pid)==14 or abs(pid)==16:
            if v1 is None:
                v1 = rt.TLorentzVector()
                v1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            elif v2 is None:
                v2 = rt.TLorentzVector()
                if pp4.Pt()>v1.Pt():
                    v2.SetPtEtaPhiM(v1.Pt(),v1.Eta(),v1.Phi(),v1.M())
                    v1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
                else:
                    v2.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            else:
                print('Huh? (v)')
        if abs(pid)<6 and pid!=0:
            if q1 is None:
                q1 = rt.TLorentzVector()
                q1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            elif q2 is None:
                q2 = rt.TLorentzVector()
                if pp4.Pt()>q1.Pt():
                    q2.SetPtEtaPhiM(q1.Pt(),q1.Eta(),q1.Phi(),q1.M())
                    q1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
                else:
                    q2.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            else:
                print('Huh? (q)')
        if abs(pid)==11:
            if e1 is None:
                e1 = rt.TLorentzVector()
                e1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            elif e2 is None:
                e2 = rt.TLorentzVector()
                if pp4.Pt()>e1.Pt():
                    e2.SetPtEtaPhiM(e1.Pt(),e1.Eta(),e1.Phi(),e1.M())
                    e1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
                else:
                    e2.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            elif e3 is None:
                e3 = rt.TLorentzVector()
                if pp4.Pt()>e1.Pt():
                    e3.SetPtEtaPhiM(e2.Pt(),e2.Eta(),e2.Phi(),e2.M())
                    e2.SetPtEtaPhiM(e1.Pt(),e1.Eta(),e1.Phi(),e1.M())
                    e1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
                elif pp4.Pt()>e2.Pt():
                    e3.SetPtEtaPhiM(e2.Pt(),e2.Eta(),e2.Phi(),e2.M())
                    e2.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
                else:
                    e3.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            else:
                print('Huh? (e)')
        if abs(pid)==13:
            if m1 is None:
                m1 = rt.TLorentzVector()
                m1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            elif m2 is None:
                m2 = rt.TLorentzVector()
                if pp4.Pt()>m1.Pt():
                    m2.SetPtEtaPhiM(m1.Pt(),m1.Eta(),m1.Phi(),m1.M())
                    m1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
                else:
                    m2.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            elif m3 is None:
                m3 = rt.TLorentzVector()
                if pp4.Pt()>m1.Pt():
                    m3.SetPtEtaPhiM(m2.Pt(),m2.Eta(),m2.Phi(),m2.M())
                    m2.SetPtEtaPhiM(m1.Pt(),m1.Eta(),m1.Phi(),m1.M())
                    m1.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
                elif pp4.Pt()>m2.Pt():
                    m3.SetPtEtaPhiM(m2.Pt(),m2.Eta(),m2.Phi(),m2.M())
                    m2.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
                else:
                    m3.SetPtEtaPhiM(pp4.Pt(),pp4.Eta(),pp4.Phi(),pp4.M())
            else:
                print('Huh? (m)')


for h in hists:
    can = rt.TCanvas("can","can",900,600)
    stack = rt.THStack("stack","stack")
    sighists = []
    binmax = 0.
    legend = rt.TLegend(0.4,0.7,0.9,0.9)
    legend.SetNColumns(2)
    for s in sample_order:
        hists[h][s].Scale(samples[s][4][com]*1000.*(1. if samples[s][0]==0. else float(samples[s][0]))*ilumi/hists[h][s].Integral())
        print(hists[h][s].Integral())
        hists[h][s].SetFillColor(samples[s][2])
        hists[h][s].SetLineColor(rt.kBlack)
        hists[h][s].SetLineWidth(2)
        if (samples[s][0]!=0.): 
            hists[h][s].SetFillColor(rt.kWhite)
            hists[h][s].SetLineColor(samples[s][2])
            hists[h][s].SetLineStyle(2)
            sighists.append(hists[h][s])
        else:
            stack.Add(hists[h][s])
        if hists[h][s].GetMaximum()>binmax:
            binmax = hists[h][s].GetMaximum()
        legend.AddEntry(hists[h][s],samples[s][1]+("" if samples[s][0]==0. else ' (x %i)'%samples[s][0]),"f")
    stack.Draw("hist")
    for s in sighists:
        s.Draw("hist same")
    legend.Draw()
    #for s in samples:
    #    hists[h][s].SetLineColor(rt.kBlack)
    #    hists[h][s].SetLineStyle(2)
    #    hists[h][s].SetLineWidth(2)
    #    if (samples[s][0]==1): hists[h][s].Draw("hist same")
        
    stack.SetTitle("")
    stack.GetXaxis().SetTitle(histnames[h])
    stack.GetYaxis().SetTitle("Events")
    stack.SetMaximum(binmax*1.4)
    tpt = rt.TPaveText(0.1,0.75,0.3,0.9,"nbNDC")
    tpt.AddText("\int L dt = %i fb^{-1}"%int(ilumi))
    tpt.AddText("\sqrt{s} = %s GeV"%com)
    tpt.SetFillStyle(0)
    tpt.SetLineStyle(0)
    tpt.SetLineWidth(0)
    tpt.SetTextFont(53)
    tpt.SetTextSize(18)
    tpt.Draw()
    can.SaveAs('%s_%s.pdf'%(h,com))
    can.SetLogy()
    stack.SetMaximum(binmax*20.)
    stack.SetMinimum(0.1)
    stack.Draw("hist")
    for s in sighists:
        s.Draw("hist same")
    legend.Draw()
    tpt.Draw()
    can.SaveAs('%s_%s_logy.pdf'%(h,com))
