import ROOT, sys, tdrstyle, os, array

#f = ROOT.TFile("{}".format(sys.argv[1]))
fname = "/eos/user/d/daebi/GEMCSCTrigger/2024E_ZMu/Muon0/2024E_ZMu_18Sep2024/240917_185009/2024E_ZMu_18Sep2024.root"

f = ROOT.TFile(fname)

event = f.Get("GEMCSCBendingAngleTester/AllLCTs") #ToC needs All LCTs
#event = f.Get("GEMCSCBendingAngleTester/HighestEMTFpT")
ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

rdf = ROOT.RDataFrame(event)

if not os.path.exists("plots/"):
  os.makedirs("plots/")


plotdir = "plots/TurnOnCurve_ZMu_zcut/"
if not os.path.exists(plotdir):
  os.makedirs(plotdir)


H_ref = 600
W_ref = 800
W = W_ref
H = H_ref

T = 0.12*H
B = 0.16*H
L = 0.16*W
R = 0.12*W

canvas = ROOT.TCanvas("c1", "c1", W, H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)

ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()


yscale = 1.3

evenodd_list = ["odd", "even"]
dRmatch = 0.4

base_cut = "has_emtf_track_match & has_LCT_match1"
residual_cut = "(abs(LCT_match_GE1_residual) < {res_cut})"
reco_cut = f"has_reco_l1_match & (abs(reco_l1_match_dR) < {dRmatch})"
mode_cut = f"((emtftrack_mode == 11) | (emtftrack_mode == 13) | (emtftrack_mode == 14) | (emtftrack_mode == 15))"
zmu_cut = f"(z_cand) & (abs(z_mass-91) < 10)"

xbins = 50
xlow = 0
xhigh = 50



for evenodd in evenodd_list:
    legend_xmin = 0.5
    legend_xmax = 0.85
    legend_ymin = 0.2
    legend_ymax = 0.6
    cut = f"{base_cut} & {reco_cut} & {mode_cut} & {zmu_cut}" #Add residual cut later for even/odd difference
    res_cut = 100
    emtfcutlist = [14]
    bacutlist = [1000, 28, 24, 16, 12, 10, 8]
    if evenodd == "even":
        cut = cut + " & (LCT_CSC_chamber%2 == 0)"
        emtfcutlist = [22]
        bacutlist = [1000, 14, 12, 8, 6, 5, 4]
        res_cut  = 20
    if evenodd == "odd":
        cut = cut + " & (LCT_CSC_chamber%2 == 1)"
        emtfcutlist = [22]
        bacutlist = [1000, 28, 24, 16, 12, 10, 8]
        res_cut = 40

    cut = cut + " & " + residual_cut.format(res_cut = res_cut)


    #cut = cut + " & emtftrack_charge == -1"
    #cut = cut + " & LCT_GE1_region == 1"
    cut = cut + " & LCT_CSC_ME1b == 1"

    for emtf_cut in emtfcutlist:
        hlist = []
        colorlist = [ROOT.kBlack, ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]
        legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
        legend.SetHeader(f"EMTFCut >= {emtf_cut}", "C")
        for i, ba_cut in enumerate(bacutlist):
            hname = f"TurnOnCurve{i}"
            hlist.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))


            xAxis_title = "Offline RECO pT"
            yAxis_title = f"Marginal Efficiency / 1GeV"

            xAxis = hlist[i].GetXaxis()
            xAxis.SetTitleOffset(2)
            xAxis.SetTitleSize(0.04)
            xAxis.SetTitle(f"{xAxis_title}")
            

            cuts_num = []
            cuts_den = []
            counts_num = []
            counts_den = []

            hnum_name = f"TurnOnCurve{i}Num"
            hnum = ROOT.TH1D(hnum_name, hnum_name, xbins, xlow, xhigh)

            cutnum = cut + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut} & emtftrack_pt >= {emtf_cut}"
            event.Project(hnum_name, "reco_l1_match_pt", cutnum)

            hden_name = f"TurnOnCurve{i}Den"
            hden = ROOT.TH1D(hden_name, hden_name, xbins, xlow, xhigh)

            cutden = cut# + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut}"
            event.Project(hden_name, "reco_l1_match_pt", cutden)

            hnum.Divide(hden)

            hlist[i] = hnum

            hlist[i].SetLineColor(colorlist[i])

            if i == 0:
                hlist[i].Draw("h")
                yAxis = hlist[i].GetYaxis()
                yAxis.SetTitleOffset(2)
                yAxis.SetTitleSize(0.04)
                yAxis.SetTitle(f"{yAxis_title}")
                yAxis.SetRangeUser(0.0, 1.0)
                xAxis = hlist[i].GetXaxis()
                xAxis.SetTitleOffset(2)
                xAxis.SetTitleSize(0.04)
                xAxis.SetTitle(f"{xAxis_title}")
            else:
                hlist[i].Draw("h same")

            legend.AddEntry(hlist[i], f"BendingAngleCut <= {ba_cut}")
            
        legend.SetTextSize(0.)
        legend.SetBorderSize(0)
        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)

        latex.SetTextAlign(12)
        latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextAlign(32)
        latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
        latex.DrawLatex(1-ROOT.gPad.GetRightMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

        frame = canvas.GetFrame()
        frame.Draw()

        canvas.SetLogy(0)

        canvas.SaveAs(plotdir+f"TurnOnCurve_emtfpt{emtf_cut}_notintegrated_{evenodd}.pdf")



    for emtf_cut in emtfcutlist:
        hlist = []
        colorlist = [ROOT.kBlack, ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]
        legend_xmin = 0.2
        legend_xmax = 0.5
        legend_ymin = 0.4
        legend_ymax = 0.8
        legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
        legend.SetHeader(f"EMTFCut >= {emtf_cut}", "C")
        for i, ba_cut in enumerate(bacutlist):
            hname = f"TurnOnCurve{i}"
            hlist.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))


            xAxis_title = "Offline RECO pT"
            yAxis_title = f"Entries"

            xAxis = hlist[i].GetXaxis()
            xAxis.SetTitleOffset(2)
            xAxis.SetTitleSize(0.04)
            xAxis.SetTitle(f"{xAxis_title}")
            

            cuts_num = []
            cuts_den = []
            counts_num = []
            counts_den = []

            hnum_name = f"TurnOnCurve{i}Num"
            hnum = ROOT.TH1D(hnum_name, hnum_name, xbins, xlow, xhigh)

            cutnum = cut + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut} & emtftrack_pt >= {emtf_cut}"
            event.Project(hnum_name, "reco_l1_match_pt", cutnum)

            hden_name = f"TurnOnCurve{i}Den"
            hden = ROOT.TH1D(hden_name, hden_name, xbins, xlow, xhigh)

            cutden = cut# + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut}"
            event.Project(hden_name, "reco_l1_match_pt", cutden)

            #hnum.Divide(hden)

            hlist[i] = hnum

            hlist[i].SetLineColor(colorlist[i])

            if i == 0:
                hlist[i].Draw("h")
                yAxis = hlist[i].GetYaxis()
                yAxis.SetTitleOffset(2)
                yAxis.SetTitleSize(0.04)
                yAxis.SetTitle(f"{yAxis_title}")
                #yAxis.SetRangeUser(0.0, 1.0)
                xAxis = hlist[i].GetXaxis()
                xAxis.SetTitleOffset(2)
                xAxis.SetTitleSize(0.04)
                xAxis.SetTitle(f"{xAxis_title}")
            else:
                hlist[i].Draw("h same")


            legend.AddEntry(hlist[i], f"BendingAngleCut <= {ba_cut}")
            
        legend.SetTextSize(0.)
        legend.SetBorderSize(0)
        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)

        latex.SetTextAlign(12)
        latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextAlign(32)
        latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
        latex.DrawLatex(1-ROOT.gPad.GetRightMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

        frame = canvas.GetFrame()
        frame.Draw()

        canvas.SetLogy(0)

        canvas.SaveAs(plotdir+f"RecoPt_Dist_emtfpt{emtf_cut}_{evenodd}.pdf")








legend_xmin = 0.5
legend_xmax = 0.85
legend_ymin = 0.2
legend_ymax = 0.6
legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
legend.SetHeader(f"EMTFCut >= 22", "C")



hlist = []
for j, evenodd in enumerate(evenodd_list):

    cut = f"{base_cut} & {reco_cut} & {mode_cut} & {zmu_cut}" #Add residual cut later for even/odd difference
    res_cut = 100
    emtfcutlist = [14]
    bacutlist = [1000]
    if evenodd == "even":
        cut = cut + " & (LCT_CSC_chamber%2 == 0)"
        emtfcutlist = [22]
        bacutlist = [1000, 8]
        res_cut  = 20
        colorlist = [ROOT.kBlue, ROOT.kCyan+1]
    if evenodd == "odd":
        cut = cut + " & (LCT_CSC_chamber%2 == 1)"
        emtfcutlist = [22]
        bacutlist = [1000, 16]
        res_cut = 40
        colorlist = [ROOT.kRed, ROOT.kMagenta]

    cut = cut + " & " + residual_cut.format(res_cut = res_cut)

    #cut = cut + " & emtftrack_charge == -1"
    #cut = cut + " & LCT_GE1_region == 1"
    cut = cut + " & LCT_CSC_ME1b == 1"


    
    hlist.append([])
    for emtf_cut in emtfcutlist:
        for i, ba_cut in enumerate(bacutlist):
            hname = f"TurnOnCurve{i}"
            hlist[j].append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))


            xAxis_title = "Offline RECO pT"
            yAxis_title = f"Marginal Efficiency / 1GeV"

            xAxis = hlist[j][i].GetXaxis()
            xAxis.SetTitleOffset(2)
            xAxis.SetTitleSize(0.04)
            xAxis.SetTitle(f"{xAxis_title}")
            

            cuts_num = []
            cuts_den = []
            counts_num = []
            counts_den = []

            hnum_name = f"TurnOnCurve{i}Num"
            hnum = ROOT.TH1D(hnum_name, hnum_name, xbins, xlow, xhigh)

            cutnum = cut + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut} & emtftrack_pt >= {emtf_cut}"
            event.Project(hnum_name, "reco_l1_match_pt", cutnum)

            hden_name = f"TurnOnCurve{i}Den"
            hden = ROOT.TH1D(hden_name, hden_name, xbins, xlow, xhigh)

            cutden = cut# + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut}"
            event.Project(hden_name, "reco_l1_match_pt", cutden)

            hnum.Divide(hden)

            hlist[j][i] = hnum

            hlist[j][i].SetLineColor(colorlist[i])

            if i == 0:
                if j == 0:
                    hlist[j][i].Draw("h")
                else:
                    hlist[j][i].Draw("h same")
                yAxis = hlist[j][i].GetYaxis()
                yAxis.SetTitleOffset(2)
                yAxis.SetTitleSize(0.04)
                yAxis.SetTitle(f"{yAxis_title}")
                yAxis.SetRangeUser(0.0, 1.0)
                xAxis = hlist[j][i].GetXaxis()
                xAxis.SetTitleOffset(2)
                xAxis.SetTitleSize(0.04)
                xAxis.SetTitle(f"{xAxis_title}")
            else:
                hlist[j][i].Draw("h same")

            legend.AddEntry(hlist[j][i], f"BendingAngleCut <= {ba_cut} {evenodd}")
            
    legend.SetTextSize(0.)
    legend.SetBorderSize(0)
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)

    latex.SetTextAlign(12)
    latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
    latex.SetTextFont(61)
    latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
    latex.SetTextFont(52)
    latex.SetTextAlign(32)
    latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
    latex.DrawLatex(1-ROOT.gPad.GetRightMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

    frame = canvas.GetFrame()
    frame.Draw()

    canvas.SetLogy(0)

    canvas.SaveAs(plotdir+f"TurnOnCurve_EvenOddComp.pdf")







legend_xmin = 0.5
legend_xmax = 0.85
legend_ymin = 0.2
legend_ymax = 0.6
cut = f"{base_cut} & {reco_cut} & {mode_cut} & {zmu_cut}" #Add residual cut later for even/odd difference
cut = cut + f" & ( ((LCT_CSC_chamber%2 == 0) & ({residual_cut.format(res_cut = 20)})) | ((LCT_CSC_chamber%2 == 1) & ({residual_cut.format(res_cut = 40)})) )"
bacutlist = [ [1000, 1000], [16, 8], [15, 7], [13, 6], [11, 5], [9, 4] ]
emtfcutlist = [22]

cut = cut + " & " + residual_cut.format(res_cut = res_cut)


#cut = cut + " & emtftrack_charge == -1"
#cut = cut + " & LCT_GE1_region == 1"
cut = cut + " & LCT_CSC_ME1b == 1"

for emtf_cut in emtfcutlist:
    hlist = []
    colorlist = [ROOT.kBlack, ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]
    legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
    legend.SetHeader(f"EMTFCut >= {emtf_cut}", "C")
    for i, ba_cut in enumerate(bacutlist):
        ba_odd = ba_cut[0]
        ba_even = ba_cut[1]
        hname = f"TurnOnCurve{i}"
        hlist.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))


        xAxis_title = "Offline RECO pT"
        yAxis_title = f"Marginal Efficiency / 1GeV"

        xAxis = hlist[i].GetXaxis()
        xAxis.SetTitleOffset(2)
        xAxis.SetTitleSize(0.04)
        xAxis.SetTitle(f"{xAxis_title}")
        

        cuts_num = []
        cuts_den = []
        counts_num = []
        counts_den = []

        hnum_name = f"TurnOnCurve{i}Num"
        hnum = ROOT.TH1D(hnum_name, hnum_name, xbins, xlow, xhigh)

        cut_ba_odd = f"((LCT_CSC_chamber%2 == 1) & (abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_odd}))"
        cut_ba_even = f"((LCT_CSC_chamber%2 == 0) & (abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_even}))"
        cutnum = cut + f" & ({cut_ba_odd} | {cut_ba_even}) & emtftrack_pt >= {emtf_cut}"
        event.Project(hnum_name, "reco_l1_match_pt", cutnum)

        hden_name = f"TurnOnCurve{i}Den"
        hden = ROOT.TH1D(hden_name, hden_name, xbins, xlow, xhigh)

        cutden = cut# + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut}"
        event.Project(hden_name, "reco_l1_match_pt", cutden)

        hnum.Divide(hden)

        hlist[i] = hnum

        hlist[i].SetLineColor(colorlist[i])

        if i == 0:
            hlist[i].Draw("h")
            yAxis = hlist[i].GetYaxis()
            yAxis.SetTitleOffset(2)
            yAxis.SetTitleSize(0.04)
            yAxis.SetTitle(f"{yAxis_title}")
            yAxis.SetRangeUser(0.0, 1.0)
            xAxis = hlist[i].GetXaxis()
            xAxis.SetTitleOffset(2)
            xAxis.SetTitleSize(0.04)
            xAxis.SetTitle(f"{xAxis_title}")
        else:
            hlist[i].Draw("h same")

        legend.AddEntry(hlist[i], f"Odd {ba_odd} Even {ba_even}")
        
    legend.SetTextSize(0.)
    legend.SetBorderSize(0)
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)

    latex.SetTextAlign(12)
    latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
    latex.SetTextFont(61)
    latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
    latex.SetTextFont(52)
    latex.SetTextAlign(32)
    latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
    latex.DrawLatex(1-ROOT.gPad.GetRightMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

    frame = canvas.GetFrame()
    frame.Draw()

    canvas.SetLogy(0)

    canvas.SaveAs(plotdir+f"TurnOnCurve_emtfpt{emtf_cut}_Combined.pdf")


