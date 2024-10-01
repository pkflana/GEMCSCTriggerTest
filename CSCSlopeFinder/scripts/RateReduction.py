import ROOT, sys, tdrstyle, os, array

fname = "/eos/user/d/daebi/GEMCSCTrigger/2024E_ZeroBias/ZeroBias_EFG_18Sep2024.root"

f = ROOT.TFile(fname)

#event = f.Get("GEMCSCBendingAngleTester/LCTL1MuonMatcher")
event = f.Get("GEMCSCBendingAngleTester/HighestEMTFpT") #Use only highest EMTFpT in event for rate
ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

rdf = ROOT.RDataFrame(event)
ROOT.RDF.Experimental.AddProgressBar(rdf)

if not os.path.exists("plots/"):
  os.makedirs("plots/")


plotdir = "plots/RateReduction_ZeroBias/"
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


legend_xmin = 0.55
legend_xmax = 0.85
legend_ymin = 0.55
legend_ymax = 0.85

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
    cut = f"{base_cut} & {reco_cut} & {mode_cut}" #Add residual cut later for even/odd difference
    res_cut = 100
    bacutlist = [1000, 28, 24, 16, 12, 10, 8]
    if evenodd == "even":
        cut = cut + " & (LCT_CSC_chamber%2 == 0)"
        bacutlist = [1000, 14, 12, 8, 6, 5, 4]
        res_cut  = 20
    if evenodd == "odd":
        cut = cut + " & (LCT_CSC_chamber%2 == 1)"
        bacutlist = [1000, 28, 24, 16, 12, 10, 8]
        res_cut = 40

    cut = cut + " & " + residual_cut.format(res_cut = res_cut)



    #cut = cut + " & emtftrack_charge == -1"
    #cut = cut + " & LCT_GE1_region == 1"
    cut = cut + " & LCT_CSC_ME1b == 1"

    rdf_tmp = rdf.Filter(cut)

    hlist = []
    h_nums = []
    colorlist = [ROOT.kBlack, ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]
    for i, ba_cut in enumerate(bacutlist):
        plot_branch = "emtftrack_pt"
        xAxis_title = "EMTFTrack pT Threshold"

        hname = f"BACut{ba_cut}"
        hlist.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))
        xAxis = hlist[i].GetXaxis()
        xAxis.SetTitleOffset(2)
        xAxis.SetTitleSize(0.04)
        xAxis.SetTitle(f"{xAxis_title}")

        cut1 = cut + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut}"
        hlist[i].SetLineColor(colorlist[i])

        #Prepare the count numbers in RDF for all cuts
        h_nums.append([])

        for xbin in range(xbins):
            cut1_num = cut1 + f" & emtftrack_pt >= {xbin}"
            h_nums[i].append(rdf_tmp.Filter(cut1_num).Count())


    nEntries = event.GetEntries(cut)

    h_nums_count = []

    legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)

    for i in range(len(h_nums)):
        h_nums_count.append([x.GetValue() for x in h_nums[i]])
        for xbin in range(xbins):
            hlist[i].SetBinContent(xbin+1, h_nums_count[i][xbin]/nEntries)
        legend.AddEntry(hlist[i], f"CSC+GE1 Bending Angle <= {bacutlist[i]}")
        xAxis = hlist[i].GetXaxis()
        xAxis.SetTitleOffset(2)
        xAxis.SetTitleSize(0.04)
        xAxis.SetTitle(f"{xAxis_title}")
        yAxis_title = "Trigger Rate (A.U.)"
        yAxis = hlist[i].GetYaxis()
        yAxis.SetTitleOffset(2)
        yAxis.SetTitleSize(0.04)
        yAxis.SetTitle(f"{yAxis_title}")   
        if i == 0:
            hlist[i].Draw("h")
        else:
            hlist[i].Draw("h same")


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

    canvas.SaveAs(plotdir+f"RateScan_{evenodd}.pdf")

    canvas.SetLogy()

    hlist[0].SetMaximum(1e0)
    hlist[0].SetMinimum(1e-5)
    canvas.Update()

    canvas.SaveAs(plotdir+f"RateScan_{evenodd}_logy.pdf")

    canvas.SetLogy(0)

    hlist[0].SetMaximum(0.1)

    canvas.SaveAs(plotdir+f"RateScan_{evenodd}_ZOOMED1.pdf")


    hlist[0].SetMaximum(0.01)

    canvas.SaveAs(plotdir+f"RateScan_{evenodd}_ZOOMED2.pdf")



    hlist[0].SetMaximum(0.002)

    canvas.SaveAs(plotdir+f"RateScan_{evenodd}_ZOOMED3.pdf")








#And last the combined dataset
cut = f"{base_cut} & {reco_cut} & {mode_cut}" #Add residual cut later for even/odd difference
cut = cut + f" & ( ((LCT_CSC_chamber%2 == 0) & ({residual_cut.format(res_cut = 20)})) | ((LCT_CSC_chamber%2 == 1) & ({residual_cut.format(res_cut = 40)})) )"
ba_odd = 16
ba_even = 8
bacutlist = [ [1000, 1000], [16, 8], [15, 7], [12, 6], [11, 5], [9, 4] ]

#cut = cut + " & emtftrack_charge == -1"
#cut = cut + " & LCT_GE1_region == 1"
cut = cut + " & LCT_CSC_ME1b == 1"

rdf_tmp = rdf.Filter(cut)

hlist = []
h_nums = []
colorlist = [ROOT.kBlack, ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]
for i, ba_cuts in enumerate(bacutlist):
    ba_odd = ba_cuts[0]
    ba_even = ba_cuts[1]
    plot_branch = "emtftrack_pt"
    xAxis_title = "EMTFTrack pT Threshold"

    hname = f"BACut{ba_cut}"
    hlist.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))
    xAxis = hlist[i].GetXaxis()
    xAxis.SetTitleOffset(2)
    xAxis.SetTitleSize(0.04)
    xAxis.SetTitle(f"{xAxis_title}")

    cut_ba_odd = f"((LCT_CSC_chamber%2 == 1) & (abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_odd}))"
    cut_ba_even = f"((LCT_CSC_chamber%2 == 1) & (abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_even}))"
    cut1 = cut + f" & ({cut_ba_odd} | {cut_ba_even})"
    hlist[i].SetLineColor(colorlist[i])

    #Prepare the count numbers in RDF for all cuts
    h_nums.append([])

    for xbin in range(xbins):
        cut1_num = cut1 + f" & emtftrack_pt >= {xbin}"
        h_nums[i].append(rdf_tmp.Filter(cut1_num).Count())


nEntries = event.GetEntries(cut)

h_nums_count = []

legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)

for i in range(len(h_nums)):
    h_nums_count.append([x.GetValue() for x in h_nums[i]])
    for xbin in range(xbins):
        hlist[i].SetBinContent(xbin+1, h_nums_count[i][xbin]/nEntries)
    legend.AddEntry(hlist[i], f"BA Odd {bacutlist[i][0]} Even {bacutlist[i][1]}")
    xAxis = hlist[i].GetXaxis()
    xAxis.SetTitleOffset(2)
    xAxis.SetTitleSize(0.04)
    xAxis.SetTitle(f"{xAxis_title}")
    yAxis_title = "Trigger Rate (A.U.)"
    yAxis = hlist[i].GetYaxis()
    yAxis.SetTitleOffset(2)
    yAxis.SetTitleSize(0.04)
    yAxis.SetTitle(f"{yAxis_title}")   
    if i == 0:
        hlist[i].Draw("h")
    else:
        hlist[i].Draw("h same")


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

canvas.SaveAs(plotdir+f"RateScan_Combined.pdf")

canvas.SetLogy()

hlist[0].SetMaximum(1e0)
hlist[0].SetMinimum(1e-5)
canvas.Update()

canvas.SaveAs(plotdir+f"RateScan_Combined_logy.pdf")

canvas.SetLogy(0)

hlist[0].SetMaximum(0.1)

canvas.SaveAs(plotdir+f"RateScan_Combined_ZOOMED1.pdf")


hlist[0].SetMaximum(0.01)

canvas.SaveAs(plotdir+f"RateScan_Combined_ZOOMED2.pdf")



hlist[0].SetMaximum(0.002)

canvas.SaveAs(plotdir+f"RateScan_Combined_ZOOMED3.pdf")




