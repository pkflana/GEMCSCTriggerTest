import ROOT
import os
import argparse

def make_plotdir(plotdir):
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

def plot_rate_scan(rdf, event, plotdir, evenodd, reglist, chargelist, bacutlist, ba_min, do_min_test, doRegionLevel, doChargeLevel, xbins=50, xlow=0, xhigh=50):
    colorlist = [ROOT.kBlack, ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]
    legend_xmin = 0.45
    legend_xmax = 0.85
    legend_ymin = 0.45
    legend_ymax = 0.85

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

    dRmatch = 0.4
    base_cut = "has_emtf_track_match & (has_LCT_match1)"
    residual_cut = "((abs(LCT_match_GE1_residual) < {res_cut}))"
    reco_cut = f"has_reco_l1_match & (abs(reco_l1_match_dR) < {dRmatch})"
    mode_cut = f"((emtftrack_mode == 11) | (emtftrack_mode == 13) | (emtftrack_mode == 14) | (emtftrack_mode == 15))"
    zmu_cut = f"(z_cand) & (abs(z_mass-91) < 10) & (z_opposite_charge)"

    for evenodd_val in evenodd:
        cut = f"{base_cut} & {mode_cut}"
        res_cut = 100
        ba_min_val = ba_min
        if evenodd_val == "even":
            cut += " & (LCT_CSC_chamber%2 == 0)"
            bacutlist_val = [1000, 12, 10, 8, 6, 5, 4]
            res_cut  = 20
            ba_min_val = -2
        else:
            cut += " & (LCT_CSC_chamber%2 == 1)"
            bacutlist_val = [1000, 28, 24, 16, 12, 10, 8]
            res_cut = 40
            ba_min_val = -4

        cut += " & " + residual_cut.format(res_cut = res_cut)
        cut += " & LCT_CSC_ME1b == 1"
        rdf_tmp = rdf.Filter(cut)
        for reg in reglist:
            if not doRegionLevel and reg in [-1, 1]:
                continue
            cut_reg = cut + f" & ((LCT_GE1_region == -1) | (LCT_GE1_region == 1))" if reg == 0 else cut + f" & (LCT_GE1_region == {reg})"
            for charge in chargelist:
                if not doChargeLevel and charge in [-1, 1]:
                    continue
                cut_charge = cut_reg + f" & ((emtftrack_charge == -1) | (emtftrack_charge == 1))" if charge == 0 else cut_reg + f" & (emtftrack_charge == {charge})"
                hlist = []
                hlist_notnorm = []
                h_nums = []
                for i, ba_cut in enumerate(bacutlist_val):
                    plot_branch = "emtftrack_pt"
                    xAxis_title = "EMTFTrack pT Threshold"
                    hname = f"BACut{ba_cut}"
                    hlist.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))
                    hlist_notnorm.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))
                    xAxis = hlist[i].GetXaxis()
                    xAxis.SetTitleOffset(2)
                    xAxis.SetTitleSize(0.04)
                    xAxis.SetTitle(f"{xAxis_title}")
                    ba_min_tmp = -1000 if (i == 0) else ba_min_val
                    if do_min_test:
                        cut1 = cut_charge + f" & ((emtftrack_charge*LCT_GE1_region*LCT_BendingAngle_GE1*(-1)) < {ba_cut}) & ((emtftrack_charge*LCT_GE1_region*LCT_BendingAngle_GE1*(-1)) >= {ba_min_tmp})"
                    else:
                        cut1 = cut_charge + f" & ((emtftrack_charge*LCT_GE1_region*LCT_BendingAngle_GE1*(-1)) < {ba_cut})"
                    hlist[i].SetLineColor(colorlist[i])
                    hlist_notnorm[i].SetLineColor(colorlist[i])
                    h_nums.append([])
                    for xbin in range(xbins):
                        cut1_num = cut1 + f" & emtftrack_pt >= {xbin}"
                        h_nums[i].append(rdf_tmp.Filter(cut1_num).Count())
                nEntries = event.GetEntries(cut)
                h_nums_count = []
                legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
                legend.SetHeader(f"{evenodd_val} -- R{reg} Charge{charge}", "C")
                for i in range(len(h_nums)):
                    h_nums_count.append([x.GetValue() for x in h_nums[i]])
                    for xbin in range(xbins):
                        hlist[i].SetBinContent(xbin+1, h_nums_count[i][xbin]/nEntries)
                        hlist_notnorm[i].SetBinContent(xbin+1, h_nums_count[i][xbin])
                    hlist[i].SetMarkerSize(0)
                    ba_min_tmp = -1000 if (i == 0) else ba_min_val
                    if do_min_test:
                        legend.AddEntry(hlist[i], f"{ba_min_tmp} <= BendingAngleCut <= {bacutlist_val[i]}")
                    else:
                        legend.AddEntry(hlist[i], f"BendingAngleCut <= {bacutlist_val[i]}")
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
                plotname = os.path.join(plotdir, f"RateScan_{evenodd_val}_R{reg}_Charge{charge}.pdf")
                canvas.SaveAs(plotname)
                canvas.SetLogy()
                hlist[0].SetMaximum(1e0)
                hlist[0].SetMinimum(1e-5)
                canvas.Update()
                plotname = os.path.join(plotdir, f"RateScan_{evenodd_val}_R{reg}_Charge{charge}_logy.pdf")
                canvas.SaveAs(plotname)
                canvas.SetLogy(0)
                hlist[0].SetMaximum(0.1)
                plotname = os.path.join(plotdir, f"RateScan_{evenodd_val}_R{reg}_Charge{charge}_ZOOMED1.pdf")
                canvas.SaveAs(plotname)
                hlist[0].SetMaximum(0.01)
                plotname = os.path.join(plotdir, f"RateScan_{evenodd_val}_R{reg}_Charge{charge}_ZOOMED2.pdf")
                canvas.SaveAs(plotname)
                hlist[0].SetMaximum(0.003)
                plotname = os.path.join(plotdir, f"RateScan_{evenodd_val}_R{reg}_Charge{charge}_ZOOMED3.pdf")
                canvas.SaveAs(plotname)
                hlist[0].SetMaximum(0.0005)
                plotname = os.path.join(plotdir, f"RateScan_{evenodd_val}_R{reg}_Charge{charge}_ZOOMED4.pdf")
                canvas.SaveAs(plotname)
                legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
                legend.SetHeader(f"{evenodd_val} -- R{reg} Charge{charge}", "C")
                for i in range(len(h_nums)):
                    ba_min_tmp = -1000 if (i == 0) else ba_min_val
                    hlist_notnorm[i].SetMarkerSize(0)
                    if do_min_test:
                        legend.AddEntry(hlist_notnorm[i], f"{ba_min_tmp} <= BendingAngleCut <= {bacutlist_val[i]}")
                    else:
                        legend.AddEntry(hlist_notnorm[i], f"BendingAngleCut <= {bacutlist_val[i]}")
                    xAxis = hlist_notnorm[i].GetXaxis()
                    xAxis.SetTitleOffset(2)
                    xAxis.SetTitleSize(0.04)
                    xAxis.SetTitle(f"{xAxis_title}")
                    yAxis_title = "Trigger Rate (A.U.)"
                    yAxis = hlist_notnorm[i].GetYaxis()
                    yAxis.SetTitleOffset(2)
                    yAxis.SetTitleSize(0.04)
                    yAxis.SetTitle(f"{yAxis_title}")
                    if i == 0:
                        hlist_notnorm[i].Draw("h")
                    else:
                        hlist_notnorm[i].Draw("h same")
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
                canvas.SetLogy(1)
                hlist_notnorm[0].SetMinimum(1e1)
                hlist_notnorm[0].SetMaximum(1e6)
                plotname = os.path.join(plotdir, f"RateScan_{evenodd_val}_R{reg}_Charge{charge}_notnorm.pdf")
                canvas.SaveAs(plotname)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True, type=str, help="Input ROOT file")
    parser.add_argument('--plotdir', default="plots/RateReduction/", type=str, help="Output plot directory")
    parser.add_argument('--xbins', type=int, default=50)
    parser.add_argument('--xlow', type=float, default=0)
    parser.add_argument('--xhigh', type=float, default=50)
    parser.add_argument('--doRegionLevel', action='store_true')
    parser.add_argument('--doChargeLevel', action='store_true')
    parser.add_argument('--do_min_test', action='store_true')
    args = parser.parse_args()

    make_plotdir(args.plotdir)
    ROOT.gROOT.SetBatch(1)
    try:
        import tdrstyle
        tdrstyle.setTDRStyle()
    except ImportError:
        print("tdrstyle not found, continuing without it.")

    ROOT.EnableImplicitMT(8)
    f = ROOT.TFile(args.infile)
    event = f.Get("GEMCSCBendingAngleTester/HighestEMTFpT")
    rdf = ROOT.RDataFrame(event)
    ROOT.RDF.Experimental.AddProgressBar(rdf)

    evenodd_list = ["odd", "even"]
    reglist = [-1, 1]
    chargelist = [-1, 1]
    bacutlist = [1000, 28, 24, 16, 12, 10, 8]  # default, will be overridden per even/odd
    ba_min = -5

    plot_rate_scan(
        rdf, event, args.plotdir, evenodd_list, reglist, chargelist, bacutlist, ba_min,
        args.do_min_test, args.doRegionLevel, args.doChargeLevel,
        xbins=args.xbins, xlow=args.xlow, xhigh=args.xhigh
    )

if __name__ == "__main__":
    main()