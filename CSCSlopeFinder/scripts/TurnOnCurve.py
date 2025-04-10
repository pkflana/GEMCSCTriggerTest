import ROOT, sys, tdrstyle, os, array

#f = ROOT.TFile("{}".format(sys.argv[1]))
fname = "/eos/user/d/daebi/GEMCSCTrigger/2024G_ZMu/2024G_ZMu_22Oct2024.root"
fname = "/eos/user/d/daebi/GEMCSCTrigger/2024G_ZMu/2024G_ZMu_27Oct2024.root"


f = ROOT.TFile(fname)

#event = f.Get("GEMCSCBendingAngleTester/AllLCTs") #ToC needs All LCTs
#event = f.Get("GEMCSCBendingAngleTester/HighestEMTFpT")
event = f.Get("GEMCSCBendingAngleTester/OnlyRecos") #All LCTs but only ones matched to a reco track (should be faster for eff study)

ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

ROOT.EnableImplicitMT(8)
rdf = ROOT.RDataFrame(event)
ROOT.RDF.Experimental.AddProgressBar(rdf)


if not os.path.exists("plots/"):
  os.makedirs("plots/")


plotdir = "plots/TurnOnCurve_ZMu30Oct/"
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
region_list = [-1, 1] #[-1, 0, 1]
charge_list = [-1, 1] #[-1, 0, 1]
dRmatch = 0.4

#base_cut = "has_emtf_track_match & (has_LCT_match1 | has_LCT_match2)"
#residual_cut = "((abs(LCT_match_GE1_residual) < {res_cut}) & (has_LCT_match1)) | ((abs(LCT_match_GE2_residual) < {res_cut}) & (has_LCT_match2))"

#Only do Layer 1 for now
base_cut = "has_emtf_track_match & (has_LCT_match1)"
residual_cut = "((abs(LCT_match_GE1_residual) < {res_cut}) & (has_LCT_match1))"
reco_cut = f"has_reco_l1_match & (abs(reco_l1_match_dR) < {dRmatch})"
mode_cut = f"((emtftrack_mode == 11) | (emtftrack_mode == 13) | (emtftrack_mode == 14) | (emtftrack_mode == 15))"
zmu_cut = f"(z_cand) & (abs(z_mass-91) < 10) & (z_opposite_charge)"

xbins = 50
xlow = 0
xhigh = 50


legend_xmin = 0.2
legend_xmax = 0.85
legend_ymin = 0.65
legend_ymax = 0.85

doRegionLevel = True
doChargeLevel = True

do_min_test = False

for evenodd in evenodd_list:

    cut = f"{base_cut} & {reco_cut} & {mode_cut} & {zmu_cut}" #Add residual cut later for even/odd difference
    res_cut = 100
    emtfcutlist = [14]
    bacutlist = [1000, 28, 24, 16, 12, 10, 8]
    ba_min = -5
    if evenodd == "even":
        cut = cut + " & (LCT_CSC_chamber%2 == 0)"
        emtfcutlist = [5, 6, 7, 8]
        bacutlist = [1000, 12, 10, 8, 6, 5, 4]
        res_cut  = 20
        ba_min = -8
    if evenodd == "odd":
        cut = cut + " & (LCT_CSC_chamber%2 == 1)"
        emtfcutlist = [5, 6, 7, 8]
        bacutlist = [1000, 28, 24, 16, 12, 10, 8]
        res_cut = 40
        ba_min = -16

    cut = cut + " & " + residual_cut.format(res_cut = res_cut)


    #cut = cut + " & emtftrack_charge == -1" #reco_l1_match_charge
    #cut = cut + " & LCT_GE1_region == 1"
    cut = cut + " & (LCT_CSC_ME1b == 1)"

    rdf_tmp = rdf.Filter(cut)

    for reg in region_list:
        if not doRegionLevel:
            if reg in [-1, 1]: continue
        cut_reg = cut + f" & ((LCT_GE1_region == -1) | (LCT_GE1_region == 1))" if reg == 0 else cut + f" & (LCT_GE1_region == {reg})"
        for charge in charge_list:
            if not doChargeLevel:
                if charge in [-1, 1]: continue
            cut_charge = cut_reg + f" & ((reco_l1_match_charge == -1) | (reco_l1_match_charge == 1))" if charge == 0 else cut_reg + f" & (reco_l1_match_charge == {charge})"

            for emtf_cut in emtfcutlist:
                hlist = []
                hnums = []
                hdens = []
                rdf_num_list = []
                rdf_den_list = []
                colorlist = [ROOT.kBlack, ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]
                legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
                legend.SetNColumns(2)
                legend.SetHeader(f"EMTFCut >= {emtf_cut} -- {evenodd} -- R{reg} Charge{charge}", "C")
                for i, ba_cut in enumerate(bacutlist):
                    hnums.append([])
                    hdens.append([])
                    hname = f"TurnOnCurve_EO{evenodd}_R{reg}_Charge{charge}_{i}"
                    hlist.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))
                    hnums.append(ROOT.TH1D(hname+"Num", hname, xbins, xlow, xhigh))
                    hdens.append(ROOT.TH1D(hname+"Den", hname, xbins, xlow, xhigh))

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

                    ba_min_tmp = -1000 if (i == 0) else ba_min

                    #cutnum = cut_charge + f" & abs(LCT_BendingAngle_GE1) <= {ba_cut} & emtftrack_pt >= {emtf_cut}"
                    if do_min_test:
                        cutnum = cut_charge + f" & ((emtftrack_charge*LCT_GE1_region*LCT_BendingAngle_GE1*(-1)) < {ba_cut}) & ((emtftrack_charge*LCT_GE1_region*LCT_BendingAngle_GE1*(-1)) >= {ba_min_tmp}) & emtftrack_pt >= {emtf_cut}"
                    else:
                        cutnum = cut_charge + f" & ((emtftrack_charge*LCT_GE1_region*LCT_BendingAngle_GE1*(-1)) < {ba_cut}) & emtftrack_pt >= {emtf_cut}"

                    cutden = cut_charge

                    rdf_num_list.append(rdf_tmp.Filter(cutnum))
                    rdf_den_list.append(rdf_tmp.Filter(cutden))


                print("Finished the loop of setting filters")
                #Now fill hists with the RDF values?
                for i, ba_cut in enumerate(bacutlist):
                    print(f"At BA cut {ba_cut}")
                    hnums[i] = rdf_num_list[i].Histo1D((f"Nums{i}", f"Nums{i}", xbins, xlow, xhigh), "reco_l1_match_pt")
                    hdens[i] = rdf_den_list[i].Histo1D((f"Dens{i}", f"Dens{i}", xbins, xlow, xhigh), "reco_l1_match_pt")

                for i, ba_cut in enumerate(bacutlist):
                    hnums[i] = hnums[i].GetValue()
                    hdens[i] = hdens[i].GetValue()
                    hlist[i] = hnums[i]
                
                for i, ba_cut in enumerate(bacutlist):
                    hlist[i].Divide(hdens[i])
                    hlist[i].SetLineColor(colorlist[i])

                    if i == 0:
                        hlist[i].Draw("h")
                        yAxis = hlist[i].GetYaxis()
                        yAxis.SetTitleOffset(2)
                        yAxis.SetTitleSize(0.04)
                        yAxis.SetTitle(f"{yAxis_title}")
                        yAxis.SetRangeUser(0.0, 1.5)
                        xAxis = hlist[i].GetXaxis()
                        xAxis.SetTitleOffset(2)
                        xAxis.SetTitleSize(0.04)
                        xAxis.SetTitle(f"{xAxis_title}")
                    else:
                        hlist[i].Draw("h same")

                    ba_min_tmp = -1000 if (i == 0) else ba_min

                    hlist[i].SetMarkerSize(0)

                    if do_min_test:
                        legend.AddEntry(hlist[i], f"{ba_min_tmp} <= BendingAngleCut <= {ba_cut}")
                    else:
                        if i == 0:
                            legend.AddEntry(hlist[i], f"Any GEM Match")
                        else:
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

                if do_min_test:
                    canvas.SaveAs(plotdir+f"TurnOnCurve_emtfpt{emtf_cut}_{evenodd}_R{reg}_Charge{charge}_bamin{ba_min}.pdf")
                else:
                    canvas.SaveAs(plotdir+f"TurnOnCurve_emtfpt{emtf_cut}_{evenodd}_R{reg}_Charge{charge}.pdf")



                legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
                legend.SetNColumns(2)
                legend.SetHeader(f"EMTFCut >= {emtf_cut} -- {evenodd} -- R{reg} Charge{charge}", "C")

                for i, ba_cut in enumerate(bacutlist):
                    if i == 0: continue
                    htmp = hlist[i]
                    htmp.Divide(hlist[0])

                    if i == 1:
                        hlist[i].Draw("h")
                        yAxis = hlist[i].GetYaxis()
                        yAxis.SetTitleOffset(2)
                        yAxis.SetTitleSize(0.04)
                        yAxis.SetTitle(f"BA Cut Efficiency on EMTF")
                        yAxis.SetRangeUser(0.5, 1.25)
                        xAxis = hlist[i].GetXaxis()
                        xAxis.SetTitleOffset(2)
                        xAxis.SetTitleSize(0.04)
                        xAxis.SetTitle(f"{xAxis_title}")
                    else:
                        hlist[i].Draw("h same")

                    ba_min_tmp = -1000 if (i == 0) else ba_min

                    hlist[i].SetMarkerSize(0)

                    if do_min_test:
                        legend.AddEntry(hlist[i], f"{ba_min_tmp} <= BendingAngleCut <= {ba_cut}")
                    else:
                        if i == 0:
                            legend.AddEntry(hlist[i], f"Any GEM Match")
                        else:
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

                if do_min_test:
                    canvas.SaveAs(plotdir+f"Efficiency_emtfpt{emtf_cut}_{evenodd}_R{reg}_Charge{charge}_bamin{ba_min}.pdf")
                else:
                    canvas.SaveAs(plotdir+f"Efficiency_emtfpt{emtf_cut}_{evenodd}_R{reg}_Charge{charge}.pdf")



        """

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
                    hname = f"RecoPt_EO{evenodd}_R{reg}_Charge{charge}_{i}"
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

                    cutnum = cut_charge + f" & (abs(BendingAngleRealistic) <= {ba_cut}) & emtftrack_pt >= {emtf_cut}"
                    event.Project(hnum_name, "reco_l1_match_pt", cutnum)

                    hden_name = f"TurnOnCurve{i}Den"
                    hden = ROOT.TH1D(hden_name, hden_name, xbins, xlow, xhigh)

                    cutden = cut_charge# + f" & abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_cut}"
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

                canvas.SaveAs(plotdir+f"RecoPt_Dist_emtfpt{emtf_cut}_{evenodd}_R{reg}_Charge{charge}.pdf")


        """













#Remove any GEM Reqs from denominator
base_cut = "has_emtf_track_match & (has_LCT_match1)"
residual_cut = "((abs(LCT_match_GE1_residual) < {res_cut}) & (has_LCT_match1))"
reco_cut = f"has_reco_l1_match & (abs(reco_l1_match_dR) < {dRmatch})"
mode_cut = f"((emtftrack_mode == 11) | (emtftrack_mode == 13) | (emtftrack_mode == 14) | (emtftrack_mode == 15))"
zmu_cut = f"(z_cand) & (abs(z_mass-91) < 10) & (z_opposite_charge)"

xbins = 50
xlow = 0
xhigh = 50



if False:
    doRegionLevel = True
    doChargeLevel = True

    for evenodd in evenodd_list:
        cut = f"{reco_cut} & {mode_cut} & {zmu_cut}" #Add residual cut later for even/odd difference
        res_cut = 100
        emtfcutlist = [14]
        bacutlist = [1000, 28, 24, 16, 12, 10, 8]
        if evenodd == "even":
            cut = cut + " & (LCT_CSC_chamber%2 == 0)"
            emtfcutlist = [14]
            bacutlist = [1000, 12, 10, 8, 6, 5, 4]
            bacutlist = [1000, 12]
            res_cut  = 20
        if evenodd == "odd":
            cut = cut + " & (LCT_CSC_chamber%2 == 1)"
            emtfcutlist = [14]
            bacutlist = [1000, 28, 24, 16, 12, 10, 8]
            bacutlist = [1000, 28]
            res_cut = 40

        #cut = cut + " & emtftrack_charge == -1" #Alexei says we should use RECO charge, for now I do not have reco charge, can use l1muon_match_charge instead
        #cut = cut + " & LCT_GE1_region == 1"
        cut = cut + " & (LCT_CSC_ME1b == 1)"

        rdf_tmp = rdf.Filter(cut)

        for reg in [-1, 0, 1]:
            if not doRegionLevel:
                if reg in [-1, 1]: continue
            cut_reg = cut + f" & ((LCT_GE1_region == -1) | (LCT_GE1_region == 1))" if reg == 0 else cut + f" & (LCT_GE1_region == {reg})"
            for charge in [-1, 0, 1]:
                if not doChargeLevel:
                    if charge in [-1, 1]: continue
                cut_charge = cut_reg + f" & ((reco_l1_match_charge == -1) | (reco_l1_match_charge == 1))" if charge == 0 else cut_reg + f" & (reco_l1_match_charge == {charge})"

                for emtf_cut in emtfcutlist:
                    hlist = []
                    hnums = []
                    hdens = []
                    rdf_num_list = []
                    rdf_den_list = []
                    colorlist = [ROOT.kBlack, ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]
                    legend = ROOT.TLegend(legend_xmin, legend_ymin, legend_xmax, legend_ymax)
                    legend.SetNColumns(2)
                    legend.SetHeader(f"EMTFCut >= {emtf_cut} -- {evenodd} -- R{reg} Charge{charge}", "C")
                    for i, ba_cut in enumerate(bacutlist):
                        hnums.append([])
                        hdens.append([])
                        hname = f"TurnOnCurve_EO{evenodd}_R{reg}_Charge{charge}_{i}"
                        hlist.append(ROOT.TH1D(hname, hname, xbins, xlow, xhigh))
                        hnums.append(ROOT.TH1D(hname+"Num", hname, xbins, xlow, xhigh))
                        hdens.append(ROOT.TH1D(hname+"Den", hname, xbins, xlow, xhigh))

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

                        gem_cut = cut + f" & {base_cut}"
                        gem_cut = gem_cut + " & " + residual_cut.format(res_cut = res_cut)
                        gem_cut = gem_cut + f" & abs(LCT_BendingAngle_GE1) <= {ba_cut}"

                        if i == 0:
                            cutnum = cut_charge + f" & emtftrack_pt >= {emtf_cut}"
                        else:
                            cutnum = cut_charge + f" & {gem_cut} & emtftrack_pt >= {emtf_cut}"
                        cutden = cut_charge

                        rdf_num_list.append(rdf_tmp.Filter(cutnum))
                        rdf_den_list.append(rdf_tmp.Filter(cutden))


                    print("Finished the loop of setting filters")
                    #Now fill hists with the RDF values?
                    for i, ba_cut in enumerate(bacutlist):
                        print(f"At BA cut {ba_cut}")
                        hnums[i] = rdf_num_list[i].Histo1D((f"Nums{i}", f"Nums{i}", xbins, xlow, xhigh), "reco_l1_match_pt")
                        hdens[i] = rdf_den_list[i].Histo1D((f"Dens{i}", f"Dens{i}", xbins, xlow, xhigh), "reco_l1_match_pt")

                    for i, ba_cut in enumerate(bacutlist):
                        hnums[i] = hnums[i].GetValue()
                        hdens[i] = hdens[i].GetValue()
                        hlist[i] = hnums[i]
                    
                    for i, ba_cut in enumerate(bacutlist):
                        hlist[i].Divide(hdens[i])
                        hlist[i].SetLineColor(colorlist[i])

                        if i == 0:
                            hlist[i].Draw("h")
                            yAxis = hlist[i].GetYaxis()
                            yAxis.SetTitleOffset(2)
                            yAxis.SetTitleSize(0.04)
                            yAxis.SetTitle(f"{yAxis_title}")
                            yAxis.SetRangeUser(0.0, 1.5)
                            xAxis = hlist[i].GetXaxis()
                            xAxis.SetTitleOffset(2)
                            xAxis.SetTitleSize(0.04)
                            xAxis.SetTitle(f"{xAxis_title}")
                        else:
                            hlist[i].Draw("h same")

                        hlist[i].SetMarkerSize(0)

                        if i == 0:
                            legend.AddEntry(hlist[i], f"No GEM Requirement")
                        else:
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

                    canvas.SaveAs(plotdir+f"TurnOnCurve_emtfpt{emtf_cut}_{evenodd}_R{reg}_Charge{charge}_NoGEMOption.pdf")

