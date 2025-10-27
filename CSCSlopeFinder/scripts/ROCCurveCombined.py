import ROOT
import os
import argparse
import numpy as np

def make_plotdir(plotdir):
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

def get_rdfs(sig_file, bkg_file, sig_tree, bkg_tree):
    f_sig = ROOT.TFile(sig_file)
    f_bkg = ROOT.TFile(bkg_file)
    event_sig = f_sig.Get(sig_tree)
    event_bkg = f_bkg.Get(bkg_tree)
    rdf_sig = ROOT.RDataFrame(event_sig)
    ROOT.RDF.Experimental.AddProgressBar(rdf_sig)
    rdf_bkg = ROOT.RDataFrame(event_bkg)
    ROOT.RDF.Experimental.AddProgressBar(rdf_bkg)
    return rdf_sig, rdf_bkg

def get_values(counts):
    return np.array([c.GetValue() for c in counts])

def compute_roc_points(rdf_sig, rdf_bkg, cut_sig, cut_bkg, var, var_range, with_ge1=False):
    sig_nums, sig_dens, bkg_nums, bkg_dens = [], [], [], []
    for cut_val in var_range:
        v = f"abs({var}) <= {cut_val}" if not with_ge1 else f"abs({var}_with_GE1) <= {cut_val}"
        sig_num = rdf_sig.Filter(f"{cut_sig} & {v}").Count()
        sig_den = rdf_sig.Filter(cut_sig).Count()
        bkg_num = rdf_bkg.Filter(f"{cut_bkg} & {v}").Count()
        bkg_den = rdf_bkg.Filter(cut_bkg).Count()
        sig_nums.append(sig_num)
        sig_dens.append(sig_den)
        bkg_nums.append(bkg_num)
        bkg_dens.append(bkg_den)
    return get_values(sig_nums), get_values(sig_dens), get_values(bkg_nums), get_values(bkg_dens)

def compute_ba_points(rdf_sig, rdf_bkg, cut_sig, cut_bkg, bacutlist):
    sig_nums, sig_dens, bkg_nums, bkg_dens = [], [], [], []
    for ba_odd, ba_even in bacutlist:
        cut_ba_odd = f"((LCT_CSC_chamber%2 == 1) & (abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_odd}))"
        cut_ba_even = f"((LCT_CSC_chamber%2 == 0) & (abs((-(LCT_eighthstrip-LCT_match_GE1_padES_aligned))) <= {ba_even}))"
        ba_cut = f"({cut_ba_odd} | {cut_ba_even})"
        sig_num = rdf_sig.Filter(f"{cut_sig} & {ba_cut}").Count()
        sig_den = rdf_sig.Filter(cut_sig).Count()
        bkg_num = rdf_bkg.Filter(f"{cut_bkg} & {ba_cut}").Count()
        bkg_den = rdf_bkg.Filter(cut_bkg).Count()
        sig_nums.append(sig_num)
        sig_dens.append(sig_den)
        bkg_nums.append(bkg_num)
        bkg_dens.append(bkg_den)
    return get_values(sig_nums), get_values(sig_dens), get_values(bkg_nums), get_values(bkg_dens)

def plot_rocs(rocs, labels, outname, title, xaxis, yaxis):
    c = ROOT.TCanvas("c", "c", 800, 600)
    mg = ROOT.TMultiGraph()
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kMagenta, ROOT.kGreen+2, ROOT.kOrange+7]
    for i, (x, y) in enumerate(rocs):
        g = ROOT.TGraph(len(x), x, y)
        g.SetLineColor(colors[i % len(colors)])
        g.SetMarkerStyle(ROOT.kFullCircle)
        mg.Add(g, "lp")
    mg.SetTitle(title)
    mg.Draw("A")
    mg.GetXaxis().SetTitle(xaxis)
    mg.GetYaxis().SetTitle(yaxis)
    mg.SetMinimum(0.5)
    legend = ROOT.TLegend(0.6, 0.2, 0.85, 0.5)
    for i, label in enumerate(labels):
        legend.AddEntry(mg.GetListOfGraphs()[i], label)
    legend.SetBorderSize(0)
    legend.Draw()
    c.SaveAs(outname)
    print(f"Saved {outname}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sigFile', required=True, type=str)
    parser.add_argument('--bkgFile', required=True, type=str)
    parser.add_argument('--sigTree', default="GEMCSCBendingAngleTester/AllLCTs", type=str)
    parser.add_argument('--bkgTree', default="GEMCSCBendingAngleTester/HighestEMTFpT", type=str)
    parser.add_argument('--plotdir', default="plots/ROCCurveCombined/", type=str)
    parser.add_argument('--recopt', type=int, default=24)
    parser.add_argument('--emtfpt', type=int, default=22)
    parser.add_argument('--slope_max', type=int, default=15)
    parser.add_argument('--ba_var', default="LCT_slope", type=str)
    parser.add_argument('--ba_var_with_ge1', default="LCT_slope_with_GE1", type=str)
    parser.add_argument('--bacutlist', nargs='+', type=int, default=[1000, 1000, 16, 8, 15, 7, 12, 6, 11, 5, 9, 4])
    args = parser.parse_args()

    make_plotdir(args.plotdir)
    ROOT.gROOT.SetBatch(1)
    # tdrstyle.setTDRStyle() # Uncomment if you have tdrstyle

    ROOT.EnableImplicitMT(16)

    dRmatch = 0.4
    base_cut = "has_emtf_track_match & has_LCT_match1"
    residual_cut = "(abs(LCT_match_GE1_residual) < 20)"  # even/odd handled in cut below
    reco_cut = f"has_reco_l1_match & (abs(reco_l1_match_dR) < {dRmatch})"
    mode_cut = "((emtftrack_mode == 11) | (emtftrack_mode == 13) | (emtftrack_mode == 14) | (emtftrack_mode == 15))"
    zmu_cut = "(z_cand) & (abs(z_mass-91) < 10)"
    cut = f"{base_cut} & {reco_cut} & {mode_cut}"
    cut += " & ( ((LCT_CSC_chamber%2 == 0) & (abs(LCT_match_GE1_residual) < 20)) | ((LCT_CSC_chamber%2 == 1) & (abs(LCT_match_GE1_residual) < 40)) )"
    cut += " & LCT_CSC_ME1b == 1"

    cut_sig = cut + f" & {reco_cut} & reco_l1_match_pt >= {args.recopt} & (reco_l1_match_pt <= ({args.recopt} + 1)) & emtftrack_pt >= {args.emtfpt} & {zmu_cut}"
    cut_bkg = cut + f" & emtftrack_pt >= {args.emtfpt}"

    rdf_sig, rdf_bkg = get_rdfs(args.sigFile, args.bkgFile, args.sigTree, args.bkgTree)
    slope_range = range(args.slope_max)
    bacutlist = [args.bacutlist[i:i+2] for i in range(0, len(args.bacutlist), 2)]

    # ROC: CSC Alone (Slope)
    sig_nums, sig_dens, bkg_nums, bkg_dens = compute_roc_points(rdf_sig, rdf_bkg, cut_sig, cut_bkg, args.ba_var, slope_range, with_ge1=False)
    x_csc_alone = bkg_nums / bkg_dens
    y_csc_alone = sig_nums / sig_dens

    # ROC: CSC+GE1 Slope
    sig_nums, sig_dens, bkg_nums, bkg_dens = compute_roc_points(rdf_sig, rdf_bkg, cut_sig, cut_bkg, args.ba_var_with_ge1, slope_range, with_ge1=True)
    x_csc_ge1_slope = bkg_nums / bkg_dens
    y_csc_ge1_slope = sig_nums / sig_dens

    # ROC: CSC+GE1 Bending Angle
    sig_nums, sig_dens, bkg_nums, bkg_dens = compute_ba_points(rdf_sig, rdf_bkg, cut_sig, cut_bkg, bacutlist)
    x_csc_ge1 = bkg_nums / bkg_dens
    y_csc_ge1 = sig_nums / sig_dens

    # Plot all on the same canvas
    rocs = [
        (x_csc_alone, y_csc_alone),
        (x_csc_ge1, y_csc_ge1),
        (x_csc_ge1_slope, y_csc_ge1_slope),
    ]
    labels = [
        "CSC Alone, <= Slope",
        "CSC+GE1, <= BendingAngle",
        "CSC+GE1, <= Slope",
    ]
    outname = os.path.join(args.plotdir, f"ROC_recopt{args.recopt}_emtfpt{args.emtfpt}.pdf")
    plot_rocs(rocs, labels, outname,
              title=f"ROC Curves (RecoPt={args.recopt}, EMTFPt={args.emtfpt})",
              xaxis=f"Trigger Rate Remaining (EMTFPt > {args.emtfpt})",
              yaxis=f"Marginal Efficiency / 1GeV (RecoPt = {args.recopt})"
    )

if __name__ == "__main__":
    main()