import ROOT
import os
import argparse
import numpy as np
from scipy.stats import beta

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])
from GEM-CSC-trg-dev.scripts.plot_tool import *

def get_rdfs(sig_file, bkg_file, sig_tree, bkg_tree):
    f_sig = ROOT.TFile(sig_file)
    f_bkg = ROOT.TFile(bkg_file)
    event_sig = f_sig.Get(sig_tree)
    event_bkg = f_bkg.Get(bkg_tree)
    rdf_sig = ROOT.RDataFrame(event_sig)
    add_progress_bar(rdf_sig)
    rdf_bkg = ROOT.RDataFrame(event_bkg)
    add_progress_bar(rdf_bkg)
    return rdf_sig, rdf_bkg

def compute_roc_points(rdf_sig, rdf_bkg, cut_sig, cut_bkg, ba_var, ba_range):
    sig_nums = []
    sig_dens = []
    bkg_nums = []
    bkg_dens = []
    for ba_cut in ba_range:
        sig_num = rdf_sig.Filter(f"{cut_sig} & abs({ba_var}) <= {ba_cut}").Count()
        sig_den = rdf_sig.Filter(cut_sig).Count()
        bkg_num = rdf_bkg.Filter(f"{cut_bkg} & abs({ba_var}) <= {ba_cut}").Count()
        bkg_den = rdf_bkg.Filter(cut_bkg).Count()
        sig_nums.append(sig_num)
        sig_dens.append(sig_den)
        bkg_nums.append(bkg_num)
        bkg_dens.append(bkg_den)
    return sig_nums, sig_dens, bkg_nums, bkg_dens

def get_values(counts):
    return np.array([c.GetValue() for c in counts])

def plot_roc(x, y, xerr=None, yerr=None, labels=None, outname="roc.pdf", title="ROC Curve"):
    c = ROOT.TCanvas("c", "c", 800, 600)
    mg = ROOT.TMultiGraph()
    g = ROOT.TGraph(len(x), x, y)
    g.SetLineColor(ROOT.kBlue)
    g.SetMarkerStyle(ROOT.kCircle)
    mg.Add(g, "lp")
    if xerr is not None and yerr is not None:
        gerr = ROOT.TGraphErrors(len(x), x, y, xerr, yerr)
        gerr.SetMarkerStyle(ROOT.kCircle)
        mg.Add(gerr, "p")
    mg.SetTitle(title)
    mg.Draw("A*")
    c.SaveAs(outname)
    print(f"Saved {outname}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sigFile', required=True, type=str, help="Signal ROOT file")
    parser.add_argument('--bkgFile', required=True, type=str, help="Background ROOT file")
    parser.add_argument('--sigTree', default="GEMCSCBendingAngleTester/OnlyRecos", type=str)
    parser.add_argument('--bkgTree', default="GEMCSCBendingAngleTester/HighestEMTFpT", type=str)
    parser.add_argument('--plotdir', default="plots/ROCCurve/", type=str)
    parser.add_argument('--layer', type=int, default=1)
    parser.add_argument('--evenodd', choices=["odd", "even"], default="odd")
    parser.add_argument('--ba_var', default="LCT_BendingAngle_GE1", type=str)
    parser.add_argument('--ba_min', type=int, default=-6)
    parser.add_argument('--ba_max', type=int, default=50)
    parser.add_argument('--doBaMin', action='store_true')
    parser.add_argument('--plotdir', default="plots/ROCCurve/", type=str)
    args = parser.parse_args()

    make_plotdir(args.plotdir)
    set_batch()
    enable_mt(16)
    # Cuts and variables
    dRmatch = 0.4
    base_cut = f"has_emtf_track_match & has_LCT_match{args.layer}"
    mode_cut = "((emtftrack_mode == 11) | (emtftrack_mode == 13) | (emtftrack_mode == 14) | (emtftrack_mode == 15))"
    residual_cut = f"(abs(LCT_match_GE{args.layer}_residual) < 20)" if args.evenodd == "even" else f"(abs(LCT_match_GE{args.layer}_residual) < 40)"
    reco_cut = f"has_reco_l1_match & (abs(reco_l1_match_dR) < {dRmatch})"
    zmu_cut = "(z_cand) & (abs(z_mass-91) < 10) & (z_opposite_charge)"
    eta_cut = "(LCT_match_GE1_roll >= 0)"
    evenodd_cut = "(LCT_CSC_chamber%2 == 0)" if args.evenodd == "even" else "(LCT_CSC_chamber%2 == 1)"

    cut_sig = f"{base_cut} & {mode_cut} & {evenodd_cut} & {residual_cut} & {reco_cut} & {zmu_cut} & LCT_CSC_ME1b == 1 & {eta_cut}"
    cut_bkg = f"{base_cut} & {mode_cut} & {evenodd_cut} & {residual_cut} & LCT_CSC_ME1b == 1 & {eta_cut}"

    ba_range = range(args.ba_min, args.ba_max+1)

    rdf_sig, rdf_bkg = get_rdfs(args.sigFile, args.bkgFile, args.sigTree, args.bkgTree)
    sig_nums, sig_dens, bkg_nums, bkg_dens = compute_roc_points(rdf_sig, rdf_bkg, cut_sig, cut_bkg, args.ba_var, ba_range)

    x = get_values(bkg_nums) / get_values(bkg_dens)
    y = get_values(sig_nums) / get_values(sig_dens)

    # Optional: error bars (binomial errors)
    xerr = np.sqrt(x * (1 - x) / get_values(bkg_dens))
    yerr = np.sqrt(y * (1 - y) / get_values(sig_dens))

    outname = os.path.join(args.plotdir, f"ROC_{args.evenodd}_L{args.layer}.pdf")
    plot_roc(x, y, xerr, yerr, outname=outname, title=f"ROC Curve {args.evenodd} L{args.layer}")

if __name__ == "__main__":
    main()