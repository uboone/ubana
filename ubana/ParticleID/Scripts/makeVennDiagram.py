## n.b. this can't be run on the gpvms because of python libraries

## also, I know it's a mess --- still learning python!

from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
import ROOT
from ROOT import TFile, TTree, TH1D, TCanvas, gROOT, TPad, TGaxis, TColor, TLegend
import pylab

file1_path = "/home/adam/Documents/data_files/pid_bnbcos.root"
file1 = TFile(file1_path, 'READ')
if (file1.IsOpen()):
    print 'File ', file1_path, '  is open'
if (file1.IsOpen() == False):
    quit()

tree = file1.Get("pidvalid/pidTree")

true_protons = 0
true_muons = 0

cut_values = [
        2.0,   # track_Lmipoverp
        0.68,  # track_Lmumip0to1_nopionkaon
        88.0,  # track_chi2p
        -90.0, # track_chi2muminusp
        12.5,  # track_PIDA
        0.0    # track_depE_minus_rangeE_mu
        ]

overlap_lmipoverp_chi2muminusp_pida_p01 = 0
overlap_lmipoverp_chi2muminusp_pida_p12 = 0
overlap_lmipoverp_chi2muminusp_pida_p02 = 0
overlap_lmipoverp_chi2muminusp_pida_p012 = 0
nooverlap_lmipoverp_chi2muminusp_pida_p0 = 0
nooverlap_lmipoverp_chi2muminusp_pida_p1 = 0
nooverlap_lmipoverp_chi2muminusp_pida_p2 = 0
overlap_lmipoverp_chi2muminusp_pida_mu01 = 0
overlap_lmipoverp_chi2muminusp_pida_mu12 = 0
overlap_lmipoverp_chi2muminusp_pida_mu02 = 0
overlap_lmipoverp_chi2muminusp_pida_mu012 = 0
nooverlap_lmipoverp_chi2muminusp_pida_mu0 = 0
nooverlap_lmipoverp_chi2muminusp_pida_mu1 = 0
nooverlap_lmipoverp_chi2muminusp_pida_mu2 = 0

overlap_lmipoverp_chi2muminusp_depErangeE_p01 = 0
overlap_lmipoverp_chi2muminusp_depErangeE_p12 = 0
overlap_lmipoverp_chi2muminusp_depErangeE_p02 = 0
overlap_lmipoverp_chi2muminusp_depErangeE_p012 = 0
nooverlap_lmipoverp_chi2muminusp_depErangeE_p0 = 0
nooverlap_lmipoverp_chi2muminusp_depErangeE_p1 = 0
nooverlap_lmipoverp_chi2muminusp_depErangeE_p2 = 0
overlap_lmipoverp_chi2muminusp_depErangeE_mu01 = 0
overlap_lmipoverp_chi2muminusp_depErangeE_mu12 = 0
overlap_lmipoverp_chi2muminusp_depErangeE_mu02 = 0
overlap_lmipoverp_chi2muminusp_depErangeE_mu012 = 0
nooverlap_lmipoverp_chi2muminusp_depErangeE_mu0 = 0
nooverlap_lmipoverp_chi2muminusp_depErangeE_mu1 = 0
nooverlap_lmipoverp_chi2muminusp_depErangeE_mu2 = 0

overlap_pida_depErangeE_p01 = 0
nooverlap_pida_depErangeE_p0 = 0
nooverlap_pida_depErangeE_p1 = 0
overlap_pida_depErangeE_mu01 = 0
nooverlap_pida_depErangeE_mu0 = 0
nooverlap_pida_depErangeE_mu1 = 0

overlap_lmipoverp_lmumip0to1nopionkaon_p01 = 0
nooverlap_lmipoverp_lmumip0to1nopionkaon_p0 = 0
nooverlap_lmipoverp_lmumip0to1nopionkaon_p1 = 0
overlap_lmipoverp_lmumip0to1nopionkaon_mu01 = 0
nooverlap_lmipoverp_lmumip0to1nopionkaon_mu0 = 0
nooverlap_lmipoverp_lmumip0to1nopionkaon_mu1 = 0

overlap_chi2p_chi2muminusp_p01 = 0
nooverlap_chi2p_chi2muminusp_p0 = 0
nooverlap_chi2p_chi2muminusp_p1 = 0
overlap_chi2p_chi2muminusp_mu01 = 0
nooverlap_chi2p_chi2muminusp_mu0 = 0
nooverlap_chi2p_chi2muminusp_mu1 = 0

for entry in tree:

    # get variables we're interested in
    track_likelihood_p = max(entry.track_likelihood_fwd_p[2], entry.track_likelihood_bwd_p[2])
    track_likelihood_mu = max(entry.track_likelihood_fwd_mu[2], entry.track_likelihood_bwd_mu[2])

    if track_likelihood_p == 0 or track_likelihood_p == -999:
        Lmipoverp = 999 
    else:
        Lmipoverp = entry.track_likelihood_fwd_mip[2]/track_likelihood_p

    Lmumip_0to1_nopionkaon = (track_likelihood_mu + entry.track_likelihood_fwd_mip[2]) + (track_likelihood_mu + entry.track_likelihood_fwd_mip[2] + track_likelihood_p)
    Chi2muminusp = entry.track_Chi2Muon[2] - entry.track_Chi2Proton[2] 
    depErangeEMu = entry.track_depE[2] - entry.track_rangeE_mu

    if entry.true_PDG == 2212:
        true_protons = true_protons+1

        if Lmipoverp <= cut_values[0] and Chi2muminusp >= cut_values[3] and entry.track_PIDA_mean[2] >= cut_values[4]:
            overlap_lmipoverp_chi2muminusp_pida_p012 = overlap_lmipoverp_chi2muminusp_pida_p012+1
        if Lmipoverp > cut_values[0] and Chi2muminusp >= cut_values[3] and entry.track_PIDA_mean[2] >= cut_values[4]:
            overlap_lmipoverp_chi2muminusp_pida_p12 = overlap_lmipoverp_chi2muminusp_pida_p12+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp < cut_values[3] and entry.track_PIDA_mean[2] >= cut_values[4]:
            overlap_lmipoverp_chi2muminusp_pida_p02 = overlap_lmipoverp_chi2muminusp_pida_p02+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp >= cut_values[3] and entry.track_PIDA_mean[2] < cut_values[4]:
            overlap_lmipoverp_chi2muminusp_pida_p01 = overlap_lmipoverp_chi2muminusp_pida_p01+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp < cut_values[3] and entry.track_PIDA_mean[2] < cut_values[4]:
            nooverlap_lmipoverp_chi2muminusp_pida_p0 = nooverlap_lmipoverp_chi2muminusp_pida_p0+1
        if Lmipoverp > cut_values[0] and Chi2muminusp >= cut_values[3] and entry.track_PIDA_mean[2] < cut_values[4]:
            nooverlap_lmipoverp_chi2muminusp_pida_p1 = nooverlap_lmipoverp_chi2muminusp_pida_p1+1
        if Lmipoverp > cut_values[0] and Chi2muminusp < cut_values[3] and entry.track_PIDA_mean[2] >= cut_values[4]:
            nooverlap_lmipoverp_chi2muminusp_pida_p2 = nooverlap_lmipoverp_chi2muminusp_pida_p2+1

        if Lmipoverp <= cut_values[0] and Chi2muminusp >= cut_values[3] and depErangeEMu >= cut_values[5]:
            overlap_lmipoverp_chi2muminusp_depErangeE_p012 = overlap_lmipoverp_chi2muminusp_depErangeE_p012+1
        if Lmipoverp > cut_values[0] and Chi2muminusp >= cut_values[3] and depErangeEMu >= cut_values[5]:
            overlap_lmipoverp_chi2muminusp_depErangeE_p12 = overlap_lmipoverp_chi2muminusp_depErangeE_p12+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp < cut_values[3] and depErangeEMu >= cut_values[5]:
            overlap_lmipoverp_chi2muminusp_depErangeE_p02 = overlap_lmipoverp_chi2muminusp_depErangeE_p02+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp >= cut_values[3] and depErangeEMu < cut_values[5]:
            overlap_lmipoverp_chi2muminusp_depErangeE_p01 = overlap_lmipoverp_chi2muminusp_depErangeE_p01+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp < cut_values[3] and depErangeEMu < cut_values[5]:
            nooverlap_lmipoverp_chi2muminusp_depErangeE_p0 = nooverlap_lmipoverp_chi2muminusp_depErangeE_p0+1
        if Lmipoverp > cut_values[0] and Chi2muminusp >= cut_values[3] and depErangeEMu < cut_values[5]:
            nooverlap_lmipoverp_chi2muminusp_depErangeE_p1 = nooverlap_lmipoverp_chi2muminusp_depErangeE_p1+1
        if Lmipoverp > cut_values[0] and Chi2muminusp < cut_values[3] and depErangeEMu >= cut_values[5]:
            nooverlap_lmipoverp_chi2muminusp_depErangeE_p2 = nooverlap_lmipoverp_chi2muminusp_depErangeE_p2+1

        if entry.track_PIDA_mean[2] >= cut_values[4] and depErangeEMu >= cut_values[5]:
            overlap_pida_depErangeE_p01 = overlap_pida_depErangeE_p01+1
        if entry.track_PIDA_mean[2] >= cut_values[4] and depErangeEMu < cut_values[5]:
            nooverlap_pida_depErangeE_p0 = nooverlap_pida_depErangeE_p0+1
        if entry.track_PIDA_mean[2] < cut_values[4] and depErangeEMu >= cut_values[5]:
            nooverlap_pida_depErangeE_p1 = nooverlap_pida_depErangeE_p1+1

        if Lmipoverp <= cut_values[0] and Lmumip_0to1_nopionkaon <= cut_values[1]:
            overlap_lmipoverp_lmumip0to1nopionkaon_p01 = overlap_lmipoverp_lmumip0to1nopionkaon_p01+1
        if Lmipoverp <= cut_values[0] and Lmumip_0to1_nopionkaon > cut_values[1]:
            nooverlap_lmipoverp_lmumip0to1nopionkaon_p0 = nooverlap_lmipoverp_lmumip0to1nopionkaon_p0+1
        if Lmipoverp > cut_values[0] and Lmumip_0to1_nopionkaon <= cut_values[1]:
            nooverlap_lmipoverp_lmumip0to1nopionkaon_p1 = nooverlap_lmipoverp_lmumip0to1nopionkaon_p1+1

        if entry.track_Chi2Proton[2] <= cut_values[2] and Chi2muminusp >= cut_values[3]:
            overlap_chi2p_chi2muminusp_p01 = overlap_chi2p_chi2muminusp_p01+1
        if entry.track_Chi2Proton[2] <= cut_values[2] and Chi2muminusp < cut_values[3]:
            nooverlap_chi2p_chi2muminusp_p0 = nooverlap_chi2p_chi2muminusp_p0+1
        if entry.track_Chi2Proton[2] > cut_values[2] and Chi2muminusp >= cut_values[3]:
            nooverlap_chi2p_chi2muminusp_p1 = nooverlap_chi2p_chi2muminusp_p1+1



    if entry.true_PDG == 13:
        true_muons = true_muons+1

        if Lmipoverp > cut_values[0] and Chi2muminusp < cut_values[3] and entry.track_PIDA_mean[2] < cut_values[4]:
            overlap_lmipoverp_chi2muminusp_pida_mu012 = overlap_lmipoverp_chi2muminusp_pida_mu012+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp < cut_values[3] and entry.track_PIDA_mean[2] < cut_values[4]:
            overlap_lmipoverp_chi2muminusp_pida_mu12 = overlap_lmipoverp_chi2muminusp_pida_mu12+1
        if Lmipoverp > cut_values[0] and Chi2muminusp >= cut_values[3] and entry.track_PIDA_mean[2] < cut_values[4]:
            overlap_lmipoverp_chi2muminusp_pida_mu02 = overlap_lmipoverp_chi2muminusp_pida_mu02+1
        if Lmipoverp > cut_values[0] and Chi2muminusp < cut_values[3] and entry.track_PIDA_mean[2] >= cut_values[4]:
            overlap_lmipoverp_chi2muminusp_pida_mu01 = overlap_lmipoverp_chi2muminusp_pida_mu01+1
        if Lmipoverp > cut_values[0] and Chi2muminusp >= cut_values[3] and entry.track_PIDA_mean[2] >= cut_values[4]:
            nooverlap_lmipoverp_chi2muminusp_pida_mu0 = nooverlap_lmipoverp_chi2muminusp_pida_mu0+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp < cut_values[3] and entry.track_PIDA_mean[2] >= cut_values[4]:
            nooverlap_lmipoverp_chi2muminusp_pida_mu1 = nooverlap_lmipoverp_chi2muminusp_pida_mu1+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp >= cut_values[3] and entry.track_PIDA_mean[2] < cut_values[4]:
            nooverlap_lmipoverp_chi2muminusp_pida_mu2 = nooverlap_lmipoverp_chi2muminusp_pida_mu2+1

        if Lmipoverp > cut_values[0] and Chi2muminusp < cut_values[3] and depErangeEMu < cut_values[5]:
            overlap_lmipoverp_chi2muminusp_depErangeE_mu012 = overlap_lmipoverp_chi2muminusp_depErangeE_mu012+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp < cut_values[3] and depErangeEMu < cut_values[5]:
            overlap_lmipoverp_chi2muminusp_depErangeE_mu12 = overlap_lmipoverp_chi2muminusp_depErangeE_mu12+1
        if Lmipoverp > cut_values[0] and Chi2muminusp >= cut_values[3] and depErangeEMu < cut_values[5]:
            overlap_lmipoverp_chi2muminusp_depErangeE_mu02 = overlap_lmipoverp_chi2muminusp_depErangeE_mu02+1
        if Lmipoverp > cut_values[0] and Chi2muminusp < cut_values[3] and depErangeEMu >= cut_values[5]:
            overlap_lmipoverp_chi2muminusp_depErangeE_mu01 = overlap_lmipoverp_chi2muminusp_depErangeE_mu01+1
        if Lmipoverp > cut_values[0] and Chi2muminusp >= cut_values[3] and depErangeEMu >= cut_values[5]:
            nooverlap_lmipoverp_chi2muminusp_depErangeE_mu0 = nooverlap_lmipoverp_chi2muminusp_depErangeE_mu0+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp < cut_values[3] and depErangeEMu >= cut_values[5]:
            nooverlap_lmipoverp_chi2muminusp_depErangeE_mu1 = nooverlap_lmipoverp_chi2muminusp_depErangeE_mu1+1
        if Lmipoverp <= cut_values[0] and Chi2muminusp >= cut_values[3] and depErangeEMu < cut_values[5]:
            nooverlap_lmipoverp_chi2muminusp_depErangeE_mu2 = nooverlap_lmipoverp_chi2muminusp_depErangeE_mu2+1

        if entry.track_PIDA_mean[2] < cut_values[4] and depErangeEMu < cut_values[5]:
            overlap_pida_depErangeE_mu01 = overlap_pida_depErangeE_mu01+1
        if entry.track_PIDA_mean[2] < cut_values[4] and depErangeEMu >= cut_values[5]:
            nooverlap_pida_depErangeE_mu0 = nooverlap_pida_depErangeE_mu0+1
        if entry.track_PIDA_mean[2] >= cut_values[4] and depErangeEMu < cut_values[5]:
            nooverlap_pida_depErangeE_mu1 = nooverlap_pida_depErangeE_mu1+1
        
        if Lmipoverp > cut_values[0] and Lmumip_0to1_nopionkaon > cut_values[1]:
            overlap_lmipoverp_lmumip0to1nopionkaon_mu01 = overlap_lmipoverp_lmumip0to1nopionkaon_mu01+1
        if Lmipoverp > cut_values[0] and Lmumip_0to1_nopionkaon <= cut_values[1]:
            nooverlap_lmipoverp_lmumip0to1nopionkaon_mu0 = nooverlap_lmipoverp_lmumip0to1nopionkaon_mu0+1
        if Lmipoverp <= cut_values[0] and Lmumip_0to1_nopionkaon > cut_values[1]:
            nooverlap_lmipoverp_lmumip0to1nopionkaon_mu1 = nooverlap_lmipoverp_lmumip0to1nopionkaon_mu1+1

        if entry.track_Chi2Proton[2] > cut_values[2] and Chi2muminusp < cut_values[3]:
            overlap_chi2p_chi2muminusp_mu01 = overlap_chi2p_chi2muminusp_mu01+1
        if entry.track_Chi2Proton[2] > cut_values[2] and Chi2muminusp >= cut_values[3]:
            nooverlap_chi2p_chi2muminusp_mu0 = nooverlap_chi2p_chi2muminusp_mu0+1
        if entry.track_Chi2Proton[2] <= cut_values[2] and Chi2muminusp < cut_values[3]:
            nooverlap_chi2p_chi2muminusp_mu1 = nooverlap_chi2p_chi2muminusp_mu1+1


# make venn diagrams

# three dimensional one -- with PIDA

figure_lmipovrp_chi2muminusp_pida, axes_lmipovrp_chi2muminusp_pida = plt.subplots(1,2)

subsets_lmipovrp_chi2muminusp_pida_p = (nooverlap_lmipoverp_chi2muminusp_pida_p0, nooverlap_lmipoverp_chi2muminusp_pida_p1, overlap_lmipoverp_chi2muminusp_pida_p01, nooverlap_lmipoverp_chi2muminusp_pida_p2, overlap_lmipoverp_chi2muminusp_pida_p02, overlap_lmipoverp_chi2muminusp_pida_p12, overlap_lmipoverp_chi2muminusp_pida_p012)
areas_lmipovrp_chi2muminusp_pida_p = [1,1,1,1,1,1,1]

for i in range (0,7):
    if subsets_lmipovrp_chi2muminusp_pida_p[i] < 0.03*true_protons:
        areas_lmipovrp_chi2muminusp_pida_p[i] = 0.03*true_protons
    else:
        areas_lmipovrp_chi2muminusp_pida_p[i] = subsets_lmipovrp_chi2muminusp_pida_p[i]

venn3_unweighted(subsets_lmipovrp_chi2muminusp_pida_p, set_labels = (r'$L_{MIP}/L_{p}$', r'$\chi^{2}_{\mu-p}$', 'PIDA'), ax=axes_lmipovrp_chi2muminusp_pida[0], subset_areas=areas_lmipovrp_chi2muminusp_pida_p)

subsets_lmipovrp_chi2muminusp_pida_mu = (nooverlap_lmipoverp_chi2muminusp_pida_mu0, nooverlap_lmipoverp_chi2muminusp_pida_mu1, overlap_lmipoverp_chi2muminusp_pida_mu01, nooverlap_lmipoverp_chi2muminusp_pida_mu2, overlap_lmipoverp_chi2muminusp_pida_mu02, overlap_lmipoverp_chi2muminusp_pida_mu12, overlap_lmipoverp_chi2muminusp_pida_mu012)
areas_lmipovrp_chi2muminusp_pida_mu = [1,1,1,1,1,1,1]

for i in range (0,7):
    if subsets_lmipovrp_chi2muminusp_pida_mu[i] < 0.03*true_muons:
        areas_lmipovrp_chi2muminusp_pida_mu[i] = 0.03*true_muons
    else:
        areas_lmipovrp_chi2muminusp_pida_mu[i] = subsets_lmipovrp_chi2muminusp_pida_mu[i]



venn3_unweighted(subsets_lmipovrp_chi2muminusp_pida_mu, set_labels = (r'$L_{MIP}/L_{p}$', r'$\chi^{2}_{\mu-p}$', 'PIDA'), ax=axes_lmipovrp_chi2muminusp_pida[1], subset_areas=areas_lmipovrp_chi2muminusp_pida_mu)

pylab.text(0.5, 1.1, 'True Selected Protons', transform=axes_lmipovrp_chi2muminusp_pida[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, 1.1, 'True Selected Muons', transform=axes_lmipovrp_chi2muminusp_pida[1].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.1, '%s protons in sample'%(true_protons), transform=axes_lmipovrp_chi2muminusp_pida[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.1, '%s muons in sample'%(true_muons), transform=axes_lmipovrp_chi2muminusp_pida[1].transAxes, horizontalalignment='center', fontsize = 11)

plt.show()
figure_lmipovrp_chi2muminusp_pida.savefig("lmipoverp_chi2muminusp_pida_overlap.png")

# three dimensional one -- with depERangeE

figure_lmipovrp_chi2muminusp_depErangeE, axes_lmipovrp_chi2muminusp_depErangeE = plt.subplots(1,2)

subsets_lmipovrp_chi2muminusp_depErangeE_p = (nooverlap_lmipoverp_chi2muminusp_depErangeE_p0, nooverlap_lmipoverp_chi2muminusp_depErangeE_p1, overlap_lmipoverp_chi2muminusp_depErangeE_p01, nooverlap_lmipoverp_chi2muminusp_depErangeE_p2, overlap_lmipoverp_chi2muminusp_depErangeE_p02, overlap_lmipoverp_chi2muminusp_depErangeE_p12, overlap_lmipoverp_chi2muminusp_depErangeE_p012)
areas_lmipovrp_chi2muminusp_depErangeE_p = [1,1,1,1,1,1,1]

for i in range (0,7):
    if subsets_lmipovrp_chi2muminusp_depErangeE_p[i] < 0.03*true_protons:
        areas_lmipovrp_chi2muminusp_depErangeE_p[i] = 0.03*true_protons
    else:
        areas_lmipovrp_chi2muminusp_depErangeE_p[i] = subsets_lmipovrp_chi2muminusp_depErangeE_p[i]

venn3_unweighted(subsets_lmipovrp_chi2muminusp_depErangeE_p, set_labels = (r'$L_{MIP}/L_{p}$', r'$\chi^{2}_{\mu-p}$', r'$E_{DEP}-E_{RANGE}$'), ax=axes_lmipovrp_chi2muminusp_depErangeE[0], subset_areas=areas_lmipovrp_chi2muminusp_depErangeE_p)

subsets_lmipovrp_chi2muminusp_depErangeE_mu = (nooverlap_lmipoverp_chi2muminusp_depErangeE_mu0, nooverlap_lmipoverp_chi2muminusp_depErangeE_mu1, overlap_lmipoverp_chi2muminusp_depErangeE_mu01, nooverlap_lmipoverp_chi2muminusp_depErangeE_mu2, overlap_lmipoverp_chi2muminusp_depErangeE_mu02, overlap_lmipoverp_chi2muminusp_depErangeE_mu12, overlap_lmipoverp_chi2muminusp_depErangeE_mu012)
areas_lmipovrp_chi2muminusp_depErangeE_mu = [1,1,1,1,1,1,1]

for i in range (0,7):
    if subsets_lmipovrp_chi2muminusp_depErangeE_mu[i] < 0.03*true_muons:
        areas_lmipovrp_chi2muminusp_depErangeE_mu[i] = 0.03*true_muons
    else:
        areas_lmipovrp_chi2muminusp_depErangeE_mu[i] = subsets_lmipovrp_chi2muminusp_depErangeE_mu[i]



venn3_unweighted(subsets_lmipovrp_chi2muminusp_depErangeE_mu, set_labels = (r'$L_{MIP}/L_{p}$', r'$\chi^{2}_{\mu-p}$', r'$E_{DEP}-E_{RANGE}$'), ax=axes_lmipovrp_chi2muminusp_depErangeE[1], subset_areas=areas_lmipovrp_chi2muminusp_depErangeE_mu)

pylab.text(0.5, 1.1, 'True Selected Protons', transform=axes_lmipovrp_chi2muminusp_depErangeE[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, 1.1, 'True Selected Muons', transform=axes_lmipovrp_chi2muminusp_depErangeE[1].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.1, '%s protons in sample'%(true_protons), transform=axes_lmipovrp_chi2muminusp_depErangeE[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.1, '%s muons in sample'%(true_muons), transform=axes_lmipovrp_chi2muminusp_depErangeE[1].transAxes, horizontalalignment='center', fontsize = 11)

plt.show()
figure_lmipovrp_chi2muminusp_depErangeE.savefig("lmipoverp_chi2muminusp_depErangeE_overlap.png")


# two dimensional one - PIDa depErangeE

figure_pida_depErangE, axes_pida_depErangE = plt.subplots(1,2)

subsets_pida_depErangE_p = (nooverlap_pida_depErangeE_p0, nooverlap_pida_depErangeE_p1, overlap_pida_depErangeE_p01)
areas_pida_depErangE_p = [1,1,1]

for i in range (0,3):
    if subsets_pida_depErangE_p[i] < 0.03*true_protons:
        areas_pida_depErangE_p[i] = 0.03*true_protons
    else:
        areas_pida_depErangE_p[i] = subsets_pida_depErangE_p[i]

venn2_unweighted(subsets_pida_depErangE_p, set_labels = ('PIDA', r'$E_{DEP} - E_{RANGE}$'), ax=axes_pida_depErangE[0], subset_areas=areas_pida_depErangE_p)

subsets_pida_depErangE_mu = (nooverlap_pida_depErangeE_mu0, nooverlap_pida_depErangeE_mu1, overlap_pida_depErangeE_mu01)
areas_pida_depErangE_mu = [1,1,1]

for i in range (0,3):
    if subsets_pida_depErangE_mu[i] < 0.03*true_muons:
        areas_pida_depErangE_mu[i] = 0.03*true_muons
    else:
        areas_pida_depErangE_mu[i] = subsets_pida_depErangE_mu[i]



venn2_unweighted(subsets_pida_depErangE_mu, set_labels = ('PIDA', r'$E_{DEP} - E_{RANGE}$'), ax=axes_pida_depErangE[1], subset_areas=areas_pida_depErangE_mu)

pylab.text(0.5, 1.1, 'True Selected Protons', transform=axes_pida_depErangE[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, 1.1, 'True Selected Muons', transform=axes_pida_depErangE[1].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.2, '%s protons in sample'%(true_protons), transform=axes_pida_depErangE[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.2, '%s muons in sample'%(true_muons), transform=axes_pida_depErangE[1].transAxes, horizontalalignment='center', fontsize = 11)

plt.show()
figure_pida_depErangE.savefig("pida_depErangeE_overlap.png")

# two dimensional one - likelihood comparisons

figure_lmipovrp_lmumip0to1nopionkaon, axes_lmipovrp_lmumip0to1nopionkaon = plt.subplots(1,2)

subsets_lmipovrp_lmumip0to1nopionkaon_p = (nooverlap_lmipoverp_lmumip0to1nopionkaon_p0, nooverlap_lmipoverp_lmumip0to1nopionkaon_p1, overlap_lmipoverp_lmumip0to1nopionkaon_p01)
areas_lmipovrp_lmumip0to1nopionkaon_p = [1,1,1]

for i in range (0,3):
    if subsets_lmipovrp_lmumip0to1nopionkaon_p[i] < 0.03*true_protons:
        areas_lmipovrp_lmumip0to1nopionkaon_p[i] = 0.03*true_protons
    else:
        areas_lmipovrp_lmumip0to1nopionkaon_p[i] = subsets_lmipovrp_lmumip0to1nopionkaon_p[i]

venn2_unweighted(subsets_lmipovrp_lmumip0to1nopionkaon_p, set_labels = (r'$L_{MIP}/L_{p}$', r'($L_{MIP}+L_{\mu})/(L_{p}+L_{MIP}+L_{\mu})$'), ax=axes_lmipovrp_lmumip0to1nopionkaon[0], subset_areas=areas_lmipovrp_lmumip0to1nopionkaon_p)

subsets_lmipovrp_lmumip0to1nopionkaon_mu = (nooverlap_lmipoverp_lmumip0to1nopionkaon_mu0, nooverlap_lmipoverp_lmumip0to1nopionkaon_mu1, overlap_lmipoverp_lmumip0to1nopionkaon_mu01)
areas_lmipovrp_lmumip0to1nopionkaon_mu = [1,1,1]

for i in range (0,3):
    if subsets_lmipovrp_lmumip0to1nopionkaon_mu[i] < 0.03*true_muons:
        areas_lmipovrp_lmumip0to1nopionkaon_mu[i] = 0.03*true_muons
    else:
        areas_lmipovrp_lmumip0to1nopionkaon_mu[i] = subsets_lmipovrp_lmumip0to1nopionkaon_mu[i]



venn2_unweighted(subsets_lmipovrp_lmumip0to1nopionkaon_mu, set_labels = (r'$L_{MIP}/L_{p}$', r'($L_{MIP}+L_{\mu})/(L_{p}+L_{MIP}+L_{\mu})$'), ax=axes_lmipovrp_lmumip0to1nopionkaon[1], subset_areas=areas_lmipovrp_lmumip0to1nopionkaon_mu)

pylab.text(0.5, 1.1, 'True Selected Protons', transform=axes_lmipovrp_lmumip0to1nopionkaon[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, 1.1, 'True Selected Muons', transform=axes_lmipovrp_lmumip0to1nopionkaon[1].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.2, '%s protons in sample'%(true_protons), transform=axes_lmipovrp_lmumip0to1nopionkaon[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.2, '%s muons in sample'%(true_muons), transform=axes_lmipovrp_lmumip0to1nopionkaon[1].transAxes, horizontalalignment='center', fontsize = 11)

plt.show()
figure_lmipovrp_lmumip0to1nopionkaon.savefig("lmipoverp_lmumip0to1nopionkaon_overlap.png")

# two dimensional one - likelihood comparisons

figure_chi2p_chi2muminusp, axes_chi2p_chi2muminusp = plt.subplots(1,2)

subsets_chi2p_chi2muminusp_p = (nooverlap_chi2p_chi2muminusp_p0, nooverlap_chi2p_chi2muminusp_p1, overlap_chi2p_chi2muminusp_p01)
areas_chi2p_chi2muminusp_p = [1,1,1]

for i in range (0,3):
    if subsets_chi2p_chi2muminusp_p[i] < 0.03*true_protons:
        areas_chi2p_chi2muminusp_p[i] = 0.03*true_protons
    else:
        areas_chi2p_chi2muminusp_p[i] = subsets_chi2p_chi2muminusp_p[i]

venn2_unweighted(subsets_chi2p_chi2muminusp_p, set_labels = (r'$\chi^{2}_P$', r'$\chi^{2}_{\mu-p}$'), ax=axes_chi2p_chi2muminusp[0], subset_areas=areas_chi2p_chi2muminusp_p)

subsets_chi2p_chi2muminusp_mu = (nooverlap_chi2p_chi2muminusp_mu0, nooverlap_chi2p_chi2muminusp_mu1, overlap_chi2p_chi2muminusp_mu01)
areas_chi2p_chi2muminusp_mu = [1,1,1]

for i in range (0,3):
    if subsets_chi2p_chi2muminusp_mu[i] < 0.03*true_muons:
        areas_chi2p_chi2muminusp_mu[i] = 0.03*true_muons
    else:
        areas_chi2p_chi2muminusp_mu[i] = subsets_chi2p_chi2muminusp_mu[i]



venn2_unweighted(subsets_chi2p_chi2muminusp_mu, set_labels = (r'$\chi^{2}_P$', r'$\chi^{2}_{\mu-p}$'), ax=axes_chi2p_chi2muminusp[1], subset_areas=areas_chi2p_chi2muminusp_mu)

pylab.text(0.5, 1.1, 'True Selected Protons', transform=axes_chi2p_chi2muminusp[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, 1.1, 'True Selected Muons', transform=axes_chi2p_chi2muminusp[1].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.2, '%s protons in sample'%(true_protons), transform=axes_chi2p_chi2muminusp[0].transAxes, horizontalalignment='center', fontsize = 11)
pylab.text(0.5, -0.2, '%s muons in sample'%(true_muons), transform=axes_chi2p_chi2muminusp[1].transAxes, horizontalalignment='center', fontsize = 11)

plt.show()
figure_chi2p_chi2muminusp.savefig("chi2p_chi2muminusp_overlap.png")
