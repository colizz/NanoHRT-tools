from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from .HeavyFlavBaseProducer import HeavyFlavBaseProducer
from ..helpers.utils import deltaR, closest, polarP4, sumP4, get_subjets, corrected_svmass, configLogger
from ..helpers.triggerHelper import passTrigger
import ROOT
import math

class MuTaggedSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(MuTaggedSampleProducer, self).__init__(channel='mutagged', **kwargs)
        self._fill_sv = False  # not filling SV vars with standard way
        self._doJetCleaning = False  # no cleaning with leptons

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(MuTaggedSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        for idx in [1, 2]:
            prefix = 'fj_%d_' % idx

            # mu-tagged variables
            self.out.branch(prefix + "nmu", "I")
            self.out.branch(prefix + "dr_fjmatched_mu_sj1", "F")
            self.out.branch(prefix + "dr_fjmatched_mu_sj2", "F")
            self.out.branch(prefix + "dimuon_pt", "F")
            self.out.branch(prefix + "dimuon_mass", "F")
            self.out.branch(prefix + "is_sj_mutagged", "O")
            self.out.branch(prefix + "sj1_mu_pt", "F")
            self.out.branch(prefix + "sj2_mu_pt", "F")
            self.out.branch(prefix + "dr_mu_sj1", "F")
            self.out.branch(prefix + "dr_mu_sj2", "F")

            # SV variables
            self.out.branch(prefix + "nsv", "F")
            self.out.branch(prefix + "sv1_masscor", "F")
            self.out.branch(prefix + "sv1_maxdxysig_masscor", "F")
            self.out.branch(prefix + "sv_sum_masscor", "F")
            self.out.branch(prefix + "sv_vecsum_masscor", "F")
            self.out.branch(prefix + "sv_flightvecsum_masscor", "F")


    def selectSoftMuons(self, event):
        event.softMuons = []
        muons = Collection(event, "Muon")
        for mu in muons:
            # use inverted isolation cut for soft muons!
            if mu.pt > 5 and abs(mu.eta) < 2.4 and mu.tightId and mu.pfRelIso04_all > 0.15:
                event.softMuons.append(mu)

        event.softMuons.sort(key=lambda x: x.pt, reverse=True)


    def matchSoftMuonsToFatJets(self, event, fatjets):
        for fj in fatjets:
            # match soft muons to fatjet
            # drcut = min(0.4, 0.5 * deltaR(*fj.subjets)) if len(fj.subjets) == 2 else 0.4
            drcut = 0.4
            for sj in fj.subjets:
                sj.mu_list = []
                for mu in event.softMuons:
                    if deltaR(mu, sj) < drcut:
                        sj.mu_list.append(mu)  
                sj.nmu = len(sj.mu_list)

            fj.mu_list = []
            for mu in event.softMuons:
                if deltaR(mu, fj) < self._jetConeSize:
                    fj.mu_list.append(mu)  
            fj.nmu = len(fj.mu_list)

            if fj.nmu >= 2: # mutagged criteria
                fj.dimuon = sumP4(fj.mu_list[0], fj.mu_list[1])

            # check if both subjets are mutagged (fatjet itself is first mutagged)
            fj.is_sj_mutagged = len(fj.subjets) == 2 and fj.subjets[0].nmu >= 1 and fj.subjets[1].nmu >= 1 and len(set(fj.subjets[0].mu_list) | set(fj.subjets[1].mu_list)) >= 2

            # if fj.nmu >= 2 and len(fj.subjets) == 2:
            #     print(fj.nmu, fj.is_sj_mutagged, fj.subjets[0].nmu, fj.subjets[1].nmu, len(set(fj.subjets[0].mu_list) | set(fj.subjets[1].mu_list)))


    def matchSVToOnlyFatJets(self, event, fatjets):
        # match SV to fatjets
        for fj in fatjets:
            fj.sv_list = []
            for sv in event.secondary_vertices:
                if deltaR(sv, fj) < self._jetConeSize:
                    fj.sv_list.append(sv)
            fj.sv_list_dxysig_sorted = sorted(fj.sv_list, key=lambda x : x.dxySig, reverse=True)


    def fillMuTaggedFatJetInfo(self, event, fatjets):
        # fill mu-tagged jet info
        for idx in [1, 2]:
            prefix = 'fj_%d_' % idx

            if len(fatjets) <= idx - 1 or not fatjets[idx - 1].is_qualified:
                # already filled zeros
                continue

            fj = fatjets[idx - 1]

            self.out.fillBranch(prefix + "nmu", fj.nmu)
            self.out.fillBranch(prefix + "dr_fjmatched_mu_sj1", min([deltaR(mu, fj.subjets[0]) for mu in fj.mu_list]))
            self.out.fillBranch(prefix + "dr_fjmatched_mu_sj2", min([deltaR(mu, fj.subjets[1]) for mu in fj.mu_list]))
            self.out.fillBranch(prefix + "dimuon_pt", fj.dimuon.pt() if fj.nmu >= 2 else -99.)
            self.out.fillBranch(prefix + "dimuon_mass", fj.dimuon.M() if fj.nmu >= 2 else -99.)
            self.out.fillBranch(prefix + "is_sj_mutagged", fj.is_sj_mutagged)
            if fj.is_sj_mutagged:
                self.out.fillBranch(prefix + "sj1_mu_pt", fj.subjets[0].mu_list[0].pt)
                self.out.fillBranch(prefix + "sj2_mu_pt", fj.subjets[1].mu_list[0].pt)
                self.out.fillBranch(prefix + "dr_mu_sj1", min([deltaR(mu, fj.subjets[0]) for mu in fj.subjets[0].mu_list]))
                self.out.fillBranch(prefix + "dr_mu_sj2", min([deltaR(mu, fj.subjets[1]) for mu in fj.subjets[1].mu_list]))
            else:
                self.out.fillBranch(prefix + "sj1_mu_pt", -99.)
                self.out.fillBranch(prefix + "sj2_mu_pt", -99.)
                self.out.fillBranch(prefix + "dr_mu_sj1", -99.)
                self.out.fillBranch(prefix + "dr_mu_sj2", -99.)

            self.out.fillBranch(prefix + "nsv", len(fj.sv_list))
            self.out.fillBranch(prefix + "sv1_masscor", corrected_svmass(fj.sv_list[0]) if len(fj.sv_list) > 0 else -99.)
            self.out.fillBranch(prefix + "sv1_maxdxysig_masscor", corrected_svmass(fj.sv_list_dxysig_sorted[0]) if len(fj.sv_list) > 0 else -99.)
            self.out.fillBranch(prefix + "sv_sum_masscor", sum(corrected_svmass(sv) for sv in fj.sv_list) if len(fj.sv_list) > 0 else -99.)

            # vector sum mass (using p4 or flight vector as the direction)
            fj.sv_corr_veclist, fj.sv_corr_flightveclist = [], []
            for sv in fj.sv_list:
                svp4 = polarP4(sv)
                svmass_corr = corrected_svmass(sv)
                svp4.SetM(svmass_corr)
                fj.sv_corr_veclist.append(svp4)

                p, costhe, sinthe = svp4.P(), math.sin(sv.pAngle), math.cos(sv.pAngle)
                pcor = svmass_corr * p * costhe / math.sqrt(sv.mass ** 2 + (p * sinthe) ** 2)
                dist = math.sqrt(sv.x ** 2 + sv.y ** 2 + sv.z ** 2)
                fj.sv_corr_flightveclist.append(
                    ROOT.Math.XYZTVector(pcor * sv.x / dist, pcor * sv.y / dist, pcor * sv.z / dist, math.sqrt(pcor ** 2 + svmass_corr ** 2))
                )

            self.out.fillBranch(prefix + "sv_vecsum_masscor", sum(fj.sv_corr_veclist, ROOT.Math.PtEtaPhiMVector()).M() if len(fj.sv_list) > 0 else -99.)
            self.out.fillBranch(prefix + "sv_flightvecsum_masscor", sum(fj.sv_corr_flightveclist, ROOT.Math.XYZTVector()).M() if len(fj.sv_list) > 0 else -99.)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        # trigger selection
        if self.year == "2018":
            passBTagMuTrig = passTrigger(event, ['HLT_BTagMu_AK8Jet300_Mu5', 'HLT_BTagMu_AK8Jet300_Mu5_noalgo', 'HLT_BTagMu_AK4Jet300_Mu5', 'HLT_BTagMu_AK4Jet300_Mu5_noalgo'])
        elif self.year == "2017":
            passBTagMuTrig = passTrigger(event, ['HLT_BTagMu_AK8Jet300_Mu5', 'HLT_BTagMu_AK4Jet300_Mu5'])
        elif self.year in ["2016APV", "2016"]:
            passBTagMuTrig = passTrigger(event, ['HLT_BTagMu_AK8Jet300_Mu5', 'HLT_BTagMu_AK4Jet300_Mu5', 'HLT_BTagMu_Jet300_Mu5'])

        # accept events only passing the trigger
        if not passBTagMuTrig:
            return False

        self.selectLeptons(event)
        self.selectSoftMuons(event)
        self.correctJetsAndMET(event)

        # accept events with >=1 fatjet
        if len(event.fatjets) < 1:
            return False
        elif len(event.fatjets) == 1:
            probed_fatjets = event.fatjets[:1]
        else:
            probed_fatjets = event.fatjets[:2]

        # mu-tagged selection
        self.matchSoftMuonsToFatJets(event, probed_fatjets)

        # check qualification (>=1 matched muon)
        for fj in probed_fatjets:
            fj.is_qualified = (fj.nmu >= 1 and fj.pt >= 350 and abs(fj.eta) <= 2.4 and fj.msoftdrop >= 40)

        if len(probed_fatjets) == 1 and not probed_fatjets[0].is_qualified:
            return False
        if len(probed_fatjets) == 2 and not any([probed_fatjets[0].is_qualified, probed_fatjets[1].is_qualified]):
            return False

        self.loadGenHistory(event, probed_fatjets)
        self.evalTagger(event, probed_fatjets)
        self.evalMassRegression(event, probed_fatjets)

        # match to SV
        self.selectSV(event)
        self.matchSVToOnlyFatJets(event, probed_fatjets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probed_fatjets)
        self.fillMuTaggedFatJetInfo(event, probed_fatjets)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def MuTaggedTree_2016(): return MuTaggedSampleProducer(year=2016)
def MuTaggedTree_2017(): return MuTaggedSampleProducer(year=2017)
def MuTaggedTree_2018(): return MuTaggedSampleProducer(year=2018)
