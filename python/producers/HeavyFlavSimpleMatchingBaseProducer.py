import os
import itertools
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from ..helpers.utils import deltaR, closest, polarP4, sumP4, get_subjets, corrected_svmass, configLogger
from ..helpers.xgbHelper import XGBEnsemble
from ..helpers.nnHelper import convert_prob, ensemble
from ..helpers.jetmetCorrector import JetMETCorrector, rndSeed

import logging
logger = logging.getLogger('nano')
configLogger('nano', loglevel=logging.INFO)

lumi_dict = {"2016APV": 19.52, "2016": 16.81, "2017": 41.48, "2018": 59.83}


class _NullObject:
    '''An null object which does not store anything, and does not raise exception.'''

    def __bool__(self):
        return False

    def __nonzero__(self):
        return False

    def __getattr__(self, name):
        pass

    def __setattr__(self, name, value):
        pass


class METObject(Object):

    def p4(self):
        return polarP4(self, eta=None, mass=None)


class HeavyFlavSimpleMatchingBaseProducer(Module, object):

    def __init__(self, channel, **kwargs):
        self._channel = channel  # 'qcd', 'photon', 'inclusive', 'muon'
        self.year = int(kwargs['year'])
        self.jetType = kwargs.get('jetType', 'ak8').lower()

        if self.jetType == 'ak8':
            self._jetConeSize = 0.8
            self._fj_name = 'FatJet'
            self._sj_name = 'SubJet'
            self._fj_gen_name = 'GenJetAK8'
            self._sj_gen_name = 'SubGenJetAK8'

        elif self.jetType == 'ak15':
            self._jetConeSize = 1.5
            self._fj_name = 'AK15Puppi'
            self._sj_name = 'AK15PuppiSubJet'
            self._fj_gen_name = 'GenJetAK15'
            self._sj_gen_name = 'GenSubJetAK15'

        else:
            raise RuntimeError('Jet type %s is not recognized!' % self.jetType)

    def beginJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.out = wrappedOutputTree

        # NOTE: branch names must start with a lower case letter
        # check keep_and_drop_output.txt
        self.out.branch("year", "I")
        self.out.branch("lumiwgt", "F")
        self.out.branch("jetR", "F")
        # self.out.branch("passmetfilters", "O")
        self.out.branch("l1PreFiringWeight", "F")
        self.out.branch("l1PreFiringWeightUp", "F")
        self.out.branch("l1PreFiringWeightDown", "F")
        # self.out.branch("nlep", "I")
        self.out.branch("ht", "F")
        self.out.branch("met", "F")
        self.out.branch("metphi", "F")

        # Large-R jets
        for idx in [1, 2, 3]:
            prefix = 'fj_%d_' % idx

            # fatjet kinematics
            self.out.branch(prefix + "is_qualified", "O")
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "rawmass", "F")
            self.out.branch(prefix + "sdmass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "tau32", "F")

            # subjets
            self.out.branch(prefix + "deltaR_sj12", "F")
            self.out.branch(prefix + "sj1_pt", "F")
            self.out.branch(prefix + "sj1_eta", "F")
            self.out.branch(prefix + "sj1_phi", "F")
            self.out.branch(prefix + "sj1_rawmass", "F")
            self.out.branch(prefix + "sj2_pt", "F")
            self.out.branch(prefix + "sj2_eta", "F")
            self.out.branch(prefix + "sj2_phi", "F")
            self.out.branch(prefix + "sj2_rawmass", "F")

            # taggers
            self.out.branch(prefix + "globalParT3_QCD", "F")
            self.out.branch(prefix + "globalParT3_TopbWev", "F")
            self.out.branch(prefix + "globalParT3_TopbWmv", "F")
            self.out.branch(prefix + "globalParT3_TopbWq", "F")
            self.out.branch(prefix + "globalParT3_TopbWqq", "F")
            self.out.branch(prefix + "globalParT3_TopbWtauhv", "F")
            self.out.branch(prefix + "globalParT3_WvsQCD", "F")
            self.out.branch(prefix + "globalParT3_XWW3q", "F")
            self.out.branch(prefix + "globalParT3_XWW4q", "F")
            self.out.branch(prefix + "globalParT3_XWWqqev", "F")
            self.out.branch(prefix + "globalParT3_XWWqqmv", "F")
            self.out.branch(prefix + "globalParT3_Xbb", "F")
            self.out.branch(prefix + "globalParT3_Xcc", "F")
            self.out.branch(prefix + "globalParT3_Xcs", "F")
            self.out.branch(prefix + "globalParT3_Xqq", "F")
            self.out.branch(prefix + "globalParT3_Xtauhtaue", "F")
            self.out.branch(prefix + "globalParT3_Xtauhtauh", "F")
            self.out.branch(prefix + "globalParT3_Xtauhtaum", "F")
            self.out.branch(prefix + "globalParT3_massCorrGeneric", "F")
            self.out.branch(prefix + "globalParT3_massCorrX2p", "F")
            self.out.branch(prefix + "globalParT3_withMassTopvsQCD", "F")
            self.out.branch(prefix + "globalParT3_withMassWvsQCD", "F")
            self.out.branch(prefix + "globalParT3_withMassZvsQCD", "F")
            self.out.branch(prefix + "particleNetLegacy_QCD", "F")
            self.out.branch(prefix + "particleNetLegacy_Xbb", "F")
            self.out.branch(prefix + "particleNetLegacy_Xcc", "F")
            self.out.branch(prefix + "particleNetLegacy_Xqq", "F")
            self.out.branch(prefix + "particleNetLegacy_mass", "F")
            self.out.branch(prefix + "particleNetWithMass_H4qvsQCD", "F")
            self.out.branch(prefix + "particleNetWithMass_HbbvsQCD", "F")
            self.out.branch(prefix + "particleNetWithMass_HccvsQCD", "F")
            self.out.branch(prefix + "particleNetWithMass_QCD", "F")
            self.out.branch(prefix + "particleNetWithMass_TvsQCD", "F")
            self.out.branch(prefix + "particleNetWithMass_WvsQCD", "F")
            self.out.branch(prefix + "particleNetWithMass_ZvsQCD", "F")
            self.out.branch(prefix + "particleNet_QCD", "F")
            self.out.branch(prefix + "particleNet_QCD0HF", "F")
            self.out.branch(prefix + "particleNet_QCD1HF", "F")
            self.out.branch(prefix + "particleNet_QCD2HF", "F")
            self.out.branch(prefix + "particleNet_WVsQCD", "F")
            self.out.branch(prefix + "particleNet_XbbVsQCD", "F")
            self.out.branch(prefix + "particleNet_XccVsQCD", "F")
            self.out.branch(prefix + "particleNet_XggVsQCD", "F")
            self.out.branch(prefix + "particleNet_XqqVsQCD", "F")
            self.out.branch(prefix + "particleNet_XteVsQCD", "F")
            self.out.branch(prefix + "particleNet_XtmVsQCD", "F")
            self.out.branch(prefix + "particleNet_XttVsQCD", "F")
            self.out.branch(prefix + "particleNet_massCorr", "F")

            # matching variables
            if self.isMC:
                self.out.branch(prefix + "partonflavour", "I")
                self.out.branch(prefix + "sj1_nbhadrons", "I")
                self.out.branch(prefix + "sj1_nchadrons", "I")
                self.out.branch(prefix + "sj1_partonflavour", "I")
                self.out.branch(prefix + "sj2_nbhadrons", "I")
                self.out.branch(prefix + "sj2_nchadrons", "I")
                self.out.branch(prefix + "sj2_partonflavour", "I")

                # info of the closest hadGenH
                self.out.branch(prefix + "dr_H", "F")
                self.out.branch(prefix + "dr_H_daus", "F")
                self.out.branch(prefix + "H_pt", "F")
                self.out.branch(prefix + "H_decay", "I")

                # info of the closest hadGenZ
                self.out.branch(prefix + "dr_Z", "F")
                self.out.branch(prefix + "dr_Z_daus", "F")
                self.out.branch(prefix + "Z_pt", "F")
                self.out.branch(prefix + "Z_decay", "I")

                # info of the closest hadGenW
                self.out.branch(prefix + "dr_W", "F")
                self.out.branch(prefix + "dr_W_daus", "F")
                self.out.branch(prefix + "W_pt", "F")
                self.out.branch(prefix + "W_decay", "I")

                # info of the closest hadGenTop
                self.out.branch(prefix + "dr_T", "F")
                self.out.branch(prefix + "dr_T_b", "F")
                self.out.branch(prefix + "dr_T_Wq_max", "F")
                self.out.branch(prefix + "dr_T_Wq_min", "F")
                self.out.branch(prefix + "T_Wq_max_pdgId", "I")
                self.out.branch(prefix + "T_Wq_min_pdgId", "I")
                self.out.branch(prefix + "T_pt", "F")

                # info of the closest lepGenTop
                self.out.branch(prefix + "dr_LepT", "F")
                self.out.branch(prefix + "dr_LepT_b", "F")
                self.out.branch(prefix + "dr_LepT_l", "F")
                self.out.branch(prefix + "dr_LepT_l_pdgId", "I")
                self.out.branch(prefix + "LepT_pt", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def correctJetsAndMET(self, event):
        # correct Jets and MET
        event.idx = event._entry if event._tree._entrylist is None else event._tree._entrylist.GetEntry(event._entry)
        event._allJets = Collection(event, "Jet")
        event.met = METObject(event, "PFMET")
        event._allFatJets = Collection(event, self._fj_name)
        event.subjets = Collection(event, self._sj_name)  # do not sort subjets after updating!!

        # link fatjet to subjets and recompute softdrop mass
        for idx, fj in enumerate(event._allFatJets):
            fj.idx = idx
            fj.subjets = get_subjets(fj, event.subjets, ('subJetIdx1', 'subJetIdx2'))

        # no jetId requirements in nano v15
        event.fatjets = [fj for fj in event._allFatJets if fj.pt > 200 and abs(fj.eta) < 2.4]
        event.ak4jets = [j for j in event._allJets if j.pt > 25 and abs(j.eta) < 2.4]

        event.ht = sum([j.pt for j in event.ak4jets])

    def loadGenHistory(self, event, fatjets):
        # gen matching
        if not self.isMC:
            return

        try:
            genparts = event.genparts
        except RuntimeError as e:
            genparts = Collection(event, "GenPart")
            for idx, gp in enumerate(genparts):
                if 'dauIdx' not in gp.__dict__:
                    gp.dauIdx = []
                if gp.genPartIdxMother >= 0:
                    mom = genparts[gp.genPartIdxMother]
                    if 'dauIdx' not in mom.__dict__:
                        mom.dauIdx = [idx]
                    else:
                        mom.dauIdx.append(idx)
            event.genparts = genparts

        def isHadronic(gp):
            if len(gp.dauIdx) == 0:
                return False
                # raise ValueError('Particle has no daughters!')
            for idx in gp.dauIdx:
                if abs(genparts[idx].pdgId) < 6:
                    return True
            return False

        def getFinal(gp):
            for idx in gp.dauIdx:
                dau = genparts[idx]
                if dau.pdgId == gp.pdgId:
                    return getFinal(dau)
            return gp

        lepGenTops = []
        hadGenTops = []
        hadGenWs = []
        hadGenZs = []
        hadGenHs = []

        for gp in genparts:
            if gp.statusFlags & (1 << 13) == 0:
                continue
            if abs(gp.pdgId) == 6:
                for idx in gp.dauIdx:
                    dau = genparts[idx]
                    if abs(dau.pdgId) == 24:
                        genW = getFinal(dau)
                        gp.genW = genW
                        if isHadronic(genW):
                            hadGenTops.append(gp)
                        else:
                            lepGenTops.append(gp)
                    elif abs(dau.pdgId) in (1, 3, 5):
                        gp.genB = dau
            elif abs(gp.pdgId) == 24:
                if isHadronic(gp):
                    hadGenWs.append(gp)
            elif abs(gp.pdgId) == 23:
                if isHadronic(gp):
                    hadGenZs.append(gp)
            elif abs(gp.pdgId) == 25:
                if isHadronic(gp):
                    hadGenHs.append(gp)

        for parton in itertools.chain(lepGenTops, hadGenTops):
            parton.daus = (parton.genB, genparts[parton.genW.dauIdx[0]], genparts[parton.genW.dauIdx[1]])
            parton.genW.daus = parton.daus[1:]
        for parton in itertools.chain(hadGenWs, hadGenZs, hadGenHs):
            parton.daus = (genparts[parton.dauIdx[0]], genparts[parton.dauIdx[1]])

        for fj in fatjets:
            fj.genH, fj.dr_H = closest(fj, hadGenHs)
            fj.genZ, fj.dr_Z = closest(fj, hadGenZs)
            fj.genW, fj.dr_W = closest(fj, hadGenWs)
            fj.genT, fj.dr_T = closest(fj, hadGenTops)
            fj.genLepT, fj.dr_LepT = closest(fj, lepGenTops)


    def fillBaseEventInfo(self, event):
        self.out.fillBranch("jetR", self._jetConeSize)
        self.out.fillBranch("year", self.year)
        self.out.fillBranch("lumiwgt", lumi_dict[self.year])

        # met_filters = bool(
        #     event.Flag_goodVertices and
        #     event.Flag_globalSuperTightHalo2016Filter and
        #     event.Flag_HBHENoiseFilter and
        #     event.Flag_HBHENoiseIsoFilter and
        #     event.Flag_EcalDeadCellTriggerPrimitiveFilter and
        #     event.Flag_BadPFMuonFilter and
        #     event.Flag_BadPFMuonDzFilter and
        #     event.Flag_eeBadScFilter
        # )
        # if self.year in ("2017", "2018"):
        #     met_filters = met_filters and event.Flag_ecalBadCalibFilter
        # self.out.fillBranch("passmetfilters", met_filters)

        # self.out.fillBranch("nlep", len(event.looseLeptons))
        self.out.fillBranch("ht", event.ht)
        self.out.fillBranch("met", event.met.pt)
        self.out.fillBranch("metphi", event.met.phi)

    def _get_filler(self, obj):

        def filler(branch, value, default=0):
            self.out.fillBranch(branch, value if obj else default)

        return filler

    def fillFatJetInfo(self, event, fatjets):
        for idx in [1, 2, 3]:
            prefix = 'fj_%d_' % idx

            if len(fatjets) <= idx - 1:
                # fill zeros if fatjet fails probe selection
                for b in self.out._branches.keys():
                    if b.startswith(prefix):
                        self.out.fillBranch(b, 0)
                continue

            fj = fatjets[idx - 1]

            # fatjet kinematics
            self.out.fillBranch(prefix + "pt", fj.pt)
            self.out.fillBranch(prefix + "eta", fj.eta)
            self.out.fillBranch(prefix + "phi", fj.phi)
            self.out.fillBranch(prefix + "rawmass", fj.mass)
            self.out.fillBranch(prefix + "sdmass", fj.msoftdrop)
            self.out.fillBranch(prefix + "tau21", fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99)
            self.out.fillBranch(prefix + "tau32", fj.tau3 / fj.tau2 if fj.tau2 > 0 else 99)

            # subjets
            self.out.fillBranch(prefix + "deltaR_sj12", deltaR(*fj.subjets) if len(fj.subjets) == 2 else 99)
            for idx_sj, sj in enumerate(fj.subjets):
                prefix_sj = prefix + 'sj%d_' % (idx_sj + 1)
                self.out.fillBranch(prefix_sj + "pt", sj.pt)
                self.out.fillBranch(prefix_sj + "eta", sj.eta)
                self.out.fillBranch(prefix_sj + "phi", sj.phi)
                self.out.fillBranch(prefix_sj + "rawmass", sj.mass)

            # taggers (nano v15 branches)
            self.out.fillBranch(prefix + "globalParT3_QCD", fj.globalParT3_QCD)
            self.out.fillBranch(prefix + "globalParT3_TopbWev", fj.globalParT3_TopbWev)
            self.out.fillBranch(prefix + "globalParT3_TopbWmv", fj.globalParT3_TopbWmv)
            self.out.fillBranch(prefix + "globalParT3_TopbWq", fj.globalParT3_TopbWq)
            self.out.fillBranch(prefix + "globalParT3_TopbWqq", fj.globalParT3_TopbWqq)
            self.out.fillBranch(prefix + "globalParT3_TopbWtauhv", fj.globalParT3_TopbWtauhv)
            self.out.fillBranch(prefix + "globalParT3_WvsQCD", fj.globalParT3_WvsQCD)
            self.out.fillBranch(prefix + "globalParT3_XWW3q", fj.globalParT3_XWW3q)
            self.out.fillBranch(prefix + "globalParT3_XWW4q", fj.globalParT3_XWW4q)
            self.out.fillBranch(prefix + "globalParT3_XWWqqev", fj.globalParT3_XWWqqev)
            self.out.fillBranch(prefix + "globalParT3_XWWqqmv", fj.globalParT3_XWWqqmv)
            self.out.fillBranch(prefix + "globalParT3_Xbb", fj.globalParT3_Xbb)
            self.out.fillBranch(prefix + "globalParT3_Xcc", fj.globalParT3_Xcc)
            self.out.fillBranch(prefix + "globalParT3_Xcs", fj.globalParT3_Xcs)
            self.out.fillBranch(prefix + "globalParT3_Xqq", fj.globalParT3_Xqq)
            self.out.fillBranch(prefix + "globalParT3_Xtauhtaue", fj.globalParT3_Xtauhtaue)
            self.out.fillBranch(prefix + "globalParT3_Xtauhtauh", fj.globalParT3_Xtauhtauh)
            self.out.fillBranch(prefix + "globalParT3_Xtauhtaum", fj.globalParT3_Xtauhtaum)
            self.out.fillBranch(prefix + "globalParT3_massCorrGeneric", fj.globalParT3_massCorrGeneric)
            self.out.fillBranch(prefix + "globalParT3_massCorrX2p", fj.globalParT3_massCorrX2p)
            self.out.fillBranch(prefix + "globalParT3_withMassTopvsQCD", fj.globalParT3_withMassTopvsQCD)
            self.out.fillBranch(prefix + "globalParT3_withMassWvsQCD", fj.globalParT3_withMassWvsQCD)
            self.out.fillBranch(prefix + "globalParT3_withMassZvsQCD", fj.globalParT3_withMassZvsQCD)
            self.out.fillBranch(prefix + "particleNetLegacy_QCD", fj.particleNetLegacy_QCD)
            self.out.fillBranch(prefix + "particleNetLegacy_Xbb", fj.particleNetLegacy_Xbb)
            self.out.fillBranch(prefix + "particleNetLegacy_Xcc", fj.particleNetLegacy_Xcc)
            self.out.fillBranch(prefix + "particleNetLegacy_Xqq", fj.particleNetLegacy_Xqq)
            self.out.fillBranch(prefix + "particleNetLegacy_mass", fj.particleNetLegacy_mass)
            self.out.fillBranch(prefix + "particleNetWithMass_H4qvsQCD", fj.particleNetWithMass_H4qvsQCD)
            self.out.fillBranch(prefix + "particleNetWithMass_HbbvsQCD", fj.particleNetWithMass_HbbvsQCD)
            self.out.fillBranch(prefix + "particleNetWithMass_HccvsQCD", fj.particleNetWithMass_HccvsQCD)
            self.out.fillBranch(prefix + "particleNetWithMass_QCD", fj.particleNetWithMass_QCD)
            self.out.fillBranch(prefix + "particleNetWithMass_TvsQCD", fj.particleNetWithMass_TvsQCD)
            self.out.fillBranch(prefix + "particleNetWithMass_WvsQCD", fj.particleNetWithMass_WvsQCD)
            self.out.fillBranch(prefix + "particleNetWithMass_ZvsQCD", fj.particleNetWithMass_ZvsQCD)
            self.out.fillBranch(prefix + "particleNet_QCD", fj.particleNet_QCD)
            self.out.fillBranch(prefix + "particleNet_QCD0HF", fj.particleNet_QCD0HF)
            self.out.fillBranch(prefix + "particleNet_QCD1HF", fj.particleNet_QCD1HF)
            self.out.fillBranch(prefix + "particleNet_QCD2HF", fj.particleNet_QCD2HF)
            self.out.fillBranch(prefix + "particleNet_WVsQCD", fj.particleNet_WVsQCD)
            self.out.fillBranch(prefix + "particleNet_XbbVsQCD", fj.particleNet_XbbVsQCD)
            self.out.fillBranch(prefix + "particleNet_XccVsQCD", fj.particleNet_XccVsQCD)
            self.out.fillBranch(prefix + "particleNet_XggVsQCD", fj.particleNet_XggVsQCD)
            self.out.fillBranch(prefix + "particleNet_XqqVsQCD", fj.particleNet_XqqVsQCD)
            self.out.fillBranch(prefix + "particleNet_XteVsQCD", fj.particleNet_XteVsQCD)
            self.out.fillBranch(prefix + "particleNet_XtmVsQCD", fj.particleNet_XtmVsQCD)
            self.out.fillBranch(prefix + "particleNet_XttVsQCD", fj.particleNet_XttVsQCD)
            self.out.fillBranch(prefix + "particleNet_massCorr", fj.particleNet_massCorr)

            # matching variables
            if self.isMC:
                try:
                    sj1 = fj.subjets[0]
                except IndexError:
                    sj1 = None
                try:
                    sj2 = fj.subjets[1]
                except IndexError:
                    sj2 = None

                self.out.fillBranch(prefix + "sj1_nbhadrons", sj1.nBHadrons if sj1 else -1)
                self.out.fillBranch(prefix + "sj1_nchadrons", sj1.nCHadrons if sj1 else -1)
                self.out.fillBranch(prefix + "sj2_nbhadrons", sj2.nBHadrons if sj2 else -1)
                self.out.fillBranch(prefix + "sj2_nchadrons", sj2.nCHadrons if sj2 else -1)
                try:
                    self.out.fillBranch(prefix + "partonflavour", fj.partonFlavour)
                    self.out.fillBranch(prefix + "sj1_partonflavour", sj1.partonFlavour if sj1 else -1)
                    self.out.fillBranch(prefix + "sj2_partonflavour", sj2.partonFlavour if sj2 else -1)
                except RuntimeError:
                    self.out.fillBranch(prefix + "partonflavour", -1)
                    self.out.fillBranch(prefix + "sj1_partonflavour", -1)
                    self.out.fillBranch(prefix + "sj2_partonflavour", -1)

                # info of the closest hadGenH
                self.out.fillBranch(prefix + "dr_H", fj.dr_H)
                self.out.fillBranch(prefix + "dr_H_daus",
                                    max([deltaR(fj, dau) for dau in fj.genH.daus]) if fj.genH else 99)
                self.out.fillBranch(prefix + "H_pt", fj.genH.pt if fj.genH else -1)
                self.out.fillBranch(prefix + "H_decay", abs(fj.genH.daus[0].pdgId) if fj.genH else 0)

                # info of the closest hadGenZ
                self.out.fillBranch(prefix + "dr_Z", fj.dr_Z)
                self.out.fillBranch(prefix + "dr_Z_daus",
                                    max([deltaR(fj, dau) for dau in fj.genZ.daus]) if fj.genZ else 99)
                self.out.fillBranch(prefix + "Z_pt", fj.genZ.pt if fj.genZ else -1)
                self.out.fillBranch(prefix + "Z_decay", abs(fj.genZ.daus[0].pdgId) if fj.genZ else 0)

                # info of the closest hadGenW
                self.out.fillBranch(prefix + "dr_W", fj.dr_W)
                self.out.fillBranch(prefix + "dr_W_daus",
                                    max([deltaR(fj, dau) for dau in fj.genW.daus]) if fj.genW else 99)
                self.out.fillBranch(prefix + "W_pt", fj.genW.pt if fj.genW else -1)
                self.out.fillBranch(prefix + "W_decay", max([abs(d.pdgId) for d in fj.genW.daus]) if fj.genW else 0)

                # info of the closest hadGenTop
                drwq1, drwq2 = [deltaR(fj, dau) for dau in fj.genT.genW.daus] if fj.genT else [99, 99]
                wq1_pdgId, wq2_pdgId = [dau.pdgId for dau in fj.genT.genW.daus] if fj.genT else [0, 0]
                if drwq1 < drwq2:
                    drwq1, drwq2 = drwq2, drwq1
                    wq1_pdgId, wq2_pdgId = wq2_pdgId, wq1_pdgId
                self.out.fillBranch(prefix + "dr_T", fj.dr_T)
                self.out.fillBranch(prefix + "dr_T_b", deltaR(fj, fj.genT.genB) if fj.genT else 99)
                self.out.fillBranch(prefix + "dr_T_Wq_max", drwq1)
                self.out.fillBranch(prefix + "dr_T_Wq_min", drwq2)
                self.out.fillBranch(prefix + "T_Wq_max_pdgId", wq1_pdgId)
                self.out.fillBranch(prefix + "T_Wq_min_pdgId", wq2_pdgId)
                self.out.fillBranch(prefix + "T_pt", fj.genT.pt if fj.genT else -1)

                # info of the closest lepGenTop
                self.out.fillBranch(prefix + "dr_LepT", fj.dr_LepT)
                self.out.fillBranch(prefix + "dr_LepT_b", deltaR(fj, fj.genLepT.genB) if fj.genLepT else 99)
                self.out.fillBranch(prefix + "dr_LepT_l", deltaR(fj, fj.genLepT.genW.daus[0]) if fj.genLepT and abs(fj.genLepT.genW.daus[0].pdgId) in [11, 13, 15] \
                                    else deltaR(fj, fj.genLepT.genW.daus[1]) if fj.genLepT else 99)
                self.out.fillBranch(prefix + "dr_LepT_l_pdgId", fj.genLepT.genW.daus[0].pdgId if fj.genLepT and abs(fj.genLepT.genW.daus[0].pdgId) in [11, 13, 15] \
                                    else fj.genLepT.genW.daus[1].pdgId if fj.genLepT else 99)
                self.out.fillBranch(prefix + "LepT_pt", fj.genLepT.pt if fj.genLepT else -1)
