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

lumi_dict = {
    2015: 19.52,
    2016: 16.81, 
    2017: 41.48, 
    2018: 59.83, 
    2022: 7971.4 / 1000,  # 5.0707 + 3.0063
    "2022EE": 26337.0 / 1000,
    2023: 17981.0 / 1000,  # runC
    "2023BPix": 9516.0 / 1000,  # runD
}


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


class HeavyFlavBaseProducer(Module, object):

    def __init__(self, channel, **kwargs):
        self._channel = channel  # 'qcd', 'photon', 'inclusive', 'muon'
        try:
            self.year = int(kwargs['year'])
        except ValueError:
            self.year = int(kwargs['year'][:5])
        self.is_run3 = (self.year >= 2022)
        self.year_str = kwargs['year']
        self.jetType = kwargs.get('jetType', 'ak8').lower()
        self._jmeSysts = {'jec': False, 'jes': None, 'jes_source': '', 'jes_uncertainty_file_prefix': '',
                          'jer': None, 'jmr': None, 'met_unclustered': None, 'smearMET': True, 'applyHEMUnc': False,
                          'jesr_extra_br': True}
        self._opts = {'sfbdt_threshold': -99,
                      'run_tagger': False, 'tagger_versions': ['V02b', 'V02c', 'V02d'],
                      'run_mass_regression': False, 'mass_regression_versions': ['V01a', 'V01b', 'V01c'],
                      'WRITE_CACHE_FILE': False}
        for k in kwargs:
            if k in self._jmeSysts:
                self._jmeSysts[k] = kwargs[k]
            else:
                self._opts[k] = kwargs[k]
        self._needsJMECorr = any([self._jmeSysts['jec'], self._jmeSysts['jes'],
                                  self._jmeSysts['jer'], self._jmeSysts['jmr'],
                                  self._jmeSysts['met_unclustered'], self._jmeSysts['applyHEMUnc']])
        self._doJetCleaning = True

        logger.info('Running %s channel for %s jets with JME systematics %s, other options %s',
                    self._channel, self.jetType, str(self._jmeSysts), str(self._opts))

        if self.jetType == 'ak8':
            self._jetConeSize = 0.8
            self._fj_name = 'FatJet'
            self._sj_name = 'SubJet'
            self._fj_gen_name = 'GenJetAK8'
            self._sj_gen_name = 'SubGenJetAK8'
            self._sfbdt_files = [
                os.path.expandvars(
                    '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/sfBDT/ak8_ul/xgb_train_qcd.model.%d' % idx)
                for idx in range(10)]  # FIXME: update to AK8 training
            self._sfbdt_vars = ['fj_2_tau21', 'fj_2_sj1_rawmass', 'fj_2_sj2_rawmass',
                                'fj_2_ntracks_sv12', 'fj_2_sj1_sv1_pt', 'fj_2_sj2_sv1_pt']
        elif self.jetType == 'ak15':
            self._jetConeSize = 1.5
            self._fj_name = 'AK15Puppi'
            self._sj_name = 'AK15PuppiSubJet'
            self._fj_gen_name = 'GenJetAK15'
            self._sj_gen_name = 'GenSubJetAK15'
            self._sfbdt_files = [
                os.path.expandvars(
                    '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/sfBDT/ak8_ul/xgb_train_qcd.model.%d' % idx)
                for idx in range(10)]
            self._sfbdt_vars = ['fj_2_tau21', 'fj_2_sj1_rawmass', 'fj_2_sj2_rawmass',
                                'fj_2_ntracks_sv12', 'fj_2_sj1_sv1_pt', 'fj_2_sj2_sv1_pt']
        else:
            raise RuntimeError('Jet type %s is not recognized!' % self.jetType)

        self._fill_sv = self._channel in ('qcd', 'photon', 'higgs', 'inclusive') and self._opts['sfbdt_threshold'] > -99

        if self._needsJMECorr:
            if not self.is_run3:
                self.jetmetCorr = JetMETCorrector(year=self.year, jetType="AK4PFchs", **self._jmeSysts)
            else:
                self.jetmetCorr = JetMETCorrector(year=self.year, jetType="AK4PFPuppi", **self._jmeSysts)
            self.fatjetCorr = JetMETCorrector(year=self.year, jetType="AK8PFPuppi", **self._jmeSysts)
            self.subjetCorr = JetMETCorrector(year=self.year, jetType="AK4PFPuppi", **self._jmeSysts)

        if self._opts['run_tagger'] or self._opts['run_mass_regression']:
            from ..helpers.makeInputs import ParticleNetTagInfoMaker
            from ..helpers.runPrediction import ParticleNetJetTagsProducer
            self.tagInfoMaker = ParticleNetTagInfoMaker(
                fatjet_branch=self._fj_name, pfcand_branch='PFCands', sv_branch='SV', jetR=self._jetConeSize)
            prefix = os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data')
            if self._opts['run_tagger']:
                self.pnTaggers = [ParticleNetJetTagsProducer(
                    '%s/ParticleNet-MD/%s/{version}/particle-net.onnx' % (prefix, self.jetType),
                    '%s/ParticleNet-MD/%s/{version}/preprocess.json' % (prefix, self.jetType),
                    version=ver, cache_suffix='tagger') for ver in self._opts['tagger_versions']]
            if self._opts['run_mass_regression']:
                self.pnMassRegressions = [ParticleNetJetTagsProducer(
                    '%s/MassRegression/%s/{version}/particle_net_regression.onnx' % (prefix, self.jetType),
                    '%s/MassRegression/%s/{version}/preprocess.json' % (prefix, self.jetType),
                    version=ver, cache_suffix='mass') for ver in self._opts['mass_regression_versions']]

        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
        self.DeepJet_WP_L = {2015: 0.0508, 2016: 0.0480, 2017: 0.0532, 2018: 0.0490, 2022: 0.0490}[self.year]
        self.DeepJet_WP_M = {2015: 0.2598, 2016: 0.2489, 2017: 0.3040, 2018: 0.2783, 2022: 0.2783}[self.year]
        self.DeepJet_WP_T = {2015: 0.6502, 2016: 0.6377, 2017: 0.7476, 2018: 0.7100, 2022: 0.7100}[self.year]

    def beginJob(self):
        if self._needsJMECorr:
            # if not self.is_run3:
            self.jetmetCorr.beginJob()
            self.fatjetCorr.beginJob()
            self.subjetCorr.beginJob()
        if self._opts['sfbdt_threshold'] > -99:
            self.xgb = XGBEnsemble(self._sfbdt_files, self._sfbdt_vars)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.hasParticleNetProb = bool(inputTree.GetBranch(self._fj_name + '_ParticleNetMD_probXbb'))

        # remove all possible h5 cache files
        for f in os.listdir('.'):
            if f.endswith('.h5'):
                os.remove(f)

        if self._opts['run_tagger']:
            for p in self.pnTaggers:
                p.load_cache(inputFile)

        if self._opts['run_mass_regression']:
            for p in self.pnMassRegressions:
                p.load_cache(inputFile)

        if self._opts['run_tagger'] or self._opts['run_mass_regression']:
            self.tagInfoMaker.init_file(inputFile, fetch_step=1000)

        self.out = wrappedOutputTree

        # NOTE: branch names must start with a lower case letter
        # check keep_and_drop_output.txt
        self.out.branch("year", "I")
        self.out.branch("lumiwgt", "F")
        self.out.branch("jetR", "F")
        self.out.branch("passmetfilters", "O")
        self.out.branch("l1PreFiringWeight", "F")
        self.out.branch("l1PreFiringWeightUp", "F")
        self.out.branch("l1PreFiringWeightDown", "F")
        self.out.branch("nlep", "I")
        self.out.branch("ht", "F")
        self.out.branch("met", "F")
        self.out.branch("metphi", "F")
        if self.isMC and self._jmeSysts['jesr_extra_br']:
            # HT with JES/JER correction
            self.out.branch("ht_jesUncFactorUp", "F")
            self.out.branch("ht_jesUncFactorDn", "F")
            self.out.branch("ht_jerSmearFactorUp", "F")
            self.out.branch("ht_jerSmearFactorDn", "F")

        # Large-R jets
        for idx in ([1, 2] if self._channel in ['qcd', 'mutagged'] else [1]):
            prefix = 'fj_%d_' % idx

            # fatjet kinematics
            self.out.branch(prefix + "is_qualified", "O")
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "rawmass", "F")
            self.out.branch(prefix + "sdmass", "F")
            self.out.branch(prefix + "regressed_mass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "tau32", "F")
            if not self.is_run3:
                self.out.branch(prefix + "btagcsvv2", "F")
            self.out.branch(prefix + "btagjp", "F")

            # subjets
            self.out.branch(prefix + "deltaR_sj12", "F")
            self.out.branch(prefix + "sj1_pt", "F")
            self.out.branch(prefix + "sj1_eta", "F")
            self.out.branch(prefix + "sj1_phi", "F")
            self.out.branch(prefix + "sj1_rawmass", "F")
            self.out.branch(prefix + "sj1_btagdeepcsv", "F")
            self.out.branch(prefix + "sj2_pt", "F")
            self.out.branch(prefix + "sj2_eta", "F")
            self.out.branch(prefix + "sj2_phi", "F")
            self.out.branch(prefix + "sj2_rawmass", "F")
            self.out.branch(prefix + "sj2_btagdeepcsv", "F")

            # taggers
            self.out.branch(prefix + "DeepAK8_TvsQCD", "F")
            self.out.branch(prefix + "DeepAK8_WvsQCD", "F")
            self.out.branch(prefix + "DeepAK8_ZvsQCD", "F")
            self.out.branch(prefix + "DeepAK8_ZHbbvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_TvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_WvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZHbbvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZHccvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_bbVsLight", "F")
            self.out.branch(prefix + "DeepAK8MD_bbVsTop", "F")

            self.out.branch(prefix + "ParticleNet_TvsQCD", "F")
            self.out.branch(prefix + "ParticleNet_WvsQCD", "F")
            self.out.branch(prefix + "ParticleNet_ZvsQCD", "F")
            if not self.is_run3:
                self.out.branch(prefix + "ParticleNetMD_Xbb", "F")
                self.out.branch(prefix + "ParticleNetMD_Xcc", "F")
                self.out.branch(prefix + "ParticleNetMD_Xqq", "F")
                self.out.branch(prefix + "ParticleNetMD_QCD", "F")
                self.out.branch(prefix + "ParticleNetMD_XbbVsQCD", "F")
                self.out.branch(prefix + "ParticleNetMD_XccVsQCD", "F")
                self.out.branch(prefix + "ParticleNetMD_XccOrXqqVsQCD", "F")
            else:
                self.out.branch(prefix + "ParticleNetLegacy_Xbb", "F")
                self.out.branch(prefix + "ParticleNetLegacy_Xcc", "F")
                self.out.branch(prefix + "ParticleNetLegacy_Xqq", "F")
                self.out.branch(prefix + "ParticleNetLegacy_QCD", "F")
                self.out.branch(prefix + "ParticleNet_XbbVsQCD", "F")
                self.out.branch(prefix + "ParticleNet_XccVsQCD", "F")
                self.out.branch(prefix + "ParticleNet_XccOrXqqVsQCD", "F")
            # Additional tagger scores from NanoAODv9
            self.out.branch(prefix + "DeepAK8MD_HbbvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_H4qvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ccVsLight", "F")
            self.out.branch(prefix + "ParticleNet_HbbvsQCD", "F")
            self.out.branch(prefix + "ParticleNet_HccvsQCD", "F")
            self.out.branch(prefix + "ParticleNet_H4qvsQCD", "F")
            self.out.branch(prefix + "ParticleNet_mass", "F")
            self.out.branch(prefix + "btagDDBvLV2", "F")
            self.out.branch(prefix + "btagDDCvBV2", "F")
            self.out.branch(prefix + "btagDDCvLV2", "F")
            self.out.branch(prefix + "btagDeepB", "F")
            self.out.branch(prefix + "btagHbb", "F")

            if self._opts['run_tagger']:
                self.out.branch(prefix + "origParticleNetMD_XccVsQCD", "F")
                self.out.branch(prefix + "origParticleNetMD_XbbVsQCD", "F")

            # matching variables
            if self.isMC:
                self.out.branch(prefix + "nbhadrons", "I")
                self.out.branch(prefix + "nchadrons", "I")
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

                # factors of JES/JER correction
                if self._jmeSysts['jesr_extra_br']:
                    self.out.branch(prefix + "jesUncFactorUp", "F")
                    self.out.branch(prefix + "jesUncFactorDn", "F")
                    self.out.branch(prefix + "jerSmearFactorUp", "F")
                    self.out.branch(prefix + "jerSmearFactorDn", "F")

            if self._fill_sv:
                # SV variables
                self.out.branch(prefix + "nsv", "I")
                self.out.branch(prefix + "nsv_ptgt25", "I")
                self.out.branch(prefix + "nsv_ptgt50", "I")
                self.out.branch(prefix + "ntracks", "I")
                self.out.branch(prefix + "ntracks_sv12", "I")

                self.out.branch(prefix + "sj1_ntracks", "I")
                self.out.branch(prefix + "sj1_nsv", "I")
                self.out.branch(prefix + "sj1_sv1_pt", "F")
                self.out.branch(prefix + "sj1_sv1_mass", "F")
                self.out.branch(prefix + "sj1_sv1_masscor", "F")
                self.out.branch(prefix + "sj1_sv1_ntracks", "I")
                self.out.branch(prefix + "sj1_sv1_dxy", "F")
                self.out.branch(prefix + "sj1_sv1_dxysig", "F")
                self.out.branch(prefix + "sj1_sv1_dlen", "F")
                self.out.branch(prefix + "sj1_sv1_dlensig", "F")
                self.out.branch(prefix + "sj1_sv1_chi2ndof", "F")
                self.out.branch(prefix + "sj1_sv1_pangle", "F")

                self.out.branch(prefix + "sj2_ntracks", "I")
                self.out.branch(prefix + "sj2_nsv", "I")
                self.out.branch(prefix + "sj2_sv1_pt", "F")
                self.out.branch(prefix + "sj2_sv1_mass", "F")
                self.out.branch(prefix + "sj2_sv1_masscor", "F")
                self.out.branch(prefix + "sj2_sv1_ntracks", "I")
                self.out.branch(prefix + "sj2_sv1_dxy", "F")
                self.out.branch(prefix + "sj2_sv1_dxysig", "F")
                self.out.branch(prefix + "sj2_sv1_dlen", "F")
                self.out.branch(prefix + "sj2_sv1_dlensig", "F")
                self.out.branch(prefix + "sj2_sv1_chi2ndof", "F")
                self.out.branch(prefix + "sj2_sv1_pangle", "F")

                self.out.branch(prefix + "sj12_masscor_dxysig", "F")

                # sfBDT
                self.out.branch(prefix + "sfBDT", "F")

                # bb/cc gen hadrons
                if self.isMC and idx==(2 if self._channel == 'qcd' else 1) and self._channel != 'higgs':
                    for hadtype in ['b', 'c']:
                        for hadidx in [1, 2]:
                            self.out.branch(prefix + "gen{}hadron{}_pt".format(hadtype, hadidx), "F")
                            self.out.branch(prefix + "gen{}hadron{}_eta".format(hadtype, hadidx), "F")
                            self.out.branch(prefix + "gen{}hadron{}_phi".format(hadtype, hadidx), "F")
                            self.out.branch(prefix + "gen{}hadron{}_mass".format(hadtype, hadidx), "F")
                            self.out.branch(prefix + "gen{}hadron{}_pdgId".format(hadtype, hadidx), "I")

                # last parton list
                if self.isMC and self._channel != 'higgs':
                    for ptsuf in ['', '50']:
                        self.out.branch(prefix + "npart{}".format(ptsuf), "I")
                        self.out.branch(prefix + "nbpart{}".format(ptsuf), "I")
                        self.out.branch(prefix + "ncpart{}".format(ptsuf), "I")
                        self.out.branch(prefix + "ngpart{}".format(ptsuf), "I")
                        self.out.branch(prefix + "part{}_sumpt".format(ptsuf), "F")
                        self.out.branch(prefix + "bpart{}_sumpt".format(ptsuf), "F")
                        self.out.branch(prefix + "cpart{}_sumpt".format(ptsuf), "F")
                        self.out.branch(prefix + "gpart{}_sumpt".format(ptsuf), "F")
         

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self._opts['run_tagger'] and self._opts['WRITE_CACHE_FILE']:
            for p in self.pnTaggers:
                p.update_cache()

        if self._opts['run_mass_regression'] and self._opts['WRITE_CACHE_FILE']:
            for p in self.pnMassRegressions:
                p.update_cache()

        # remove all h5 cache files
        if self._opts['run_tagger'] or self._opts['run_mass_regression']:
            for f in os.listdir('.'):
                if f.endswith('.h5'):
                    os.remove(f)

    def selectLeptons(self, event):
        # do lepton selection
        event.looseLeptons = []  # used for jet lepton cleaning & lepton counting

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if not self.is_run3:
                electron_mvaIso_WP90 = el.mvaFall17V2noIso_WP90
            else:
                # Run 3
                electron_mvaIso_WP90 = el.mvaIso_WP90

            if el.pt > 10 and abs(el.eta) < 2.5 and abs(el.dxy) < 0.05 and abs(el.dz) < 0.2 \
                and electron_mvaIso_WP90 and el.miniPFRelIso_all < 0.4:
                    event.looseLeptons.append(el)

        muons = Collection(event, "Muon")
        for mu in muons:
            if mu.pt > 10 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.05 and abs(mu.dz) < 0.2 \
                    and mu.looseId and mu.miniPFRelIso_all < 0.4:
                event.looseLeptons.append(mu)

        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)

    def correctJetsAndMET(self, event):
        # correct Jets and MET
        event.idx = event._entry if event._tree._entrylist is None else event._tree._entrylist.GetEntry(event._entry)
        event._allJets = Collection(event, "Jet")
        event.met = METObject(event, "MET")
        event._allFatJets = Collection(event, self._fj_name)
        event.subjets = Collection(event, self._sj_name)  # do not sort subjets after updating!!

        # ## do some hack here... use uncorrected jet pT!
        # for idx, j in enumerate(event._allFatJets):
        #     j.rawP4 = polarP4(j) * (1. - j.rawFactor)
        #     j.pt = j.rawP4.pt()
        #     j.mass = j.rawP4.mass()
        if self._needsJMECorr:
            if not self.is_run3:
                rho = event.fixedGridRhoFastjetAll
            else:
                # Run 3
                # print(f"event: {event._tree}")
                # branches = event._tree.GetListOfBranches()
                # branch_names = [branch.GetName() for branch in branches]
                # print(f"event._tree Branch names: {branch_names}")
                rho = event.Rho_fixedGridRhoFastjetAll
            # correct AK4 jets and MET
            self.jetmetCorr.setSeed(rndSeed(event, event._allJets))
            self.jetmetCorr.correctJetAndMET(jets=event._allJets, lowPtJets=Collection(event, "CorrT1METJet"),
                                            met=event.met, rawMET=METObject(event, "RawMET"),
                                            defaultMET=METObject(event, "MET"),
                                            rho=rho, genjets=Collection(event, 'GenJet') if self.isMC else None,
                                            isMC=self.isMC, runNumber=event.run)
            event._allJets = sorted(event._allJets, key=lambda x: x.pt, reverse=True)  # sort by pt after updating

            # correct fatjets
            self.fatjetCorr.setSeed(rndSeed(event, event._allFatJets))
            self.fatjetCorr.correctJetAndMET(jets=event._allFatJets, met=None, rho=rho,
                                             genjets=Collection(event, self._fj_gen_name) if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            # correct subjets
            self.subjetCorr.setSeed(rndSeed(event, event.subjets))
            self.subjetCorr.correctJetAndMET(jets=event.subjets, met=None, rho=rho,
                                             genjets=Collection(event, self._sj_gen_name) if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)

        # jet mass resolution smearing
        if self.isMC and self._jmeSysts['jmr']:
            raise NotImplementedError

        # link fatjet to subjets and recompute softdrop mass
        for idx, fj in enumerate(event._allFatJets):
            fj.idx = idx
            fj.is_qualified = True
            fj.subjets = get_subjets(fj, event.subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = sumP4(*fj.subjets).M()
        event._allFatJets = sorted(event._allFatJets, key=lambda x: x.pt, reverse=True)  # sort by pt

        # select lepton-cleaned jets
        if self._doJetCleaning:
            event.fatjets = [fj for fj in event._allFatJets if fj.pt > 200 and abs(fj.eta) < 2.4 and (
                fj.jetId & 2) and closest(fj, event.looseLeptons)[1] >= self._jetConeSize]
            event.ak4jets = [j for j in event._allJets if j.pt > 25 and abs(j.eta) < 2.4 and (
                j.jetId & 4) and closest(j, event.looseLeptons)[1] >= 0.4]
        else:
            event.fatjets = [fj for fj in event._allFatJets if fj.pt > 200 and abs(fj.eta) < 2.4 and (
                fj.jetId & 2)]
            event.ak4jets = [j for j in event._allJets if j.pt > 25 and abs(j.eta) < 2.4 and (
                j.jetId & 4)]
        event.ht = sum([j.pt for j in event.ak4jets])
        if self.isMC and self._jmeSysts['jesr_extra_br']:
            event.ht_jesUncFactorUp = sum([j.pt * j.jesUncFactorUp for j in event.ak4jets])
            event.ht_jesUncFactorDn = sum([j.pt * j.jesUncFactorDn for j in event.ak4jets])
            event.ht_jerSmearFactorUp = sum([j.pt * j.jerSmearFactorUp for j in event.ak4jets])
            event.ht_jerSmearFactorDn = sum([j.pt * j.jerSmearFactorDn for j in event.ak4jets])

    def selectSV(self, event):
        event._allSV = Collection(event, "SV")
        event.secondary_vertices = []
        for sv in event._allSV:
            # if sv.ntracks > 2 and abs(sv.dxy) < 3. and sv.dlenSig > 4:
            # if sv.dlenSig > 4:
            event.secondary_vertices.append(sv)
        event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x: x.pt, reverse=True)  # sort by pt
        # event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x : x.dxySig, reverse=True)  # sort by dxysig

    def matchSVToFatJets(self, event, fatjets):
        # match SV to fatjets
        for fj in fatjets:
            fj.sv_list = []
            for sv in event.secondary_vertices:
                if deltaR(sv, fj) < self._jetConeSize:
                    fj.sv_list.append(sv)
            # match SV to subjets
            drcut = min(0.4, 0.5 * deltaR(*fj.subjets)) if len(fj.subjets) == 2 else 0.4
            for sj in fj.subjets:
                sj.sv_list = []
                for sv in event.secondary_vertices:
                    if deltaR(sv, sj) < drcut:
                        sj.sv_list.append(sv)

            fj.nsv_ptgt25 = 0
            fj.nsv_ptgt50 = 0
            fj.ntracks = 0
            fj.ntracks_sv12 = 0
            for isv, sv in enumerate(fj.sv_list):
                fj.ntracks += sv.ntracks
                if isv < 2:
                    fj.ntracks_sv12 += sv.ntracks
                if sv.pt > 25:
                    fj.nsv_ptgt25 += 1
                if sv.pt > 50:
                    fj.nsv_ptgt50 += 1

            # sfBDT & sj12_masscor_dxysig
            fj.sfBDT = -1
            fj.sj12_masscor_dxysig = 0
            if len(fj.subjets) == 2:
                sj1, sj2 = fj.subjets
                if len(sj1.sv_list) > 0 and len(sj2.sv_list) > 0:
                    sj1_sv, sj2_sv = sj1.sv_list[0], sj2.sv_list[0]
                    sfbdt_inputs = {
                        'fj_2_tau21': fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99,
                        'fj_2_sj1_rawmass': sj1.mass,
                        'fj_2_sj2_rawmass': sj2.mass,
                        'fj_2_ntracks_sv12': fj.ntracks_sv12,
                        'fj_2_sj1_sv1_pt': sj1_sv.pt,
                        'fj_2_sj2_sv1_pt': sj2_sv.pt,
                    }
                    if hasattr(self, 'xgb'):
                        fj.sfBDT = self.xgb.eval(sfbdt_inputs, model_idx=(event.event % 10))
                    fj.sj12_masscor_dxysig = corrected_svmass(sj1_sv if sj1_sv.dxySig > sj2_sv.dxySig else sj2_sv)

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

        if self._fill_sv and self._channel != 'higgs':
            # bb/cc matching
            # FIXME: only available for qcd & ggh(cc/bb) sample
            probe_fj = event.fatjets[1 if self._channel == 'qcd' else 0]
            probe_fj.genBhadron, probe_fj.genChadron = [], []
            for gp in genparts:
                if gp.pdgId in [5, -5] and gp.genPartIdxMother>=0 and genparts[gp.genPartIdxMother].pdgId in [21, 25] and deltaR(gp, probe_fj)<=self._jetConeSize:
                    if len(probe_fj.genBhadron)==0 or (len(probe_fj.genBhadron)>0 and gp.genPartIdxMother==probe_fj.genBhadron[0].genPartIdxMother):
                        probe_fj.genBhadron.append(gp)
                if gp.pdgId in [4, -4] and gp.genPartIdxMother>=0 and genparts[gp.genPartIdxMother].pdgId in [21, 25] and deltaR(gp, probe_fj)<=self._jetConeSize:
                    if len(probe_fj.genChadron)==0 or (len(probe_fj.genChadron)>0 and gp.genPartIdxMother==probe_fj.genChadron[0].genPartIdxMother):
                        probe_fj.genChadron.append(gp)
            probe_fj.genBhadron.sort(key=lambda x: x.pt, reverse=True)
            probe_fj.genChadron.sort(key=lambda x: x.pt, reverse=True)
            # null padding
            probe_fj.genBhadron += [_NullObject() for _ in range(2-len(probe_fj.genBhadron))]
            probe_fj.genChadron += [_NullObject() for _ in range(2-len(probe_fj.genChadron))]

            # last parton information
            for ifj in range(2 if self._channel == 'qcd' else 1):
                fj = event.fatjets[ifj]
                fj.npart, fj.nbpart, fj.ncpart, fj.ngpart, fj.part_sumpt, fj.bpart_sumpt, fj.cpart_sumpt, fj.gpart_sumpt = 0, 0, 0, 0, 0, 0, 0, 0
                fj.npart50, fj.nbpart50, fj.ncpart50, fj.ngpart50, fj.part50_sumpt, fj.bpart50_sumpt, fj.cpart50_sumpt, fj.gpart50_sumpt = 0, 0, 0, 0, 0, 0, 0, 0
                for gp in genparts:
                    if gp.status>70 and gp.status<80 and (gp.statusFlags & (1 << 13)) and abs(gp.pdgId) in [1,2,3,4,5,6,21] and gp.pt>=5 and deltaR(gp, fj)<=self._jetConeSize:
                        fj.npart += 1; fj.part_sumpt += gp.pt
                        if gp.pdgId in [5, -5]:
                            fj.nbpart += 1; fj.bpart_sumpt += gp.pt
                        elif gp.pdgId in [4, -4]:
                            fj.ncpart += 1; fj.cpart_sumpt += gp.pt
                        elif gp.pdgId == 21:
                            fj.ngpart += 1; fj.gpart_sumpt += gp.pt
                        if gp.pt>=50:
                            fj.npart50 += 1; fj.part50_sumpt += gp.pt
                            if gp.pdgId in [5, -5]:
                                fj.nbpart50 += 1; fj.bpart50_sumpt += gp.pt
                            elif gp.pdgId in [4, -4]:
                                fj.ncpart50 += 1; fj.cpart50_sumpt += gp.pt
                            elif gp.pdgId == 21:
                                fj.ngpart50 += 1; fj.gpart50_sumpt += gp.pt


    def evalTagger(self, event, jets):
        for j in jets:
            if self._opts['run_tagger']:
                outputs = [p.predict_with_cache(self.tagInfoMaker, event.idx, j.idx, j) for p in self.pnTaggers]
                outputs = ensemble(outputs, np.mean)
                try:
                    j.pn_Xbb = outputs['probXbb']
                    j.pn_Xcc = outputs['probXcc']
                    j.pn_Xqq = outputs['probXqq']
                    j.pn_QCD = convert_prob(outputs, None, prefix='prob')
                except RuntimeError:
                    # TODO: may need to update this in Run 3
                    raise NotImplementedError("Need to update for Run 3")
            else:
                if self.hasParticleNetProb:
                    try:
                        j.pn_Xbb = j.ParticleNetMD_probXbb
                        j.pn_Xcc = j.ParticleNetMD_probXcc
                        j.pn_Xqq = j.ParticleNetMD_probXqq
                        j.pn_QCD = convert_prob(j, None, prefix='ParticleNetMD_prob')
                    except RuntimeError:
                        # TODO: may need to update this in Run 3
                        raise NotImplementedError("Need to update for Run 3")
                else:
                    if not self.is_run3:
                        j.pn_Xbb = j.particleNetMD_Xbb
                        j.pn_Xcc = j.particleNetMD_Xcc
                        j.pn_Xqq = j.particleNetMD_Xqq
                        j.pn_QCD = j.particleNetMD_QCD
                    else:
                        # Run 3
                        j.pn_legacy_Xbb = j.particleNetLegacy_Xbb
                        j.pn_legacy_Xcc = j.particleNetLegacy_Xcc
                        j.pn_legacy_Xqq = j.particleNetLegacy_Xqq
                        j.pn_legacy_QCD = j.particleNetLegacy_QCD
                        
            if not self.is_run3:            
                j.pn_XbbVsQCD = convert_prob(j, ['Xbb'], ['QCD'], prefix='pn_')
                j.pn_XccVsQCD = convert_prob(j, ['Xcc'], ['QCD'], prefix='pn_')
                j.pn_XccOrXqqVsQCD = convert_prob(j, ['Xcc', 'Xqq'], ['QCD'], prefix='pn_')
            else:
                # Run 3
                j.pn_XbbVsQCD = j.particleNet_XbbVsQCD
                j.pn_XccVsQCD = j.particleNet_XccVsQCD
                pn_XqqVsQCD = j.particleNet_XqqVsQCD
                # j.pn_QCD = j.particleNet_QCD
                j.pn_XccOrXqqVsQCD = j.pn_XccVsQCD + pn_XqqVsQCD

    def evalMassRegression(self, event, jets):
        for j in jets:
            if self._opts['run_mass_regression']:
                outputs = [p.predict_with_cache(self.tagInfoMaker, event.idx, j.idx, j) for p in self.pnMassRegressions]
                j.regressed_mass = ensemble(outputs, np.median)['mass']
            else:
                try:
                    j.regressed_mass = j.particleNet_mass
                except RuntimeError:
                    j.regressed_mass = 0

    def fillBaseEventInfo(self, event):
        self.out.fillBranch("jetR", self._jetConeSize)
        self.out.fillBranch("year", self.year)
        self.out.fillBranch("lumiwgt", lumi_dict[self.year_str])

        met_filters = bool(
            event.Flag_goodVertices and
            event.Flag_globalSuperTightHalo2016Filter and
            event.Flag_HBHENoiseFilter and
            event.Flag_HBHENoiseIsoFilter and
            event.Flag_EcalDeadCellTriggerPrimitiveFilter and
            event.Flag_BadPFMuonFilter and
            event.Flag_BadPFMuonDzFilter and
            event.Flag_eeBadScFilter
        )
        if self.year in (2017, 2018):
            met_filters = met_filters and event.Flag_ecalBadCalibFilter
        self.out.fillBranch("passmetfilters", met_filters)

        # L1 prefire weights
        if self.year <= 2017:
            self.out.fillBranch("l1PreFiringWeight", event.L1PreFiringWeight_Nom)
            self.out.fillBranch("l1PreFiringWeightUp", event.L1PreFiringWeight_Up)
            self.out.fillBranch("l1PreFiringWeightDown", event.L1PreFiringWeight_Dn)
        else:
            self.out.fillBranch("l1PreFiringWeight", 1.0)
            self.out.fillBranch("l1PreFiringWeightUp", 1.0)
            self.out.fillBranch("l1PreFiringWeightDown", 1.0)

        self.out.fillBranch("nlep", len(event.looseLeptons))
        self.out.fillBranch("ht", event.ht)
        if self.isMC and self._jmeSysts['jesr_extra_br']:
            self.out.fillBranch("ht_jesUncFactorUp", event.ht_jesUncFactorUp)
            self.out.fillBranch("ht_jesUncFactorDn", event.ht_jesUncFactorDn)
            self.out.fillBranch("ht_jerSmearFactorUp", event.ht_jerSmearFactorUp)
            self.out.fillBranch("ht_jerSmearFactorDn", event.ht_jerSmearFactorDn)
        self.out.fillBranch("met", event.met.pt)
        self.out.fillBranch("metphi", event.met.phi)

    def _get_filler(self, obj):

        def filler(branch, value, default=0):
            self.out.fillBranch(branch, value if obj else default)

        return filler

    def fillFatJetInfo(self, event, fatjets):
        for idx in ([1, 2] if self._channel in ['qcd', 'mutagged'] else [1]):
            prefix = 'fj_%d_' % idx

            if len(fatjets) <= idx - 1 or not fatjets[idx - 1].is_qualified:
                # fill zeros if fatjet fails probe selection
                for b in self.out._branches.keys():
                    if b.startswith(prefix):
                        self.out.fillBranch(b, 0)
                continue

            fj = fatjets[idx - 1]

            # fatjet kinematics
            self.out.fillBranch(prefix + "is_qualified", fj.is_qualified)
            self.out.fillBranch(prefix + "pt", fj.pt)
            self.out.fillBranch(prefix + "eta", fj.eta)
            self.out.fillBranch(prefix + "phi", fj.phi)
            self.out.fillBranch(prefix + "rawmass", fj.mass)
            self.out.fillBranch(prefix + "sdmass", fj.msoftdrop)
            self.out.fillBranch(prefix + "regressed_mass", fj.regressed_mass)
            self.out.fillBranch(prefix + "tau21", fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99)
            self.out.fillBranch(prefix + "tau32", fj.tau3 / fj.tau2 if fj.tau2 > 0 else 99)
            if not self.is_run3:
                self.out.fillBranch(prefix + "btagcsvv2", fj.btagCSVV2)
            else:
                # Run 3: not used
                pass
            try:
                self.out.fillBranch(prefix + "btagjp", fj.btagJP)
            except RuntimeError:
                self.out.fillBranch(prefix + "btagjp", -1)

            # subjets
            self.out.fillBranch(prefix + "deltaR_sj12", deltaR(*fj.subjets) if len(fj.subjets) == 2 else 99)
            for idx_sj, sj in enumerate(fj.subjets):
                prefix_sj = prefix + 'sj%d_' % (idx_sj + 1)
                self.out.fillBranch(prefix_sj + "pt", sj.pt)
                self.out.fillBranch(prefix_sj + "eta", sj.eta)
                self.out.fillBranch(prefix_sj + "phi", sj.phi)
                self.out.fillBranch(prefix_sj + "rawmass", sj.mass)
                try:
                    self.out.fillBranch(prefix_sj + "btagdeepcsv", sj.btagDeepB)
                except RuntimeError:
                    self.out.fillBranch(prefix_sj + "btagdeepcsv", -1)

            # taggers
            try:
                # Full
                self.out.fillBranch(prefix + "DeepAK8_TvsQCD", fj.deepTag_TvsQCD)
                self.out.fillBranch(prefix + "DeepAK8_WvsQCD", fj.deepTag_WvsQCD)
                self.out.fillBranch(prefix + "DeepAK8_ZvsQCD", fj.deepTag_ZvsQCD)
                # MD
                self.out.fillBranch(prefix + "DeepAK8MD_TvsQCD", fj.deepTagMD_TvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_WvsQCD", fj.deepTagMD_WvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ZvsQCD", fj.deepTagMD_ZvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHbbvsQCD", fj.deepTagMD_ZHbbvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHccvsQCD", fj.deepTagMD_ZHccvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsLight", fj.deepTagMD_bbvsLight)
                try:
                    bbVsTop = (1 / (1 + (fj.deepTagMD_TvsQCD / fj.deepTagMD_HbbvsQCD) * (1 - fj.deepTagMD_HbbvsQCD) / (1 - fj.deepTagMD_TvsQCD)))  # noqa
                except ZeroDivisionError:
                    bbVsTop = 0
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsTop", bbVsTop)
            except RuntimeError:
                # if no DeepAK8 branches
                self.out.fillBranch(prefix + "DeepAK8_TvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8_WvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8_ZvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_TvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_WvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ZvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHbbvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHccvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsLight", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsTop", -1)

            try:
                self.out.fillBranch(prefix + "DeepAK8_ZHbbvsQCD",
                                    convert_prob(fj, ['Zbb', 'Hbb'], prefix='deepTag_prob'))
            except RuntimeError:
                # if no DeepAK8 raw probs
                self.out.fillBranch(prefix + "DeepAK8_ZHbbvsQCD", -1)

            # ParticleNet
            if self.hasParticleNetProb:
                self.out.fillBranch(prefix + "ParticleNet_TvsQCD",
                                    convert_prob(fj, ['Tbcq', 'Tbqq'], prefix='ParticleNet_prob'))
                self.out.fillBranch(prefix + "ParticleNet_WvsQCD",
                                    convert_prob(fj, ['Wcq', 'Wqq'], prefix='ParticleNet_prob'))
                self.out.fillBranch(prefix + "ParticleNet_ZvsQCD",
                                    convert_prob(fj, ['Zbb', 'Zcc', 'Zqq'], prefix='ParticleNet_prob'))
            else:
                try:
                    # nominal ParticleNet from official NanoAOD
                    self.out.fillBranch(prefix + "ParticleNet_TvsQCD", fj.particleNet_TvsQCD)
                    self.out.fillBranch(prefix + "ParticleNet_WvsQCD", fj.particleNet_WvsQCD)
                    self.out.fillBranch(prefix + "ParticleNet_ZvsQCD", fj.particleNet_ZvsQCD)
                except RuntimeError:
                    # if no nominal ParticleNet
                    self.out.fillBranch(prefix + "ParticleNet_TvsQCD", -1)
                    self.out.fillBranch(prefix + "ParticleNet_WvsQCD", -1)
                    self.out.fillBranch(prefix + "ParticleNet_ZvsQCD", -1)

            # ParticleNet-MD
            if not self.is_run3:
                self.out.fillBranch(prefix + "ParticleNetMD_Xbb", fj.pn_Xbb)
                self.out.fillBranch(prefix + "ParticleNetMD_Xcc", fj.pn_Xcc)
                self.out.fillBranch(prefix + "ParticleNetMD_Xqq", fj.pn_Xqq)
                self.out.fillBranch(prefix + "ParticleNetMD_QCD", fj.pn_QCD)
                self.out.fillBranch(prefix + "ParticleNetMD_XbbVsQCD", fj.pn_XbbVsQCD)
                self.out.fillBranch(prefix + "ParticleNetMD_XccVsQCD", fj.pn_XccVsQCD)
                self.out.fillBranch(prefix + "ParticleNetMD_XccOrXqqVsQCD", fj.pn_XccOrXqqVsQCD)
            else:
                self.out.fillBranch(prefix + "ParticleNetLegacy_Xbb", fj.pn_legacy_Xbb)
                self.out.fillBranch(prefix + "ParticleNetLegacy_Xcc", fj.pn_legacy_Xcc)
                self.out.fillBranch(prefix + "ParticleNetLegacy_Xqq", fj.pn_legacy_Xqq)
                self.out.fillBranch(prefix + "ParticleNetLegacy_QCD", fj.pn_legacy_QCD)
                self.out.fillBranch(prefix + "ParticleNet_XbbVsQCD", fj.pn_XbbVsQCD)
                self.out.fillBranch(prefix + "ParticleNet_XccVsQCD", fj.pn_XccVsQCD)
                self.out.fillBranch(prefix + "ParticleNet_XccOrXqqVsQCD", fj.pn_XccOrXqqVsQCD)

            if self._opts['run_tagger']:
                self.out.fillBranch(prefix + "origParticleNetMD_XccVsQCD",
                                    convert_prob(fj, ['Xcc'], None, prefix='ParticleNetMD_prob'))
                self.out.fillBranch(prefix + "origParticleNetMD_XbbVsQCD",
                                    convert_prob(fj, ['Xbb'], None, prefix='ParticleNetMD_prob'))

            # Additional tagger scores from NanoAODv9
            try:
                self.out.fillBranch(prefix + "DeepAK8MD_HbbvsQCD", fj.deepTagMD_HbbvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_H4qvsQCD", fj.deepTagMD_H4qvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ccVsLight", fj.deepTagMD_ccvsLight)
            except RuntimeError:
                self.out.fillBranch(prefix + "DeepAK8MD_HbbvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_H4qvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ccVsLight", -1)
            try:
                self.out.fillBranch(prefix + "ParticleNet_HbbvsQCD", fj.particleNet_HbbvsQCD)
                self.out.fillBranch(prefix + "ParticleNet_HccvsQCD", fj.particleNet_HccvsQCD)
                self.out.fillBranch(prefix + "ParticleNet_H4qvsQCD", fj.particleNet_H4qvsQCD)
            except RuntimeError:
                self.out.fillBranch(prefix + "ParticleNet_HbbvsQCD", -1)
                self.out.fillBranch(prefix + "ParticleNet_HccvsQCD", -1)
                self.out.fillBranch(prefix + "ParticleNet_H4qvsQCD", -1)
            try:
                self.out.fillBranch(prefix + "ParticleNet_mass", fj.particleNet_mass)
            except RuntimeError:
                self.out.fillBranch(prefix + "ParticleNet_mass", -1)
            try:
                self.out.fillBranch(prefix + "btagDDBvLV2", fj.btagDDBvLV2)
                self.out.fillBranch(prefix + "btagDDCvBV2", fj.btagDDCvBV2)
                self.out.fillBranch(prefix + "btagDDCvLV2", fj.btagDDCvLV2)
                self.out.fillBranch(prefix + "btagDeepB", fj.btagDeepB)
                self.out.fillBranch(prefix + "btagHbb", fj.btagHbb)
            except RuntimeError:
                self.out.fillBranch(prefix + "btagDDBvLV2", -1)
                self.out.fillBranch(prefix + "btagDDCvBV2", -1)
                self.out.fillBranch(prefix + "btagDDCvLV2", -1)
                self.out.fillBranch(prefix + "btagDeepB", -1)
                self.out.fillBranch(prefix + "btagHbb", -1)

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

                self.out.fillBranch(prefix + "nbhadrons", fj.nBHadrons)
                self.out.fillBranch(prefix + "nchadrons", fj.nCHadrons)
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

                if self._jmeSysts['jesr_extra_br']:
                    self.out.fillBranch(prefix + "jesUncFactorUp", fj.jesUncFactorUp)
                    self.out.fillBranch(prefix + "jesUncFactorDn", fj.jesUncFactorDn)
                    self.out.fillBranch(prefix + "jerSmearFactorUp", fj.jerSmearFactorUp)
                    self.out.fillBranch(prefix + "jerSmearFactorDn", fj.jerSmearFactorDn)

            if self._fill_sv:
                # SV variables
                self.out.fillBranch(prefix + "nsv", len(fj.sv_list))
                self.out.fillBranch(prefix + "nsv_ptgt25", fj.nsv_ptgt25)
                self.out.fillBranch(prefix + "nsv_ptgt50", fj.nsv_ptgt50)
                self.out.fillBranch(prefix + "ntracks", fj.ntracks)
                self.out.fillBranch(prefix + "ntracks_sv12", fj.ntracks_sv12)

                for idx_sj in (0, 1):
                    prefix_sj = prefix + 'sj%d_' % (idx_sj + 1)
                    try:
                        sj = fj.subjets[idx_sj]
                    except IndexError:
                        # fill zeros if not enough subjets
                        for b in self.out._branches.keys():
                            if b.startswith(prefix_sj):
                                self.out.fillBranch(b, 0)
                        continue

                    self.out.fillBranch(prefix_sj + "ntracks", sum([sv.ntracks for sv in sj.sv_list]))
                    self.out.fillBranch(prefix_sj + "nsv", len(sj.sv_list))
                    sv = sj.sv_list[0] if len(sj.sv_list) else _NullObject()
                    fill_sv = self._get_filler(sv)  # wrapper, fill default value if sv=None
                    fill_sv(prefix_sj + "sv1_pt", sv.pt)
                    fill_sv(prefix_sj + "sv1_mass", sv.mass)
                    fill_sv(prefix_sj + "sv1_masscor", corrected_svmass(sv) if sv else 0)
                    fill_sv(prefix_sj + "sv1_ntracks", sv.ntracks)
                    fill_sv(prefix_sj + "sv1_dxy", sv.dxy)
                    fill_sv(prefix_sj + "sv1_dxysig", sv.dxySig)
                    fill_sv(prefix_sj + "sv1_dlen", sv.dlen)
                    fill_sv(prefix_sj + "sv1_dlensig", sv.dlenSig)
                    fill_sv(prefix_sj + "sv1_chi2ndof", sv.chi2)
                    fill_sv(prefix_sj + "sv1_pangle", sv.pAngle)
                self.out.fillBranch(prefix + "sj12_masscor_dxysig", fj.sj12_masscor_dxysig)

                # sfBDT
                self.out.fillBranch(prefix + "sfBDT", fj.sfBDT)

                if self.isMC and idx==(2 if self._channel == 'qcd' else 1) and self._channel != 'higgs':
                    for hadtype in ['b', 'c']:
                        for hadidx in [1, 2]:
                            gp = fj.genBhadron[hadidx - 1] if hadtype=='b' else fj.genChadron[hadidx - 1]
                            fill_gp = self._get_filler(gp)  # wrapper, fill default value if sv=None
                            fill_gp(prefix + "gen{}hadron{}_pt".format(hadtype, hadidx), gp.pt)
                            fill_gp(prefix + "gen{}hadron{}_eta".format(hadtype, hadidx), gp.eta)
                            fill_gp(prefix + "gen{}hadron{}_phi".format(hadtype, hadidx), gp.phi)
                            fill_gp(prefix + "gen{}hadron{}_mass".format(hadtype, hadidx), gp.mass)
                            fill_gp(prefix + "gen{}hadron{}_pdgId".format(hadtype, hadidx), gp.pdgId)

                if self.isMC and self._channel != 'higgs':
                    self.out.fillBranch(prefix + "npart", fj.npart)
                    self.out.fillBranch(prefix + "nbpart", fj.nbpart)
                    self.out.fillBranch(prefix + "ncpart", fj.ncpart)
                    self.out.fillBranch(prefix + "ngpart", fj.ngpart)
                    self.out.fillBranch(prefix + "part_sumpt", fj.part_sumpt)
                    self.out.fillBranch(prefix + "bpart_sumpt", fj.bpart_sumpt)
                    self.out.fillBranch(prefix + "cpart_sumpt", fj.cpart_sumpt)
                    self.out.fillBranch(prefix + "gpart_sumpt", fj.gpart_sumpt)
                    self.out.fillBranch(prefix + "npart50", fj.npart50)
                    self.out.fillBranch(prefix + "nbpart50", fj.nbpart50)
                    self.out.fillBranch(prefix + "ncpart50", fj.ncpart50)
                    self.out.fillBranch(prefix + "ngpart50", fj.ngpart50)
                    self.out.fillBranch(prefix + "part50_sumpt", fj.part50_sumpt)
                    self.out.fillBranch(prefix + "bpart50_sumpt", fj.bpart50_sumpt)
                    self.out.fillBranch(prefix + "cpart50_sumpt", fj.cpart50_sumpt)
                    self.out.fillBranch(prefix + "gpart50_sumpt", fj.gpart50_sumpt)
