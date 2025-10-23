from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from .HeavyFlavBaseProducer import HeavyFlavBaseProducer
from ..helpers.utils import deltaR, polarP4, configLogger

import logging
logger = logging.getLogger('qcd')
configLogger('qcd', loglevel=logging.INFO)

class QCDSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        self._require_sv_cut = kwargs.pop('require_sv_cut', True)
        if self._require_sv_cut == False:
            logger.info('require_sv_cut is disabled! SV selection will not be applied.')

        self._run_gen_hadron_nsubs = kwargs.pop('run_gen_hadron_nsubs', False)
        if self._run_gen_hadron_nsubs == True:
            logger.info('run_gen_hadron_nsubs is enabled! Gen hadron N-subjettiness variables will be filled. ' + \
                'MAKE SURE that the input NanoAOD file contains full GenParticles!')

        super(QCDSampleProducer, self).__init__(channel='qcd', **kwargs)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(QCDSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        # trigger variables
        self.out.branch("passHTTrig", "O")

        # gen hadron N-subjettiness variables
        if self.isMC and self._run_gen_hadron_nsubs:
            for idx in [1, 2]:
                prefix = 'fj_%d_' % idx
                self.out.branch(prefix + "gen_hadron_tau1", "F")
                self.out.branch(prefix + "gen_hadron_tau2", "F")
                self.out.branch(prefix + "gen_hadron_tau3", "F")
                self.out.branch(prefix + "gen_hadron_tau4", "F")
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.selectLeptons(event)
        self.correctJetsAndMET(event)

        if len(event.fatjets) < 2:
            return False
        probe_jets = event.fatjets[:2]

        self.selectSV(event)
        if self._require_sv_cut and len(event.secondary_vertices) < 2:
            return False

        self.matchSVToFatJets(event, probe_jets)

        # check if any of two leading fatjets is qualified
        for fj in probe_jets:
            fj.is_qualified = False
            if len(fj.subjets) == 2 and fj.msoftdrop > 50 and fj.msoftdrop < 200:
                if not self._opts['run_sfbdt']:
                    fj.is_qualified = True
                else:
                    if fj.sfBDT > self._opts['sfbdt_threshold']:
                        fj.is_qualified = True

        if all([fj.is_qualified == False for fj in probe_jets]):
            # no fatjet is qualified, reject event
            return False

        self.loadGenHistory(event, probe_jets)
        if self._run_gen_hadron_nsubs:
            self.calculateGenHadronNSubjettiness(event, probe_jets)

        self.evalTagger(event, probe_jets)
        self.evalMassRegression(event, probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)
        if self._run_gen_hadron_nsubs:
            self.fillGenHadronNSubjettiness(event, probe_jets)

        if self.year in ["2016APV", "2016"]:
            self.out.fillBranch("passHTTrig", event.HLT_PFHT900)
        else:
            self.out.fillBranch("passHTTrig", event.HLT_PFHT1050)

        return True


    def calculateGenHadronNSubjettiness(self, event, fatjets):
        """Cluster first-generation GEN hadrons and compute N-subjettiness variables.
        Note: should only run for special MC samples with full GenParticles!
        """
        if not self.isMC:
            logger.warning('calculateGenHadronNSubjettiness is only supported for MC samples')
            return

        import fastjet
        import numpy as np
        
        for fj in fatjets:
            # collect all initial hadrons matched to the fatjet
            genhadrons = [gp for gp in Collection(event, "GenPart") 
                         if abs(gp.pdgId) > 100 and gp.genPartIdxMother >= 0 
                         and abs(event.genparts[gp.genPartIdxMother].pdgId) in [1, 2, 3, 4, 5, 6]
                         and deltaR(gp, fj) <= self._jetConeSize]
            
            if len(genhadrons) == 0:
                for nsub in [1, 2, 3, 4]:
                    setattr(fj, f'gen_hadron_tau{nsub}', -1)
                continue
            
            # prepare particle array for fastjet
            particles = []
            for gp in genhadrons:
                gp_p4 = polarP4(gp)
                p = fastjet.PseudoJet(gp_p4.px(), gp_p4.py(), gp_p4.pz(), gp_p4.energy())
                particles.append(p)
            
            # cluster particles using kt algorithm
            jet_def = fastjet.JetDefinition(fastjet.kt_algorithm, self._jetConeSize)
            cs = fastjet.ClusterSequence(particles, jet_def)

            # calculate tau1 to tau4
            for nsub in [1, 2, 3, 4]:
                if len(genhadrons) < nsub:
                    setattr(fj, f'gen_hadron_tau{nsub}', -1)
                else:
                    # get N exclusive subjets
                    subjets = cs.exclusive_jets(nsub)
                    
                    # calculate tau_N = (1/R_0) * Σ_k p_T,k * min(ΔR_k,j1, ΔR_k,j2, ..., ΔR_k,jN)
                    pt_sum = 0
                    tau_sum = 0
                    for gp in genhadrons:
                        # find minimum ΔR to any subjet
                        min_dr = min([deltaR(gp.eta, gp.phi, j.eta(), j.phi()) for j in subjets])
                        pt_sum += gp.pt
                        tau_sum += gp.pt * min_dr
                    setattr(fj, f'gen_hadron_tau{nsub}', tau_sum / (pt_sum * self._jetConeSize) if pt_sum > 0 else 0)

    def fillGenHadronNSubjettiness(self, event, fatjets):
        """Fill gen hadron N-subjettiness variables into output branches.
        """
        for idx in [1, 2]:
            prefix = 'fj_%d_' % idx
            if idx - 1 < len(fatjets) and fatjets[idx - 1].is_qualified:
                fj = fatjets[idx - 1]
                self.out.fillBranch(prefix + "gen_hadron_tau1", fj.gen_hadron_tau1)
                self.out.fillBranch(prefix + "gen_hadron_tau2", fj.gen_hadron_tau2)
                self.out.fillBranch(prefix + "gen_hadron_tau3", fj.gen_hadron_tau3)
                self.out.fillBranch(prefix + "gen_hadron_tau4", fj.gen_hadron_tau4)
