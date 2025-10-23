from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from .HeavyFlavSimpleMatchingBaseProducer import HeavyFlavSimpleMatchingBaseProducer


class SimpleMatchingProducer(HeavyFlavSimpleMatchingBaseProducer):

    def __init__(self, **kwargs):
        super(SimpleMatchingProducer, self).__init__(channel='qcd', **kwargs)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(SimpleMatchingProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.correctJetsAndMET(event)
        if len(event.fatjets) == 0:
            return False

        self.loadGenHistory(event, event.fatjets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, event.fatjets)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def SimpleMatchingTree_2016(): return SimpleMatchingProducer(year=2016)
def SimpleMatchingTree_2017(): return SimpleMatchingProducer(year=2017)
def SimpleMatchingTree_2018(): return SimpleMatchingProducer(year=2018)
