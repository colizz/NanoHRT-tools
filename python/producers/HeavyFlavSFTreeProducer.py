from .HeavyFlavPhotonSampleProducer import PhotonSampleProducer
from .HeavyFlavQCDSampleProducer import QCDSampleProducer
from .HeavyFlavMuonSampleProducer import MuonSampleProducer
from .HeavyFlavDibosonSampleProducer import DibosonSampleProducer
from .HeavyFlavInclusiveSampleProducer import InclusiveSampleProducer
from .HeavyFlavHiggsSampleProducer import HiggsSampleProducer
from .HeavyFlavMuTaggedSampleProducer import MuTaggedSampleProducer
from .HeavyFlavSimpleMatchingProducer import SimpleMatchingProducer

def heavyFlavSFTreeFromConfig():
    import yaml
    with open('heavyFlavSFTree_cfg.json') as f:
        cfg = yaml.safe_load(f)
        channel = cfg['channel']
        del cfg['channel']
    if channel == 'photon':
        return PhotonSampleProducer(**cfg)
    elif channel == 'qcd':
        return QCDSampleProducer(**cfg)
    elif channel == 'muon':
        return MuonSampleProducer(**cfg)
    elif channel == 'diboson':
        return DibosonSampleProducer(**cfg)
    elif channel == 'inclusive':
        return InclusiveSampleProducer(**cfg)
    elif channel == 'higgs':
        return HiggsSampleProducer(**cfg)
    elif channel == 'mutagged':
        return MuTaggedSampleProducer(**cfg)
    elif channel.startswith('simple-matching-'):
        return SimpleMatchingProducer(**cfg)
    else:
        return RuntimeError('Unsupported channel %s' % channel)
