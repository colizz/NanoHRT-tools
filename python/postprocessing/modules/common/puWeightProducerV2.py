from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import os
import numpy as np
import correctionlib
correctionlib.register_pyroot_binding()


class puWeightProducerV2(Module):
    """
    PU weights from correctionlib JSON (.json.gz)
    Creates:
      puWeight, puWeightUp, puWeightDown
    """

    def __init__(self,
                 jsonfile,
                 key,
                 name="puWeight",
                 nvtx_var="Pileup_nTrueInt",
                 doSysVar=True,
                 clip_npu=(0.0, 99.0),
                 clip_w=(0.0, 10.0),
                 label_nom="nominal",
                 label_up="up",
                 label_down="down"):
        self.jsonfile = jsonfile
        self.key = key
        self.name = name
        self.nvtxVar = nvtx_var
        self.doSysVar = doSysVar
        self.clip_npu = clip_npu
        self.clip_w = clip_w
        self.label_nom = label_nom
        self.label_up = label_up
        self.label_down = label_down
        self._corr = None
        self.out = None

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # load correction once
        if self._corr is None:
            if self.verbose:
                print(f"[PU] Loading JSON: {self.jsonfile}")
                print(f"[PU] Key: {self.key}")
            cset = correctionlib.CorrectionSet.from_file(self.jsonfile)
            self._corr = cset[self.key]

        self.out = wrappedOutputTree
        self.out.branch(self.name, "F")
        if self.doSysVar:
            self.out.branch(self.name + "Up", "F")
            self.out.branch(self.name + "Down", "F")

    def analyze(self, event):
        w = 1.0
        w_up = 1.0
        w_down = 1.0

        if hasattr(event, self.nvtxVar):
            npu = float(getattr(event, self.nvtxVar))
            npu = float(np.clip(npu, self.clip_npu[0], self.clip_npu[1]))

            w = float(self._corr.evaluate(npu, self.label_nom))

            if self.doSysVar:
                w_up = float(self._corr.evaluate(npu, self.label_up))
                w_down = float(self._corr.evaluate(npu, self.label_down))

            # safety clipping 
            w = float(np.clip(w, self.clip_w[0], self.clip_w[1]))
            if self.doSysVar:
                w_up = float(np.clip(w_up, self.clip_w[0], self.clip_w[1]))
                w_down = float(np.clip(w_down, self.clip_w[0], self.clip_w[1]))

        self.out.fillBranch(self.name, w)
        if self.doSysVar:
            self.out.fillBranch(self.name + "Up", w_up)
            self.out.fillBranch(self.name + "Down", w_down)

        return True


# Standard lambda functions for different eras

puWeight_2016APV_V9 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run2-2016preVFP-UL-NanoAODv9/latest/puWeights.json.gz"
    key="Collisions16_UltraLegacy_goldenJSON",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)

puWeight_2016_V9 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run2-2016postVFP-UL-NanoAODv9/latest/puWeights.json.gz"
    key="Collisions16_UltraLegacy_goldenJSON",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)

puWeight_2017_V9 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run2-2017-UL-NanoAODv9/latest/puWeights.json.gz"
    key="Collisions17_UltraLegacy_goldenJSON",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)

puWeight_2018_V9 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run2-2018-UL-NanoAODv9/latest/puWeights.json.gz"
    key="Collisions18_UltraLegacy_goldenJSON",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)

puWeight_2022_V12 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-22CDSep23-Summer22-NanoAODv12/latest/puWeights.json.gz"
    key="Collisions2022_355100_357900_eraBCD_GoldenJson",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)

puWeight_2022EE_V12 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/puWeights.json.gz"
    key="Collisions2022_359022_362760_eraEFG_GoldenJson",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)

puWeight_2023_V12 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-23CSep23-Summer23-NanoAODv12/latest/puWeights.json.gz"
    key="Collisions2023_366403_369802_eraBC_GoldenJson",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)

puWeight_2023BPix_V12 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/puWeights.json.gz"
    key="Collisions2023_369803_370790_eraD_GoldenJson",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)

puWeight_2024_V15 = lambda: puWeightProducerV2(
    jsonfile="/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/puWeights_BCDEFGHI.json.gz"
    key="Collisions24_BCDEFGHI_goldenJSON",
    name="puWeight",
    nvtx_var="Pileup_nTrueInt",
    doSysVar=True,
    label_nom="nominal",
    label_up="up",
    label_down="down"
)
