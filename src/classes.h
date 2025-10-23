#include "PhysicsTools/NanoHRTTools/interface/PyJetResolutionWrapper.h"
#include "PhysicsTools/NanoHRTTools/interface/PyJetResolutionScaleFactorWrapper.h"
#include "PhysicsTools/NanoHRTTools/interface/PyJetParametersWrapper.h"
#include "PhysicsTools/NanoHRTTools/interface/WeightCalculatorFromHistogram.h"
#include "PhysicsTools/NanoHRTTools/interface/ReduceMantissa.h"

PyJetResolutionWrapper jetRes;
PyJetResolutionScaleFactorWrapper jetResScaleFactor;
PyJetParametersWrapper jetParams;
WeightCalculatorFromHistogram wcalc;
ReduceMantissaToNbitsRounding red(12);
