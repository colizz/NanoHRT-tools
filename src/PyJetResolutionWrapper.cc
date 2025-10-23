#include "PhysicsTools/NanoHRTTools/interface/PyJetResolutionWrapper.h"

PyJetResolutionWrapper::PyJetResolutionWrapper(const std::string& filename)
   : JME::JetResolution(filename)
{}
 
PyJetResolutionWrapper::PyJetResolutionWrapper(const JME::JetResolutionObject& object)
  : JME::JetResolution(object)
{}

PyJetResolutionWrapper::PyJetResolutionWrapper()
  : JME::JetResolution()
{}
