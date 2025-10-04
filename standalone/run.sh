#!/bin/bash

# HTCondor job script for NanoAOD to n-tuple conversion
# Arguments: $1=input_file, $2=eos_path

set -e  # Exit on error

INPUT_FILE=$1
OUTPUT_EOS_PATH=$2

echo "========================================="
echo "Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: $(pwd)"
echo "========================================="
echo "Input file: $INPUT_FILE"
echo "Output EOS path: $OUTPUT_EOS_PATH"
echo "========================================="

# Setup ROOT environment (modify according to your setup)
if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
    source /cvmfs/cms.cern.ch/cmsset_default.sh
fi

# Setup CMSSW environment
export RELEASE=CMSSW_15_0_10
if [ -r $RELEASE/src ] ; then
  echo release $RELEASE already exists
else
  scram p CMSSW $RELEASE
fi
cd $RELEASE/src
eval `scram runtime -sh`
cd ../..

# Check if ROOT is available
if ! command -v root &> /dev/null; then
    echo "ERROR: ROOT is not available!"
    exit 1
fi

echo "ROOT version: $(root-config --version)"

# Create local output filename
LOCAL_OUTPUT="output.root"

# Run the ROOT macro
echo "========================================="
echo "Starting ROOT macro execution..."
echo "========================================="

root -l -b -q "nanoAOD_to_ntuple.C++(\"$INPUT_FILE\", \"$LOCAL_OUTPUT\")"

RETURN_CODE=$?

if [ $RETURN_CODE -ne 0 ]; then
    echo "ERROR: ROOT macro failed with return code $RETURN_CODE"
    exit $RETURN_CODE
fi

# Check if output file was created
if [ ! -f "$LOCAL_OUTPUT" ]; then
    echo "ERROR: Output file $LOCAL_OUTPUT was not created!"
    exit 1
fi

echo "========================================="
echo "Output file created successfully"
echo "File size: $(ls -lh $LOCAL_OUTPUT | awk '{print $5}')"
echo "========================================="

# Copy output to EOS
echo "Copying output to EOS..."
xrdcp -f $LOCAL_OUTPUT $OUTPUT_EOS_PATH

if [ $? -eq 0 ]; then
    echo "Successfully copied to: $OUTPUT_EOS_PATH"
else
    echo "ERROR: Failed to copy output to EOS!"
    echo "Local output file remains: $LOCAL_OUTPUT"
    exit 1
fi

echo "========================================="
echo "Job completed successfully at: $(date)"
echo "========================================="

touch dummy.cc
exit 0
