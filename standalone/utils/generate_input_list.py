#!/usr/bin/env python3

import yaml
import subprocess
import sys

def get_files_from_das(dataset):
    """Query DAS to get list of files for a dataset."""
    cmd = f'dasgoclient -query="file dataset={dataset}"'
    
    try:
        result = subprocess.run(
            cmd, 
            shell=True, 
            capture_output=True, 
            text=True, 
            timeout=300
        )
        
        if result.returncode != 0:
            print(f"ERROR querying DAS for {dataset}: {result.stderr}", file=sys.stderr)
            return []
        
        files = [line.strip() for line in result.stdout.strip().split('\n') if line.strip()]
        return files
    
    except subprocess.TimeoutExpired:
        print(f"ERROR: DAS query timeout for {dataset}", file=sys.stderr)
        return []
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return []

def main():
    yaml_file = "/afs/cern.ch/work/c/coli/hcc/nano/ul/CMSSW_11_1_0_pre5_PY3/src/PhysicsTools/NanoHRTTools/run/samples/simple-matching-allcamp_2018_MC.yaml"
    out_string = ""
    
    # Read YAML file
    try:
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"ERROR: {yaml_file} not found", file=sys.stderr)
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"ERROR parsing YAML: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Process each dataset
    for key, value in data.items():
        if value is None or not isinstance(value, list) or len(value) == 0:
            continue
        print(f'hadd {key}_merged.root {key}_file*.root')
        
    #     dataset = value[0]  # Get the dataset path
    #     print(f"# Processing {key}: {dataset}")
        
    #     # Query DAS for files
    #     files = get_files_from_das(dataset)
        
    #     if not files:
    #         print(f"# WARNING: No files found for {key}", file=sys.stderr)
    #         continue
        
    #     if key.startswith("qcd"):
    #         files = files[:20]  # Limit to first 20 files for QCD samples
    #     else:
    #         files = files[:5]   # Limit to first 5 files for other samples

    #     # Print files in the requested format
    #     for idx, filepath in enumerate(files):
    #         output_name = f"{key}_file{idx}"
    #         out_string += f"{filepath} {output_name}\n"
        
    #     print(f"# Found {len(files)} files for {key}")

    # # Print the final output string
    # with open("/afs/cern.ch/work/c/coli/hcc/nano/ul/CMSSW_11_1_0_pre5_PY3/src/PhysicsTools/NanoHRTTools/standalone/samples/test.txt", "w") as out_file:
    #     out_file.write(out_string)
    # print(out_string)

if __name__ == "__main__":
    main()

