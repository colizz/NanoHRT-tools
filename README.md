# NanoHRT-tools

## Set up CMSSW (on EL9 machines)

```bash
cmsrel CMSSW_15_0_10
cd CMSSW_15_0_10/src
cmsenv
```
> Note: no need to set up official NanoAOD-tools as it has been integrated into CMSSW.

## Get customized NanoAOD tools for HeavyResTagging (NanoHRT-tools)

```bash
git clone https://github.com/colizz/NanoHRT-tools.git PhysicsTools/NanoHRTTools -b dev/nanov15
```

## Compile

```bash
scram b -j8
```

## Update Note

**October 2025:**
This update adapts the framework based on the Run 2 UL setup [1] and two subsequent improvements for early Run 3 (2022/2023, for processing NanoAOD v12) [2,3], and makes it compatible with all currently used NanoAOD versions (NanoAOD v9, v12, v15). When running on NanoAOD v9/v12 samples, the framework gives consistent results to [2,3]

[1] https://github.com/colizz/NanoHRT-tools/tree/dev-UL-0201

[2] https://github.com/lpaizano/NanoHRT-tools/tree/dev/run3

[3] https://github.com/zichunhao/NanoHRT-tools/tree/wz-calibration

Changes:

- Moved files in `src/interface/python/data` from original NanoAOD-tools to NanoHRT-tools if they are not migrated to CMSSW's NanoAOD-tools.
- Specialized support for different NanoAOD versions, including: fatjet taggers, jet b-tag WPs, usage of MET branches, jet corrections.
- Alignment with the latest data campaigns: luminosity values, golden JSON, PU reweighting files (FIXME), lepton ID/isolation, JEC/JER.
- Refactoring of the `qcd` channel
- Updates to JetID logic: in nanoAOD v12, `Jet_jetId` is preserved but [re-computation is recommended](https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV#nanoAOD_Flags).
- Added `jet_veto_maps` following the logic in [2], updated to use the [latest minimal jet selection criteria](https://cms-jerc.web.cern.ch/Recommendations/#run-3).

<details>

<summary>**Cross validation with early NanoHRT-tools branches**</summary>

**1. Validation with Run 2 UL setup for the `qcd` channel (deriving sfBDT SFs) [1]**

[1] https://github.com/colizz/NanoHRT-tools/tree/dev-UL-0201

Configure the [`runHeavyFlavTrees.py`](run/runHeavyFlavTrees.py) script by updating the `default_config` dictionary
```python
default_config.update({
    'nano_version': 'V9',
    'fill_sv': True,
})
```

Then run the production
```bash
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/val_nanov9 --jet-type ak8 --channel qcd --sample-dir samples_nanov9 --year 2018 -n 1
```

**2. Validation with Run 2 UL setup for the `muon` channel (deriving top/W SFs) [1,1a]**

[1] https://github.com/colizz/NanoHRT-tools/tree/dev-UL-0201

[1a] https://github.com/hqucms/NanoHRT-tools/tree/dev/UL

Configure the [`runHeavyFlavTrees.py`](run/runHeavyFlavTrees.py) script by updating the `default_config` dictionary
```python
default_config.update({
    'nano_version': 'V9',
})
```

Then run the production
```bash
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/val_nanov9 --jet-type ak8 --channel muon --sample-dir samples_nanov9 --year 2018 -n 1
```

**3. Validation with early Run 3 setup for the `muon` channel (deriving top/W SFs) [2]**

[2] https://github.com/lpaizano/NanoHRT-tools/tree/dev/run3

Configure the [`runHeavyFlavTrees.py`](run/runHeavyFlavTrees.py) script by updating the `default_config` dictionary
```python
default_config.update({
    'nano_version': 'V12',
    'use_existing_jet_ids': True, # a jetId bug has been identified. Latest recommendation is to re-compute jetId via jet branches (set it to False) but here we use the existing jetId for cross validation
    'jec': True, # should re-compute JECs for NanoAOD v12
})
```

Then run the production
```bash
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/val_nanov12 --jet-type ak8 --channel muon --sample-dir samples_nanov12 --year 2022EE -n 1
```

**4. Validation with early Run 3 setup for the `qcd` channel (deriving sfBDT SFs) [3]**

[3] https://github.com/zichunhao/NanoHRT-tools/tree/wz-calibration

Configure the [`runHeavyFlavTrees.py`](run/runHeavyFlavTrees.py) script by updating the `default_config` dictionary
```python
default_config.update({
    'nano_version': 'V12',
    'fill_sv': True,
    'custom_tagger_list': ["globalParT_QCD0HF", "globalParT_QCD1HF", "globalParT_QCD2HF", "globalParT_Xbb", "globalParT_Xcc", "globalParT_XbbVsQCD", "globalParT_massRes", "globalParT_massVis"], # presented in DAZSLE custom NanoAOD v12 samples
    'use_existing_jet_ids': True, # a jetId bug has been identified. Latest recommendation is to re-compute jetId via jet branches (set it to False) but here we use the existing jetId for cross validation
    'jec': True, # should re-compute JECs for NanoAOD v12
})
```

Then run the production
```bash
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/val_nanov12 --jet-type ak8 --channel qcd --sample-dir samples_nanov12 --year 2022EE -n 1
```

</details>

## Production

<!-- ```bash
cd PhysicsTools/NanoHRTTools/run
```

##### Make trees to produce ntuples for heavy flavour tagging (bb/cc) measurement

```bash
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20230926_ULNanoV9 --jet-type ak8 --channel qcd --year 2018 --sfbdt 0 -n 1
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20230926_ULNanoV9 --jet-type ak8 --channel qcd --year 2018 --run-data --sfbdt 0 -n 1

python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20230926_ULNanoV9 --jet-type ak8 --channel qcd --year 2017 --sfbdt 0 -n 1
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20230926_ULNanoV9 --jet-type ak8 --channel qcd --year 2017 --run-data --sfbdt 0 -n 1

python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20230926_ULNanoV9 --jet-type ak8 --channel qcd --year 2016 --sfbdt 0 -n 1
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20230926_ULNanoV9 --jet-type ak8 --channel qcd --year 2016 --run-data --sfbdt 0 -n 1

python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20230926_ULNanoV9 --jet-type ak8 --channel qcd --year 2015 --sfbdt 0 -n 1
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20230926_ULNanoV9 --jet-type ak8 --channel qcd --year 2015 --run-data --sfbdt 0 -n 1
```

where, `/eos/<some-eos-path-on-lxplus>/` is some path on EOS you have write access to.

Follow the instruction on screen to submit condor jobs. After all condor jobs finish, run the same command appended with ` --post`, to merge the trees.
-->

<details>

<summary>**Production recipes for bookkeeping (keep updating)**</summary>

**For `qcd` channel:**

**A. For generating gen hadron N-subjettiness variables for sfBDT training.**

Updating the `default_config` dictionary:
```python
default_config.update({
    'nano_version': 'V15',
    'fill_sv': True,
    'require_sv_cut': False, 'run_gen_hadron_nsubs': True, # for qcd channel -> dedicated for generating gen hadron N-subjettiness variables for sfBDT training
    'jec': True,
})
```
Running the production (after properly configuring the samples to run in e.g. [`run/samples_nanov15/qcd_2024_MC.yaml`](run/samples_nanov15/qcd_2024_MC.yaml)):
```bash
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20251024_ULNanoV15_gen_hadron_nsubs --jet-type ak8 --channel qcd --sample-dir samples_nanov15 --year 2024 -n 1
```

**B. For nominal `qcd` channel production.**

Updating the `default_config` dictionary:
```python
default_config.update({
    'nano_version': 'V15',
    'fill_sv': True,
    'jec': True,
})
```
Running the production (after properly configuring the samples to run in e.g. [`run/samples_nanov15/qcd_2024_MC.yaml`](run/samples_nanov15/qcd_2024_MC.yaml)):
```bash
python runHeavyFlavTrees.py -o /eos/<some-eos-path-on-lxplus>/20251024_ULNanoV15 --jet-type ak8 --channel qcd --sample-dir samples_nanov15 --year 2024 -n 1
```

</details>

```bash
python runPostProcessing.py [-i /path/of/input] -o /path/to/output -d datasets.yaml --friend 
-I PhysicsTools.NanoHRTTools.producers.hrtMCTreeProducer hrtMCTree -n 1
```

To merge the trees, run the same command but add `--post -w ''` (i.e., set `-w` to an empty string (`''`) -- we do not add the cross sections, but simply reweight signals to match the QCD spectrum afterwards).


##### Make trees for heavy flavour tagging (bb/cc) or top/W data/MC comparison and scale factor measurement:

```bash
python runHeavyFlavTrees.py -i /eos/uscms/store/user/lpcjme/noreplica/NanoHRT/path/to/input -o /path/to/output 
(--sample-dir custom_samples) --jet-type [ak8,ak15] --channel [photon|qcd|muon|inclusive|higgs|mutagged|simple-matching] --year [2016APV|2016|2017|2018] -n 10 
(--batch) (--run-data) (--run-syst)
(--condor-extras '+AccountingGroup = "group_u_CMST3.all"')
```

Command line options:

  - the preselection and basic configurations for each channel is coded in `runHRTTrees.py`. Remember to set them to correct values before submitting jobs (the values will go to `metadata.json` after a job is created).
  - add `--run-data` to make data trees
  - add `--run-syst` to make the systematic trees
  - can run data & MC for multiple years together w/ e.g., `--year 2016APV,2016,2017,2018`. The `--run-data` option will be ignored in this case. Add also `--run-syst` to make the systematic trees.
  - use `--sample-dir` to specify the directory containing the sample lists. Currently we maintain two sets of sample lists: the default one is under samples_* (e.g. `--sample-dir [samples_nanov9](run/samples_nanov9)`) which is used for running over official NanoAOD datasets remotely, and the other one is [custom_samples](run/custom_samples) which is used for running over privately produced NanoAOD datasets locally. To run over the private produced samples, ones needs to add `--sample-dir custom_samples` to the command line.
  - the `--batch` option will submit jobs to condor automatically without confirmation
  - remove `-i` to run over remote files (e.g., official NanoAOD, or private NanoAOD published on DAS); consider adding `--prefetch` to copy files first before running
  - **[NEW]** use `--condor-extras` to pass extra options to condor job description file.
     
More options of `runPostProcessing.py` or `runHRTTrees.py` (a wrapper of `runPostProcessing.py`) can be found with `python runPostProcessing.py -h` or `python runHRTTrees.py -h`, e.g.,

 - To resubmit failed jobs, run the same command but add `--resubmit`.

 - To add cross section weights and merge output trees according to the config file, run the same command but add `--post`. The cross section file to use can be set with the `-w` option. 


### Truth-matching criteria

For maximal flexibility, a number of truth-matching varibles are defined in [HeavyFlavBaseProducer](python/producers/HeavyFlavBaseProducer.py) for hadronically decaying top quarks and W, Z, Higgs bosons. For W/Z/H we define:

 - `fj_idx_dr_X`: deltaR of the fatjet to the nearest **hadronically decaying** X particle. **If found, this top quark `X` is then used to define all the following variables.** Default to 99 if no hadronically decaying X in the event.
 - `fj_idx_dr_X_daus`: max deltaR between the fatjet and the two quarks from X decay.
 - `fj_idx_X_pt`: pt of X
 - `fj_idx_X_decay`: max abs(pdgId) of the two quarks from X decay. For H/Z, this means 5: bb, 4: cc, <4: qq. For W, this means 4: cx, <4: qq. Default to 0 if no hadronically decaying X in the event.

Top quark is treated a bit differently:

 - `fj_idx_dr_T`: deltaR of the fatjet to the nearest **hadronically decaying** top quark. **If found, this top quark `T` is then used to define all the following variables.** Default to 99 if no hadronically decaying top in the event.
 - `fj_idx_dr_T_b`: deltaR between the fatjet and the b quark from the hadronic `T` decay.
 - `fj_idx_dr_T_Wq_(max|min)`: max|min deltaR between the fatjet and the two quarks from the W decay.
 - `fj_idx_T_Wq_(max|min)_pdgId`: pdgId (**w/o** taking the absolute value) of the corresponding two quarks from W decay.
 - `fj_idx_T_pt`: pt of `T`

#### Truth-matching criteria for top/W tagging scale factors


 - top-matched: all three quarks contained in the fatjet
   - `fj_1_dr_T_b<jetR && fj_1_dr_T_Wq_max<jetR`
 - W-matched: only the two W quarks contained, the b quark is outside the jet cone (if the W is from top quark decay)
   - `((fj_1_T_Wq_max_pdgId==0 && fj_1_dr_W_daus<jetR) || (fj_1_T_Wq_max_pdgId!=0 && fj_1_dr_T_b>=jetR && fj_1_dr_T_Wq_max<jetR))`
   - **[Note]** the first part is mainly intended for tW events where the top quark decays leptonically, and the W boson decays hadronically. This can be a sizeable contribution to the W-matched events and needs to be taken into account properly. The trick here makes use of the fact that `fj_1_T_Wq_max_pdgId` is non-zero only if there is a hadronic top in the event.
 - unmatched: defined as `(NOT top-matched) and (NOT W-matched)`, i.e.,
   - `!(fj_1_dr_T_b<jetR && fj_1_dr_T_Wq_max<jetR) && !((fj_1_T_Wq_max_pdgId==0 && fj_1_dr_W_daus<jetR) || (fj_1_T_Wq_max_pdgId!=0 && fj_1_dr_T_b>=jetR && fj_1_dr_T_Wq_max<jetR))`

[Extra] For selecting specifically W->cx decays from the W-matched jets:

 - W(cx)-matched:
   - `((fj_1_T_Wq_max_pdgId==0 && fj_1_dr_W_daus<jetR && fj_1_W_decay==4) || (fj_1_T_Wq_max_pdgId!=0 && fj_1_dr_T_b>=jetR && fj_1_dr_T_Wq_max<jetR && (abs(fj_1_T_Wq_max_pdgId)==4 || abs(fj_1_T_Wq_min_pdgId)==4)))`

### Checklist when updating to new data-taking years / production campaigns

- [ ] triggers
- [ ] lumi values
- [ ] golden JSON
- [ ] PU rewgt
- [ ] lepton ID/ISO
- [ ] b-tag WP
- [ ] JEC/JER
- [ ] MET filters
- [ ] MET recipes (if any)
- [ ] samples (check also those in PRODUCTION status)
