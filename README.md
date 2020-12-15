# NanoTuples

Custom NanoAOD ntuple producers with additional boosted jet taggers and their PF candidates.

<!-- TOC -->

- [NanoTuples](#nanotuples)
    - [Version](#version)
    - [Setup](#setup)
        - [Set up CMSSW](#set-up-cmssw)
        - [Merge CMSSW branch](#merge-cmssw-branch)
        - [Get customized NanoAOD producers](#get-customized-nanoaod-producers)
        - [Install a faster version of ONNXRuntime](#install-a-faster-version-of-onnxruntime)
        - [Compile](#compile)
        - [Test](#test)
    - [Production](#production)

<!-- /TOC -->

------

## Version

The current version is based on [NanoAODv7](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv7).

Customizations:

- [ParticleNet-MD](https://cds.cern.ch/record/2707946?ln=en) (V00) for AK8 jets
- PF candidates for AK8 jets
- [*not enabled by default*] AK15 jets w/ ParticleNet-MD (V01)

------

## Setup

### Set up CMSSW

```bash
cmsrel CMSSW_11_1_0_pre5
cd CMSSW_11_1_0_pre5/src
cmsenv
```

### Merge CMSSW branch

```bash
git cms-merge-topic -u hqucms:particle-net-onnx-variable-len
```

### Get customized NanoAOD producers

```bash
git clone https://github.com/mkilpatr/NanoTuples.git PhysicsTools/NanoTuples -b HiggsToTauTau
```

### Install a faster version of ONNXRuntime

```bash
$CMSSW_BASE/src/PhysicsTools/NanoTuples/scripts/install_onnxruntime.sh
```

### Compile

```bash
scram b -j16
```

### Test

MC (2016, 94X, MiniAODv3):

```bash
cmsDriver.py test_nanoTuples_mc2016 -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_mcRun2_asymptotic_v8 --step NANO --nThreads 1 --era Run2_2016,run2_nanoAOD_94X2016 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC --filein /store/mc/RunIISummer16MiniAODv3/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/80000/FEC61D42-F5F5-E811-8435-001E67A4061D.root --fileout file:nano_mc2016.root --customise_commands "process.options.wantSummary = cms.untracked.bool(True)" >& test_mc2016.log &

less +F test_mc2016.log
```

Data (2016, 94X, MiniAODv3):

```bash
cmsDriver.py test_nanoTuples_data2016 -n 1000 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v13 --step NANO --nThreads 1 --era Run2_2016,run2_nanoAOD_94X2016 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeData --filein /store/data/Run2016H/MET/MINIAOD/17Jul2018-v2/00000/0A0B71F7-75B8-E811-BAB7-0425C5DE7BE4.root --fileout file:nano_data2016.root --customise_commands "process.options.wantSummary = cms.untracked.bool(True)" >& test_data2016.log &

less +F test_data2016.log
```


MC (2017, 94X, MiniAODv2):

```bash
cmsDriver.py test_nanoTuples_mc2017 -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_mc2017_realistic_v8 --step NANO --nThreads 1 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC --filein /store/mc/RunIIFall17MiniAODv2/DY1JetsToLL_M-50_LHEZpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/60000/F492A0D0-3F56-E811-9387-FA163EB32A35.root --fileout file:nano_mc2017.root --customise_commands "process.options.wantSummary = cms.untracked.bool(True)" >& test_mc2017.log &

less +F test_mc2017.log
```

Data (2017, 94X, MiniAODv2):

```bash
cmsDriver.py test_nanoTuples_data2017 -n 1000 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v13 --step NANO --nThreads 1 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeData --filein /store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/EC099452-C938-E811-9922-0CC47A7C354C.root --fileout file:nano_data2017.root --customise_commands "process.options.wantSummary = cms.untracked.bool(True)" >& test_data2017.log &

less +F test_data2017.log
```


MC (2018, 102X):

```bash
cmsDriver.py test_nanoTuples_mc2018 -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v21 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC --filein /store/mc/RunIIAutumn18MiniAOD/WJetsToLNu_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/FFF61350-8D94-6D49-8FBD-C53BBBA7A1E9.root --fileout file:nano_mc2018.root --customise_commands "process.options.wantSummary = cms.untracked.bool(True)" >& test_mc2018.log &

less +F test_mc2018.log
```

Data (2018ABC, 102X):

```bash
cmsDriver.py test_nanoTuples_data2018abc -n 1000 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v13 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeData --filein /store/data/Run2018B/JetHT/MINIAOD/17Sep2018-v1/60000/FE3C69F0-A0BC-8941-92AF-B0DA1A6270BF.root --fileout file:nano_data2018b.root --customise_commands "process.options.wantSummary = cms.untracked.bool(True)" >& test_data2018b.log &

less +F test_data2018b.log
```

Data (2018D, 102X):

```bash
cmsDriver.py test_nanoTuples_data2018d -n 1000 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v16 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeData --filein /store/data/Run2018D/SingleMuon/MINIAOD/22Jan2019-v2/60002/FB8A58FA-1EC5-2F4E-A1C1-E02EBC1D0DF0.root --fileout file:nano_data2018d.root --customise_commands "process.options.wantSummary = cms.untracked.bool(True)" >& test_data2018d.log &

less +F test_data2018b.log
```

------

## Production

**Step 0**: switch to the crab production directory and set up grid proxy, CRAB environment, etc.

```bash
cd $CMSSW_BASE/src/PhysicsTools/NanoTuples/crab
# set up grid proxy
voms-proxy-init -rfc -voms cms --valid 168:00
# set up CRAB env (must be done after cmsenv)
source /cvmfs/cms.cern.ch/common/crab-setup.sh
```

**Step 1**: generate the python config file with `cmsDriver.py` with the following commands:


MC (2016, 94X, MiniAODv3):

```bash
cmsDriver.py mc2016 -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_mcRun2_asymptotic_v8 --step NANO --nThreads 1 --era Run2_2016,run2_nanoAOD_94X2016 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC --filein file:step-1.root --fileout file:nano.root --no_exec
```

Data (2016, 94X, MiniAODv3):

```bash
cmsDriver.py data2016 -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v13 --step NANO --nThreads 1 --era Run2_2016,run2_nanoAOD_94X2016 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeData --filein file:step-1.root --fileout file:nano.root --no_exec
```


MC (2017, 94X, MiniAODv2):

```bash
cmsDriver.py mc2017 -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_mc2017_realistic_v8 --step NANO --nThreads 1 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC --filein file:step-1.root --fileout file:nano.root --no_exec
```

Data (2017, 94X, MiniAODv2):

```bash
cmsDriver.py data2017 -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v13 --step NANO --nThreads 1 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeData --filein file:step-1.root --fileout file:nano.root --no_exec
```

MC (2018, 102X):

```bash
cmsDriver.py mc2018 -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v21 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC --filein file:step-1.root --fileout file:nano.root --no_exec
```

Data (2018ABC, 102X):

```bash
cmsDriver.py data2018abc -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v13 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeData --filein file:step-1.root --fileout file:nano.root --no_exec
```

Data (2018D, 102X):

```bash
cmsDriver.py data2018d -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v16 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeData --filein file:step-1.root --fileout file:nano.root --no_exec
```


**Step 2**: use the `crab.py` script to submit the CRAB jobs:

For MC:

`python crab.py -p mc_NANO.py --site T2_CH_CERN -o /store/user/$USER/outputdir -t NanoTuples-[version] -i mc.txt --num-cores 1 --send-external -s FileBased -n 2 --work-area crab_projects_mc --dryrun`

For data:

`python crab.py -p data_NANO.py --site T2_CH_CERN -o /store/user/$USER/outputdir -t NanoTuples-[version] -i data.txt --num-cores 1 --send-external -s EventAwareLumiBased -n 100000 -j [json_file] --work-area crab_projects_data --dryrun`


A JSON file can be applied for data samples with the `-j` options.

Golden JSON, 2016:

```
https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
```

Golden JSON, 2017:

```
https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
```

Golden JSON, 2018:

```
https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
```

These command will perform a "dryrun" to print out the CRAB configuration files. Please check everything is correct (e.g., the output path, version number, requested number of cores, etc.) before submitting the actual jobs. To actually submit the jobs to CRAB, just remove the `--dryrun` option at the end.

**Step 3**: check job status

The status of the CRAB jobs can be checked with:

```bash
./crab.py --status --work-area crab_projects_*  --options "maxjobruntime=2500 maxmemory=2500" && ./crab.py --summary
```

Note that this will also **resubmit** failed jobs automatically.

The crab dashboard can also be used to get a quick overview of the job status:

- [https://monit-grafana.cern.ch/d/cmsTMGlobal/cms-tasks-monitoring-globalview?orgId=11](https://monit-grafana.cern.ch/d/cmsTMGlobal/cms-tasks-monitoring-globalview?orgId=11)

More options of this `crab.py` script can be found with:

```bash
./crab.py -h
```
