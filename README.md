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
git clone git@github.com:mkilpatr/NanoHiggsTuples.git PhysicsTools/NanoHiggsTuples
```

### Install a faster version of ONNXRuntime

```bash
$CMSSW_BASE/src/PhysicsTools/NanoHiggsTuples/scripts/install_onnxruntime.sh
```

### Compile

```bash
scram b -j16
```

### Test

MC (2018, 102X):

```bash
cmsDriver.py test_nanoTuples_mc2018 -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v21 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC --filein /store/mc/RunIIAutumn18MiniAOD/WJetsToLNu_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/FFF61350-8D94-6D49-8FBD-C53BBBA7A1E9.root --fileout file:nano_mc2018.root --customise_commands "process.options.wantSummary = cms.untracked.bool(True)" >& test_mc2018.log &

less +F test_mc2018.log
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

MC (2018, 102X):

```bash
cmsDriver.py mc2018 -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v21 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC --filein file:step-1.root --fileout file:nano.root --no_exec
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
