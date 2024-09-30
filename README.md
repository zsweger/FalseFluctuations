# False Fluctuations Paper

This repository documents the code used in the analysis of my paper on how to avoid measuring false fluctuations signals. The paper used UrQMD simulations of Au+Au collisions at $\sqrt{s_{NN}}$ = 3.9 GeV, and two toy models of detector scenarios.

This code was written to be run on Brookhaven National Lab's RCF computing clusters.

## Setup
Before making any paper figures, you first need to generate UrQMD samples.

```bash
setup 64bits
cd generateUrQMD/
```

Change the number of jobs in `conf.mk`, line 2. Each job simulates 8000 collisions.

Set the collision center-of-mass energy per nucleon in `urqmd.cc` line 54.

```bash
rm nohup.out
make
make clone
nohup ./queue_daemon > nohup.out &
```

Once all the jobs have finished, run
```bash
bash pull_data.sh
```
And check that the root files are in the directory `roots/`

Our UrQMD files are now generated and we can move on to processing them.

## Analysis

This step creates root files containing event trees of proton number given pT and y cuts, refMult3, nPart, and impact parameter.

```bash
setup 32bits
cd generateTrees/
```

Update the energy in `EnergyConfig.h`.

Clean the directory.
```bash
./clean.sh
```
And compile:
```bash
make
```

Now generate the file list of UrQMD files you generated:
```
realpath ../generateUrQMD/roots/*.root > UrQMD.list
```

Check that the file list name in each of the `Analysis*.xml` files matches this name.

Update the paths in `Analysis*.xml`.

## Instructions to Reproduce Figures 1, 6, 7, 9, and 10

All of these figures use a 0.2 percent (fraction=0.002) pileup rate. Make sure you ran `make` as instructed above, and do:
```bash
./submitSinglePileupRate.sh
```
This generates files in `outdir/`. We need to combine these, excluding TTrees:
```bash
hadd -T PileupFigures.root outdir/*0.002pileupRate.root
``` 
This file contains several histograms.

1. The histogram to reproduce Figure 1 is named `refMult3AndProtons_all`
2. The histograms for Figure 6 are named `refMult3_0` (no pileup) and `refMult3_1` (with pileup).
3. The baseline TH2 in Figure 7 is named `refMult3AndProtons_all`. The AB panel is supplemented with `refMult3AndProtons_singledouble`, the BA panel is supplemented with `refMult3AndProtons_doublesingle`, and the BB panel is supplemented with `refMult3AndProtons_doubledouble`.
4. The histograms for the left side of Figure 9 are named `npAA`, `npAB`, `npBA`, and `npBB`.
5. The analysis window shown on the right side of Figure 10 is contained in the histogram titled `hpT_y`. 

To calculate the cumulants for the right side of Figure 9, and the left side of Figure 10, do the following:
```bash
cd ../calculateCumulants/
setup 64bits
cp calc_toymodel1.cc calc.cc
make
``` 

Create a filelist with the root files generated in the previous step:
```bash
realpath ../generateTrees/outdir/*0.002pileupRate.root > file.list
```

To generate the cumulants that went into the right side of Figure 9, do
- `./run 0` for panel AA
- `./run 1` for panel AB
- `./run 2` for panel BA
- `./run 3` for panel BB

To generate the cumulants used in the left side of Figure 10, do 
- `./run 0`  for the truth
- `./run 8`  for p>1.5GeV
- `./run 9`  for p>1.6GeV
- `./run 10` for p>1.7GeV
- `./run 11` for p>1.8GeV
- `./run 12` for p>1.9GeV 


## Instructions to Reproduce Figures 8, 12, and 13

All of these figures use a 1 percent (fraction=0.01) detector failure rate. Make sure you ran `make` as instructed above, and do:
```bash
./submitSingleFailureRate.sh
```
This generates files in `outdir/`. We need to combine these, excluding TTrees:
```bash
hadd -T FailureFigures.root outdir/*0.01failureRate.root
```
This file contains several histograms.

1. The baseline TH2 in Figure 7 is named `refMult3AndProtons_all`. The AB panel is supplemented with `refMult3AndProtons_godbad`, the BA panel is supplemented with `refMult3AndProtons_badgod`, and the BB panel is supplemented with `refMult3AndProtons_badbad`.
2. The histograms for the left side of Figure 12 are named `npAA`, `npAB`, `npBA`, and `npBB`.
3. The analysis window shown on the right side of Figure 13 is contained in the histogram titled `hpT_y`.

To calculate the cumulants for the right side of Figure 12, and the left side of Figure 13, do the following:
```bash
cd ../calculateCumulants/
setup 64bits
cp calc_toymodel2.cc calc.cc
make
```

Create a filelist with the root files generated in the previous step:
```bash
realpath ../generateTrees/outdir/*0.01failureRate.root > file.list
```

To generate the cumulants that went into the right side of Figure 12, do
- `./run 0` for panel AA
- `./run 1` for panel AB
- `./run 2` for panel BA
- `./run 3` for panel BB

To generate the cumulants used in the left side of Figure 13, do
- `./run 0`  for the truth
- `./run 11` for p>1.7GeV
- `./run 12` for p>1.8GeV
- `./run 13` for p>1.9GeV
- `./run 14` for p>2.0GeV
- `./run 15` for p>2.1GeV 



## Instructions to Reproduce Figure 2
To generate the histograms used in Figure 2 of the paper, make sure you ran `make` as instructed above, and do:
```bash
./Calc_Fig2
```
Check for the output root file: `outdir/Figure2Sample.root`.
The histograms in Figure 2 had about 500k samples each. The normal histogram is named `gausProtons`, and the leptokurtotic one is named `kurtoticProtons`. Draw these on the same canvas to reproduce the figure.

To calculate the cumulants, set up with:
```bash
cd ../calculateCumulants/
setup 64bits
cp calc_fig2.cc calc.cc
make
```

Create a filelist with just the root file from the previous step:
```bash
realpath ../generateTrees/outdir/Figure2Sample.* > file.list
```

To calculate the cumulants for the gaussian histogram do:
```bash
./run 0
```
To calculate the cumulants for the leptokurtotic histogram do:
```bash
./run 1
```

## Instructions to Reproduce Figure 11

The pileup fractions used in Figure 11 were:
- 0.00001
- 0.00002
- 0.00005
- 0.0001
- 0.0002
- 0.0005
- 0.001
- 0.002
- 0.005
- 0.01



Make sure you ran `make` as instructed above, and do:
```bash
./submitPileupScan.sh
```
Once all the jobs have finished, calculate the cumulants by doing the following:
```bash
cd ../calculateCumulants/
setup 64bits
cp calc_toymodel1.cc calc.cc
make
```

No need to create a new filelist or delete an old one. Just do
```bash
nohup ./submitPileupScan.sh > nohup_pileupscan.out &
```

The output files are in the `Output` directory.


## Instructions to Reproduce Figure 14

The failure rates used in Figure 14 were:
- 0.00001
- 0.00002
- 0.00005
- 0.0001
- 0.0002
- 0.0005
- 0.001
- 0.002
- 0.005
- 0.01



Make sure you ran `make` as instructed above, and do:
```bash
./submitFailureScan.sh
```

Once all the jobs have finished, calculate the cumulants by doing the following:
```bash
cd ../calculateCumulants/
setup 64bits
cp calc_toymodel2.cc calc.cc
make
```

No need to create a new filelist or delete an old one. Just do
```bash
nohup ./submitFailureScan.sh > nohup_failurescan.out &
```

The output files are in the `Output` directory.

