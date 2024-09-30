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

## Instructions to Reproduce Figures 1, 6, 7, 9, and 10

Make sure you ran `make` as instructed above, and do:
```bash
star-submit Analysis_ToyModel1.xml
```
This generates files in `outdir/`. We need to combine these, excluding TTrees:
```bash
hadd -T PileupFigures.root outdir/*.root
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
``` 

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
ls ../generateTrees/outdir/Figure2Sample.* > file.list
```

To calculate the cumulants for the gaussian histogram do:
```bash
./run 0
```
To calculate the cumulants for the leptokurtotic histogram do:
```bash
./run 1
```

## Instructions to Reproduce
