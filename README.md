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
rm -rf log
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
ls ../roots/*.root > UrQMD.list
```

Check that the file list name in each of the `Analysis*.xml` files matches this name.

### Instructions to Reproduce Figure 2
To generate the histograms used in Figure 2 of the paper, make sure you ran `make` as instructed above, and do:
```bash
./Calc_Fig2
```
Check for the output root file: `outdir/Figure2Sample.root`.
The histograms in Figure 2 had about 500k samples each. The normal histogram is named `gausProtons`, and the leptokurtotic one is named `kurtoticProtons`. Draw these on the same canvas to reproduce the figure.

To calculate the cumulants:
```bash
cd ../calculateCumulants/
setup 64bits
cp calc_fig2.cc calc.cc
make
```

To calculate the cumulants for the gaussian histogram do:
```bash
./run 0
```
To calculate the cumulants for the leptokurtotic histogram do:
```bash
./run 1
```
