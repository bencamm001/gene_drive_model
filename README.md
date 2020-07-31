# Readme V1

## Requirements
Julia v1.3.1 or higher


## Quickrun
Run the model_18.jl file in a Julia REPL then run testing_19_polished.jl.
This will run a gene drive scenario and output a plot showing the frequency of various alleles.

As is the case with julia, the first time running the model in the REPL will take much longer than subsequent runs, especially if packages are missing.

## Default Parameters
The preset parameters in testing_19_polished.jl define a gene drive scenario as follows:
Two populations each totalling 1,000,000 individuals with three alleles present; a wild type (WT) allele, a gene drive resistant allele (R) and a gene drive allele (D). The population is fully panmictic and has low levels (10<sup>-4</sup>) of migration between the two populations. Selection pressures are present in both populations with the exposure level in the first population at 0.7 and the second at 0.2. The fitness cost associated with the gene drive is 0.45 and is dominant. The gene drive converts WT alleles at a rate of 0.8 and R alleles at 0.02. The starting frequency of the gene drive in the first population is 0.001 and 0.0 in the second. The R allele starts a frequency of 0.01 in the first population and 0.3 in the second.

Running the model as is should result in the image below.

![default_run](https://user-images.githubusercontent.com/27834989/89003909-5aa1a000-d344-11ea-80df-a9c6f24b1b0a.png "Default Simulation")

The solid blue line represents the frequency of the gene drive allele, the red dashed line represents the frequency of the drive resistant allele and the green dashed line represents the freqeuncy of the wild type allele. With these parameters, the gene drive is able to fixate in the target population but not the neighbouring population.

## How to change parameters
All changes to runs are made in the testing file, the model file just has all the necessary functions.

### Population parameters:

The number of populations modelled is defined by *popsi* (line 10). The number of populations changes the number of values needed for inbreeding, pressure, exposure, migration rates, starting frequencies.

Number of alleles in the population is defined by *indiv* (line 10).

Population size is defined by *popi* (line 10) and the carrying capacity of each population is *max_pop* (line 20).

The number of generations the simulation runs for is defined by *gens* (line 17).

The inbreeding coefficient is defined by *Fis* (line 23) with each element in the array being for each separate population.

The reproductive output of each individual is defined by *clutch* (line 26).

The starting frequency of the gene drive is defined by *starting_freq* (line 235).

The starting frequency of the drive resistant allele is defined by *res_freq* (line 236).

The migration rate between populations is defined by *mig_rates* (line 268).



### Selection parameters
*Pressure* (line 46) dictates whether there is any selection pressures in each population.

*Exposure* (line 50) dictates the extent to which a selection pressure affects a population, e.g. how extensively a herbicide is applied to a region.

The genetic interaction of fitness costs is defined by *dominance* (line 59) with each element in the array showing the interaction for each fitness cost.

Fitness of individual genotypes can be perturbed by *env* (line 64) that adds a random environmental effect to fitness.

The fitness of each allele is defined by *traits1* (line 69). These values can be determined by inferring from the sequence against an arbitrary most fit sequence, or assigned manually. The first row of values is the baseline fitness from the allele, the second row is the fitness cost of the gene drive and the third row is the fitness cost of resistance to the gene drive.

### Sequence parameters
The mutations rate is defined by *mu_rate* (line 93).

The distribution of mutations is controlled by *gamma* (line 96)

### Drive parameters

There are 5 methods to determine the conversion efficiency of the drive: set, normal, beta, bipolar and sequence.

**Set:**
The conversion efficiency of the drive is defined by *conv_set* (line 115) and a drive resistant allele can be set with *res* (line 116).

**Normal:**
The conversion efficiency for each allele is determined through a normal distribution as defined by the mean, *nmean* (line 125), the standard deviation, *nstd* (line 126).

**Beta:**
The conversion efficiency for each allele is determined through a beta distribution as defined by shape 1, *shape1* (line 120), and shape 2, *shape2* (line 121).

**Bipolar:**
The conversion efficiency for each allele is determined through two normal distributions that are sampled at from at a rate defined by a binomial distribution. The frequency that distribution one is sampled from is defined by *percent_peak1* (line 105) with the mean and standard deviation for the two normal distribtuions being *peak1, se1* and *peak2, se2* respectively.

**Sequence:**
The conversion efficiecny of each allele can be inferred from a gRNA sequence. A gRNA is defined by *grna* (line 132) and the binding efficiecny of the gRNA can be modulated with *priming* (line 135). Presently the conversion efficiency is a percent similarity between the gRNA sequence and the allele sequence, however, in future this will be adjusted to reflect the importance of gRNA positional mutations.


