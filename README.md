# AutoParts

AutoParts implements a non-parametric Bayesian model that simultaneously estimates the parameters of a nucleotide phylogenetic model and the partition schemes across predefined components of a sequence alignment. These components, which are specified by the user, are clustered using a [Chinese restaurant process](https://en.wikipedia.org/wiki/Chinese_restaurant_process) (CRP), a representation of the Dirichlet process prior. This prior groups components of the alignment that appear to have evolved under the same evolutionary process while allowing the complexity of the partition scheme to adapt to the data.

Independent partition schemes are inferred for several components of the phylogenetic model, including exchangeability rates, base frequencies, tree length, and among-site rate variation. Thus, different subsets of the alignment may share parameter values for one component of the model while differing for another. Tree topology and branch-length proportions are shared across the predefined data partitions.

# Installation

To compile this software, you must have CMake (>= 3.10) and a C++17-compatible compiler. AutoParts also requires Eigen3 as a dependency, but will attempt to fetch it automatically if it is not already installed.

To compile, run:

```bash
git clone https://github.com/jhuelsenbeck/AnnaAutoparts.git
cd AnnaAutoparts
cmake -S . -B build
cmake --build build
```

The resulting executable, named `autoparts`, will be located in the `build` directory.

# Example

AutoParts takes a specially formatted input file of the form:

```text
NUM_TAXA NUM_CHAR
Taxon1 ATCG...
Taxon2 AGCT...
Taxon3 TAGC...
charset Partition1 = 1-1000
charset Partition2 = 1001-2000/3
charset Partition3 = 1002-2000/3
charset Partition4 = 1003-2000/3
```

The first line contains two numbers separated by a space: the number of taxa and the number of characters (nucleotides). This is followed by the taxon names and their corresponding sequences, with each taxon and sequence on the same line. Finally, the file specifies the user-defined components, or data partitions, of the alignment. In the example above, `Partition1` contains the first 1000 nucleotides of the sequence. The next three partitions span nucleotides 1001–2000 but are divided so that each partition contains every third nucleotide. For a protein-coding sequence, `Partition2`, `Partition3`, and `Partition4` might therefore represent the first, second, and third codon positions, respectively.

Examples of input files can be found in the `Analyses` directory of this repository.

Example usage of AutoParts may look like:

```bash
autoparts <Source_Directory>/Analyses/Empirical/cynmix_analyses/cynmix.in -o cynmix.out -s 100 -n 2000000 -p 1000 -b 20000 -g 4 -mT 1.2 -c no
```

This command runs an MCMC analysis on the `cynmix` dataset. A sample from the chain is written every 100 iterations (`-s 100`), and the chain is run for 2,000,000 iterations (`-n 2000000`). Progress is printed to the screen every 1,000 iterations (`-p 1000`). The analysis also includes an additional 20,000 burn-in iterations (`-b 20000`), during which the MCMC proposal mechanisms may be tuned. Among-site rate heterogeneity is modeled using four discrete-gamma categories (`-g 4`). The mean of the concentration-parameter hyperprior is specified so that the prior expected number of process partitions is 1.2 (`-mT 1.2`), and the concentration parameter is jointly inferred rather than fixed (`-c no`).

# Usage
| Argument                          | Option name                            | Description                                                                                           |
| --------------------------------- | -------------------------------------- | ----------------------------------------------------------------------------------------------------- |
| **Data input/output**             |                                        |                                                                                                       |
| `-i`                              | Input file                             | Input file name                                                                                       |
| `-t`                              | Tree file                              | Tree file name; constrains the analysis to a fixed tree                                               |
| `-o`                              | Output file                            | Output file name                                                                                      |
| **Model parameters**              |                                        |                                                                                                       |
| `-lenMean`                        | Tree length mean                       | Prior mean for tree length                                                                            |
| `-lenSD`                          | Tree length standard deviation         | Prior standard deviation for tree length                                                              |
| `-lenLam`                         | Tree length exponential parameter      | Exponential parameter for tree length                                                                 |
| `-e`                              | ASRV shape exponential parameter       | Exponential parameter for the gamma shape parameter describing among-site rate variation (ASRV)       |
| `-g`                              | Gamma categories                       | Number of discrete gamma categories                                                                   |
| **Dirichlet process prior (DPP)** |                                        |                                                                                                       |
| `-c`                              | Fix concentration parameter            | Specifies whether the concentration parameter is fixed (`yes`) or treated as a random variable (`no`) |
| `-k`                              | Prior mean number of categories        | Prior mean number of categories when the concentration parameter is fixed                             |
| `-m`                              | Concentration parameter prior mean     | Prior mean of the concentration parameter when it is treated as a random variable                     |
| `-mT`                             | Expected number of tables              | Prior mean (parameterized in terms of expected categories) of the concentration parameter when is treated as a random variable            |
| `-v`                              | Concentration parameter prior variance | Prior variance of the concentration parameter when it is treated as a random variable                 |
| `-eLength`                        | Expected tree-length tables            | Expected number of tables for tree length when the concentration parameter is fixed                   |
| `-eShape`                         | Expected gamma-shape tables            | Expected number of tables for the gamma shape parameter when the concentration parameter is fixed     |
| `-ePi`                            | Expected base-frequency tables         | Expected number of tables for base frequencies when the concentration parameter is fixed              |
| `-eTheta`                         | Expected exchangeability-rate tables   | Expected number of tables for exchangeability rates when the concentration parameter is fixed         |
| **MCMC**                          |                                        |                                                                                                       |
| `-n`                              | MCMC cycles                            | Number of MCMC cycles                                                                                 |
| `-p`                              | Print frequency                        | Frequency at which progress is printed to the screen                                                  |
| `-s`                              | Sampling frequency                     | Frequency at which samples are written to file                                                        |
| `-b`                              | Burn-in                                | Number of burn-in iterations                                                                          |
| `-tu`                             | Tuning frequency                       | Frequency at which tunable parameters are adjusted during burn-in                                     |
| `-tLocal`                         | Tree topology tuning                   | MCMC tuning parameter for the tree-topology proposal                                                  |
| `-tTBR`                           | TBR heat tuning                        | MCMC tuning parameter for the heat parameter used in biased TBR proposals                             |
| `-tBrlen`                         | Branch-proportion tuning               | MCMC tuning parameter for the branch-proportions update                                               |
| `-tShape`                         | Gamma-shape tuning                     | MCMC tuning parameter for the gamma shape parameter                                                   |
| `-tFreqs`                         | Base-frequency simplex tuning          | MCMC tuning parameter for base frequencies using the Dirichlet simplex proposal                       |
| `-tRates`                         | Substitution-rate simplex tuning       | MCMC tuning parameter for substitution rates using the Dirichlet simplex proposal                     |
| `-tFreqsSingle`                   | Single base-frequency tuning           | MCMC tuning parameter for base frequencies using the Beta simplex proposal                            |
| `-tRatesSingle`                   | Single substitution-rate tuning        | MCMC tuning parameter for substitution rates using the Beta simplex proposal                          |
| `-tLength`                        | Tree-length tuning                     | MCMC tuning parameter for the tree-length parameter                                                   |

