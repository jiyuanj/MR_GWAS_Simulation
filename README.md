# MR-GWAS-MR-Simulation

## Introduction

This project simulates GWAS data with varying levels of pleiotropy. Specifically, we define the following parameters:

- $`\beta_{jX}`$: GWAS summary statistics for trait $`X`$ (true effect size of SNP on exposure)
- $`\sigma_{jX}`$: Standard error of $`\beta_{jX}`$
- $`N`$: Sample size for the GWAS
- $`\gamma`$: True causal effect
- $`f_j`$: Effect allele frequency for SNP $`j`$
- $`u_j`$: Pleiotropy from each SNP $`j`$

The simulation produces the following summary statistics:

- **Simulated summary statistics**: $`\hat{\beta}_{jY}`$, $`\sigma_{jY}`$, where
  - $`\hat{\beta}_{jY}`$: SNP effects on trait $`Y`$ incorporating both causal and pleiotropic effects
  - $`\sigma_{jY}`$: Calculated standard error of $`\hat{\beta}_{jY}`$

Different forms of pleiotropy $`u_j`$ are simulated to generate $`Y`$-related data. At last large scale simulation with different pleiotropy tests are conducted to check the robustness of them. 

## Setup

We begin by copying **461 true instruments’ effects** and corresponding standard errors from GWAS data on systolic blood pressure(SBP), which are $`\beta_{jX}`$ and $`\sigma_{jX}`$.

### Simulation Procedure:

1. **Simulated SNP effects on exposure**:  
   $`
   \hat{\beta}_{jX} = \beta_{jX} + \mathcal{N}(0, \sigma_{jX}^2)
   `$
2. **Standard error of SNP effects on outcome**:  
   $`
   \sigma_{jY} = \sqrt{\frac{1}{2 f_j (1 - f_j) N}}
   `$

3. **Simulated SNP effects on trait $ Y $**:  
   $`
   \beta_{jY} = \gamma \cdot \hat{\beta}_{jX} + u_j
   `$
   where $`u_j`$ follows two types of pleiotropy:

   - **Q-style Pleiotropy**:  
     $`
     u_j \sim \mathcal{N}(0, \sigma)
     `$
   - **PRESSO-style Pleiotropy**: A subset of SNPs (5, 10, or 20 out of 461) where  
     $`
     u_j \sim \mathcal{N}(0, \sigma)
     `$

4. **Estimated SNP effects on trait $` Y `$**:  
   $`
   \hat{\beta}_{jY} = \beta_{jY} + \mathcal{N}(0, \sigma_{jY}^2)
   `$

5. **Summary statistics estimation methods**:  
   The estimated SNP effects on trait $` Y `$ incorporate both causal and pleiotropic effects and are analyzed using methods like:
   - Inverse Variance Weighted (IVW)
   - MR-RAPS
   - PRESSO
   - Additional robustness tests.

## Guide for Running Code
To run the code on computing cluster, download this repoistory as a folder. Run this code in terminal
```bash
scp -r "path to this folder" uniquename@greatlakes-xfer.arc-ts.umich.edu:
```
This will upload this downloaded folder in your pc to home directory on your computing cluster. 
Then log into the computing cluster by running 
```bash
ssh uniquename@greatlakes.arc-ts.umich.edu
```
Please make sure that packages **TwoSampleMR**, **mr.divw**, **MRPRESSO** have loaded in your R module on computing cluster. 
Tehn create a batch file by running
```bash
nano your_work.slurm
```
Make sure that add a line of code "cd ..." to access the folder so that script can be ran. Then submit work by running
```bash
sbatch your_work.slurm
```
You can check status of work by running
```bash
## Check Error
cat error.txt
## Check Output
cat output.txt
```
The returned result of simulations will be saved as csv files in "**results**" folder after all code are ran.

## Coding Explanation
Three functions are used to simulate GWAS data and run tests on them.


This function intakes bunch of parameters for GWAS data and output a data.frame of GWAS data with Q-style of pleiotropy.
```r
q_gwas_data_fun(SNP, f_j, beta.exposure, se.exposure, se.outcome, gamma=0.1, 
                            ple.mean, ple.se, match_sign = TRUE)
```
Intake arguments:
- SNP: SNP ID
- f_j: effect allele frequency
- beta.exposure: $`\hat{\beta}_{jX}`$, effect of SNP on exposure with Gaussian noise added.
- se.exposure: $`\sigma_{jX}`$, standard error of $`\beta_{jX}`$
- se.outcome: $`\sigma_{jY}`$, standard error of $`\beta_{jY}`$ by calculation
- gamma: $`\gamma`$, true causal effect
- ple.mean: Mean of simulated pleiotropy
- ple.se: Standard error of simulated pleiotropy
- match_sign: Make sign of simulated pleiotropy match with direction of beta.exposure or not

This function intakes bunch of parameters for GWAS data and output a data.frame of GWAS data with PRESSO-style of pleiotropy.
```r
presso_gwas_data_fun(SNP, f_j, beta.exposure, se.exposure, se.outcome,
gamma, random_indice, mean, var, match_sign)
```
Intake arguments:
- SNP: SNP ID
- f_j: effect allele frequency
- beta.exposure: $`\hat{\beta}_{jX}`$, effect of SNP on exposure with Gaussian noise added.
- se.exposure: $`\sigma_{jX}`$, standard error of $`\beta_{jX}`$
- se.outcome: $`\sigma_{jY}`$, standard error of $`\beta_{jY}`$ by calculation
- gamma: $`\gamma`$, true causal effect
- random_indice: SNPs' that are going to be added pleiotropy to be an outlier
- mean: Mean of simulated pleiotropy
- var: Standard error of simulated pleiotropy
- match_sign: Make sign of simulated pleiotropy match with direction of beta.exposure or not

The output of ```q_gwas_data_fun``` and ```presso_gwas_data_fun``` are both GWAS data.frame with columns compatible with ```mr()``` functions and of dimension: $`\text{Number of SNPS} \times 11`$. Eleven columns are:
- SNP: SNP ID
- f_j
- beta.exposure
- se.exposure
- beta.outcome
- se.outcome
- id.exposure
- exposure
- id.outcome
- outcome
- mr.keep

This function intakes bunch of parameters for GWAS data and output a data.frame of detailed estimated and test result in each cycle of simulation. This function will call previous two GWAS data generating functions while running. 
```r
power(num_sim, SNP, f_j, beta.exposure, se.exposure, se.outcome, gamma=0.1,
ple.mean, ple.se, match_sign = TRUE, type)
```
Intake arguments:
- numsim: Number of cyclyes of simualtions
- SNP: SNP ID
- f_j: effect allele frequency
- beta.exposure: $`\hat{\beta}_{jX}`$, effect of SNP on exposure with Gaussian noise added.
- se.exposure: $`\sigma_{jX}`$, standard error of $`\beta_{jX}`$
- se.outcome: $`\sigma_{jY}`$, standard error of $`\beta_{jY}`$ by calculation
- gamma: $`\gamma`$, true causal effect
- ple.mean: Mean of simulated pleiotropy
- ple.var: Standard error of simulated pleiotropy
- match_sign: Make sign of simulated pleiotropy match with direction of beta.exposure or not
- type: To simulate Q-style or PRESSO-style pleiotropy GWAS data

The output is a data.frame of dimension: $`\text{Number of Simulations} \times 9`$. Nine columns are:
- Q test p-value using IVW
- Egger Intercept test p-value
- Q test p-value using Radial IVW
- IVW estimate of $`\hat{\gamma}`$
- Radial IVW estimate of $`\hat{\gamma}`$
- Weighted median estimate of $`\hat{\gamma}`$
- Egger estimate of $`\hat{\gamma}`$
- Weighted mode estimate of $`\hat{\gamma}`$
- RAPS estimate of $`\hat{\gamma}`$

The output also display messages about power of Q test IVW, Q test Radial IVW and Egger intercept test. 

## Comment
When running ```mr_heterogeneity()```, the Q statistics differ a lot when using IVW and IVW-Radial. This is because when using IVW-Radial, this function automatically uses modified Q test, where the Q statistics is computed by $`Q = \sum_{j=1}^{p} \frac{(\hat{\Gamma}_j - \hat{\beta} \hat{\gamma}_j)^2}{\sigma^2_{\hat{Y}_j} + \beta^2 \sigma^2_{X_j}}`$. When using standard IVW, $`Q = \sum_{j=1}^{p} \frac{(\hat{\Gamma}_j - \hat{\beta} \hat{\gamma}_j)^2}{\sigma^2_{\hat{Y}_j}}`$.

