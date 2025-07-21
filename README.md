# FastST
Microbiomes play crucial roles in human health, disease development, and global ecosystem functioning. Understanding the origins, movements, and compositions of microbial communities is essential for unraveling the principles governing microbial ecology. Microbial source tracking (MST) approaches have emerged as valuable tools for quantifying the proportions of different microbial sources within target communities, enabling researchers to track transmissions between hosts and environments, identify similarities between microbiome samples, and determine sources of contamination in various settings. Current MST methods like SourceTracker and FEAST have advanced the field by employing Bayesian and expectation-maximization approaches, respectively, but are limited by computational inefficiency with high-dimensional data and inability to infer directionality in source-sink relationships.

This study presents a novel computational framework for microbial source tracking called SourceTracer. SourceTracer infers the relative contributions of source environments to sink microbiomes while also determining directionality when source-sink relationships are not predefined. Through extensive simulation studies with varying numbers of sources and complexity, SourceTracer demonstrates superior performance in both accuracy and computational efficiency compared to FEAST and SourceTracker, maintaining consistent execution times even as the number of source environments increases. Furthermore, the proposed method achieved over 90\% accuracy in directionality inference across all tested scenarios, even when multiple major sources are present, broadening its applicability in practical microbiome research and environmental monitoring. 

This is a source code for SourceTracer and data simulation in **"SourceTracer: An efficient tool for inferring decomposition and directionality of microbial community"** paper.

## Requirements
* R (>= 4.1.2)
* npreg (v.1.1.0)
* gtools (v.3.9.5)

## Usage
Clone the repository or download source code files, and install the required packages if needed.

## Data Simulation
* Run the below command to generate sink data by using mixed-proportion model, for simulation purpose :
```
source("funST.R")
gendata = function(N, K, C, nmajor = 2, pmajor = 0.9, lb.amajor = 0.1)
```
* **Input**
  * `N`: scalar, number of taxa groups in source/sink
  * `K`: scalar, number of KNOWN sources
  * `C`: scalar, total taxa count in sink
  * `nmajor`: scalar, number of major sources
  * `pmajor`: scalar, sum of major proportions
  * `lb.amajor`: scalar, lower bound of major proportions

* **Output** - a list with the following components:
  * `x.vec`: vector of length N, sink data 
  * `beta.vec`: vector of length N, abundance of sink 
  * `gamma.mat`: matrix of size (K+1) by N, abundance of KNOWN+UNKNOWN sources (note: sum up each row to 1, The 1st row is reserved for the unobserved source)
  * `alpha.vec`: vector of length K+1, proportion of KNOWN+UNKNOWN sources (note: sum up to 1, Unobserved source is always put in the 1st element)

## Estimate the source proportions
* Run the below command to read the source codes :
```
source("funST.R")
alpha.est.vec = glsest(x.vec, gamma.mat, useGLS = T)
```
* **Input**
  * `x.vec`: vector of length N, sink data (note: N is the number of taxa groups in source/sink)
  * `gamma.mat`: matrix of size K by N, abundance of KNOWN sources (note: K is the number of known sources and N is the number of taxa groups in source/sink)
  * `useGLS`: scalar, whether use GLS (TRUE) or OLS (FALSE)
 
* **Output**
  * a vector of estimated proportions for ALL K + 1 sources

## Inference of directionality
* Run the below command to infer direction from source-sink data :
```
source("funST.R")
direc.est.vec = DIRinfer(x.mat, gammaSS.mat)
```
* **Input**
  * `x.mat`: matrix of size K+1 by N, source-sink data (note: K is the number of known sources and N is the number of taxa groups in source/sink)
  * `gammaSS.mat`: matrix of size K+1 by N, abundance of KNOWN source and sink

* **Output**
  * scalar, index of sink from 1 : (K+1)

## R scripts for data simulation and testing SourceTracer in the manuscript
* simu1.R : ESTIMATION OF SOURCE PROPORTIONS
* simu2.R : INFERENCE OF DIRECTIONALITY

## Contact
If you have any questions or problems, please contact to **joungmin AT vt.edu**.
