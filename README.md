# FastST
Microbiomes play crucial roles in human health, disease development, and global ecosystem functioning. Understanding the origins, movements, and compositions of microbial communities is essential for unraveling the principles governing microbial ecology. Microbial source tracking (MST) approaches have emerged as valuable tools for quantifying the proportions of different microbial sources within target communities, enabling researchers to track transmissions between hosts and environments, identify similarities between microbiome samples, and determine sources of contamination in various settings. Current MST methods like SourceTracker and FEAST have advanced the field by employing Bayesian and expectation-maximization approaches, respectively, but are limited by computational inefficiency with high-dimensional data and inability to infer directionality in source-sink relationships.

This study presents a novel computational framework for microbial source tracking called SourceTracer. SourceTracer infers the relative contributions of source environments to sink microbiomes while also determining directionality when source-sink relationships are not predefined. Through extensive simulation studies with varying numbers of sources and complexity, SourceTracer demonstrates superior performance in both accuracy and computational efficiency compared to FEAST and SourceTracker, maintaining consistent execution times even as the number of source environments increases. Furthermore, the proposed method achieved over 90\% accuracy in directionality inference across all tested scenarios, even when multiple major sources are present, broadening its applicability in practical microbiome research and environmental monitoring. 

This is a source code for FastST, data simulation, and three scenario studies in **"FastST: An efficient tool for inferring decomposition and directionality of microbial community"** paper.

## Requirements
* R (>= 4.1.2)
* npreg (v.1.1.0)
* gtools (v.3.9.5)

## Usage
Clone the repository or download source code files, and install the required packages if needed. **We wrote a detailed software tutorial to run our FastST in "README.docx", and highly recommend users to read this document before running.**


## Simple Usage
* Make sure "FastST" is the current working directory. If you are running our tools in Linux or Mac, please change the permission for the standalone files running :
```
chmod +x ./Code/FastST
```

* Users can simply run following to analyze their own microbiome data. The file “userdata.txt” contains an N by K matrix where the 1st column is the observed sink counts on N taxa, and the rest columns are the counts of the observed K sources.
```
./Code/FastST userdata.txt
```

* After running FastST, a `Result` directory will be created in the current working directory, containing the output file `./Result /results.txt`. In this file, the first column represents the estimated proportions computed by FastST, and the second column shows the standard errors. The first row corresponds to the unknown source. 

* Again, we highly recommend users to read **README.docx** before running.
* `example_comparison_tools` directory contains an example of our running code for the case of K = 10 in Scenario 1 (fully simulated microbiome data), demonstrating how we executed the comparison tools.
* `simulation_dataset` directory contains the fully simulated microbiome data in Scenario 1.



## Contact
If you have any questions or problems, please contact to **joungmin AT vt.edu**.
