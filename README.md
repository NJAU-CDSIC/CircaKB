# CircaKB Data Analysis

This repository deposited the codes for analyzing the circadian patterns of gene expression.

## Code files

- `initData.R`: Import the data file (.csv)  and obtain the time-course gene expression matrix.
- `batchProcessing.R`: Based on the experimental conditions,  7 models for circadian oscillation detection or 5 models for differential rhythmicity analysis will be implemented through this interface.
- `Methods.R`: This script file includes the source code for all 12 computational models.
- `params_x.json`: Parameters asociated with uploading dataset (import to initData.R)
- `Example_data`: This folder hosts three example datasets.
- `RUN.sh`: Command lines to run the scripts on the shell.

### Configuration

Data analysis has been tested in Ubuntu 20.04 using R-4.2.1.
Ensure your environment has the following dependencies:
| Package       | Version  |
| ------------- | -------- |
| jsonlite      | v1.8.8   |
| stringr       | v1.5.1   |
| dplyr         | v1.1.4   |
| tidyr         | v1.3.0   |
| fdrtool       | v1.2.17  |
| cosinor2      | v0.2.1   |
| rain          | v1.30.0  |
| MetaCycle     | v1.2.0   |
| GeneCycle     | v1.1.5   |
| DODR          | v0.99.2  |
| circacompare  | v0.2.0   |
| diffCircadian | v0.0.0   |
| nloptr        | v2.0.3   |
| limorhyde     | v1.0.1   |
| limma         | v3.52.4  |

## How to Use

--Command line:

./RUN.sh

--Model output: will generate a file called xxx.csv

--For example:
1. Open a terminal and navigate to the project directory.
2. Run the following command to execute the `initData.R` file and read in the parameters from `params.json`: `./RUN.sh`

