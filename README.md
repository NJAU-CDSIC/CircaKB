# CircaKB Data Analysis

This repository deposited the codes for analyzing the circadian patterns of gene expression.

## Code files

- `initData.R`: Import the data file (.csv)  and obtain the time-course gene expression matrix.
- `batchProcessing.R`: Based on the experimental conditions,  7 models for circadian oscillation detection or 5 models for differential rhythmicity analysis will be implemented through this interface.
- `Methods.R`: This script file includes the source code for all 12 computational models. The code for each model can be freely obtained from the original authors via their publications.
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


--Running the CircaKB analysis requires configuring the `params_x.json` file.

--Command line:

To execute the script, you have two options:

1.Using `./RUN.sh`

2.Using `nohup Rscript --vanilla initData.R params_1.json 2>&1 | tee initData.out || { echo 'initData_1.R failed' ; exit 1; } &`

--Model output:  will generate a folder containing result files, which are stored in .csv format.

--For example:
1. Configure the params_1.json file with the following parameters:
{
  "dup": 1,
  "startTime": 8,
  "samplingCount": 6,
  "interval": 4,
  "idenCounts": 2,
  "dataMode": "mode1",
  "singleGroup": ["H", "L"],
  "compareGroup": ["HVs.L"],
  "singleMethods": [1, 2, 3, 4],
  "diffMethods": [1, 2, 3, 4],
  "geneNum": 100,
  "fileDir": "Example_data/Example 1.csv",
  "saveAddr": "Result"
}


2. Run the following command:
nohup Rscript --vanilla initData.R params_1.json 2>&1 | tee initData.out || { echo 'initData_1.R failed' ; exit 1; } &


