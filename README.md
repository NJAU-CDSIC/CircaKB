# CircaKB Data Analysis

This section is for data analysis in the CircaKB backend, used to process gene expression matrices and generate corresponding results saved to .csv files.

## File Structure

- `initData.R`: Receives user input parameters, initializes necessary information, reads in the gene expression matrix CSV file, and calls `batchProcessing.R` for batch algorithm processing.
- `batchProcessing.R`: Batch calls 7 Circadian oscillation and 5 Differential rhythmicity detection algorithms.
- `Methods.R`: Specific algorithm implementation and result saving.
- `params_x.json`: Parameters required to execute the algorithms, read by `initData.R`.
- `Example`: Contains three types of example data.
- `RUN.sh`: Command lines to run the scripts on the shell are deposited in RUN.sh.

## How to Use
1. Data analysis has been tested in Ubuntu 20.04 using R-4.2.1.
2. Configure `params_x.json`.
3. Use `./RUN.sh` to start the program.

### Preparation

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
| parallel      | v4.2.1   |

### Running Method

Example:
1. Open a terminal and navigate to the project directory.
2. Run the following command to execute the `initData.R` file and read in the parameters from `params.json`: `./RUN.sh`
