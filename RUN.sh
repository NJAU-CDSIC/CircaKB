## Analyze data of type mode1
nohup Rscript --vanilla initData.R /Volumes/zz_Seagate/CDSIC/R/UGithub/params_1.json  2>&1 | tee initData.out || { echo 'initData_1.R failed' ; exit 1; } &

## Analyze data of type mode2
nohup Rscript --vanilla initData.R /Volumes/zz_Seagate/CDSIC/R/UGithub/params_2.json  2>&1 | tee initData_2.out || { echo 'initData_2.R failed' ; exit 1; } &

## Analyze data of type mode3
nohup Rscript --vanilla initData.R params_3.json  2>&1 | tee initData_3.out || { echo 'initData_3.R failed' ; exit 1; } &
