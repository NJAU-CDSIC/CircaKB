nohup Rscript --vanilla initData.R /Volumes/zz_Seagate/CDSIC/R/UGithub/params.json  2>&1 | tee initData.out || { echo 'initData.R failed' ; exit 1; } &


