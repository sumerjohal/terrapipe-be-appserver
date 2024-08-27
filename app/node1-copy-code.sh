#!/bin/bash

# Ensure the PATH is set correctly
export PATH=/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin

# Define the source and destination paths
SOURCE_PATH="admin_user@66.220.3.88:/mnt/md1/NLDAS/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/Month=*/"
DESTINATION_PATH="/mnt/md1/NLDAS/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/"

# SCP command to copy files
/usr/bin/scp -r $SOURCE_PATH $DESTINATION_PATH

# Log the copy time
LOGFILE="/home/user/terrapipe/app/logs/$(date +%Y%m%d)_output.log"
echo "Files copied at ====> $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOGFILE"
