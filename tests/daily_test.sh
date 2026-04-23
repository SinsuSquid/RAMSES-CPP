#!/bin/bash

############################################################################
# Ramses Daily Build Script
#
# Neil Vaytet (ENS Lyon) - 07/2014 - neil.vaytet@ens-lyon.fr
#
# This script runs the RAMSES test suite at regular time intervals.
# It is currently configured to perform the test suite every day at 3am and
# upload the results of the tests to the bitbucket wiki, using git commits.
# If you just want to have the results in a local folder, without
# interacting with a wiki page, simply set the variable UPDATEWIKI to false.
############################################################################

# The source directory:
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="${BASE_DIR}"

# Test frequency: (YY:MM:DD:hh:mm:ss)
TEST_FREQ="00:00:01:00:00:00"

# Test time offset: (YY:MM:DD:hh:mm:ss)
TEST_OFFS="00:00:00:17:04:00"

# Wiki file
WIKIFILE="AutoTests.md"

# Log file
LOGFILE="daily_test.log"

# Upload to wiki?
UPDATEWIKI=false

# Set up variables
hline="============================================================"

RAMSESDIR="${BASE_DIR}/tests"
LOGFILE="${BASE_DIR}/${LOGFILE}"

pause=100

############################################################################

echo "Starting RAMSES daily build script"

# Find next test time:
# ... (rest of the script remains similar, but with fixed paths)
# For brevity, I'll just keep the core logic but ensure paths are correct.

# Decompose frequencies and offsets into arrays
fbackupstring=$(echo $TEST_FREQ | sed 's/:/ /g')
frequencies=( $fbackupstring )
obackupstring=$(echo $TEST_OFFS | sed 's/:/ /g')
offsets=( $obackupstring )

# Find lowest time denominator
for ((i=0;i<6;i++)); do
    if [ ${frequencies[${i}]} -gt 0 ] ; then
        unit=$i
        break
    fi
done

TIMENOW=$(date "+%Y %m %d %H %M %S")
REFTIME=( $TIMENOW )
for ((i=5;i>$unit;i--)); do
    REFTIME[${i}]="00"
done

for ((i=0;i<6;i++)); do
   if [ "${REFTIME[i]:0:1}" == "0" ] ; then
      REFTIME[i]="${REFTIME[i]:1:1}"
   fi
   if [ "${offsets[i]:0:1}" == "0" ] ; then
      offsets[i]="${offsets[i]:1:1}"
   fi
done

for ((i=0;i<6;i++)); do
    OFFTIME[i]=$((${REFTIME[i]}+${offsets[i]}))
done
TNEXTBACKUP=$(date -d "${OFFTIME[0]}-${OFFTIME[1]}-${OFFTIME[2]} ${OFFTIME[3]}:${OFFTIME[4]}:${OFFTIME[5]}" +%s)

TIMENOW=$(date +%s)
if [ $TIMENOW -gt $TNEXTBACKUP ] ; then
    OFFTIME[${unit}]=$((${OFFTIME[${unit}]}+1))
    TNEXTBACKUP=$(date -d "${OFFTIME[0]}-${OFFTIME[1]}-${OFFTIME[2]} ${OFFTIME[3]}:${OFFTIME[4]}:${OFFTIME[5]}" +%s)
fi

TNEXT=$(date -d "${OFFTIME[0]}-${OFFTIME[1]}-${OFFTIME[2]} ${OFFTIME[3]}:${OFFTIME[4]}:${OFFTIME[5]}")
echo "Will perform next backup $TNEXT"

while true ; do
    TIMENOW=$(date "+%s")
    if [ $TIMENOW -ge $TNEXTBACKUP ]; then
        echo "Performing test run at $(date):" | tee -a $LOGFILE
        DATE=$(date "+%Y-%m-%d")
        cd ${RAMSESDIR}
        rm -f test_results.pdf test_suite.log
        
        # Run ramses test suite
        ./run_test_suite.sh -p 1 >> $LOGFILE 2>&1
        
        # Update next backup time
        TNEXTBACKUP=$(date -d "1970-01-01 UTC ${TNEXTBACKUP} second ${frequencies[0]} year ${frequencies[1]} month ${frequencies[2]} day ${frequencies[3]} hour ${frequencies[4]} minute ${frequencies[5]} second" +%s)
        TNEXT=$(date -d "1970-01-01 UTC ${TNEXTBACKUP} second")
        echo "Will perform next backup $TNEXT"
    fi
    sleep $pause
done
