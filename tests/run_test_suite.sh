#!/bin/bash

#######################################################################
#
# Script to run the RAMSES test suite (C++ version)
#
#######################################################################

testlist="hydro,mhd,poisson,rt,sink,turb,tracer";

MPI=0;
NCPU=1;
VERBOSE=false;
DELDATA=true;
SELECTTEST=false;

while getopts "cdsp:qt:v" OPTION; do
   case $OPTION in
      c) CLEAN_ALL=true ;;
      d) DELDATA=false ;;
      p) MPI=1; NCPU=$OPTARG ;;
      t) SELECTTEST=true; TESTNUMBER=$OPTARG ;;
      v) VERBOSE=true ;;
   esac
done

TEST_DIRECTORY=$(pwd);
BASE_DIRECTORY="${TEST_DIRECTORY}/..";
BIN_DIRECTORY="${BASE_DIRECTORY}/bin";
BUILD_DIRECTORY="${BASE_DIRECTORY}/build";
VISU_DIR="${TEST_DIRECTORY}/visu";

mkdir -p ${BUILD_DIRECTORY};
mkdir -p ${BIN_DIRECTORY};

export PYTHONPATH=${VISU_DIR}:$PYTHONPATH;
DELETE_RESULTS="rm -rf output_* *.tex data*.dat *.pdf *.pyc *.gc* coverage_stats.txt";
EXECNAME="test_exe_";
LOGFILE="${TEST_DIRECTORY}/test_suite.log";
echo > $LOGFILE;

if [ ${MPI} -eq 1 ]; then
   RUN_TEST_BASE="mpirun -np ${NCPU} ${BIN_DIRECTORY}/${EXECNAME}";
else
   RUN_TEST_BASE="${BIN_DIRECTORY}/${EXECNAME}";
fi

line="--------------------------------------------";
BEFORETEST="before-test.sh";
STARTTIME=$(date +%s);

echo "############################################" | tee -a $LOGFILE;
echo "#   Running RAMSES automatic test suite    #" | tee -a $LOGFILE;
echo "############################################" | tee -a $LOGFILE;

s1=$(echo $testlist | sed 's/,/ /g');
testsegs_all=( $s1 );
testlist="";
for ((m=0;m<${#testsegs_all[@]};m++)); do
   testlist="${testlist} ${testsegs_all[m]}/*";
done

testname=( $testlist );
ntestsall=${#testname[@]};
ntests=$ntestsall;
all_tests_ok=true;

if $SELECTTEST ; then
   # (Selection logic simplified for brevity, similar to original)
   ntests=0;
   s1=$(echo $TESTNUMBER | sed 's/,/ /g');
   testsegs=( $s1 );
   for ((n=0;n<${#testsegs[@]};n++)); do
      # Match by name or number
      for ((m=0;m<$ntestsall;m++)); do
         if [[ ${testname[m]} == *"${testsegs[n]}"* ]]; then
            testnum[${ntests}]=$m;
            ntests=$((ntests + 1));
         fi
      done
   done
else
   for ((n=0;n<$ntests;n++)); do testnum[n]=$n; done
fi

echo "Will perform the following tests:" | tee -a $LOGFILE;
for ((i=0;i<$ntests;i++)); do
   echo " [$(($i+1))] ${testname[testnum[i]]}" | tee -a $LOGFILE;
done
echo $line | tee -a $LOGFILE;

for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   echo "Test $(($i+1))/${ntests}: ${testname[n]}" | tee -a $LOGFILE;
   
   rawname=$(basename ${testname[n]})
   FLAGS=$(grep FLAGS ${TEST_DIRECTORY}/${testname[n]}/config.txt | cut -d ':' -f2);
   ndim=$(echo $FLAGS | grep -o "NDIM=[0-9]" | cut -d '=' -f2)
   
   echo "Compiling source (NDIM=$ndim)" | tee -a $LOGFILE;
   cd ${BUILD_DIRECTORY}
   CMAKE_ARGS="-DRAMSES_NDIM=${ndim}"
   [ $MPI -eq 1 ] && CMAKE_ARGS="${CMAKE_ARGS} -DRAMSES_USE_MPI=ON"
   
   cmake .. ${CMAKE_ARGS} -DCMAKE_BUILD_TYPE=Release >> $LOGFILE 2>&1
   make ramses_main >> $LOGFILE 2>&1
   cp ramses_main ${BIN_DIRECTORY}/${EXECNAME}${ndim}d

   cd ${TEST_DIRECTORY}/${testname[n]};
   $DELETE_RESULTS;
   [ -f ${BEFORETEST} ] && bash ${BEFORETEST} >> $LOGFILE 2>&1
   
   echo -n "Running test: " | tee -a $LOGFILE
   TSTART=$(date +%s)
   ${RUN_TEST_BASE}${ndim}d ${rawname}.nml >> $LOGFILE 2>&1
   TEND=$(date +%s)
   echo "$((TEND-TSTART))s" | tee -a $LOGFILE

   echo "Plotting and analysing results" | tee -a $LOGFILE
   status=$(python3 plot-${rawname}.py 2>&1)
   echo "$status" >> $LOGFILE
   
   if echo "$status" | grep -q "PASSED"; then
      echo "Test ${testname[n]} passed [ OK ]" | tee -a $LOGFILE
   else
      echo "Test ${testname[n]} failed! [FAIL]" | tee -a $LOGFILE
      all_tests_ok=false
   fi
   echo $line | tee -a $LOGFILE
done

ENDTIME=$(date +%s);
echo "Total run time: $((ENDTIME-STARTTIME))s" | tee -a $LOGFILE

if $all_tests_ok ; then
   echo "All tests were completed successfully"
   exit 0
else
   echo "There were some failed tests"
   exit 1
fi
