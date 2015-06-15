#!/bin/bash
# Written by Sébastien Wertz
# Resubmit condor jobs which have a non-empty error output file

condor_dir=condor

jobs=`find ${condor_dir}/logs/ -name "*.err" -size +100c -execdir sh -c "basename {} | sed 's/.err//g'" \;`
echo "Resubmit jobs:"
echo ${jobs}

cmd_file=${condor_dir}/inputs/resubmit.cmd
cp ${condor_dir}/base.cmd ${cmd_file}

for job in ${jobs}
do
  echo "" >> ${cmd_file}
  echo "Executable = ${condor_dir}/inputs/${job}.sh" >> ${cmd_file}
  echo "output = ${condor_dir}/logs/${job}.out" >> ${cmd_file}
  echo "error = ${condor_dir}/logs/${job}.err" >> ${cmd_file}
  echo "log = ${condor_dir}/logs/${job}.log" >> ${cmd_file}
  echo "Queue 1" >> ${cmd_file}
done

condor_submit ${cmd_file}

