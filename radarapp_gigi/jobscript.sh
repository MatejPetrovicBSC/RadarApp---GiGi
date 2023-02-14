#!/bin/bash
#SBATCH --partition=fpga-sdv
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --constraint=release:22.05

#############################################
## Construct SDV_HOST from X86_HOST        ##
#############################################
# X86_HOST: pickle-{1,2,3}
# SDV_HOST: fpga-sdv-{1,2,3}


X86_HOST="`hostname`"
SDV_HOST="fpga-sdv-`echo ${X86_HOST} | cut -d '-' -f 2`"

printf "******************************\n"
printf "* x86 node: %s\n" "${X86_HOST}"
printf "* SDV node: %s\n" "${SDV_HOST}"
printf "******************************\n\n"

#############################################
## Commands to be executed on the SDV node ##
#############################################
# Make sure SSH keys are configured properly!
# Add "-o StrictHostKeyChecking=no" to disable key check prompt
#
# Make sure you `cd` to the working directory path!
#     You can use the $SLURM_SUBMIT_DIR variable to do so
#     You can also use a global path for <myscript>

module load llvm/EPI-0.7-development-cross

for trigger_ID in 1 2 3 4 5
do 
	ssh -o StrictHostKeyChecking=no ${X86_HOST} "sudo /apps/riscv/fgpa-sdv/ila2prv/sdv-ila.sh ${trigger_ID} 0 ~/radarapp/radarapp/ilatraces/radar_trigger${trigger_ID}" & 
	sleep 10
	ssh -o StrictHostKeyChecking=no ${SDV_HOST} "cd radarapp/radarapp/; ./covariancev2_main-ILA-vector"
	wait
       	ssh -o StrictHostKeyChecking=no ${X86_HOST} "/apps/riscv/fgpa-sdv/ila2prv/ila2prv ~/radarapp/radarapp/ilatraces/radar_trigger${trigger_ID}.ila ~/radarapp/radarapp/ilatraces extend"
done
