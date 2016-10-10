#!/bin/bash
#SBATCH -J Test        # Job Name
#SBATCH -o Test.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -p gpudev  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 04:00:00     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -n 16

./boruvka input/NY.txt > NYcode_output
#./b1 input/BAY.txt > BAYcode_output
#./b1 input/COL.txt > COLcode_output
#./b1 input/FLA.txt > FLAcode_output
#./b1 input/NW.txt > NWcode_output
#./b1 input/NE.txt > NEcode_output
#./b1 input/CAL.txt > CALcode_output
#./b1 input/LKS.txt > LKScode_output
#./b1 input/E.txt > Ecode_output

./hybrid input/NY.txt > NYcodeh_output
#./ndj input/BAY.txt > BAYcodeh_output
#./ndj input/COL.txt > COLcodeh_output
#./ndj input/FLA.txt > FLAcodeh_output
#./ndj input/NW.txt > NWcodeh_output
#./ndj input/NE.txt > NEcodeh_output
#./ndj input/CAL.txt > CALcodeh_output
#./ndj input/LKS.txt > LKScodeh_output
#./ndj input/E.txt > Ecodeh_output

./hybrid_omp input/NY.txt > NYcodehy_output
#./hy input/BAY.txt > BAYcodehy_output
#./hy input/COL.txt > COLcodehy_output
#./hy input/FLA.txt > FLAcodehy_output
#./hy input/NW.txt > NWcodehy_output
#./hy input/NE.txt > NEcodehy_output
#./hy input/CAL.txt > CALcodehy_output
#./hy input/LKS.txt > LKScodehy_output
#./hy input/E.txt > Ecodehy_output

