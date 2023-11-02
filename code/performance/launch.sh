#!/bin/bash
#SBATCH --job-name=abl-perf
#SBATCH --time=02:00:00
set -euo pipefail

arch="$(echo "$SLURM_JOB_PARTITION" | cut -d'-' -f1)"
np=$((SLURM_NTASKS_PER_NODE*SLURM_JOB_NUM_NODES))
echo "$(date): starting benchmark for $np processes on $arch"

# benchmark suite: measure computation of individual terms
# run largest case first to crash early if out of memory
for nh in 256 128; do
    nvpps="4 2 1" # nv per core
    [[ nh -le 256 ]] && nvpps="16 8 $nvpps"
    for nvpp in $nvpps; do
        nv=$((np*nvpp))
        echo "$(date): running benchmark configuration $np-$nh-$nv"
        ibrun -n $np $JULIACMD ../../code/BoundaryLayerDynamics.jl/benchmark/suite.jl $nh $nv suite-$arch-$np-$nh-$nv.json
    done
done

# integration benchmark: measure computation of full time steps
[[ "$arch" = icx ]] && nv=1280 || nv=768
echo "$(date): running integration benchmark with nv=$nv"
ibrun -n $np $JULIACMD ../../code/BoundaryLayerDynamics.jl/benchmark/integration.jl 256 $nv "integration-$arch-$np.json"

# fortran code: measure computation of full time steps
export LES_CURRENT_LAUNCH=1 LES_TOTAL_LAUNCHES=1
cd .cache/fortran-les/dns
echo "$(date): running fortran dns benchmark"
ibrun -n $np ../bin/les-mpi-$np > ../../../integration-$arch-$np-fortran-dns.log
cd ../les
echo "$(date): running fortran les benchmark"
ibrun -n $np ../bin/les-mpi-$np > ../../../integration-$arch-$np-fortran-les.log
echo "$(date): benchmark complete"
