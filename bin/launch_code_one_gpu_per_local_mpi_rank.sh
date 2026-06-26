#!/bin/bash

# ==============================================
# Launcher for runs that use 1 GPU per MPI rank. 
# For use with applications that do not select the device 
# per local MPI rank within the code, such as 
# standard langauge parallel codes on GPUs.
# ==============================================

# === LOCAL RANK SELECTION ===
# === Uncomment the correct LOCAL_RANK 
# === source for your environment:

LOCAL_RANK=${SLURM_LOCALID}                 # Slurm
# LOCAL_RANK=${PBS_VNODENUM}                # PBS / Torque
# LOCAL_RANK=${OMPI_COMM_WORLD_LOCAL_RANK}  # OpenMPI
# LOCAL_RANK=${PMI_LOCAL_RANK}              # MPICH / Intel MPI / Hydra
# LOCAL_RANK=${MPI_LOCALRANKID}             # Intel OneAPI MPI >=25 

# === If nothing was set, default to 0:
LOCAL_RANK=${LOCAL_RANK:-0}

# === Make the local MPI rank only see the GPU it should use ===
# === Pick combination to uncomment based on your environment:
export CUDA_VISIBLE_DEVICES=${LOCAL_RANK}                 # NVIDIA
#export ROCR_VISIBLE_DEVICES=${LOCAL_RANK}                # AMD ROCm
#export HIP_VISIBLE_DEVICES=${LOCAL_RANK}                 # AMD HIP
#export ONEAPI_DEVICE_SELECTOR=level_zero:${LOCAL_RANK}   # Intel oneAPI / SYCL
#export ZE_AFFINITY_MASK=${LOCAL_RANK}                    # Intel Level Zero
#export I_MPI_OFFLOAD_PIN=0                               # Intel (needed if using one of the above)

# === Fancy launch: use numactl to bind the cpu and memory 
# === location (not usable for all configs)
#exec numactl --cpunodebind=${LOCAL_RANK} \
#             --membind=${LOCAL_RANK} \
#             "$@"

# === Basic launch ===
exec "$@"



