#! bin/bash
ssh slurm-login << 'EOF'
conda activate cellpose
srun \
  --job-name=dosefitr_cellpose \
  --partition=gpuq \
  --gres=gpu:A100:1 \
  --ntasks=1 \
  --cpus-per-task=96 \
  --time 24:00:00 \
  --mem=448GB \
  python dosefitr_cellpose.py data/composite \
    --channels 1 \
    --cellprob -1 \
    --flow 0.5 \
    --niter 2000 \
    --diameter 100 \
    --quiet
EOF
