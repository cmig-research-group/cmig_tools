#!/bin/bash
set -euo pipefail

WORKSPACE_HOST="169.228.60.16"
REMOTE_DIR="~/fema"
SSH_USER="sjones"  # Change if needed
SSH_OPTS="-i ~/.ssh/id_ecdsa -o StrictHostKeyChecking=no"  # Adjust for your key

# List of compiled directories from Makefile
ARTIFACT_DIRS=(
  FEMA_DEAP_gencache
  FEMA_DEAP_worker
  FEMA_DEAP_wrapper
  ../FEMA/FEMA_wrapper
)

echo "==> Checking that artifacts exist..."
for dir in "${ARTIFACT_DIRS[@]}"; do
  if [ ! -d "$dir" ]; then
    echo "❌ Missing build output: $dir"
    exit 1
  fi
done

echo "==> Syncing artifacts to ${WORKSPACE_HOST}:${REMOTE_DIR} ..."
rsync -avz -e "ssh ${SSH_OPTS}" "${ARTIFACT_DIRS[@]}" "${SSH_USER}@${WORKSPACE_HOST}:${REMOTE_DIR}/"

echo "✅ Deployment complete."

