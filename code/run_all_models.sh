#!/bin/bash
# Batch run all GAMS FBA models for Wild_Type and O2_mutant
set -e

BASE_DIR="/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY"

module load gams/24.7.4

# Wild Type models
for dap in B6 B8 B10 B12 B15 B18 B22 B30; do
    DIR="$BASE_DIR/Wild_Type/$dap"
    echo "========================================="
    echo "Running Wild_Type $dap ..."
    cd "$DIR"
    gams Endosperm_Model_FBA.gms
    if [ -f results_FBA.txt ]; then
        echo "  -> SUCCESS: results_FBA.txt generated"
    else
        echo "  -> WARNING: results_FBA.txt NOT found"
    fi
done

# O2 mutant models
for dap in O6 O8 O10 O12 O15 O18 O22 O30; do
    DIR="$BASE_DIR/O2_mutant/$dap"
    echo "========================================="
    echo "Running O2_mutant $dap ..."
    cd "$DIR"
    gams Endosperm_Model_FBA.gms
    if [ -f results_FBA.txt ]; then
        echo "  -> SUCCESS: results_FBA.txt generated"
    else
        echo "  -> WARNING: results_FBA.txt NOT found"
    fi
done

echo "========================================="
echo "All models completed!"
