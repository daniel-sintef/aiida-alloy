#!/usr/bin/env sh
../launch_workflow_alalloy_scf.py \
    --code_node "pw-v6.3" \
    --structure_group_name "Al6xxxDB_structures" \
    --workchain_group_name "Al6xxxDB_structures_calc" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefv1.1_alalloy" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 21600 \
    --max_active_calculations 300 \
    --sleep_interval 600 \
    --run_debug
