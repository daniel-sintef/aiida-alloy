#!/usr/bin/env sh
../launch_workflow_alalloy_scf.py \
    --code_node "pw-v6.3" \
    --structure_group_name "BetaPrime_structures" \
    --workchain_group_name "BetaPrime_structures_scf" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 28800 \
    --max_active_calculations 300 \
    --sleep_interval 600 \
#    --run_debug
