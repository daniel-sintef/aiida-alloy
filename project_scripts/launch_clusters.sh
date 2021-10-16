#!/usr/bin/env sh
../aiida_launch_workflow_alalloy.py \
    --code_node "pw6.3_helvetios-bfgspatch" \
    --structure_group_name "Cluster-CuinAl_structures" \
    --workchain_group_name "Cluster-CuinAl_relax" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 64800 \
    --max_active_calculations 300 \
    --sleep_interval 600 \
    -cm relax \
    --number_of_nodes 16
