#!/usr/bin/env sh
#    --code_node "pw6.3_daint-bfgspatch" \
../aiida_launch_workflow_alalloy.py \
    --code_node "pw6.3_helvetios-bfgspatch" \
    --structure_group_name "TMP-Zn-AntiSite_structures" \
    --workchain_group_name "OQMD-AntiSite_relax" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 85000 \
    --max_active_calculations 150 \
    --sleep_interval 600 \
    --calc_method relax
