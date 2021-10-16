#!/usr/bin/env sh
echo "using FIDIS"
../aiida_launch_workflow_alalloy.py \
    --code_node "pw6.3_fidis_albert" \
    --structure_group_name "SimpleAlDeformations" \
    --workchain_group_name "SimpleAlDeformations_scf" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 21600 \
    --max_active_calculations 300 \
    --sleep_interval 600 \
    -mac 400 
