#!/usr/bin/env sh
#    --code_node "pw6.3_eiger" \
../aiida_launch_workflow_alalloy.py \
    --code_node "pw6.3_helvetios-bfgspatch" \
    --structure_group_label "ThetaP-GSF_structures" \
    --workchain_group_label "ThetaP-GSF_vc-relax-zonly" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 85000 \
    --max_active_calculations 300 \
    --sleep_interval 600 \
    --z_movement_only \
    --z_cellrelax_only \
    -cm vc-relax
