#!/usr/bin/env sh
#TEST
../aiida_launch_workflow_alalloy.py \
    --code_node 276635 \
    --structure_group_name "test_primitivesf_alonly" \
    --workchain_group_name "test_primitivesf_alonly_vc-zonly" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 10 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 21600 \
    --max_active_calculations 300 \
    --sleep_interval 600 \
    -cm vc-relax \
    -mac 500 \
    -nk 1 \
    --press_conv_thr 1 \
    --z_cellrelax_only \
    --keep_workdir \
    --number_of_nodes 1 

