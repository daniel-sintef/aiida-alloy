#!/usr/bin/env sh
echo "WARNING not using production parameter settings!"
echo "WARNING not using production kpoint settings!"
echo "WARNING using local code"
echo "WARNING using very high mac"
#fb5cf8be-87dd-4757-92da-97d59e1e31b9 #very-fast setting
#9b370584-3f56-471c-a724-dbaadf022ec5 #production setting

#    --kptper_recipang 80 \
../aiida_launch_workflow_alalloy.py \
    --code_node pw6.3_daint-bfgspatch \
    --structure_group_label "TEST-Alonly_structures" \
    --workchain_group_label "TEST-Alonly_elastic" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 60000 \
    --max_active_calculations 300 \
    --sleep_interval 600 \
    -cm elastic  \
    --use_conventional_structure \
    -mac 1000 \
    -nk 1 \
    --number_of_nodes 1
