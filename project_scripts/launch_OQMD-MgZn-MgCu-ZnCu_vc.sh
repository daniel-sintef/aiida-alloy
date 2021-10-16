#!/usr/bin/env sh
#    --code_node "pw6.3_daint-bfgspatch" \
#    --code_node "pw6.3_fidis_albert" \

# CAREFUL! CHECK NUMBER NODES & NK
../aiida_launch_workflow_alalloy.py \
    --code_node "pw6.3_daint-bfgspatch" \
    --structure_group_label "OQMD-MgZn-MgCu-ZnCu_structures" \
    --workchain_group_label "OQMD-MgZn-MgCu-ZnCu_vc-relax" \
    --base_parameter_node "166a7865-5796-42f3-86c7-7327c4de7d84" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 85000 \
    --max_active_calculations 300 \
    --sleep_interval 600  \
    --calc_method vc-relax \
    --number_of_nodes 1 \
    -nk 1
