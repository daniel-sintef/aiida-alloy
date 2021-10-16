#!/usr/bin/env sh
#    --max_atoms_submit 30 \
#    --max_active_calculations 4 \
echo "WARNING elastic may overlaunch calcs!"
#    --code_node "pw6.3_daint-bfgspatch" \
#    --code_node "pw6.3_helvetios-bfgspatch" \
../aiida_launch_workflow_alalloy.py \
    --code_node "pw6.3_daint-bfgspatch" \
    --structure_group_name "OQMD_structures" \
    --workchain_group_name "OQMD_elastic" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSPefV1.1_Zn-dnlPAW_3" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 60000 \
    --sleep_interval 600 \
    -cm elastic  \
    --strain_magnitudes '-0.01,-0.005,0.005,0.01' \
    --use_conventional_structure \
    --max_active_elastic 10 \
    --max_atoms_submit 30 
