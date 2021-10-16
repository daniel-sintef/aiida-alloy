#../aiida_create_antisite_structures.py -se Al,Cu,Vac -sg OQMD-AntiSite_structures --input_group OQMD-RelaxSymmStable_structures --target_supercellsize 108
#../aiida_create_antisite_structures.py -se Cu,Mg,Zn -sg OQMD-AntiSite_structures --input_group TMP-ZnElastic_structures --target_supercellsize 108
#../aiida_create_antisite_structures.py -se Cu,Mg,Zn -sg OQMD-AntiSite_structures  --input_structures 601632  --target_supercellsize 108
../aiida_create_antisite_structures.py -se Cu,Mg,Zn,Al -sg etaPrime  --input_structures eeaae96e-7024-4e0b-9c5c-3a1dbdc1ea66   --target_supercellsize 18

