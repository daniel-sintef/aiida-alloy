#../aiida_create_antisite_structures.py -se Al,Cu,Vac -sg OQMD-AntiSite_structures --input_group OQMD-RelaxSymmStable_structures --target_supercellsize 108
#../aiida_create_antisite_structures.py -se Cu,Mg,Zn -sg OQMD-AntiSite_structures --input_group TMP-ZnElastic_structures --target_supercellsize 108
#../aiida_create_antisite_structures.py -se Cu,Mg,Zn -sg OQMD-AntiSite_structures  --input_structures 601632  --target_supercellsize 108
#ELEMENTS=Li,Na,Ca,Sc,Ti,Cr,Mn,Co,Ni,Ge,Sr,Y,Zr,Pd,Ag,In,Sn,Sb,La,Pt,Au,Ce,Nd,Gd,Dy
ELEMENTS=Cu,Mg,Si,Zn
ELEMENTS=Al
../aiida_create_antisite_structures_test.py -se $ELEMENTS  -sg Tphase-antisite-many_structures --input_structures 433d9d35-e533-4885-8531-0aab53015ab0  --target_supercellsize 108


