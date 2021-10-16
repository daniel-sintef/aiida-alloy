ELEMENTS=Cu,Mg,Si,Zn,Al,Li,Na,Ca,Sc,Ti,Cr,Mn,Co,Ni,Ge,Sr,Y,Zr,Pd,Ag,In,Sn,Sb,La,Pt,Au,Ce,Nd,Gd,Dy
#df246522-37de-4d18-a671-d2effa4d715c etaP_I cell param
#34c8c1d5-3173-4503-8239-08a14d2f6393 vc-relax
../aiida_create_antisite_structures.py -se $ELEMENTS  -sg EtaPIV-antisite-many-vcrelax_structures --input_structures 34c8c1d5-3173-4503-8239-08a14d2f6393 --target_supercellsize 108


