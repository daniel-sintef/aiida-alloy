import numpy as np

concentrations = []
elements = ['Al','Mg','Cu']
v = 0.00

i = 0
for x in np.linspace(0,1,20):
    for y in np.linspace(0,1,20):
        v_div = 3.
        conc_a = x
        conc_a = max(conc_a-v/v_div,0)
        conc_b = (1-x)*y
        conc_b = max(conc_b-v/v_div,0)
        conc_c = (1-x)*(1-y)
        conc_c = max(conc_c-v/v_div,0)
        if conc_a == 0 or conc_b == 0 or conc_c == 0:
            continue
        print("{:.5f},{:.5f},{:.5f},{:.5f}".format(conc_a, conc_b, conc_c, v))
        i += 1
