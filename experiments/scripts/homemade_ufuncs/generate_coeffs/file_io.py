import sys
import os

def write_to_coeffs(func_name,func_def,coeff_arrays,
                    lower_bounds,upper_bounds):
    
    labels = ['A','B','C','D','E','F','G']

    N = len(coeff_arrays)
    if N > len(labels):
        print('error in writing to file coeffs! too many different arrays')
        return -1

    
    
    with open(f'../{func_name}_coeffs.tmp','w') as f:
    
        f.write('\n\n/* Starting with')
        f.write(f' *\t {func_name}(x) = {func_def}')
        f.write(' * function. */\n\n')

        
        for i in range(N):
            lab = labels[i]
            low = lower_bounds[i]
            upp = upper_bounds[i]

            coeff = coeff_arrays[i]
            f.write(f'/* coeffs for {func_name}(x) with {low} < x < {upp} */\n')
            f.write(f'static double {func_name}_{lab}[] = \n')
            f.write('  {\n')
            
            for c in coeff:
                
                f.write(f'    {c:.15e},\n')
                
            f.write('  };\n')


