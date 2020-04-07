import sys, os

from numpy.distutils.exec_command import exec_command




def compile_fortran(fname, module_name, extra_flags=''):

    folder = os.path.dirname(os.path.realpath(__file__))


    args = (f" -c -m {module_name} {fname} "
            f"--f90flags={extra_flags}")
    command = (f'cd "{folder}" && "{sys.executable}" -c '
               f'"import numpy.f2py as f2py;f2py.main()" {args}')
    print(command)
    status, output = exec_command(command)
    return status, output, command


if __name__=="__main__":

    module_name = 'fortrantools'
    fname = 'fortrantools.f95'

    compile_fortran(fname,module_name,
                    extra_flags='-fopenmp -lgomp')


