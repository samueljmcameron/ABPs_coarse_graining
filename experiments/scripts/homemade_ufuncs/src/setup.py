def configuration(parent_package='', top_path=None):

    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('.',
                           parent_package,
                           top_path)

    config.add_extension('nufuncs',['nufuncsmodule.c',
                                    'c_funcs.c',
                                    's_and_q_funcs/s1.c',
                                    's_and_q_funcs/s2.c',
                                    's_and_q_funcs/q1.c',
                                    's_and_q_funcs/q2.c',
                                    's_and_q_funcs/chbevl.c'])


    return config

if __name__ == "__main__":

    from numpy.distutils.core import setup

    setup(configuration=configuration,
          include_dirs = ['s_and_q_funcs/'])
