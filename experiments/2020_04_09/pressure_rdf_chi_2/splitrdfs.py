
import sys

sys.path.append('../../analysis_scripts')

from split_rdf_files import splitfile

if __name__ == "__main__":
    
    inpath = '../../2020_04_08/winkler_pressure/data/'
    outpath = '../../2020_04_08/winkler_pressure/rdf_splits/'

    rhos = ['0.05','0.1','0.2','0.3','0.4','0.5','0.6']

    fps = ['0','0.25','0.5','0.75','1','2','4','8']

    for rho in rhos:

        for fp in fps:

            splitfile(inpath,outpath,fp,rho)
