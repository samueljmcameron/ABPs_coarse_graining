
import sys

def splitfile(inpath,outpath,fp,rho,binnum=2000):


    filename = inpath + f'g_{fp}_{rho}.rdf'

    with open(filename,'r') as f:
        count = 0
        while True:
            line = f.readline()
            try:
                if line[0] == "#":
                    continue
            except:
                print('end of file reached!')
                break
            if line.count(' ') == 1:
                a = line.split(' ')
                newf = outpath + f'g_{fp}_{rho}_{count}.rdf'
                count += 1
                lines = [next(f) for i in range(binnum)]
                header = '#' + a[0] + '\n'
                lines = [header] + lines
                with open(newf,'w') as newfile:
                    newfile.writelines(lines)


            
    return

if __name__ == "__main__":

    inpath = '../2020_03_31/interactions_pressure/data/'
    outpath = '../2020_03_31/interactions_pressure/processed_correlations/'

    rhos = ['0.05','0.1','0.2','0.3','0.4','0.5','0.6','0.7']

    fps = [0,1,5,10,20,40,80]

    for rho in rhos:

        for fp in fps:

            splitfile(inpath,outpath,fp,rho)
