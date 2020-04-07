
import sys

def splitfile(path,fp,rho):


    filename = path + f'g_{fp}_{rho}.rdf'

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
                newf = path + f'g_{fp}_{rho}_{count}.rdf'
                count += 1
                lines = [next(f) for i in range(2000)]
                header = '#' + a[0] + '\n'
                lines = [header] + lines
                with open(newf,'w') as newfile:
                    newfile.writelines(lines)


            
    return

if __name__ == "__main__":

    path = 'correlations/'

    rhos = ['0.05','0.1','0.2','0.4']

    fps = [0,1,5,10,20,40,60,80,100]

    for rho in rhos:

        for fp in fps:

            splitfile(path,fp,rho)
