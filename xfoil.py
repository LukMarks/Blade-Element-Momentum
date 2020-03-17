import os
import subprocess as sp
path  ='.\\xfoil_output'

def xfoil(Airfoil,alpha,reynolds,mach,inter,np,path,L,D):
# Files to be deleted after each run....
    Coef = list()
    #os.system('del dump')
    #os.remove("xfoil_input.txt")
    #os.system('rm xfoil_input')
    #os.system('rm screen')
    os.remove("xfoil_output.txt")
    #os.system('rm xfoil_output')

    itens  = [Airfoil, alpha, alpha, '0', mach,reynolds,inter,np ]# Key names to look in the setup file

    nome = '\n'+itens[0]
    rey = '\n'+str(itens[5])
    alfa = '\n'+str(itens[1])
    inter = '\n'+str(itens[6])
    mach = '\n'+str(itens[4])
    np = '\n'+str(itens[7])
    path = 'load'+nome+'\npane\nppar\nn'+np+'\n\n\noper\nvisc'+rey+'\nmach'+mach+'\niter'+inter+'\npacc\nxfoil_output.txt\n\naseq'+alfa+alfa+'\n1\n\nquit'
 
    '''
    f=open('xfoil_input.txt','w')
    f.write('%s \n' % ('load'))
    f.write('%s \n' % (itens[0]))
    f.write('%s \n' % ('pane'))
    f.write('%s \n' % ('ppar'))
    f.write('%s \n' % ('n'))
    f.write('%s \n' % (itens[7]))
    f.write('%s \n' % ('   '))
    f.write('%s \n' % ('   '))
    f.write('%s \n' % ('oper'))
    f.write('%s \n' % ('visc'))
    f.write('%s \n' % (itens[5]))
    f.write('%s \n' % ('mach'))
    f.write('%s \n' % (itens[4]))
    f.write('%s \n' % ('iter'))
    f.write('%s \n' % (itens[6]))
    f.write('%s \n' % ('pacc'))
    f.write('%s \n' % ('xfoil_output'))
    f.write('%s \n' % ('   '))
    f.write('%s \n' % ('aseq'))
    f.write('%s \n' % (itens[1]))
    f.write('%s \n' % (itens[2]))
    f.write('%s \n' % (itens[3]))
    f.write('%s \n' % ('  '))
    f.write('%s \n' % ('quit'))
    f.close()
    '''

    n = path
    p = sp.Popen(['xfoil'],
		 stdin=sp.PIPE,
		 stdout=sp.PIPE,
		 stderr=sp.STDOUT)

    grep_stdout = p.communicate(input=n.encode())[0]


    i    = 0
    aoa  = list()
    cl   = list()
    cd   = list()
    cdp  = list()
    cm   = list()
    xutr = list()
    xltr = list()

    #if os.path.isfile(path):

    f      = open('xfoil_output.txt','r')

    for line in f:
        if (i > 11):

           aoa.append(line.strip().split()[0])
           cl.append(line.strip().split()[1])
           cd.append(line.strip().split()[2])
           cdp.append(line.strip().split()[3])
           cm.append(line.strip().split()[4])
           xutr.append(line.strip().split()[5])
           xltr.append(line.strip().split()[6])
           #print((line.strip().split()[6]))
           res = float(line.strip().split()[2])
        else:
             res = 1.0
        i += 1

    f.close()

        #colocar maximo
    j = 0

    if len(cl)==0:
        Cl =0.0000001#float(L[len(cl)-1])
        Cd =0.0000001#float(D[len(cd)-1])
        #if len(L)==0:
         #   Cl =float(L[len(cl)-2])
          #  Cd =float(D[len(cd)-2])

        Coef.append(Cl)
        Coef.append(Cd)
        #print('não foi')
        j +=1
        Coef.append(j)
    else:
        Cl =float(cl[len(cl)-1])
        Cd =float(cd[len(cd)-1])
        Coef.append(Cl)
        Coef.append(Cd)
        j=j
        Coef.append(j)

    return Coef
