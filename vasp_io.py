'''
Created on 2014.  05.  14.
 
@author: nwan
@author: hanwooh
'''
import numpy as np 
import re

def findReg(reg,file):
    f = open(file)
    lines = f.readlines()
    f.close()
    find = []
    for line in lines:
        m = re.findall(reg, line)
        find += m
    return find

def readEIGENVAL(fileName='EIGENVAL'):
    # read EIGENVAL
    f=open(fileName)
    buffer=f.readlines()
    [nKpt,nBand]=[int(i) for i in buffer[5].split()][1:]
    print [nBand,nKpt]
    bandInfo=[]
    #print 'NOW READING ENERGY PART OF EIGENVAL FILE'
    for j in range(nKpt):
        for i in range(nBand):
            eigenval = buffer[i + 8 + (nBand+2)*j].split()
            eigenval = float(eigenval[1])
            bandInfo.append([i+1,j+1,eigenval])    
    f.close()
    # [bandNum, kptNum , eigenval]
    return [nKpt,nBand,bandInfo]

def readPROCAR(fileName='PROCAR', orbital=-1):
    f = open(fileName)
    buffer = f.readlines()

    ''' # of k-points:   14         # of bands: 450         # of ions:  85 '''
    nKpt   = int(re.search('(?<=# of k-points:)\s+\d+',buffer[1]).group(0) )
    nBands = int(re.search('(?<=# of bands:)\s+\d+'   ,buffer[1]).group(0) )
    nIons  = int(re.search('(?<=# of ions:)\s+\d+'    ,buffer[1]).group(0) )
    
    nOrbits = len(buffer[7].split())-1
    Proj = np.zeros((nKpt,nBands,nIons,nOrbits))
    Kpts = np.zeros((nKpt,3))
    Eigs = np.zeros((nKpt,nBands))
    Occs = np.zeros((nKpt,nBands))

    kptInfoLength =1
    for i in [line.find(' k-point') for line in buffer[3+1:]]:
        kptInfoLength-=i
        if i==0:break


    for kpt in range(nKpt):
        # read k-th band projection to ion orbital 
        # read k-point
        kptLineNum = 2 + 1 + kptInfoLength* kpt
        kptLine = buffer[kptLineNum]

        kVec = re.search('(?<=:)\s*([-]?[0-9]*\.?[0-9]+)\s*([-]?[0-9]*\.?[0-9]+)\s*([-]?[0-9]*\.?[0-9]+)',kptLine)
        kVec = np.array( [float(kVec.group(1)), float(kVec.group(2)), float(kVec.group(3))] )
        Kpts[kpt,:] = kVec
        
        kp_weight = float( re.search('(?<=weight =)\s*[-]?[0-9]*\.?[0-9]+',kptLine).group(0) )

        bandInfoLength =1
        for i in [line.find('band') for line in buffer[kptLineNum+2+1:]]:
            bandInfoLength-=i
            if i==0:break
        # print bandInfoLength

        for band in range(nBands): 
            # bandLineNum = kptLineNum + 2 + band * (nIons *3 +6) #works for only non-S*L coupling
            bandLineNum = kptLineNum + 2 + band *bandInfoLength
            eig = float( re.search('(?<=energy)\s+[-]?[0-9]*\.?[0-9]+',buffer[bandLineNum]).group(0) )
            occ = float( re.search('(?<=occ.)\s+[-]?[0-9]*\.?[0-9]+',buffer[bandLineNum]).group(0) )
            
            Eigs[kpt,band] = eig
            Occs[kpt,band] = occ

            for ion in range(nIons):
                ionLineNum = bandLineNum +3 + ion

                orbital_proj = [float(o) for o in buffer[ionLineNum].split()[1:]]
                Proj[kpt,band,ion,:] = orbital_proj

    return nKpt, nBands, nIons, Kpts, Eigs, Proj, Occs

def readCONTCAR(fileName='CONTCAR', rtspecies = False):
    latticeVecs=[]
    atomSet=[]
    atomSetDirect=[]
    sd=False
    f=open(fileName,'r')
    # read first & second line 
    f.readline()
    latConst=float(f.readline())
    # read lattice vectors
    latVec=np.array([float(i)*latConst for i in f.readline().split()])
    latticeVecs.append(latVec)
    latVec=np.array([float(i)*latConst for i in f.readline().split()])
    latticeVecs.append(latVec)
    latVec=np.array([float(i)*latConst for i in f.readline().split()])
    latticeVecs.append(latVec)
    
    # read species
    species=f.readline().split()
    numSpecies=[int(i) for i in f.readline().split()]
    
    line = f.readline().strip()

    if line=='Selective dynamics':        
        DorC=f.readline()
    else:
        DorC=line
    
    # read coordinate 
    k=0
    for symbol in species:
        for n in range(numSpecies[k]):
            coord=np.array([float(i) for i in f.readline().split()[:3]])

            atomSetDirect.append([symbol,coord])
            if DorC[0]=='D' or  DorC[0]=='d' : # Direct
                coord = latticeVecs[0]*coord[0]+latticeVecs[1]*coord[1]+latticeVecs[2]*coord[2]
            else: 
                print "check coord! it's not direct form"

            atomSet.append([symbol,coord])
        k+=1
    f.close()

    for i,latVec in enumerate(latticeVecs):
        latticeVecs[i]=latVec/latConst
    if rtspecies==True:
        return [latConst,latticeVecs,atomSetDirect,species]
    else:
        return [latConst,latticeVecs,atomSetDirect]


def readLOCPOT(finaName='LOCPOT'):
    state = 0
    value = []
    LOCPOT = []
    sc = []

    number = 0
    grid = []

    f = open(finaName)
    buffer = f.readlines()
    f.close()
    # length = len(buffer)
    a = float(buffer[1])
    #~ print a
    for i in xrange(2,5) :
        sc.append(buffer[i].split())

    scLatVecx = np.array([float(i) for i in sc[0]])
    scLatVecy = np.array([float(i) for i in sc[1]])
    scLatVecz = np.array([float(i) for i in sc[2]])
    #~ print [scLatVecx,scLatVecy,scLatVecz]

    atom = buffer[5+1].split()
    number = sum([int(i) for i in atom])
    #print number

    grid = buffer[8+1+number].split()
    gridx = int(grid[0])
    gridy = int(grid[1])
    gridz = int(grid[2])
    
    #~ print [gridx,gridy,gridz]

    lenthPerline = len( buffer[9+1+number].split())
    for i in xrange(9+1+number,9+1+number+gridx*gridy*gridz/lenthPerline+1) :
        temp = buffer[i].split()
        for j in xrange(len(temp)) :
                value.append(float(temp[j]))
    
    value = np.array(value) 
    value = value[:gridx*gridy*gridz]
    
    LOCPOT = value.reshape(gridz,gridy,gridx).T

    return  a,[scLatVecx,scLatVecy,scLatVecz],[gridx,gridy,gridz] , LOCPOT

def readLOCPOT_lowMemory(finaName='LOCPOT'):
    ''' probably you don't want to use this'''
    state = 0
    value = []
    LOCPOT = []
    sc = []

    number = 0
    grid = []

    f = open(finaName)
    # buffer = f.readlines()
    # length = len(buffer)

    '''read lattice constant'''
    f.readline() # system
    a = float(f.readline())

    ''' read supercell lattice vector '''
    for i in xrange(2,5) :
        sc.append(f.readline().split())

    scLatVecx = np.array([float(i) for i in sc[0]])
    scLatVecy = np.array([float(i) for i in sc[1]])
    scLatVecz = np.array([float(i) for i in sc[2]])
    #~ print [scLatVecx,scLatVecy,scLatVecz]

    '''   read atom species and number of atoms   '''
    atom = f.readline().split()
    noAtom = np.array([ int(a) for a in f.readline().split()])


    '''   read atomic positions   '''
    for i in range(noAtom.sum()+2):
        f.readline()

    
    '''   read grid points   '''
    gridx,gridy,gridz = [ int(a) for a in f.readline().split()]
    # print gridx,gridy,gridz

    value = []
    for line in f:
        value += [float(a) for a in line.split()]
        # value = np.append(value,np.array( [float(a) for a in line.split()]))
    value = np.array(value).reshape(gridz,gridy,gridx).T


def readEIGENVAL(fileName='EIGENVAL'):
    # read EIGENVAL
    f=open(fileName)
    buffer=f.readlines()
    [nKpt,nBand]=[int(i) for i in buffer[5].split()][1:]
    print [nBand,nKpt]
    bandInfo = []
    kpoints =[]
    eigenvals =np.zeros((nKpt,nBand))
    #print 'NOW READING ENERGY PART OF EIGENVAL FILE'
    for j in range(nKpt):
        kpoint =np.array(buffer[-1 + 8 + (nBand+2)*j].split())[:3]
        kpoint = np.array([float(k) for k in kpoint])
        kpoints.append(kpoint)

        for i in range(nBand):
            eigenval = buffer[i + 8 + (nBand+2)*j].split()
            eigenval = float(eigenval[1])
            eigenvals[j,i] = eigenval
            #bandInfo.append([i+1,j+1,eigenval])    
    f.close()
    
    return kpoints,eigenvals

def readDOSCAR(fileName,atomNum):
    f=open(fileName)
    data=f.readlines()
    f.close()
    numAtom=data[0]
    numAtom=numAtom.split()
    numAtom=int(numAtom[0])
    #print('numAtom',numAtom)
    '''
    read PDOS
    '''
    '''
    go to data of atom selected
    '''
    eSet=[]
    sDOSSet=[]
    pDOSSet=[]
    dDOSSet=[]

    head=data[5].split()
    numRow=int(head[2])
    fermiE=float(head[3])

    eSet=[]
    tDOSSet=[]
    iDOSSet=[]
    '''read total dos'''
    for i in range(6,6+numRow):
        row=data[i].split()
        e=float(row[0])-fermiE
        tDOS=float(row[1])
        iDOS=float(row[2])

        eSet.append(e)
        tDOSSet.append(tDOS)
        iDOSSet.append(iDOS)
    tDOSSet = np.array(tDOSSet)


    ''' read average pdos or specific atom pdos '''
    if atomNum>0:
        atomSet = [atomNum]
    else:
        atomSet = range(numAtom)


    ''' read pdos  '''
    sDOSSetSum = np.zeros(numRow)
    pDOSSetSum = np.zeros(numRow)
    dDOSSetSum = np.zeros(numRow)
    for atomNum_i in range(numAtom):
        # eSet=[]
        sDOSSet=[]
        pDOSSet=[]
        dDOSSet=[]  
        sRow=5+(atomNum_i+1)*(numRow+1)
        head=data[sRow].split()
        
        
        #print('numRow',numRow)

        for i in range(sRow+1,sRow+1+numRow):
            row=data[i].split()
            # e=float(row[0])-fermiE
            sDOS=float(row[1])
            pDOS=float(row[2])
            dDOS=float(row[3])
            # eSet.append(e)
            sDOSSet.append(sDOS)
            pDOSSet.append(pDOS)
            dDOSSet.append(dDOS)
        sDOSSetSum += np.array(sDOSSet)
        pDOSSetSum += np.array(pDOSSet)
        dDOSSetSum += np.array(dDOSSet)
    # print numRow
    sDOSSet = sDOSSetSum
    pDOSSet = pDOSSetSum
    dDOSSet = dDOSSetSum
   
    eSet   =np.array(eSet   )
    return [eSet,tDOSSet,sDOSSet,pDOSSet,dDOSSet]


def writePOSCAR(output,latConst,latticeVecs,atomSetDirect,lSelective=False,lDirect=True):
    species = [atom[0] for atom in atomSetDirect]
    species1 = list(set(species))
    species1.sort(key=species.index)
    species=species1

    f = open(output,'w')
    f.write('system\n')
    f.write(str(latConst)+'\n')
    

    for latticeVec in latticeVecs:
        for i in range(3):
            f.write(str(latticeVec[i]) + ' ')
        f.write('\n')
    
    for s in species:
        f.write(s + ' ')    
    f.write('\n')

    for s in species:
        n = len([atom for atom in atomSetDirect if atom[0]==s])
        f.write(str(n)+' ' )
    f.write('\n')

    if lSelective: f.write('Selective dynamics\n')
    if lDirect: f.write('Direct\n')
    else : f.write('Cartesian\n')
    for atom in atomSetDirect:
        f.write(str(atom[1][0])+' '+str(atom[1][1])+' '+ str(atom[1][2])+' ' )
        if lSelective:
            # print lSelective
            f.write(str(atom[2][0])+' '+str(atom[2][1])+' '+ str(atom[2][2])+' ' )
        f.write('\n')

    f.close()


def getNELECT(OUTCAR):
    '   NELECT =     338.0000    total number of electrons'
    nelect = findReg('(?<=NELECT =)\s+\d+', OUTCAR)
    return int(nelect[0])

def getVBM(EIGENVAL,band_no = 0):
    if band_no == 0:
        band_no = getNELECT('OUTCAR')/2
    eigenval = findReg('(?<=^'+ str(band_no).rjust(4) +')\s+[-]?[0-9]*\.?[0-9]+',EIGENVAL)
    vbm = [float(e) for e in eigenval]
    vbm = max(vbm)
    return vbm

def getCBM(EIGENVAL,band_no = 0):
    if band_no == 0:
        band_no = getNELECT('OUTCAR') / 2 + 1
    eigenval = findReg('(?<=^'+ str(band_no).rjust(4) +')\s+[-]?[0-9]*\.?[0-9]+',EIGENVAL)
    cbm = [float(e) for e in eigenval]
    cbm = min(cbm)
    return cbm

def getEgap(OUTCAR, EIGENVAL):
    n_elect = getNELECT(OUTCAR)
    vbm = getVBM(EIGENVAL, n_elect/2)
    cbm = getCBM(EIGENVAL, n_elect/2+1)
    egap = cbm - vbm
    return egap

def getTotE(OSZICAR):
    ' 1 F= -.54010511E+03 E0= -.54010511E+03  d E =-.786584E-14'
    energy = findReg('(?<=F=)\s+[-+]?[0-9]*\.?[0-9]+[eE][-+]?[0-9]+? ', OSZICAR)
    totE = [float(e) for e in energy]
    totE = totE[len(totE)-1]
    return totE

def get_tot_E_outcar(outcar, enthalpy=True):
    # get total energy from line of outcar
    '  free  energy   TOTEN  =       -53.472728 eV'
    '''
    FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
    ---------------------------------------------------
    free  energy   TOTEN  =      -178.22800012 eV

    energy  without entropy=     -178.22769568  energy(sigma->0) =     -178.22789864
    enthalpy is  TOTEN    =      -153.45223930 eV   P V=       24.77576082
    '''
    if enthalpy:
        reg = '(?<=enthalpy is  TOTEN)\s*=\s+[-+]?[0-9]*\.?[0-9]+'
    else:
        reg = '(?<=free  energy   TOTEN)\s*=\s+[-+]?[0-9]*\.?[0-9]+'
    find = []
    for line in outcar:
        m = re.findall(reg, line)
        find += m
    tot_E = float(find[-1].replace('=',''))
    return tot_E

def readOUTCAR(OUTCAR='OUTCAR'):
    f = open(OUTCAR, 'r')
    outcar = f.readlines()
    f.close()
    return outcar

def writeOUTCAR(outcar, output_file='OUTCAR'):
    f = open(output_file, 'w')
    f.writelines(outcar)
    f.close()

def get_enthalpy(dir_name):
    outcar = readOUTCAR('{}/OUTCAR' % dir_name)
    # get total energy from line of outcar
    '  free  energy   TOTEN  =       -53.472728 eV'
    reg = '(?<=TOTEN  =)\s+[-+]?[0-9]*\.?[0-9]+'
    find = []
    for line in outcar:
        m = re.findall(reg, line)
        find += m
    tot_E = float(find[-1])
    return tot_E

if __name__ == '__main__':
    outcar = readOUTCAR('./test_B_run/00000/OUTCAR')
    print get_tot_E_outcar(outcar)