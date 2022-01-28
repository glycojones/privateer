#
#   PDB version 2 to version 3 converter. This routine also will correct DNA/RNA names
#
def usage():
    print "%s -i input -o output -a altnames " %sys.argv[0]
    print "%s --input input --output output --altnames altnames " %sys.argv[0]
    print 'input    - input pdb file: It must be defied'
    print 'output   - output pdb file: If it is not defined then it is created using input file name'
    print 'altnames - file with alternative names: If it is not defined the it is taken from %s' % os.getenv("CLIBD_MON", "")

def substring_replace(s,t,i1,j1):
    #  Replace s[i:j] by t. Take care of lengths
    i,j=i1,j1
    k=j-i
    
    if(k<0) :
        i,j=j,i
        k=-k
    if len(t) < k :
        tt = t
        for i in range(len(t),k) :
            tt = tt + " "
    elif(len(t) > k):
        tt = t[:k]
    else:
        tt=t
    s=s[:i]+tt+s[j:]
    return s

def known_pdb_corrs(allines):
    n=0
    for aa in allines:
        resname = aa[17:20].strip()
        if (aa[0:4]=="ATOM" or aa[0:4]=="HETA") :
            if(resname=="Ad") :
                resname = "DA"
            elif (resname == "Cd") :
                resname = "DC"
            elif (resname == "Td") :
                resname = "DT"
            elif (resname == "Gd") :
                resname = "DG"
            elif (resname == "Ud") :
                resname = "DU"
            elif (resname == "Ar") :
                resname == "A"
            elif (resname == "Cr") :
                resname == "C"
            elif (resname == "Tr") :
                resname = "T"
            elif (resname == "Ur") :
                resname = "U"
                allines[n] = substring_replace(aa,resname,16,19)
            n=n+1
        return allines
    #
    #    There might some more to fiddle around. They need to be put here are put in another function

def known_dnarna_corrs(alllines):
    n = 0
    for aa in allines:
        resname = aa[17:20].strip()
        if (aa[0:4]=="ATOM" or aa[0:4]=="HETA") :
            if(resname in ) :
                if( check_dnarna() == 0) :
                    substring_replace(aa)
                
def process_alts(residue,alts):
    altsthis = []
    resname = residue[1][17:20]
    print "resname %s" %resname
    n = 0
    for bb in alts:
        if(bb[0] == resname):
            n = n +1
            altsthis.append(bb)

    if(n <1) : return residue
    
    natall = 0
    nalt   = 0
    nat    = 0
    atreal1 = []
    atreal2 = []
    for ats in residue:
        natall = natall+1
        atom = ats[12:16].strip()
        print atom
        for ats_alt in altsthis:
            if atom==ats_alt[2] :
                nalt=nalt+1
                atreal1.append([ats_alt[1],ats_alt[3]])
            elif atom==ats_alt[1] :
                atreal2.append([ats_alt[1],ats_alt[3]])
                nat = nat + 1
    if (nalt == natall) :
        #
        # replace atoms
        nat = 0
        for ats in residue:
            nat = nat + 1
            atom = ats[12:4].strp()
            for ats_alt in altsthis:
                if(atom==ats_alt[2]) :
                    #
                    #  Take care: atom names may need to be written
                    #  in a fashion pdb wants
                    atom = atreal1[nat][0].strip()
                    element = atreal1[nat][1].strip()
                    if (len(element) == 1) :
                        if(len(atom) < 4): atom = " "+atom
                        element = " "+element
                        ats = substring_replace(ats,atom,12,16)
                        l=len(ats)
                        residue[nat-1] = substring_replace(ats,element,76,l)
                        break
    elif(nat==natall):
        print nat
    elif(nat !=natall) :
        print nat, natall
        print "we have a problem"
        print "At least some of the atoms in residues %s have neither alt no real atom correspondence" %resnum
        sys.exit()
    return residue

#
#   Main body of the program
#
import getopt, sys, os

#
#   Read and interpret command line options
try:
    opts, args = getopt.getopt(sys.argv[1:],"i:a:o:h",["input=","altnames=","output=","help"])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
#
pdbout   = " "
altnames = " "
pdbin    = " "
for o,a in opts:
    if o in ("-i","--input"):
        pdbin=a
    elif o in ("-a", "--altnames"):
        altnames = a
    elif o in ("-o","--output"):
        pdbout = a
    elif o in ("-h","--help"):
        usage()
        sys.exit()

#
#  Make sure that necessary files exist
if(pdbin == "" or not os.path.isfile(pdbin)):
    print "Input pdb file must be defined and it must exist"
    usage()
    sys.exit()

if(pdbout == " " or len(pdbout) < 0):
    l = pdbin.find(".pdb")
    if( l > 0) :
        pdbout=pdbin[0:l]+"_v3.pdb"
    else:
        pdbout = "v3.pdb"
if altnames==" " :
   altnames = os.getenv("CLIBD_MON", "")+"pdb_alt_names.txt"

if( not os.path.isfile(altnames)) :
    print "File with altnames - %s does not exist" %altnames

#
#   Now ready to do the job

# Read the list of monomers with their alts

altatoms = []
i=0
j=0
with open(altnames,"r") as alts:
    for aaa in alts:
        aaa=aaa.strip()
        bbb =aaa.split()
        if(bbb[0] == "resname") :
            #
            #  New residue
            bsave = bbb[1]
        elif(bbb[0] == "atom") :
            #
            # read the residue further
            altatoms.append([bsave,bbb[1],bbb[3],bbb[5]])
#
#  Now read pdb file and organise it. Remove hydrogens and aniso cards

dictPDB = ["REMARK","SCALE1","CRYST1","MTRX1","MTRX2","MTRX3","ATOM","HETATM","LINK","MODRES","ORIGX","HEADER"]
allines=[]
with open(pdbin,"r") as pdbin:
    for line1 in pdbin:
        line = line1.rstrip()
        element = " "
        if  (line[0:4] == 'ATOM' or line[0:4] == 'HETA') and len(line) > 77 :
           element = line[76:78].strip()
        for aa in dictPDB:
            if line[0:6].rstrip() == aa and element != 'H' :
                #
                #  Remove charges from atom names. They should be on the element nam
                atom = line[12:16]
                ip=atom.find("+")
                if(ip < 0) : ip=atom.find("-")
                if(ip > 0) : atom = atom[0:ip]
                l=len(atom)
                if atom[l-1:l] == "," : atom = substring_replace(atom,"'",l-1,l)
                line = substring_replace(line,atom,12,16)
                allines.append(line)

allines = known_pdb_corrs(allines)


#
#   Correct known potential problems


#
#   Start processing this pdb
resnumOld ="Something Unreasonable"
residue = []

#out = open(pdbout,"w")
toout = []
nres = 0
for aa in allines:
    if aa[0:4]=="ATOM" or aa[0:6]=="HETATM" :
        resnum = aa[21:27]
        resname = aa[17:20]
        if (resnum == resnumOld) :
            nres = nres + 1
            residue.append(aa)
        else:
            if resnumOld != "Something Unreasonable" :
                #
                #   Process this residue
                residue1 = process_alts(residue,altatoms)

                for line in residue1:
                    toout.append(line)
            nres = 1
            residue = [aa]
            resnumOld = resnum
    else:
        if nres > 0:
            residue1 = process_alts(residue,altatoms)
            for line in reside1:
                print line
                toout.append(line)            
            nres = 0
        toout.append(aa)
print nres
if nres > 0:
    residue1 = process_alts(residue,altatoms)
    for line in residue1:
        print line
        toout.append(line)            
out = open(pdbout,"w")
for aa in toout:
    out.write(aa+"\n")

out.close()
