import os , rosetta , pyrosetta
from Bio.PDB import *
pyrosetta.init()

class HomologPrep:
    aaTable = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'K', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    def __init__(self, pdbID , cwd):
        self.pdbID =  pdbID
        self.cwd = cwd
    def pdbFetch(self):
        """Fetch PDB file from Homolog Input"""
        fetcher = PDBList()
        pdbFile = fetcher.retrieve_pdb_file(self.pdbID, False, self.cwd)
        preFix = ('pdb' + self.pdbID + '.ent')
        os.rename(preFix, (self.pdbID + '.pdb'))
        self.pdbFile = self.pdbID + '.pdb'
        print 'Fetch PDB and filename fix in directory : Complete'
    def hlPose(self):
        """Create Pose Object for Homolog"""
        self.hlPose = pyrosetta.pose_from_file(self.pdbFile)
        print 'Pose Homolog : Complete'
    def hlInfo(self):
        """Gather Positional Information in PDBfile (resPos, dnaPos, h2oPos)"""
        self.hlInfo = []
        p = PDBParser()
        structure = p.get_structure(self.pdbID, self.pdbFile)
        tempResList = Selection.unfold_entities(structure, 'R') #gathers every residue in structure to list
        for r in tempResList:
            self.hlInfo.append(r.get_resname())
        print 'Gather Information about Homolog : Complete'
        HomologPrep.dnaSeqCheck(self)
        HomologPrep.h2oCheck(self)
        HomologPrep.resCheck(self)
        HomologPrep.dimerCheck(self)
        HomologPrep.hlRanges(self)
        print 'Homolog Amino Acid Sequence : ' , self.hlSeq
        print 'Bound DNA Sequence : ' , self.dnaSeq
    def dnaSeqCheck(self):
        """Identify DNA Base Positions in PDB File"""
        self.dnaPos = {} #positions and ID of the DNA in pdbFile
        dnaCheck = [ ' DC' , ' DT', ' DA', ' DG' ]
        position = 0
        for entry in self.hlInfo: #find dna in file
            if entry in dnaCheck:
                self.dnaPos[position] = entry[2]
            position += 1
        print 'Locate DNA Positions in File: Complete'
        dnaL = []
        for pos, base in self.dnaPos.items(): #assemble dnaSeq
            dnaL.append(base)
        dnaPal = ''.join(dnaL)
        self.dnaSeq = dnaPal[0 : len(dnaPal)/2]
    def h2oCheck(self):
        """Identify Water Positions in PDB File"""
        self.h2oPos = []
        position = 0
        for entry in self.hlInfo: #find h2o in file
            if entry == 'HOH':
                self.h2oPos.append(position)
            position += 1
        print 'Locate H2O Positions in File: Complete'
    def resCheck(self):
        """Identify Residue Positions in PDB file"""
        self.resPos = {}
        resL = []
        position = 0
        for res in self.hlInfo: #find res in file
            if position in self.dnaPos:
                position += 1
            elif position in self.h2oPos:
                position += 1
            else:
                self.resPos[position] = res
                position += 1
        print "Locate Residue Positions in File: Complete"
        for p,r in self.resPos.items():
            resL.append(HomologPrep.aaTable[r])
        self.hlSeq = ''.join(resL)
    def dimerCheck(self):
        """Check for Dimerization in PDB file"""
        a = self.hlSeq[0 : len(self.hlSeq)/2]
        b = self.hlSeq[(len(self.hlSeq)/2) : len(self.hlSeq)]
        if a == b:
            self.dimerCheck = 'y'
        else:
            self.dimerCheck = 'n'
        print "Check for Dimerization: Complete"
        print "Homolog is a Dimer (y/n): y"
    def hlRanges(self):
        """Identify Range of DNA , HOH , and Residue in Homolog and Store in self.rangeDict"""
        toRange = [self.dnaPos , self.h2oPos , self.resPos]
        self.rangeDict = {}
        for x in toRange:
            try:
                minP = min(x.items())[0]
                maxP = max(x.items())[0]
            except AttributeError:
                minP = min(x)
                maxP = max(x)
            rangeMake = range(minP , maxP+1)
            getRange = []
            clearL = []
            for pos in rangeMake:
                if pos in x:
                    clearL.append(pos)
                    if pos == (maxP):
                        getRange.append((min(clearL) , max(clearL)))
                else:
                    if len(clearL) == 0:
                        pass
                    else:
                        getRange.append((min(clearL) , max(clearL)))
                        clearL = []           
            if x is self.dnaPos:
                self.rangeDict['dnaRange'] = getRange
                print 'Identify DNA Range: Complete'
            elif x is self.h2oPos:
                self.rangeDict['h2oRange'] = getRange
                print 'Identify H2O Range: Complete'
            elif x is self.resPos:
                self.rangeDict['resRange'] = getRange
                print 'Identify Amino Acid Range: Complete'

            

###############
#    MAIN     #
###############
pdbID = '1a0a'
cwd = os.getcwd()
homolog = HomologPrep(pdbID , cwd)
HomologPrep.pdbFetch(homolog)
HomologPrep.hlPose(homolog)
HomologPrep.hlInfo(homolog)

