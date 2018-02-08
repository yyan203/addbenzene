from math import exp, log10
import sys, argparse, pdb
import random


class Atom:

    def __init__(self, id, type, x, y, z):
        assert isinstance(x, float)
        self.xyz_ = [x, y, z]
        self.id_ = id
        self.type_ = type


class allatoms:

    def __init__(self):
        self.atoms = {}  # store atoms information, use atom ID as key
        self.bonds = {}  # store bonds information, use atom ID as key
        self.atomnum = 0
        self.bondnum = 0
        self.maxID = -1

    # read single frame xyz format: ID elementname X Y Z"
    # read atom info,  first column is also line number
    # 1 Zn 25.9831 14.4739 14.674
    # 2 Zn 24.08 20.3523 21.7369
    # 3 Zn 26.5129 15.4012 19.8225
    def read_xyz(self, filename):
        with open(filename) as f:
            for index, line in enumerate(f):
                if index == 1:
                    self.atomnum = int(line)
                if index == 5:
                    t = line.split(' ')
                    self.atoms[int(t[0])] = Atom(int(t[0]), str(t[1]), float(t[2]), float(t[3]), float(t[4]))
                    self.maxID = max(self.maxID, int(t[0]))
        print("XYZ file reading:", self.atomnum, " atoms!")

    # read single frame bonds information: type atom1-ID  atom2-ID
    # read bond info
    #  1  1  6660
    #  3  6655  8703
    def read_bonds(self, filename):
        with open(filename) as f:
            for index, line in enumerate(f):
                if index == 3:
                    t = line.split(' ')
                    if t[1] not in self.bonds:
                        self.bonds[t[1]] = set([t[2]])
                    else:
                        self.bonds[t[1]].add(t[2])
                    if t[2] not in self.bonds:
                        self.bonds[t[2]] = set([t[1]])
                    else:
                        self.bonds[t[2]].add(t[1])
                    self.bondnum += 1
        print("Read:", self.bondnum, " bonds!")


    #  after delete atoms and add new atoms, the maximum ID of atoms might larger than self.atomnum
    #  because some ID has no atoms associated
    #  assuming there are only  Zn  H  C  N elements

    def outputxyz(self, xyzfile):
        Zn, H, C, N = {}, {}, {}, {}
        #nZn, nH, nC, nN = 1, 1, 1, 1
        j = 1
        for i in self.atoms:
            if self.atoms[i].type_ is "Zn":
                Zn[j] = self.atom[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ is "H":
                H[j] = self.atom[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ is "C":
                C[j] = self.atom[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ is "N":
                N[j] = self.atom[i].id_; j += 1

        f = open(xyzfile, 'w')
        f.write('%d\n' % self.atomnum)
        f.write('add benzine to Zif4\n')
        for i in Zn:
            xyz = self.atoms[Zn[i]].xyz_,
            f.write("%d %s %f %f %f\n" % (i, self.atoms[Zn[i]].type_, xyz[0], xyz[1],xyz[1]))
        for i in H:
            xyz = self.atoms[H[i]].xyz_,
            f.write("%d %s %f %f %f\n" % (i, self.atoms[H[i]].type_, xyz[0], xyz[1],xyz[1]))
        for i in C:
            xyz = self.atoms[C[i]].xyz_,
            f.write("%d %s %f %f %f\n" % (i, self.atoms[C[i]].type_, xyz[0], xyz[1],xyz[1]))
        for i in N:
            xyz = self.atoms[C[i]].xyz_,
            f.write("%d %s %f %f %f\n" % (i, self.atoms[C[i]].type_, xyz[0], xyz[1],xyz[1]))
        f.close()

    def outputbond(self, bondfile):
        oldID2newID = {}
        newID2oldID = {}
        j = 1
        newbonds = {}
        bondtype = {"Zn-N": 1, "N-Zn": 1, "H-C" : 2, "C-H" : 2, "C-N" : 3, "N-C" : 3, "C-C" : 4}
        for i in self.atoms:
            if self.atoms[i].type_ is "Zn":
                oldID2newID[self.atom[i].id_] = j; j += 1
                newID2oldID[j] = self.atom[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ is "H":
                oldID2newID[self.atom[i].id_] = j; j += 1
                newID2oldID[j] = self.atom[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ is "C":
                oldID2newID[self.atom[i].id_] = j; j += 1
                newID2oldID[j] = self.atom[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ is "N":
                oldID2newID[self.atom[i].id_] = j; j += 1
                newID2oldID[j] = self.atom[i].id_; j += 1
        j = 1
        while j <= self.maxID:
            if j in self.bonds:
                newJ = oldID2newID[j]
                neighbour = self.bonds[j]
                for k in neighbour:
                    newK = oldID2newID[k]
                    assert newJ is not newK, "bonds error, same atoms!"
                    if newJ < newK:
                        if newJ in newbonds:
                            newbonds[newJ].add(newK)
                        else:
                            newbonds[newJ] = set([newK])
                    else:
                        if newK in newbonds:
                            newbonds[newK].add(newJ)
                        else:
                            newbonds[newK] = set([newJ])
            j += 1

        f = open(bondfile, 'w')
        for i in newbonds:
            for j in newbonds[i]:
                type = bondtype[self.atoms[newID2oldID[i]] + "-" + self.atoms[newID2oldID[j]]]
                f.write("%d %d %d\n" % (type, i, j))
        f.close()


    # delete atoms and its associated bonds
    def delete_atom(self, atomID):
        assert atomID in self.atoms, "Atom %d does not exist!" % atomID
        temp = self.bonds.pop(atomID);
        self.bondnum -= len(temp)
        for i in self.bonds:
            if atomID in self.bonds[i]:
                self.bonds[i].pop(atomID)
        self.atoms.pop(atomID)
        self.atomnum -= 1

    def add_atom(self, type, coord):
        self.atomnum += 1
        newid = self.atomnum
        while newid in self.atoms:
            newid += 1
        self.atoms[newid] = Atom(newid, type, coord[0], coord[1], coord[2])
        self.maxID = max(self.maxID, newid)
        return newid

    def add_bond(self, id1, id2):
        self.bondnum += 1
        assert id2 not in self.bonds[id1], "Bond %d -> %d already exist!" % (id1, id2)

        if id1 not in self.bonds:
            self.bonds[id1] = set([id2])
        else:
            self.bonds[id1].add(id2)
        if id2 not in self.bonds:
            self.bonds[id2] = set([id1])
        else:
            self.bonds[id2].add(id1)



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("NUM", type=int, help="Number of benzine to add")
    parser.add_argument("XYZ", help="XYZ_file(format: ID elementname X Y Z")
    parser.add_argument("BOND", help="Original bond list")
    parser.add_argument("OUT_bond", help="Output new bond list with benzine")
    parser.add_argument("OUT_xyz",  help="Output new xyz  file with benzine")
    args = parser.parse_args()
    print(args.XYZ)

    ben = args.NUM
    mysystem = allatoms()
    mysystem.read_xyz(args.XYZ)
    mysystem.read_bonds(args.BOND)

    # add benzine
    # remember visited N atoms
    Nlist = set()
    for i in range(1, args.NUM + 1):
        # find first N
        flag = 1
        random = None
        while flag is not 0:
            random = random.randint(6145, 8192)
            if random in Nlist: continue
            else:
                Nlist.add(random)
                print("choose N ID:", random)
                flag = 0
        neighbor = mysystem.bonds[random]

        # get the C and N in N-C-N bonds
        Cindex, Nindex = None, None
        for j in neighbor:
            if mysystem.atoms[j].type_ is "C":
                for k in mysystem.bonds[j]:
                    if mysystem.atoms[k].type_ is "N" and k is not random:
                        Cindex, Nindex = j, k
                        Nlist.add(Nindex); break

        # get the two H in H-C-C-H in imidazole
        # get the two C in H-C-C-H in imidazole

        Ca, Cb, Ha, Hb = None, None, None, None
        for j in neighbor:
            if mysystem.atoms[j].type_ is "C":
                for k in mysystem.bonds[j]:
                    if mysystem.atoms[k].type_ is "N" and k is random:
                        Ca = j; break
        for j in mysystem.bonds[Ca]:
            if mysystem.atoms[j].type_ is "C":
                Cb = j; break

        for j in mysystem.bonds[Ca]:
            if mysystem.atoms[j].type_ is "H":
                Ha = j

        for j in mysystem.bonds[Cb]:
            if mysystem.atoms[j].type_ is "H":
                Hb = j

        print("N:", Nindex, " H1:", Ha, " H2:", Hb)

        # delete the two H atoms
        mysystem.delete_atom(Ha)
        mysystem.delete_atom(Hb)

        # add 4 C  and  4  H  as shown below (nC#, nH#)
        #
        #   C-C  1.39     C-N  1.34      C-H 1.02     N-Zn 1.97
        #
        #        nH1         nH2
        #           \       /
        #           nC1-——nC2
        #          /         \
        #   nH3--nC3    O    nC4—-nH4     #  O is the center
        #          \         /
        #           C2 ——— C1
        #           /       \
        #  Zn ---  N2         N1 --- Zn
        #           \       /
        #             \   /
        #               C

        N1, N2, C1, C2 = mysystem.atoms[random], mysystem.atoms[Nindex], mysystem.atoms[Ca], mysystem.atoms[Cb]
        newC1, newC2, newC3, newC4 = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
        newH1, newH2, newH3, newH4 = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]

        O = [0.0, 0.0, 0.0]
        for l in range(3):
            newC1[l] = N1[l] + 3.1 / 1.34 * (C1[l] - N1[l])
            newH1[l] = N1[l] + 4.1 / 1.34 * (C1[l] - N1[l])
            newC2[l] = N2[l] + 3.1 / 1.34 * (C2[l] - N2[l])
            newH2[l] = N2[l] + 4.1 / 1.34 * (C2[l] - N2[l])
            O[l] = N1[l] + 2.73 / 1.34 * (C1[l] - N1[l])
            newC3[l] = O[l] + C2[l] - C1[l]
            newH3[l] = O[l] + 2.41 * (C2[l] - C1[l])
            newC4[l] = O[l] + C1[l] - C2[l]
            newH4[l] = O[l] + 2.41 * (C1[l] - C2[l])
        C1id, C2id = mysystem.add_atom("C", newC1), mysystem.add_atom("C", newC2)
        C3id, C4id = mysystem.add_atom("C", newC3), mysystem.add_atom("C", newC4)
        H1id, H2id = mysystem.add_atom("C", newH1), mysystem.add_atom("C", newH2)
        H3id, H4id = mysystem.add_atom("H", newH3), mysystem.add_atom("H", newH4)
        mysystem.add_bond(C1id, H1id)
        mysystem.add_bond(C2id, H2id)
        mysystem.add_bond(C3id, H3id)
        mysystem.add_bond(C4id, H4id)
        mysystem.add_bond(C1id, C2id)
        mysystem.add_bond(C1id, C3id)
        mysystem.add_bond(C2id, C4id)
        mysystem.add_bond(C3id, Cb)
        mysystem.add_bond(C4id, Ca)

    # output XYZ and BONDS information
    mysystem.outputxyz(args.XYZ)
    mysystem.outputbond(args.BOND)

    print("Done")
