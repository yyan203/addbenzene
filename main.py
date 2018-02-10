from math import exp, log10
import sys, argparse, pdb
import random as rd


class Atom:

    def __init__(self, id, type, x, y, z):
        assert isinstance(x, float)
        self.id_ = id
        self.type_ = type
        self.xyz_ = [x, y, z]


class allatoms:

    def __init__(self):
        self.atoms = {}  # store atoms information, use atom ID as key
        self.bonds = {}  # store bonds information, use atom ID as key
        self.atomnum = 0
        self.bondnum = 0
        self.maxID = -1

    # read single frame xyz format: ID elementname X Y Z ( must contain head lines, i.e. # of atom and comment line"
    # read atom info,  first column is also line number
    # 1 Zn 25.9831 14.4739 14.674
    # 2 Zn 24.08 20.3523 21.7369
    # 3 Zn 26.5129 15.4012 19.8225
    def read_xyz(self, filename):
        with open(filename) as f:
            for index, line in enumerate(f):
                if index == 0:
                    self.atomnum = int(line)
                if index > 1:
                    t = line.split(' ')
                    #print(str(t[1]))
                    self.atoms[int(t[0])] = Atom(int(t[0]), str(t[1]), float(t[2]), float(t[3]), float(t[4]))
                    #print(self.atoms[int(t[0])].xyz_)
                    self.maxID = max(self.maxID, int(t[0]))
        #print(self.atoms[10].xyz_)
        #print(self.atoms[10].type_)
        print("XYZ file import:", self.atomnum, " atoms!")

    # read single frame bonds information: type atom1-ID  atom2-ID
    # read bond info
    #  1  1  6660
    #  3  6655  8703
    def read_bonds(self, filename):
        with open(filename) as f:
            for index, line in enumerate(f):
                t = line.rstrip().split()
                if int(t[1]) not in self.bonds:
                    self.bonds[int(t[1])] = set([int(t[2])])
                else:
                    self.bonds[int(t[1])].add(int(t[2]))
                if int(t[2]) not in self.bonds:
                    self.bonds[int(t[2])] = set([int(t[1])])
                else:
                    self.bonds[int(t[2])].add(int(t[1]))
                self.bondnum += 1
        print("Import ", self.bondnum, " bonds!")
        #print(self.bonds)


    #  after delete atoms and add new atoms, the maximum ID of atoms might larger than self.atomnum
    #  because some ID has no atoms associated
    #  assuming there are only  Zn  H  C  N elements

    def outputxyz(self, xyzfile):
        Zn, H, C, N = {}, {}, {}, {}
        #nZn, nH, nC, nN = 1, 1, 1, 1
        j = 1
        #print(self.atoms[10].type_,"yes here")
        for i in self.atoms:
            if self.atoms[i].type_ == "Zn":
                Zn[j] = self.atoms[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "H":
                H[j] = self.atoms[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "C":
                C[j] = self.atoms[i].id_; j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "N":
                N[j] = self.atoms[i].id_; j += 1

        f = open(xyzfile, 'w')
        f.write('%d\n' % self.atomnum)
        f.write('add benzine to Zif4\n')
        for i in sorted(Zn.keys()):
            xyz = self.atoms[Zn[i]].xyz_
            #print(Zn[i], xyz)
            f.write("%d %s %f %f %f\n" % (i, self.atoms[Zn[i]].type_, xyz[0], xyz[1], xyz[2]))
        for i in sorted(H.keys()):
            xyz = self.atoms[H[i]].xyz_
            #print(xyz)
            f.write("%d %s %f %f %f\n" % (i, self.atoms[H[i]].type_, xyz[0], xyz[1], xyz[2]))
        for i in sorted(C.keys()):
            xyz = self.atoms[C[i]].xyz_
            f.write("%d %s %f %f %f\n" % (i, self.atoms[C[i]].type_, xyz[0], xyz[1], xyz[2]))
        for i in sorted(N.keys()):
            xyz = self.atoms[N[i]].xyz_
            f.write("%d %s %f %f %f\n" % (i, self.atoms[N[i]].type_, xyz[0], xyz[1], xyz[2]))
        f.close()

    def outputbond(self, bondfile):
        oldID2newID = {}
        newID2oldID = {}
        j = 1
        newbonds = {}
        bondtype = {"Zn-N": 1, "N-Zn": 1, "H-C": 2, "C-H": 2, "C-N": 3, "N-C": 3, "C-C": 4}
        for i in self.atoms:
            if self.atoms[i].type_ == "Zn":
                oldID2newID[self.atoms[i].id_] = j
                newID2oldID[j] = self.atoms[i].id_
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "H":
                oldID2newID[self.atoms[i].id_] = j
                newID2oldID[j] = self.atoms[i].id_
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "C":
                oldID2newID[self.atoms[i].id_] = j
                newID2oldID[j] = self.atoms[i].id_
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "N":
                oldID2newID[self.atoms[i].id_] = j
                newID2oldID[j] = self.atoms[i].id_
                j += 1
        #print("old2new ID:",oldID2newID)
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
        #print(newbonds)
        for i in newbonds:
            for j in newbonds[i]:
                if i < j:
                    #print(newID2oldID[i], newID2oldID[j], self.atoms[newID2oldID[i]].type_, self.atoms[newID2oldID[j]].type_ )
                    typ = bondtype[self.atoms[newID2oldID[i]].type_ + "-" + self.atoms[newID2oldID[j]].type_]
                    f.write("%d %d %d\n" % (int(typ), i, j))
        f.close()


    # delete atoms and its associated bonds
    def delete_atom(self, atomID):
        assert atomID in self.atoms, "Atom %d does not exist!" % atomID
        neigh = self.bonds.pop(atomID)
        self.bondnum -= len(neigh)
        #print(self.bondnum)
        for i in neigh:
            if len(self.bonds[i]) == 1:
                self.bonds.pop(i)
            else:
                self.bonds[i].remove(atomID)
        self.atoms.pop(atomID)
        self.atomnum -= 1

    def add_atom(self, typ, coord):
        self.atomnum += 1
        newid = self.atomnum
        while newid in self.atoms:
            newid += 1
        self.atoms[newid] = Atom(newid, typ, coord[0], coord[1], coord[2])
        self.maxID = max(self.maxID, newid)
        return newid

    def add_bond(self, id1, id2):
        self.bondnum += 1
        assert id1 not in self.bonds or id2 not in self.bonds[id1], "Bond %d -> %d already exist!" % (id1, id2)

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
    parser.add_argument("OUTBOND", help="Output new bond list with benzine")
    parser.add_argument("OUTXYZ",  help="Output new xyz  file with benzine")
    parser.add_argument("Lx", type=float, help="Box length X")
    parser.add_argument("Ly", type=float, help="Box length Y")
    parser.add_argument("Lz", type=float, help="Box length Z")
    parser.add_argument("Nmin", type=int, help="minimum N index")
    parser.add_argument("Nmax", type=int, help="maximum N index")
    args = parser.parse_args()
    print(args.XYZ)

    ben = args.NUM
    mysystem = allatoms()
    mysystem.read_xyz(args.XYZ)
    mysystem.read_bonds(args.BOND)

    # add benzine
    # remember visited N atoms
    Nlist = set()
    boxlen = [args.Lx, args.Ly, args.Lz]  # for later use to bring atom position back into the box space
    for i in range(1, args.NUM + 1):
        # find first N
        flag = 1
        random = None
        while flag is not 0:
            width = args.Nmax - args.Nmin + 1
            random = rd.randint(1, width) * rd.randint(1, width) % width + args.Nmin
            #random = 267
            #random = 235
            if random in Nlist: continue
            else:
                Nlist.add(random)
                print("\n#============= Adding %dth benzine ring ==============#" % (i))
                print("choose N ID:", random)
                flag = 0
        assert random in mysystem.bonds, "Random Nitrogen ID not in the atom list!!"
        neighbor = mysystem.bonds[random]
        #print(neighbor)

        # get the C and N in N-C-N bonds
        Cindex, Nindex = None, None
        for j in neighbor:
            if mysystem.atoms[j].type_ is "C":
                #print(j)
                for k in mysystem.bonds[j]:
                    if mysystem.atoms[k].type_ is "N" and k != random:
                        Cindex, Nindex = j, k
                        Nlist.add(Nindex)
                        break

        print("Cindex,Nindex", Cindex, Nindex)
        # get the two H in H-C-C-H in imidazole
        # get the two C in H-C-C-H in imidazole

        Ca, Cb, Ha, Hb = None, None, None, None
        for j in neighbor:
            if mysystem.atoms[j].type_ is "C":
                for k in mysystem.bonds[j]:
                    if mysystem.atoms[k].type_ is "C":
                        Ca = j
                        #print("Ca",Ca)
                        break
        for j in mysystem.bonds[Ca]:
            if mysystem.atoms[j].type_ is "C":
                #print("Cb",Cb)
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
        print("Bond # after delete 2 H", mysystem.bondnum)

        # add 4 C  and  4  H  as shown below (nC#, nH#)
        #
        #   C-C  1.39     C-N  1.34      C-H 1.02     N-Zn 1.97
        #
        #        nH1    Z    nH2
        #           \       /
        #           nC1-M—nC2
        #          /         \
        #   nH3--nC3    O    nC4—-nH4     #  O is the center
        #          \         /
        #           C2 —o— C1
        #           /       \
        #  Zn ---  N2   n   N1  --- Zn
        #           \       /
        #             \   /
        #               C

        N1, N2, C1, C2 = mysystem.atoms[random].xyz_, mysystem.atoms[Nindex].xyz_, \
                         mysystem.atoms[Ca].xyz_, mysystem.atoms[Cb].xyz_
        # remember translation operation, so all atoms can get back to their original position
        N1m, N2m, C1m, C2m = [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]
        # print("N1", N1)
        # print("N2", N2)
        for l in range(3):
            if N1[l] - N2[l] > 0.5 * boxlen[l]:
                N2[l] += boxlen[l]
                N2m[l] += 1
            if N1[l] - N2[l] < - 0.5 * boxlen[l]:
                N1[l] += boxlen[l]
                N1m[l] += 1
            if C1[l] - C2[l] > 0.5 * boxlen[l]:
                C2[l] += boxlen[l]
                C2m[l] += 1
            if C1[l] - C2[l] < - 0.5 * boxlen[l]:
                C1[l] += boxlen[l]
                C1m[l] += 1
        #print(N1m,N2m,C1m,C2m)

        #print(random,Nindex,Ca,Cb)
        newC1, newC2, newC3, newC4 = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
        newH1, newH2, newH3, newH4 = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]

        O = [0.0, 0.0, 0.0]
        M = [0.0, 0.0, 0.0]
        Z = [0.0, 0.0, 0.0]
        o = [0.0, 0.0, 0.0]
        n = [0.0, 0.0, 0.0]
        for l in range(3):
            n[l] = 0.5 * (N1[l] + N2[l])
            o[l] = 0.5 * (C1[l] + C2[l])
            o[l] = n[l] + 0.5 * (o[l] - n[l])
            C1[l] = o[l] + 0.5 * (N1[l] - n[l])
            C2[l] = o[l] - 0.5 * (N1[l] - n[l])
            O[l] = n[l] + 1.5 * (o[l] - n[l]) # 1.7
            M[l] = n[l] + 2.2 * (o[l] - n[l]) # 2.4
            Z[l] = n[l] + 2.8 * (o[l] - n[l]) # 3.0

        for l in range(3):
            newC1[l] = M[l] + 0.9 / 2 * (C2[l] - C1[l])
            newH1[l] = Z[l] + 1.3 / 2 * (C2[l] - C1[l])
            newC2[l] = M[l] - 0.9 / 2 * (C2[l] - C1[l])
            newH2[l] = Z[l] - 1.3 / 2 * (C2[l] - C1[l])
            newC3[l] = O[l] + 0.6 * (C2[l] - C1[l])  # 0.9
            newH3[l] = O[l] + 1.0 * (C2[l] - C1[l])  # 1.5
            newC4[l] = O[l] - 0.6 * (C2[l] - C1[l])  # 0.9
            newH4[l] = O[l] - 1.0 * (C2[l] - C1[l])  # 1.5
        # period boundary
        for l in range(3):

            if N1m[l] > 0:
                N1[l] -= boxlen[l]
            if N2m[l] > 0:
                N2[l] -= boxlen[l]
            if C1m[l] > 0:
                C1[l] -= boxlen[l]
            if C2m[l] > 0:
                C2[l] -= boxlen[l]

            while newC1[l] < 0:
                newC1[l] += boxlen[l]
            while newC1[l] > boxlen[l]:
                newC1[l] -= boxlen[l]
            while newC2[l] < 0:
                newC2[l] += boxlen[l]
            while newC2[l] > boxlen[l]:
                newC2[l] -= boxlen[l]
            while newC3[l] < 0:
                newC3[l] += boxlen[l]
            while newC3[l] > boxlen[l]:
                newC3[l] -= boxlen[l]
            while newC4[l] < 0:
                newC4[l] += boxlen[l]
            while newC4[l] > boxlen[l]:
                newC4[l] -= boxlen[l]
            while newH1[l] < 0:
                newH1[l] += boxlen[l]
            while newH1[l] > boxlen[l]:
                newH1[l] -= boxlen[l]
            while newH2[l] < 0:
                newH2[l] += boxlen[l]
            while newH2[l] > boxlen[l]:
                newH2[l] -= boxlen[l]
            while newH3[l] < 0:
                newH3[l] += boxlen[l]
            while newH3[l] > boxlen[l]:
                newH3[l] -= boxlen[l]
            while newH4[l] < 0:
                newH4[l] += boxlen[l]
            while newH4[l] > boxlen[l]:
                newH4[l] -= boxlen[l]

        C1id, C2id = mysystem.add_atom("C", newC1), mysystem.add_atom("C", newC2)
        C3id, C4id = mysystem.add_atom("C", newC3), mysystem.add_atom("C", newC4)
        H1id, H2id = mysystem.add_atom("H", newH1), mysystem.add_atom("H", newH2)
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
        print("Bond # after add bIm", mysystem.bondnum)

#abd local
    # output XYZ and BONDS information
    mysystem.outputxyz(args.OUTXYZ)
    mysystem.outputbond(args.OUTBOND)

    print("\n\nJob completed!!!")
#yj remote
if __name__ == "__main__":
    main()
