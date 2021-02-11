import numpy
from numpy import linalg as la
import Data
import Tesselations
import matplotlib.pyplot as plt
import Laguerre
def Plattice():
    data = Data.DATA
    lattice = numpy.array(Data.DATA['Lattice'])
    P_lattice = numpy.array(Data.DATA['Plattice'])
    S, R = Tesselations.semi1(300, 100, 1)
    fixl = len(S)
    length = fixl
    matrx = {
        '0': numpy.empty((0, 2*len(S)), dtype = complex),
        '1': numpy.empty((0, 2*len(S)), dtype = complex),
    }
    ranges = {
        '1': [],
        '2': [],
        '3': [],
        '4': [],
        '5': [],
        '6': [],
        '7': [],
        '8': [],
        '9': []
    }
    for i in range(len(lattice.T)):
        if i==0:
            ranges[str(int(lattice[1, i]))].append(0)
        if i>=1:
            ranges[str(int(lattice[1, i-1]))].append(ranges[str(int(lattice[1, i-1]))][0]+int(lattice[2, i-1])-1)
            ranges[str(int(lattice[1, i]))].append(ranges[str(int(lattice[1, i-1]))][1]+1)
        for j in range(int(lattice[2, i])):
            if int(lattice[1, i]) == 1:
                S, R= Tesselations.semi1(300, 100, 1)
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/1.5/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/1.5/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/2/int(lattice[2, i]))*lattice[4, i]
            if int(lattice[1, i]) == 2:
                S, R= Tesselations.semi2(300, 100, 1)
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/1.5/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/1.5/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/2/int(lattice[2, i]))*lattice[4, i]    
            if int(lattice[1, i]) == 3:
                S, R = Tesselations.triangle(420, 100, Theta = numpy.sqrt(3))
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/1.5/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/1.5/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/300)*lattice[4, i] 
            if int(lattice[1, i])== 4:
                S, R= Tesselations.semi1(300, 100, 1.22499)
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/3/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/3/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/400)*lattice[4, i] 
            if int(lattice[1, i]) == 5:
                S, R= Tesselations.semi1(300, 100, 1.31607)
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/3/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/3/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/4/int(lattice[2, i]))*lattice[4, i]        
            if int(lattice[1, i]) == 6:
                S, R= Tesselations.semi6(300, 100, 1)
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/3/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/3/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/4/int(lattice[2, i]))*lattice[4, i] 
            if int(lattice[1, i]) == 7:
                S, R= Tesselations.semi3(300, 100, 1)
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/3/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/3/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/4/int(lattice[2, i]))*lattice[4, i] 
            if int(lattice[1, i]) == 8:
                S, R= Tesselations.semi4(300, 100, 1)
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/3/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/3/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/4/int(lattice[2, i]))*lattice[4, i] 
            if int(lattice[1, i]) == 9:
                S, R= Tesselations.semi5(300, 100, 1)
                for g in range(len(S)):
                    S[g,:] += numpy.array([numpy.random.normal(0.001, R[g]*j/3/int(lattice[2, i])), numpy.random.normal(0.001,R[g]*j/3/int(lattice[2, i]))])*lattice[3, i]
                    R[g] += numpy.random.normal(0.001, R[g]*j/4/int(lattice[2, i]))*lattice[4, i] 
            S = S[0:fixl,:]
            R = R[0:fixl]
            Si = numpy.zeros(len(S), dtype = complex)
            Si = numpy.vectorize(complex)(S[:, 0], S[:, 1])
            matrx['0']= numpy.append(matrx['0'], [numpy.append(Si,R, axis = 0)], axis = 0)
    U, s, Vt = la.svd(matrx['0'], full_matrices=False)
    if len(P_lattice.T)>0:
        fig4 = plt.figure(figsize = (10, 2.5*len(P_lattice.T)), constrained_layout=True)
        fig4.set_constrained_layout_pads(w_pad=4/72, h_pad=4/72, hspace=0, wspace=0)
        spec4 = fig4.add_gridspec(ncols=3, nrows=len(P_lattice.T))
        S1 = numpy.diag(s)
        col = ['b', 'r']
        for j in range(len(P_lattice.T)):
            s1 = 0
            for g in range(len(s)):
                s1 += s[g]/numpy.sum(s)
                if s1>=P_lattice[1, j]/100:
                    k=g+1
                    break
            matrx['1'] = numpy.dot(U[:, 0:k], numpy.dot(S1[0:k, 0:k], Vt[0:k, :]))
            for i in range(3):
                start = ranges[str(int(P_lattice[0, j]))][0]
                index = start + int(P_lattice[2+i, j])-1
                S = numpy.zeros((length, 2))
                R = numpy.zeros(length)
                ax = fig4.add_subplot(spec4[j, i])
                ax.axis('off')
                for z in range(2):
                    for g in range(len(S)):
                        if z == 0:
                            S[g, 0] = numpy.real(matrx[str(z)][index, g])
                            S[g, 1] = numpy.imag(matrx[str(z)][index, g])
                        if z == 1:
                            S[g, 0] = numpy.real(matrx[str(z)][index, g])
                            S[g, 1] = numpy.imag(matrx[str(z)][index, g])   
                    R[:] = numpy.real(matrx[str(z)][index, 1*length:2*length])
                    # R = numpy.ones(length)*100*numpy.sqrt(2/400/numpy.sqrt(3))
                    tri_list, V = Laguerre.get_power_triangulation(S, R)
                    # Compute the Voronoi cells
                    voronoi_cell_map = Laguerre.get_voronoi_cells(S, V, tri_list)
                    # Display the result
                    edge_map = Laguerre.display(S, R, tri_list, voronoi_cell_map, 100)
                    elem_con, node_coord = Laguerre.cleanV1(100, V, edge_map)
                    # ax.plot(S[:, 0], S[:, 1], 'r.', markersize = 2)
                    for i in range(len(elem_con)):
                        ax.plot([node_coord[int(elem_con[i, 0]), 0], node_coord[int(elem_con[i, 1]), 0]], [node_coord[int(elem_con[i, 0]), 1], node_coord[int(elem_con[i, 1]), 1]], col[z], linewidth = 1)
                    ax.set_xlim(numpy.mean(S[:, 0])-1.8*numpy.std(S[:, 0]), numpy.mean(S[:, 0])+1.8*numpy.std(S[:, 0]))
                    ax.set_ylim(numpy.mean(S[:, 1])-1.8*numpy.std(S[:, 1]), numpy.mean(S[:, 1])+1.8*numpy.std(S[:, 1]))
                    ax.set_aspect(1)
        plt.show()
        plt.savefig('compare.png')
    scum = numpy.zeros(len(s))
    for i in range(len(s)):
        scum[i] = numpy.sum(s[:i])
    plt.plot(numpy.arange(0, 100, 1), scum[0:100]/numpy.sum(s))
    plt.plot(numpy.arange(0, 100, 1), scum[0:100]/numpy.sum(s), 'r.')
    plt.xticks(numpy.arange(0, 101, 10))
    plt.yticks(numpy.linspace(0, 1, 11))
    plt.show()
    return U, s, Vt, matrx['0'], scum
        