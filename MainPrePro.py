
import numpy as np
import ipywidgets as widgets
from ipywidgets import Layout
import Data
def prevals():
    minp = [2, 2, 2, 3, 4]
    constraints = []
    flag = np.array([1, 0, 0, 0, 1, 0, 0, 0, 0, 0])
    props = np.array([[0, 29500, 29500, 0.3, 0.3, 29500 / (2 * (1 + 0.3))]])
    nodes2 =np.array([[1, 0, 3, 1, 1], [2, 4, 3, 1, 1], [3, 0, 0, 1, 1], [4, 4, 0, 1, 1]])
    ien2 = np.array([[1, 3, 4, 200, 400], [2, 1, 4, 200, 150], [3, 4, 2, 200, 200]])
    return minp, nodes2, ien2
class MainPrePro:
    def prevals(self):
        lattices = [('Honeycomb', 1), ('Square', 2), ('Triangle', 3), ('3.6.3.6', 4), ('3.12.12', 5), ('3.3.3.4.4', 6), ('3.4.6.4', 7), ('4.6.12', 8), ('4.8.8', 9)]
        nodes =np.array([[1, 'Honeycomb', 100, 1, 1, 1, 1, 0, 0], [2, 'Square', 100, 1, 1, 1, 1, 0, 0], [3, 'Triangle', 100, 1, 1, 1, 1, 0, 0], [4, '3.6.3.6', 100, 1, 1, 1, 1, 0, 0]])
        ien = np.array([[1, 87, 1, 50, 100], [2, 87, 1, 50, 100], [3, 87, 1, 50, 100]])
        return lattices, nodes, ien
    def Mesh(self, minp):
        self.meshtitle = widgets.Label( value = 'Mesh')
        self.meshtext = ['nsd', 'ndf', 'nen', 'nel', 'nnp']
        self.mlabel = widgets.HBox([widgets.Label(value=self.meshtext[j], layout = widgets.Layout(width='55px')) for j in range(len(self.meshtext)-2)])
        self.minp1 = widgets.HBox([widgets.IntText(value=self.minp[j], layout = widgets.Layout(width='55px'), disabled=True) for j in range(len(self.minp)-2)])
        self.msubmit = widgets.Button(description="Submit", layout= widgets.Layout(border = 'solid 1px black'))
        self.rmesh = widgets.VBox([self.mlabel, widgets.HBox([self.minp1])])
        return self.rmesh, self.msubmit
    
    def add_node(self, b):
        New_row = widgets.HBox(self.nitems[len(self.nitems)-1])
        self.noder0.children = self.noder0.children + (New_row,)
    def del_node(self, b):
        del_row = list(self.noder0.children)
        del_row = del_row[:-1]
        self.noder0.children = tuple(del_row)
    
    def nodes_widget(self, nodes, minp):
        dims = 2
        self.nodeTitle = widgets.Label(value='Lattices for SVD')
        self.nodetext = ['#','Lattice','Number', 'Perturb', 'Polydisperse',]
        self.nodes = nodes
        self.minp = minp
        self.node = [[] for i in range(len(self.nodes))]
        self.nitems = [[] for i in range(len(self.nodes))]
        ADDNODE = widgets.Button(description="Add Lattice", layout= widgets.Layout(border = 'solid 1px black'))
        DELNODE = widgets.Button(description="Remove Lattice", layout= widgets.Layout(border = 'solid 1px black'))
        self.nlabel= widgets.HBox([widgets.Label(value=self.nodetext[j], layout = widgets.Layout(width='100px')) for j in range(len(self.nodetext))])
        for i in range(4):
            if(i<len(self.nodes)):
                self.nitems[i].append(widgets.IntText(value = i+1, disabled = True, layout = widgets.Layout(width='100px')))
                self.nitems[i].append(widgets.Dropdown(options = self.minp, value = i+1, layout = widgets.Layout(width='100px')))
                for j in range(2, len(self.nodetext)):
                    self.nitems[i].append(widgets.FloatText(value=self.nodes[i, j], layout = widgets.Layout(width='100px')))
                # self.nitems[i].append(widgets.Text(value = '1, 50, 100'))
                self.node[i] = widgets.HBox(self.nitems[i])
        self.noder0 = widgets.VBox([self.node[j] for j in range(4)])
        self.brow = widgets.HBox([ADDNODE, DELNODE])
        self.rnode = widgets.VBox([self.nodeTitle, self.nlabel, self.noder0, self.brow], layout= widgets.Layout(border = 'solid 1px black'))
        ADDNODE.on_click(self.add_node)
        DELNODE.on_click(self.del_node)
        self.rnode
        return self.rnode, ADDNODE, DELNODE

    def elems_widget(self, ien):
        self.ienTitle = widgets.Label(value='Plots')
        self.ientext = ['Lattice','Accuracy(%)','Plot1','Plot2', 'Plot3']
        self.ien = ien
        self.ien1 = [[] for i in range(len(self.ien))]
        self.Citems = [[] for i in range(len(self.ien))]
        ADDien = widgets.Button(description="Add Lattice", layout= widgets.Layout(border = 'solid 1px black'))
        DELien = widgets.Button(description="Remove Lattice", layout= widgets.Layout(border = 'solid 1px black'))
        self.nlabel= widgets.HBox([widgets.Label(value=self.ientext[j], layout = widgets.Layout(width='100px')) for j in range(len(self.ientext))])
        for i in range(3):
            if(i<len(self.ien)):
                self.Citems[i].append(widgets.Dropdown(options = self.minp, value = i+1, layout = widgets.Layout(width='100px')))
                for j in range(1, len(self.ientext)):
                    self.Citems[i].append(widgets.FloatText(value=self.ien[i, j], layout = widgets.Layout(width='100px')))
                self.ien1[i] = widgets.HBox(self.Citems[i])
        self.ien1r0 = widgets.VBox([self.ien1[j] for j in range(3)])
        self.brow = widgets.HBox([ADDien, DELien])
        self.rien1 = widgets.VBox([self.ienTitle, self.nlabel, self.ien1r0, self.brow], layout= widgets.Layout(border = 'solid 1px black'))
        ADDien.on_click(self.add_ien)
        DELien.on_click(self.del_ien)
        self.rien1
        return self.rien1, ADDien

    def add_ien(self, b):
        New_row = widgets.HBox(self.Citems[len(self.Citems)-1])
        self.ien1r0.children = self.ien1r0.children + (New_row,)

    def del_ien(self, b):
        del_row = list(self.ien1r0.children)
        del_row = del_row[:-1]
        self.ien1r0.children = tuple(del_row)

    def fsubmit(self, b):
        self.xn = np.zeros((5, len(self.noder0.children)))
        # self.idb = np.zeros((2, len(self.noder0.children)))
        # self.f = np.zeros((2, len(self.noder0.children)))
        # self.g = np.zeros((self.dims, len(self.noder0.children)))
        for j in range(len(self.noder0.children)):
            for i in range(5):
                self.xn[i, j] = self.noder0.children[j].children[i].value
                # self.idb[i, j] = self.noder0.children[j].children[i+1+2].value
                # self.f[i, j] = self.noder0.children[j].children[i+1+2*2].value
                # self.g[i, j] = self.noder0.children[j].children[i+1+7].value
        self.ien = np.zeros((5, len(self.ien1r0.children)))
        # self.A = np.zeros(len(self.ien1r0.children))
        # self.E = np.zeros(len(self.ien1r0.children))
        for j in range(len(self.ien1r0.children)):
            for i in range(5):
                self.ien[i, j] = self.ien1r0.children[j].children[i].value
            # self.E[j] = self.ien1r0.children[j].children[3].value
            # self.A[j] = self.ien1r0.children[j].children[4].value
        # self.nsd = self.dims
        # self.ndf = self.minp[1]
        # self.nen = self.minp[2]

        Data.DATA['Lattice'] = self.xn
        # Data.DATA['Perturb'] = self.idb
        Data.DATA['Plattice'] = self.ien
        # Data.DATA['SVD'] = self.f
        # Data.DATA['Plot'] = self.g
        # Data.DATA['nsd'] = self.nsd
        # Data.DATA['ndf'] = self.ndf
        # Data.DATA['nen'] = self.nen
        # Data.DATA['nel'] = len(self.ien1r0.children)
        # Data.DATA['nnp'] = len(self.noder0.children)

    def run(self):
        self.minp, self.nodes, self.ien = self.prevals()
        self.rnode, self.ADDNODE, self.DELNODE = self.nodes_widget(self.nodes, self.minp)
        self.rien1, self.ADDNODE = self.elems_widget(self.ien)
        # self.rmesh, self.msubmit = self.Mesh(self.minp)
        self.Submit = widgets.Button(description = 'Submit', layout= widgets.Layout(border = 'solid 1px black'))
        self.page = widgets.VBox([self.rnode, self.rien1, self.Submit])
        self.Submit.on_click(self.fsubmit)
        return self.page


    