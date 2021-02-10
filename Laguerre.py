import itertools
import numpy
from scipy.spatial import ConvexHull
import math
from math import pi
from matplotlib.collections import LineCollection
from matplotlib import pyplot as plot
from copy import deepcopy
from functools import reduce
import pdb
from zipfile import ZipFile 
import os 
# --- Misc. geometry code -----------------------------------------------------
'''
Pick N points uniformly from the unit disc
This sampling algorithm does not use rejection sampling.
'''
def disc_uniform_pick(N):
	point = numpy.random.rand(N, 2)*100
	return point



def norm2(X):
	return numpy.sqrt(numpy.sum(X ** 2))



def normalized(X):
	return X / norm2(X)



# --- Delaunay triangulation --------------------------------------------------

def get_triangle_normal(A, B, C):
	return normalized(numpy.cross(A, B) + numpy.cross(B, C) + numpy.cross(C, A))



def get_power_circumcenter(A, B, C):
	N = get_triangle_normal(A, B, C)
	return (-.5 / N[2]) * N[:2]



def is_ccw_triangle(A, B, C):
	M = numpy.concatenate([numpy.stack([A, B, C]), numpy.ones((3, 1))], axis = 1)
	return numpy.linalg.det(M) > 0



def get_power_triangulation(S, R):
	# Compute the lifted weighted points
	S_norm = numpy.sum(S ** 2, axis = 1) - R ** 2
	S_lifted = numpy.concatenate([S, S_norm[:,None]], axis = 1)

	# Special case for 3 points
	if S.shape[0] == 3:
		if is_ccw_triangle(S[0], S[1], S[2]):
			return [[0, 1, 2]], numpy.array([get_power_circumcenter(*S_lifted)])
		else:
			return [[0, 2, 1]], numpy.array([get_power_circumcenter(*S_lifted)])

	# Compute the convex hull of the lifted weighted points
	hull = ConvexHull(S_lifted)
	
	# Extract the Delaunay triangulation from the lower hull
	tri_list = tuple([a, b, c] if is_ccw_triangle(S[a], S[b], S[c]) else [a, c, b]  for (a, b, c), eq in zip(hull.simplices, hull.equations) if eq[2] <= 0)
	
	# Compute the Voronoi points
	V = numpy.array([get_power_circumcenter(*S_lifted[tri]) for tri in tri_list])

	# Job done
	return tri_list, V


# --- Compute Voronoi cells ---------------------------------------------------

'''
Compute the segments and half-lines that delimits each Voronoi cell
  * The segments are oriented so that they are in CCW order
  * Each cell is a list of (i, j), (A, U, tmin, tmax) where
     * i, j are the indices of two ends of the segment. Segments end points are
       the circumcenters. If i or j is set to None, then it's an infinite end
     * A is the origin of the segment
     * U is the direction of the segment, as a unit vector
     * tmin is the parameter for the left end of the segment. Can be -1, for minus infinity
     * tmax is the parameter for the right end of the segment. Can be -1, for infinity
     * Therefore, the endpoints are [A + tmin * U, A + tmax * U]
'''
def get_voronoi_cells(S, V, tri_list):
	# Keep track of which circles are included in the triangulation
	vertices_set = frozenset(itertools.chain(*tri_list))

	# Keep track of which edge separate which triangles
	edge_map = { }
	for i, tri in enumerate(tri_list):
		for edge in itertools.combinations(tri, 2):
			edge = tuple(sorted(edge))
			if edge in edge_map:
				edge_map[edge].append(i)
			else:
				edge_map[edge] = [i]

	# For each triangle
	voronoi_cell_map = { i : [] for i in vertices_set }

	for i, (a, b, c) in enumerate(tri_list):
		# For each edge of the triangle
		for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
		# Finite Voronoi edge
			edge = tuple(sorted((u, v)))
			if len(edge_map[edge]) == 2:
				j, k = edge_map[edge]
				if k == i:
					j, k = k, j
				
				# Compute the segment parameters
				U = V[k] - V[j]
				U_norm = norm2(U)

				# Add the segment
				voronoi_cell_map[u].append(((j, k), (V[j], U / U_norm, 0, U_norm)))
			else: 
			# Infinite Voronoi edge
				# Compute the segment parameters
				A, B, C, D = S[u], S[v], S[w], V[i]
				U = normalized(B - A)
				I = A + numpy.dot(D - A, U) * U
				W = normalized(I - D)
				if numpy.dot(W, I - C) < 0:
					W = -W	
			
				# Add the segment
				voronoi_cell_map[u].append(((edge_map[edge][0], -1), (D,  W, 0, None)))				
				voronoi_cell_map[v].append(((-1, edge_map[edge][0]), (D, -W, None, 0)))				

	# Order the segments
	def order_segment_list(segment_list):
		# Pick the first element
		first = min((seg[0][0], i) for i, seg in enumerate(segment_list))[1]

		# In-place ordering
		segment_list[0], segment_list[first] = segment_list[first], segment_list[0]
		for i in range(len(segment_list) - 1):
			for j in range(i + 1, len(segment_list)):
				if segment_list[i][0][1] == segment_list[j][0][0]:
					segment_list[i+1], segment_list[j] = segment_list[j], segment_list[i+1]
					break

		# Job done
		return segment_list

	# Job done
	return { i : order_segment_list(segment_list) for i, segment_list in voronoi_cell_map.items() }



# --- Plot all the things -----------------------------------------------------

def display(S, R, tri_list, voronoi_cell_map, Side):
	# Setup
	# fig, ax = plot.subplots()
	# plot.axis('equal')
	# plot.axis('off')	

	# Set min/max display size, as Matplotlib does it wrong
	min_corner = numpy.amin(S, axis = 0) - numpy.max(R)
	max_corner = numpy.amax(S, axis = 0) + numpy.max(R)

	# # Plot the samples
	# for Si, Ri in zip(S, R):
	# 	ax.add_artist(plot.Circle(Si, Ri, fill = True, alpha = .4, lw = 0., color = '#8080f0', zorder = 1))

	# # Plot the power triangulation
	# edge_set = frozenset(tuple(sorted(edge)) for tri in tri_list for edge in itertools.combinations(tri, 2))
	# line_list = LineCollection([(S[i], S[j]) for i, j in edge_set], lw = 1., colors = '.9')
	# line_list.set_zorder(0)
	# ax.add_collection(line_list)

	# Plot the Voronoi cells
	edge_map = { }
	for segment_list in voronoi_cell_map.values():
		for edge, (A, U, tmin, tmax) in segment_list:
			edge = tuple(sorted(edge))
			if edge not in edge_map:
				if tmax is None:
					tmax = 10
				if tmin is None:
					tmin = -10
				edge_map[edge] = (A + tmin * U, A + tmax * U)
	
	# line_list = LineCollection(edge_map.values(), lw = 1., colors = 'k')
	# line_list.set_zorder(0)
	# ax.add_collection(line_list)
	# ax.set_xlim(xmin=0.0, xmax=Side)
	# ax.set_ylim(ymin=0.0, ymax=Side)
	# ax.plot(S[:, 0], S[:, 1], 'r.')
	# Job done
	# plot.xlim((0, Side))
	# plot.ylim((0, Side))
	# plot.show()
	return edge_map

def reorder(node_C, V, min = 3):
	node_con = []
	for node in node_C:
		node = numpy.unique(node)
		points = []
		for i in node:
			points.append(V[i])
		points = numpy.array(points)
		if(len(points)>2):
		  try:
			  hull = ConvexHull(points)
		  except Exception as e:
			  pdb.set_trace()
		  connect = []
		  for z in range(len(hull.vertices)):
		    connect.append(node[hull.vertices[z]])
		  node_con.append(numpy.array(connect))
		elif len(node)> min:
			node_con.append(node)
	node_con = numpy.array(node_con)
	return node_con

def get_cell_edges(S, R, V, lim):

	# Compute the lifted weighted points
	S_norm = numpy.sum(S ** 2, axis = 1) - R ** 2
	S_lifted = numpy.concatenate([S, S_norm[:,None]], axis = 1)
 
	# Compute the convex hull of the lifted weighted points
	hull = ConvexHull(S_lifted)
 
	# Extract the Delaunay triangulation from the lower hull
	tri_list = tuple([a, b, c] if is_ccw_triangle(S[a], S[b], S[c]) else [a, c, b]  for (a, b, c), eq in zip(hull.simplices, hull.equations) if eq[2] <= 0)
	node_connectivity = []
	for i in range(len(S)):
			node_con = []
			for tri in tri_list:
				for k in tri:
					if i == k:
						for z in range(len(V)):
							if numpy.round(get_power_circumcenter(*S_lifted[tri])[0], 4) == V[z, 0]:
								if numpy.round(get_power_circumcenter(*S_lifted[tri])[1], 4) == V[z, 1]:
									if z not in node_con:
									  node_con.append(z)
			node_con = numpy.array(node_con)
			node_connectivity.append(node_con)
		
	node_connectivity = numpy.array(node_connectivity)
	node_C = []
	for node in node_connectivity:
		nodec = []
		for i in range(len(node)):
			nodec.append(node[i])
		nodec = numpy.array(nodec)
		node_C.append(nodec)
	node_C = numpy.array(node_C)
	return node_connectivity

def elemcon(node_con):
	elem_con = []
	for node in node_con:
		for i in range(len(node)):
			if(i==len(node)-1):
				elem_con.append([node[i], node[0]]) if node[i]<node[0] else elem_con.append([node[0], node[i]])
			else:
				elem_con.append([node[i], node[i+1]]) if node[i]<node[i+1] else elem_con.append([node[i+1], node[i]])
	return numpy.array(numpy.unique(elem_con, axis = 0))
 
def cleanV(V, Side):
	nV = numpy.around(V, decimals=4)
	nV1 = nV[numpy.isfinite(nV).any(axis=1)]
	newV = nV1[numpy.unique(numpy.round(nV1, 1), axis=0, return_index=True)[1]]
	nV = []
	for i in range(len(newV)):
		if abs(newV[i, 0])<2*Side and abs(newV[i, 1]<2*Side) and newV[i, 0]>=-Side/5 and newV[i, 1]>=-Side/5:
			nV.append(newV[i])
	nV = numpy.array(nV)
	return nV

def cleanV1(Side, V, edge_map):
  Side = 100
  V1 = numpy.empty((0, 3))
  for i in range(len(V)):
    V1 = numpy.append(V1, [[i, V[i, 0], V[i, 1]]], axis = 0)
  V2 = numpy.empty((0, 3))
  for i in range(len(V1)):
    if abs(V1[i, 1])< Side*10 and abs(V1[i, 2])< Side*10:
      V2 = numpy.append(V2, [V1[i,:]], axis = 0)
  # V2 = V2[numpy.unique(numpy.round(V2[:, 1:], 1), axis=0, return_index=True)[1]]
  revmap = {}
  for i in range(len(V2)):
    key = V2[i, 0]
    revmap[key] = i
  elem_con = numpy.empty((0, 2))
  for edge in edge_map.keys():
    if edge[0] in revmap:
      if edge[1] in revmap:
        elem_con = numpy.append(elem_con, [[revmap[edge[0]], revmap[edge[1]]]], axis =0)
  return elem_con, V2[:, 1:]

def write_to_text(array, name):
  file = open(name, "w+")
  for row in array:
    content = ''
    for j in row:
      content += str(j)+' ' 
    content += '\n'
    file.write(content) 
  file.close()

def seed_pert(S, R, c, w):
  for i in range(0,len(S)):
    count = 0
    r = Side/numpy.sqrt(Seeds)
    attempt = 0
    while count == 0:
      New = S[i, :] + numpy.random.normal(0, c*r/4, (1, 2))
      newr = R[i] + R[i]*numpy.random.normal(0, w*r/4, (1, 2))
      attempt+=1
      if attempt>200:
        count = 1
      if numpy.all(numpy.sum((S-New)**2, axis = 1) - newr>0):
        S[i, :] = New
        count = 1
  return R
def weight_pert(S, R):
  for i in range(0,len(S)):
    count = 0
    r = Side/numpy.sqrt(Seeds)
    attempt = 0
    while count == 0:
      New = S[i, :] + numpy.random.normal(0, r/4, (1, 2))
      attempt+=1
      if attempt>200:
        count = 1
      if numpy.all(numpy.sum((S-New)**2, axis = 1) - 0.5>0):
        S[i, :] = New
        count = 1
  return R
