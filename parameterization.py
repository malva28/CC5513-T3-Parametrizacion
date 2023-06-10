import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import eigs, spsolve
import openmesh
import argparse
import os
import polyscope as ps
from transformations import normalize, get_bounding_box_extremes


# Esta función calcula el ángulo entre dos aristas u y v
def myangle(u, v):
    du = np.linalg.norm(u)
    dv = np.linalg.norm(v)

    du = max(du, 1e-8)
    dv = max(dv, 1e-8)

    return np.arccos(np.dot(u, v) / (du * dv))


# Esta función computa la matriz Laplaciana de una malla
# Como la matriz Laplaciana tiene muchos elementos ceros, se utiliza una matriz dispersa (lil_matrix)
def laplacian(mesh, weight_method):
    n = mesh.n_vertices()  # Num. vértices de la malla
    print(f"Num. vertices {n}")
    W = lil_matrix((n, n), dtype=float)  # Se crea una matriz dispersa de n x n
    print(W.shape)

    points = mesh.points()  # Se obtienen las coordenadas de los vértices de la malla

    # Para cada vertice de la malla
    for i, v in enumerate(mesh.vertices()):
        f_it = openmesh.VertexFaceIter(mesh, v)  # Se obtienen las caras que comparten el vértice v
        for f in f_it:  # Para cada cara f
            v_it = openmesh.FaceVertexIter(mesh, f)  # Se obtienen los vértices de la cara f
            L = []
            for vv in v_it:  # Se obtienen los vértices compartidos para esa cara
                if vv.idx() != i:
                    L.append(vv.idx())
            j = L[0]
            k = L[1]

            # Se obtienen las coordenadas de los vértices en una cara
            vi = points[i, :]
            vj = points[j, :]
            vk = points[k, :]

            res_ij, res_ik = weight_method(vi, vj, vk)

            # Se acumulan los cotangentes para las aristas ij e ik
            W[i, j] = W[i, j] + res_ij
            W[i, k] = W[i, k] + res_ik

    # Se suman todas las filas de la matriz W y se calcula la inversa de cada suma
    # rows represents Di values, sum axis=1 sums all the values j in Di.
    # These are the denominators of lambda ij, each corresponds to lambda i
    S = 1.0 / W.sum(axis=1)
    print(S.shape)

    # Se calcula la matriz Laplaciana normalizada. La diagonal es 1 y los demás elementos son -1/n. La idea es que
    # la suma de cada fila sea 0

    # spdiags(S) * W multiplies row i of W by the the value ii of the diagonal
    W = eye(n, n) - spdiags(np.squeeze(S), 0, n, n) * W
    return W


def harmonic(point_i, point_j, point_k):
    # Se calculan los ángulos alpha y beta
    alpha = myangle(point_i - point_k, point_j - point_k)
    beta = myangle(point_i - point_j, point_k - point_j)

    cot_ij = 1.0 / np.tan(alpha)
    cot_ik = 1.0 / np.tan(beta)
    return cot_ij, cot_ik


def mean_value(point_i, point_j, point_k):
    # Se calculan los ángulos gamma y delta
    gamma = delta = myangle(point_j - point_i, point_k - point_i)

    tan_ij = np.tan(gamma/2) / np.linalg.norm(point_i - point_j)
    tan_ik = np.tan(delta/2) / np.linalg.norm(point_i - point_k)
    return tan_ij, tan_ik


def uniform(point_i, point_j, point_k):
    return 1, 1


def harmonic_laplacian(mesh):
    return laplacian(mesh, harmonic)


def uniform_laplacian(mesh):
    return laplacian(mesh, uniform)


def mean_value_laplacian(mesh):
    return laplacian(mesh, mean_value)


# Esta función calcula el borde de una malla. Recibe como parámetros la matriz F de caras y el número de vértices
# Devuelve un arreglo con los índices de los puntos pertenecientes al borde.
# Puntos consecutivos están conectados por aristas
def compute_boundary(F, numVert):
    numFaces = F.shape[0]

    # Creamos una matriz de adyacencia para vértices
    A = np.zeros((numVert, numVert))

    # Se suma un valor a A[i,j] si los vértices i y j son adyacentes
    for i in range(numFaces):
        f = F[i, :]
        A[f[0], f[1]] = A[f[0], f[1]] + 1
        A[f[0], f[2]] = A[f[0], f[2]] + 1
        A[f[2], f[1]] = A[f[2], f[1]] + 1

    # Como la malla es un grafo no dirigido, también consideramos lo simétrico
    A = A + A.T

    # Una arista de frontera es aquella que tiene un valor 1 en la matriz A
    result = np.where(A == 1)
    boundary = [result[1][0], result[0][0]]  # Calculamos una primera arista de frontera

    # A partir de la primera arista de frontera, se calcula el resto del borde
    s = boundary[1]
    i = 1

    while i < numVert:
        # Buscar algún elemento con valor 1 en la fila s
        vals = []
        for j in range(numVert):
            if A[s, j] == 1:
                vals.append(j)

        assert (len(vals) == 2)
        if vals[0] == boundary[i - 1]:
            s = vals[1]
        else:
            s = vals[0]

        if s != boundary[0]:
            boundary.append(s)
        else:
            break
        i = i + 1
    return boundary


def compute_barycentric_coordinates(input_mesh, boundary_indices, bound_distance, laplacian_fun=harmonic_laplacian):
    # Se calcula la parametrización de la frontera en una circunferencia,
    # usando como guía la longitud de cada arista de la frontera
    numBound = len(boundary_indices)
    mesh_points = input_mesh.points()
    num_vert = mesh_points.shape[0]
    print("Chord length:", bound_distance)
    vb = mesh_points[boundary_indices, :]
    sel = list(range(1, numBound))
    sel.append(0)
    vb2 = vb[sel, :]

    # vb2 - vb computes the vector difference between consecutive points, then, the norm transforms this into distances
    D = np.cumsum(np.linalg.norm(vb2 - vb, axis=1))
    # cumsum sums all distnaces, so the last item in (D - D[0]) should be d, hence,
    # t is normalized so the distances lie in values [0,1]
    t = (D - D[0]) / bound_distance
    # polar coordinates. Each value in 2*pi*t corresponds to an angle
    xy_boundary = np.hstack([np.expand_dims(np.cos(2 * np.pi * t), axis=1),
                             np.expand_dims(np.sin(2 * np.pi * t), axis=1)])
    print(xy_boundary)

    # Se calcula la matriz Laplaciana
    L = laplacian_fun(input_mesh)

    # Removemos las filas de los vértices de la frontera
    for i in boundary_indices:
        L[i, :] = 0
        L[i, i] = 1

    # Creamos un vector X en donde de calcularán las coordenadas X de la parametrización
    x = np.zeros((num_vert, 1))
    # Se asignan las coordenadas X de los vértices de la frontera
    for i, b in enumerate(boundary_indices):
        x[b, 0] = xy_boundary[i, 0]

    # Se resuelve el sistema de ecuaciones
    xx = spsolve(L, x)

    # Hacemos lo mismo para coordenadas Y
    y = np.zeros((num_vert, 1))
    for i, b in enumerate(boundary_indices):
        y[b, 0] = xy_boundary[i, 1]
    yy = spsolve(L, y)

    min_bound, max_bound = get_bounding_box_extremes(mesh_points)
    z_dist = max_bound[2] - min_bound[2]
    # Se crea un arreglo de numpy con las coordenadas X e Y de la parametrización, se agregan coordenadas Z = 1 para poder dibujar
    coord3D = np.column_stack((xx, yy, np.ones(num_vert) * 1.3* z_dist))

    # Creamos un array con solamente las coordenadas X, Y de la parametrización. La usamos para dibujar en Polyscope
    coord2D = np.column_stack((xx, yy))
    print(coord2D.shape)
    return coord3D, coord2D


if __name__ == "__main__":
    # El script recibe como parámetro el nombre de un archivo .off
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, default="data/bunny.off", help="File containing the input mesh data")
    parser.add_argument("--weight", choices=["uniform", "harmonic", "mean_value"], default="harmonic",
                        help='Specifies the type of weight used when rendering one model')
    parser.add_argument("--hide_barycentric", action='store_false', help='Hides the circle with the barycentric coordinates in polyscope.')

    opt = parser.parse_args()

    # Usamos OpenMesh para leer la malla
    mesh = openmesh.read_trimesh(opt.file)
    points = mesh.points()

    normalize(points)
    bb_min, bb_max = get_bounding_box_extremes(points)
    box_factor = 2/ (bb_max[2] - bb_min[2])
    points *= box_factor
    bb_min, bb_max = get_bounding_box_extremes(points)
    print(bb_min, bb_max, bb_max[2] - bb_min[2])

    numVert = points.shape[0]

    # Calculamos la frontera
    boundary = compute_boundary(mesh.face_vertex_indices(), points.shape[0])


    # Calculamos la longitud de la frontera
    # d is the total distance
    d = 0
    lastBound = boundary[-1]

    for i in boundary:
        d = d + np.linalg.norm(points[i, :] - points[lastBound, :])
        lastBound = i

    if opt.weight == "uniform":
        lap_fun = uniform_laplacian
    elif opt.weight == "harmonic":
        lap_fun = harmonic_laplacian
    else:
        lap_fun = mean_value_laplacian

    coord, coord2 = compute_barycentric_coordinates(mesh, boundary, d, lap_fun)

    # Dibujamos en Polyscope
    ps.init()
    # Malla original con la parametrización. Internamente Polyscope aplica una textura usando la parametrización
    ps_mesh = ps.register_surface_mesh("mesh", points, mesh.face_vertex_indices())
    ps_mesh.add_parameterization_quantity("parameterization", coord2, defined_on='vertices', enabled=True)

    # Mostramos la parametrización como si fuera una superficie más
    ps_mesh2 = ps.register_surface_mesh("mesh2",
                                        coord,
                                        mesh.face_vertex_indices(),
                                        edge_width=1,
                                        edge_color=[0, 0, 1],
                                        enabled=opt.hide_barycentric,
                                        transparency=0.7)
    ps.show()
