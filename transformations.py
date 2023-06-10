# Script tomado del repositorio de Cubolink
# https://github.com/Cubolink/Partial-Axial-Symmetry
import numpy as np


def get_bounding_box_extremes(points):
    # Find box-hull diagonal extremes
    (min_x, max_x) = (np.infty, -np.infty)
    (min_y, max_y) = (np.infty, -np.infty)
    (min_z, max_z) = (np.infty, -np.infty)
    for point in points:
        min_x = point[0] if point[0] < min_x else min_x
        min_y = point[1] if point[1] < min_y else min_y
        min_z = point[2] if point[2] < min_z else min_z

        max_x = point[0] if point[0] > max_x else max_x
        max_y = point[1] if point[1] > max_y else max_y
        max_z = point[2] if point[2] > max_z else max_z

    return np.array([min_x, min_y, min_z]), np.array([max_x, max_y, max_z])


def normalize(points):
    extreme_1, extreme_2 = get_bounding_box_extremes(points)
    # re-center to (0, 0, 0)
    recenter(points, (extreme_1+extreme_2)/2)

    # Scale by 1/box_diagonal
    distance = np.linalg.norm(extreme_1 - extreme_2)
    for point in points:
        point /= distance


def recenter(points, new_center):
    for point in points:
        point -= new_center


def rotate(points, rotation_matrix):
    for point in points:
        point[:] = np.matmul(rotation_matrix, point)
