import sys
import subprocess
import numpy as np
import openmesh
import argparse
import polyscope as ps
from parameterization import compute_boundary, compute_barycentric_coordinates
from parameterization import harmonic_laplacian, uniform_laplacian, mean_value_laplacian
from transformations import normalize, get_bounding_box_extremes

if __name__ == "__main__":
    # El script recibe como parámetro el nombre de un archivo .off
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, default="data/bunny.off", help="File containing the input mesh data")
    display_group = parser.add_mutually_exclusive_group()
    display_group.add_argument('--show_three', action='store_true',
                               help='Display at once the model three times with each different laplacian method.')
    display_group.add_argument('--harmonic_vs_mean', action='store_true',
                               help='Display overlapping models with the harmonic and mean-value weights applied, for better comparison')
    display_group.add_argument('--show_only_one', action='store_true',
                               help='Displays the model once. Must specify the method used for weight assignment')

    parser.add_argument("--weight", choices=["uniform", "harmonic", "mean_value"], default="harmonic",
                        required='--show_only_one' in sys.argv,
                        help='Specifies the type of weight used when rendering one model')
    parser.add_argument("--hide_barycentric", action='store_false',
                        help='Hides the circle with the barycentric coordinates in polyscope.')

    opt = parser.parse_args()
    print(opt)

    show_three = not opt.harmonic_vs_mean and not opt.show_only_one

    if opt.show_only_one:
        weight = opt.weight
        python_exe = sys.executable
        args = [python_exe, 'parameterization.py',
                '--file', opt.file,
                '--weight', weight]
        if not opt.hide_barycentric:
            args.append("--hide_barycentric")
        process = subprocess.run(args)

        ret = process.returncode
        print("== Subprocess finished with exit code {} ==".format(ret))

    else:
        # Usamos OpenMesh para leer la malla
        mesh = openmesh.read_trimesh(opt.file)
        points = mesh.points()
        normalize(points)
        bb_min, bb_max = get_bounding_box_extremes(points)
        box_factor = 2 / (bb_max[2] - bb_min[2])
        points *= box_factor

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

        coords_3D = []
        coords_2D = []

        lap_dict = {
            "uniform": uniform_laplacian,
            "harmonic": harmonic_laplacian,
            "mean_value": mean_value_laplacian
        }

        for lap_fun in list(lap_dict.values()):
            coord, coord2 = compute_barycentric_coordinates(mesh, boundary, d, lap_fun)
            coords_3D.append(coord)
            coords_2D.append(coord2)

        if show_three:
            # Dibujamos en Polyscope
            ps.init()
            # Malla original con la parametrización. Internamente Polyscope aplica una textura usando la parametrización
            min_x = np.min([points[:, 0], coords_3D[1][:, 0]])
            max_x = np.max([points[:, 0], coords_3D[1][:, 0]])
            distance_x = np.abs(max_x - min_x)

            lap_labels = list(lap_dict.keys())

            for i in range(len(lap_dict)):
                shift = (i - 1) * 1.1 * distance_x
                shifted_points = np.copy(points)
                shifted_points[:, 0] += shift
                coord = coords_3D[i]
                coord[:, 0] += shift

                coord2 = coords_2D[i]

                ps_mesh = ps.register_surface_mesh("mesh_" + lap_labels[i], shifted_points, mesh.face_vertex_indices())
                ps_mesh.add_parameterization_quantity("parameterization_" + lap_labels[i], coord2,
                                                      defined_on='vertices', enabled=True)

                # Mostramos la parametrización como si fuera una superficie más
                ps_mesh2 = ps.register_surface_mesh("mesh2D_" + lap_labels[i],
                                                    coord,
                                                    mesh.face_vertex_indices(),
                                                    edge_width=1,
                                                    edge_color=[0, 0, 1],
                                                    enabled=opt.hide_barycentric,
                                                    transparency=0.7)
            ps.show()

        elif opt.harmonic_vs_mean:
            ps.init()
            # Malla original con la parametrización. Internamente Polyscope aplica una textura usando la parametrización

            lap_labels = list(lap_dict.keys())

            for i in range(1, len(lap_dict)):
                coord = coords_3D[i]
                coord2 = coords_2D[i]

                ps_mesh = ps.register_surface_mesh("mesh_" + lap_labels[i], points, mesh.face_vertex_indices())
                ps_mesh.add_parameterization_quantity("parameterization_" + lap_labels[i], coord2,
                                                      defined_on='vertices',
                                                      enabled=True)

                # Mostramos la parametrización como si fuera una superficie más
                ps_mesh2 = ps.register_surface_mesh("mesh2D_" + lap_labels[i],
                                                    coord,
                                                    mesh.face_vertex_indices(),
                                                    edge_width=1,
                                                    edge_color=[0, 0, 1],
                                                    enabled=opt.hide_barycentric,
                                                    transparency=0.7)
            ps.show()
