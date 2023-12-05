"""
Create segmented mirrors from parameters.
"""
import numpy as np
import copy
import optic_object_wavefronts as optcad
import merlict
import posixpath


def make_example_config_heritage():
    return {
        "focal_length": 4.889,
        "DaviesCotton_over_parabolic_mixing_factor": 0.5,
        "max_outer_aperture_radius": 2.1,
        "min_inner_aperture_radius": 0.2,
        "outer_aperture_shape_hex": 0,
        "facet_inner_hex_radius": 0.30,
        "gap_between_facets": 0.01,
    }


def add_to_frame_heritage(
    frame,
    scenery,
    config,
    facet_body_width,
    boundary_layer_facet_front="facet/front",
    boundary_layer_facet_body="facet/body",
    facet_fn=7,
    ref="segmented_mirror",
):
    a = config["DaviesCotton_over_parabolic_mixing_factor"]

    approx_num_facets_on_outer_radius = (
        config["max_outer_aperture_radius"] / config["facet_inner_hex_radius"]
    )
    fn_circle = int(np.ceil(2.0 * np.pi * approx_num_facets_on_outer_radius))

    outer_radius_facet_supports = (
        config["max_outer_aperture_radius"] - config["facet_inner_hex_radius"]
    )
    if config["outer_aperture_shape_hex"] == 1:
        aperture_outer_polygon = (
            optcad.geometry.regular_polygon.make_vertices_xy(
                outer_radius=outer_radius_facet_supports,
                fn=6,
                rot=np.pi / 6,
            )
        )
    else:
        aperture_outer_polygon = (
            optcad.geometry.regular_polygon.make_vertices_xy(
                outer_radius=outer_radius_facet_supports,
                fn=fn_circle,
                rot=0.0,
            )
        )

    inner_radius_facet_supports = (
        config["min_inner_aperture_radius"] + config["facet_inner_hex_radius"]
    )
    aperture_inner_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
        outer_radius=inner_radius_facet_supports,
        fn=fn_circle,
    )

    return add_to_frame_in_scenery(
        frame=frame,
        scenery=scenery,
        focal_length=config["focal_length"],
        aperture_outer_polygon=aperture_outer_polygon,
        aperture_inner_polygon=aperture_inner_polygon,
        facet_inner_hex_radius=config["facet_inner_hex_radius"],
        gap_between_facets=config["gap_between_facets"],
        davies_cotton_weight=a,
        parabola_weight=(1.0 - a),
        sphere_weight=0.0,
        mean_distance_of_facet_centers_to_focal_point_is_focal_length=True,
        facet_rotation="individual",
        boundary_layer_facet_front=boundary_layer_facet_front,
        boundary_layer_facet_body=boundary_layer_facet_body,
        facet_fn=facet_fn,
        facet_body_width=facet_body_width,
        ref=ref,
    )


def add_to_frame(
    frame,
    focal_length,
    aperture_outer_polygon,
    aperture_inner_polygon,
    facet_inner_hex_radius,
    gap_between_facets,
    davies_cotton_weight,
    parabola_weight,
    sphere_weight,
    mean_distance_of_facet_centers_to_focal_point_is_focal_length,
    facet_rotation,
    facet_body_width,
    boundary_layer_facet_front="facet/front",
    boundary_layer_facet_body="facet/body",
    facet_fn=7,
    ref="segmented_mirror",
):
    """
    Parameters
    ----------
    frame : dict
        A frame in the geometry.
    geometry : dict
        The geometry.
    focal_length : float
        Distance from aperture's principal plane where image forms.
    aperture_outer_polygon : dict
        Polygon marking the outer bound for facet-centers.
    aperture_inner_polygon : dict
        Polygon marking the inner bound for facet-centers.
    facet_inner_hex_radius : float
        Inner radius of hexagonal facet.
    gap_between_facets : float
        Distance of gap between neighboring facets.
    mean_distance_of_facet_centers_to_focal_point_is_focal_length : bool
        Whether to shift the facet-centers in z or not.
    davies_cotton_weight : float
        Weight 0-1,
    parabola_weight : float
        Weight 0-1,
    sphere_weight : float
        Weight 0-1,
    facet_rotation : str
        How to rotate the facets.
    facet_body_width : float
        Smallest width of facet's body from back to working-surface.
    facet_fn : int
        Density of vertices and faces in facet.
    ref : str
        A name to distinguish this mirror from others.
    """

    # facet
    # -----
    facet_outer_radius = facet_inner_hex_radius * (2.0 / np.sqrt(3.0))
    facet_curvature_radius = 2.0 * focal_length
    facet_object_key = ref + "facet"

    facet = optcad.primitives.spherical_planar_lens_hexagonal.init(
        outer_radius=facet_outer_radius,
        curvature_radius=facet_curvature_radius,
        width=facet_body_width,
        fn=facet_fn,
        ref="facet",
    )
    facet_mtl_to_boundary_layers_map = {
        "facet/front": posixpath.join(ref, boundary_layer_facet_front),
        "facet/back": posixpath.join(ref, boundary_layer_facet_body),
        "facet/side": posixpath.join(ref, boundary_layer_facet_body),
    }

    objs = {}
    objs[facet_object_key] = optcad.export.reduce_mesh_to_obj(facet)

    grid_spacing = (2 * facet_inner_hex_radius) + gap_between_facets

    facet_centers = init_facet_centers_xy(
        aperture_outer_polygon=aperture_outer_polygon,
        aperture_inner_polygon=aperture_inner_polygon,
        grid_spacing=grid_spacing,
        grid_style="hexagonal",
        grid_rotation=np.pi / 2,
        center_of_grid=[0.0, 0.0],
        ref="facet_centers",
    )

    facet_centers = set_facet_centers_z(
        facet_centers=facet_centers,
        focal_length=focal_length,
        davies_cotton_weight=davies_cotton_weight,
        parabola_weight=parabola_weight,
        sphere_weight=sphere_weight,
    )

    if mean_distance_of_facet_centers_to_focal_point_is_focal_length:
        facet_centers = shift_facet_centers_z_so_mean_distance_of_facet_centers_to_focal_point_is_focal_length(
            facet_centers=facet_centers,
            focal_length=focal_length,
            max_delta=1e-6 * focal_length,
            max_iterations=1000,
        )

    # orientation
    # -----------
    facet_id = 0
    for fkey in facet_centers:
        if facet_rotation == "individual":
            focal_point = [0, 0, focal_length]
            rot = init_facet_rotation_with_focal_point(
                facet_center=facet_centers[fkey],
                focal_point=focal_point,
            )
        elif facet_rotation == "sphere":
            normal_point = [0, 0, 2 * focal_length]
            rot = init_facet_rotation_with_normal_point(
                facet_center=facet_centers[fkey],
                normal_point=normal_point,
            )
        else:
            raise KeyError
        child = {
            "id": int(facet_id),
            "pos": facet_centers[fkey],
            "rot": rot,
            "obj": facet_object_key,
            "mtl": facet_mtl_to_boundary_layers_map,
        }
        frame["children"].append(child)
        facet_id += 1

    return objs


def init_facet_centers_xy(
    aperture_outer_polygon=optcad.geometry.regular_polygon.make_vertices_xy(
        outer_radius=1.0
    ),
    aperture_inner_polygon=optcad.geometry.regular_polygon.make_vertices_xy(
        outer_radius=0.5
    ),
    grid_spacing=0.1,
    grid_style="hexagonal",
    grid_rotation=np.pi / 2,
    center_of_grid=[0.0, 0.0],
    ref="grid",
):
    _, min_max_distances = optcad.polygon.find_min_max_distant_to_point(
        polygon=aperture_outer_polygon, point=center_of_grid
    )
    outer_radius = min_max_distances[1]
    fN = 2 * int(np.ceil(outer_radius / grid_spacing))

    if grid_style == "hexagonal":
        _grid = optcad.geometry.grid.hexagonal.init_from_spacing(
            spacing=grid_spacing, ref=ref, fN=fN
        )
    elif grid_style == "rectangular":
        _grid = optcad.geometry.grid.rectangular.init_from_spacing(
            spacing=grid_spacing, ref=ref, fN=fN
        )
    else:
        assert False, "grid style {:s} is unknown.".format(grid_style)

    _grid = optcad.polygon.rotate_z(_grid, grid_rotation)
    mask_inside_outer = optcad.polygon.mask_vertices_inside(
        vertices=_grid, polygon=aperture_outer_polygon
    )
    if aperture_inner_polygon is None:
        mask = mask_inside_outer
    else:
        mask_inside_inner = optcad.polygon.mask_vertices_inside(
            vertices=_grid, polygon=aperture_inner_polygon
        )
        mask_outside_inner = np.logical_not(mask_inside_inner)
        mask = np.logical_and(mask_inside_outer, mask_outside_inner)
    return optcad.polygon.keep_vertices_in_mask(vertices=_grid, mask=mask)


def set_facet_centers_z(
    facet_centers,
    focal_length,
    davies_cotton_weight,
    parabola_weight,
    sphere_weight,
):
    assert focal_length > 0.0

    assert 0 <= davies_cotton_weight <= 1.0
    assert 0 <= parabola_weight <= 1.0
    assert 0 <= sphere_weight <= 1.0
    assert davies_cotton_weight + parabola_weight + sphere_weight == 1.0

    for fkey in facet_centers:
        davies_cotton_z = optcad.geometry.sphere.surface_height(
            x=facet_centers[fkey][0],
            y=facet_centers[fkey][1],
            curvature_radius=focal_length,
        )
        prabola_z = optcad.geometry.parabola.surface_height(
            x=facet_centers[fkey][0],
            y=facet_centers[fkey][1],
            focal_length=focal_length,
        )
        sphere_z = optcad.geometry.sphere.surface_height(
            x=facet_centers[fkey][0],
            y=facet_centers[fkey][1],
            curvature_radius=2.0 * focal_length,
        )
        facet_centers[fkey][2] = (
            davies_cotton_z * davies_cotton_weight
            + prabola_z * parabola_weight
            + sphere_z * sphere_weight
        )

    return facet_centers


def shift_facet_centers_z_so_mean_distance_of_facet_centers_to_focal_point_is_focal_length(
    facet_centers,
    focal_length,
    max_delta,
    max_iterations=1000,
):
    focal_point = [0.0, 0.0, focal_length]
    i = 0
    while True:
        delta = focal_length - vertices_mean_distance_to_vertex(
            vertices=facet_centers,
            vertex=focal_point,
        )
        if delta < max_delta:
            break
        i += 1
        assert i < max_iterations
        shift = [0.0, 0.0, -delta / 2]
        facet_centers = vertices_add(vertices=facet_centers, vertex=shift)
    return facet_centers


def vertices_mean_distance_to_vertex(vertices, vertex):
    vertex = np.array(vertex)
    distances = []
    for fkey in vertices:
        d = np.linalg.norm(vertex - vertices[fkey])
        distances.append(d)
    return np.mean(distances)


def vertices_add(vertices, vertex):
    vertex = np.array(vertex)
    for fkey in vertices:
        vertices[fkey] += vertex
    return vertices


def init_facet_rotation_with_normal_point(
    facet_center,
    normal_point,
    direction_incoming_light=[0.0, 0.0, -1.0],
):
    axis, angle = facet_rotation_axis_and_angle(
        facet_center=facet_center,
        target_point=normal_point,
        direction_incoming_light=direction_incoming_light,
    )

    if np.abs(angle) > 0.0:
        return {
            "repr": "axis_angle",
            "axis": axis.tolist(),
            "angle_deg": float(np.rad2deg(angle)),
        }
    else:
        return {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]}


def init_facet_rotation_with_focal_point(
    facet_center,
    focal_point,
    direction_incoming_light=[0.0, 0.0, -1.0],
):
    axis, angle = facet_rotation_axis_and_angle(
        facet_center=facet_center,
        target_point=focal_point,
        direction_incoming_light=direction_incoming_light,
    )
    angle = 0.5 * angle

    if np.abs(angle) > 0.0:
        return {
            "repr": "axis_angle",
            "axis": axis.tolist(),
            "angle_deg": float(np.rad2deg(angle)),
        }
    else:
        return {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]}


def facet_rotation_axis_and_angle(
    facet_center,
    target_point,
    direction_incoming_light=[0.0, 0.0, -1.0],
):
    facet_center = np.array(facet_center)
    target_point = np.array(target_point)

    direction_incoming_light = np.array(direction_incoming_light)
    _direction_of_light_travel_norm = np.linalg.norm(direction_incoming_light)
    assert _direction_of_light_travel_norm > 0.0
    direction_incoming_light /= _direction_of_light_travel_norm
    direction_of_light_coming_from = -1.0 * direction_incoming_light

    direction_facet_to_target = -facet_center + target_point
    direction_facet_to_target_norm = np.linalg.norm(direction_facet_to_target)
    assert direction_facet_to_target_norm > 0.0
    direction_facet_to_target /= direction_facet_to_target_norm

    rotation_angle = -1.0 * np.arccos(
        np.dot(direction_facet_to_target, direction_of_light_coming_from)
    )
    if rotation_angle == 0.0:
        rotation_axis = None
    else:
        rotation_axis = np.cross(
            direction_facet_to_target, direction_of_light_coming_from
        )
    return rotation_axis, rotation_angle
