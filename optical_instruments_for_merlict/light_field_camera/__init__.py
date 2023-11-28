"""
Create a light-field camera from parameters.
"""
import numpy as np
import os
import posixpath
import optic_object_wavefronts as optcad
import merlict
from .. import segmented_mirror


def make_example_config():
    return {
        "expected_primary_optics_focal_length": 16.0,
        "expected_primary_optics_max_aperture_radius": 5.0,
        "max_FoV_diameter": np.deg2rad(4.5),
        "pixel_FoV_hex_flat2flat": np.deg2rad(0.2),
        "housing_overhead": 0,
        "num_paxel_on_pixel_diagonal": 5,
        "curved": 1,
    }


def make_geometry(
    primary_optics_focal_length,
    primary_optics_diameter,
    camera_field_of_view_polygon,
    eye_field_of_view_full_angle,
    eyes_point_towards_center_of_primary_optics,
    body_overhead,
    screen_curvature_radius=1e12,
):
    """
    primary_optics_focal_length : float
        The expected focal-length of the primary optics (probably a mirror).
    primary_optics_diameter : float
        The expected diameter of the primary optics.
    camera_field_of_view_polygon : OrderedDict[(str, [float,float,float])
        In units of rad. The outer perimiter of the field-of-view. Only eyes
        which fit into it will be added to the camera.
    eye_field_of_view_full_angle : float
        In units of rad.
    eyes_point_towards_center_of_primary_optics : bool
        If true, the camera's screen is curved and better suited for cameras
        with larger field-of-views.
    """
    assert primary_optics_focal_length > 0.0
    assert primary_optics_diameter > 0.0
    assert len(camera_field_of_view_polygon) >= 3
    assert eye_field_of_view_full_angle > 0.0
    assert body_overhead > 0.0

    c = {}
    c["primary_optics"] = {}
    c["primary_optics"]["focal_length"] = float(primary_optics_focal_length)
    c["primary_optics"]["diameter"] = float(primary_optics_diameter)
    c["primary_optics"]["radius"] = 0.5 * c["primary_optics"]["diameter"]
    c["primary_optics"]["fstop"] = (
        c["primary_optics"]["focal_length"] / c["primary_optics"]["diameter"]
    )
    c["primary_optics"]["center"] = np.array(
        [0.0, 0.0, -c["primary_optics"]["focal_length"]]
    )

    c["field_of_view"] = {}
    _, c["field_of_view"]["polygon"] = optcad.polygon.to_keys_and_numpy_array(
        polygon=camera_field_of_view_polygon
    )

    (
        c["field_of_view"]["cxlim"],
        c["field_of_view"]["cylim"],
        _,
    ) = optcad.polygon.limits(polygon=camera_field_of_view_polygon)

    _cxout = np.max(np.abs(c["field_of_view"]["cxlim"]))
    _cyout = np.max(np.abs(c["field_of_view"]["cylim"]))
    c["field_of_view"][
        "half_angle_to_safely_include_all_potential_eyes"
    ] = 1.2 * np.hypot(_cxout, _cyout)

    c["eye"] = {}
    c["eye"]["hull"] = {}
    c["eye"]["hull"]["field_of_view"] = {}
    c["eye"]["hull"]["field_of_view"]["full_angle"] = float(
        eye_field_of_view_full_angle
    )
    c["eye"]["hull"]["field_of_view"]["half_angle"] = (
        0.5 * c["eye"]["hull"]["field_of_view"]["full_angle"]
    )
    c["eye"]["hull"]["radius"] = (
        np.tan(c["eye"]["hull"]["field_of_view"]["half_angle"])
        * c["primary_optics"]["focal_length"]
    )
    c["eye"]["hull"]["diameter"] = 2.0 * c["eye"]["hull"]["radius"]

    grid_fn = int(
        np.ceil(
            c["field_of_view"][
                "half_angle_to_safely_include_all_potential_eyes"
            ]
            / c["eye"]["hull"]["field_of_view"]["full_angle"]
        )
    )
    _anticipated_viewing_directions_of_eyes = (
        optcad.geometry.grid.hexagonal.init_from_spacing(
            spacing=c["eye"]["hull"]["field_of_view"]["full_angle"],
            ref="",
            fN=grid_fn,
        )
    )
    anticipated_viewing_directions_of_eyes = (
        optcad.polygon.get_vertices_inside(
            vertices=_anticipated_viewing_directions_of_eyes,
            polygon=camera_field_of_view_polygon,
        )
    )

    c["screen"] = {}
    c["screen"]["curvature_radius"] = float(screen_curvature_radius)
    c["screen"]["eyes"] = []

    for ci, ckey in enumerate(anticipated_viewing_directions_of_eyes):
        eye = {}
        eye["continuous_eye_id"] = ci
        eye["hexagonal_grid_id"] = ckey

        eye["pos"] = (
            c["primary_optics"]["focal_length"]
            * anticipated_viewing_directions_of_eyes[ckey]
        )
        _sphere_z = optcad.geometry.sphere.surface_height(
            x=eye["pos"][0],
            y=eye["pos"][1],
            curvature_radius=c["screen"]["curvature_radius"],
        )
        eye["pos"][2] = -1.0 * _sphere_z

        if eyes_point_towards_center_of_primary_optics:
            _axis, _angle = segmented_mirror.facet_rotation_axis_and_angle(
                facet_center=eye["pos"],
                target_point=c["primary_optics"]["center"],
                direction_incoming_light=np.array([0.0, 0.0, 1.0]),
            )
            if np.abs(_angle) > 0.0:
                eye["rot"] = {
                    "repr": "axis_angle",
                    "axis": _axis,
                    "angle_deg": float(np.rad2deg(_angle)),
                }
            else:
                eye["rot"] = {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]}
        else:
            eye["rot"] = {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]}

        c["screen"]["eyes"].append(eye)

    eyes_positions = {}
    for i in range(len(c["screen"]["eyes"])):
        vkey = c["screen"]["eyes"][i]["hexagonal_grid_id"]
        eyes_positions[vkey] = c["screen"]["eyes"][i]["pos"]

    _screen_eyes_positions_voronoi_cells = (
        optcad.geometry.grid.hexagonal.init_voronoi_cells_from_centers(
            centers=eyes_positions,
            centers_spacing=c["eye"]["hull"]["diameter"],
            rot=0.0,
        )
    )

    c["screen"][
        "polygon"
    ] = optcad.geometry.grid.hexagonal.find_hull_of_voronoi_cells(
        voronoi_cells=_screen_eyes_positions_voronoi_cells,
        centers=eyes_positions,
        centers_spacing=c["eye"]["hull"]["diameter"],
    )

    fovpoly = {}
    for vkey in camera_field_of_view_polygon:
        fovpoly[vkey] = (
            c["primary_optics"]["focal_length"]
            * camera_field_of_view_polygon[vkey]
        )

    c["body"] = {}
    c["body"]["polygon"] = optcad.minkowski.minkowski_hull_xy(
        poly1=fovpoly,
        poly2=optcad.geometry.regular_polygon.make_vertices_xy(
            outer_radius=body_overhead,
            fn=max([24, len(camera_field_of_view_polygon)]),
        ),
        ref="body/outer_bound",
    )

    c["body"]["mesh"] = optcad.primitives.template_curved_surface.init(
        outer_polygon=c["body"]["polygon"],
        curvature_config={
            "curvature_radius": c["primary_optics"]["focal_length"]
        },
        curvature_height_function=optcad.geometry.sphere.surface_height,
        curvature_surface_normal_function=optcad.geometry.sphere.surface_normal,
        inner_polygon=c["screen"]["polygon"],
        fn_hex_grid=3 * grid_fn,
        ref="body",
    )

    return c


def add_to_frame_in_scenery(
    frame,
    scenery,
    light_field_camera_geometry,
):
    """
    Parameters
    ----------
    frame : dict
        This light-field camera will be made a child of this frame. The frame
        is itself a child of scenery["geometry"]["relations"].
    scenery : dict
        A merlict scenery in 'sceneryPy' representation.
    light_field_camera_geometry : dict
        The light-field camera's geometry. A dict with all the relevant
        positions and orientations pre calculated.
    """

    return scenery


"""

def add_segmented_mirror_to_frame_in_scenery(
    frame,
    scenery,

    camera_field_of_view,
    outer_medium="vacuum",
    camera_num_photosensors_on_diagonal=5,
    cameras_point_towards_center_of_primary_optics=True,
    camera_surface_mirror="perfect_mirror",
    camera_surface_body="perfect_absorber",
    camera_photosensor_surface="perfect_absorber/rgb_12_12_12",
    camera_photosensor_gap=1e-3,
    camera_lens_medium="glass",
    camera_lens_fn=7,
    ref="light_field_sensor",
):

    Parameters
    ----------
    frame : dict
        A frame in the scenery.
    scenery : dict
        The scenery.
    config : dict
        The geometry of the working-surface in Sebastian's format used since
        his Master-thesis.

    ref : str
        A name to distinguish this mirror from others.
    join = posixpath.join


    # camera
    # ------
    camera_ref = join(ref, "camera")
    camera_housing_inner_radius = (
        camera_field_of_view * primary_optics_focal_length
    )
    camera_housing_outer_radius = (2 / np.sqrt(3)) * camera_housing_inner_radius

    expected_primary_optics_fstop_number = (
        primary_optics_focal_length
        / primary_optics_diameter
    )

    camera_geometry = LightFieldCameraModule.make_geometry(
        housing_outer_radius=camera_housing_outer_radius,
        housing_wall_width=,
        housing_height=,
        lens_curvature_radius=,
        lens_fn=camera_lens_fn,
        photosensor_num_on_diagonal=,
        photosensor_gap=,
        photosensor_plane_distance=,
    )
    camera_mesh = LightFieldCameraModule.init(
        camera_geometry=camera_geometry,
        ref=camera_ref,
    )
    camera_mtl = {
        join(ref, "cam", "lens", "top"): join(ref, "lens"),
        join(ref, "cam", "lens", "bottom"): join(ref + "l"),
        join(ref, "cam", "lens", "side"): join(ref + "h"),
        join(ref, "cam", "housing", "top"): join(ref + "h"),
        join(ref, "cam", "housing", "bottom"): join(ref + "h"),
        join(ref, "cam", "housing", "outer"): join(ref + "h"),
        join(ref, "cam", "housing", "inner"): join(ref + "m"),
    }
    num_photosensors = len(
        camera_geometry["photosensor"]["grid"]["positions"]
    )
    for i in range(num_photosensors):
        mtlkey = join(ref, cam, "photosensor_{:06d}".format(i))
        camera_mtl[mtlkey] = ref + "p"

    # add objects
    # -----------
    assert camera_ref not in scenery["objects"]
    scenery["objects"][camera_ref] = camera_mesh

    # facet-supports
    # --------------
    approx_num_facets_on_outer_radius = (
        config["max_outer_aperture_radius"] / config["facet_inner_hex_radius"]
    )

    fn_circle = int(np.ceil(2.0 * np.pi * approx_num_facets_on_outer_radius))
    grid_spacing = (
        2.0 * config["facet_inner_hex_radius"] + config["gap_between_facets"]
    )

    # outer bound xy
    # --------------
    outer_radius_facet_supports = (
        config["max_outer_aperture_radius"] - config["facet_inner_hex_radius"]
    )
    inner_radius_facet_supports = (
        config["min_inner_aperture_radius"] + config["facet_inner_hex_radius"]
    )

    aperture_outer_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
        outer_radius=outer_radius_facet_supports, fn=fn_circle, rot=0.0,
    )

    facet_centers = init_facet_centers_xy(
        aperture_outer_polygon=aperture_outer_polygon,
        aperture_inner_polygon=optcad.geometry.regular_polygon.make_vertices_xy(
            outer_radius=inner_radius_facet_supports, fn=fn_circle,
        ),
        grid_spacing=grid_spacing,
        grid_style="hexagonal",
        grid_rotation=np.pi / 2,
        center_of_grid=[0.0, 0.0],
        ref="facet_centers",
    )

    facet_centers = set_facet_centers_z(
        facet_centers=facet_centers,
        focal_length=config["focal_length"],
        DaviesCotton_over_parabolic_mixing_factor=config[
            "DaviesCotton_over_parabolic_mixing_factor"
        ],
        max_delta=1e-6 * config["focal_length"],
        max_iterations=1000,
    )

    # orientation
    # -----------
    focal_point = [0, 0, config["focal_length"]]
    camera_id = 0
    for fkey in camera_centers:
        child = {
            "id": int(camera_id),
            "pos": camera_centers[fkey],
            "rot": init_camera_rotation(
                camera_centers=camera_centers[fkey], focal_point=focal_point,
            ),
            "obj": camera_ref,
            "mtl": camera_mtl_to_boundary_layers_map,
        }
        frame["children"].append(child)
        camera_id += 1

    return scenery
"""
