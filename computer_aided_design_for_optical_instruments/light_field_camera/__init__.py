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


def _init_principal_aperture(focal_length, inner_diameter):
    """
    Parameters
    ----------
    focal_length : float
        The focal-length of the principal aperture which is infornt
        (which comes before) the light-field camera. In case of the
        Cherenkov-plenoscope this is related to the mirror.
    inner_diameter : float
        The largest diameter of a ring to fit inside the aperture. In case of
        a hexagonal mirror, this is the inner diameter of the hexagon.

    Returns
    -------
    principal_aperture : dict
        Properties of the principal aperture plane w.r.t. the frame of the
        light-field camera.
    """
    assert focal_length > 0.0
    assert inner_diameter > 0.0
    c = {}
    c["focal_length"] = float(focal_length)
    c["inner_diameter"] = float(inner_diameter)
    c["inner_radius"] = 0.5 * c["inner_diameter"]
    c["fstop"] = c["focal_length"] / c["inner_diameter"]
    c["origin"] = np.array([0.0, 0.0, -c["focal_length"]])
    c["optical_axis"] = np.array([0.0, 0.0, 1.0])
    return c


def _init_screen(polygon, eyes_spacing, curvature_radius, eyes_facing_point):
    """
    Parameters
    ----------
    polygon : dict str -> 3D np.array
        Inside this polygon, eyes will be placed.
    eyes_spacing : float
        The spacing of the eyes in a hexagonal grid. This is the inner diameter
        of the hexagons in the grid.
    curvature_radius : float
        The screen's radius of curvature.
    eyes_facing_point : 3D np.array
        The position the eyes are facing to.
    """
    assert eyes_spacing > 0.0
    eyes_facing_point = np.array(eyes_facing_point)

    c = {}
    c["curvature_radius"] = float(curvature_radius)
    c["eyes_centers_limit"] = {}
    c["eyes_centers_limit"]["polygon"] = polygon

    (
        c["eyes_centers_limit"]["xlim"],
        c["eyes_centers_limit"]["ylim"],
        _,
    ) = optcad.polygon.limits(polygon=c["eyes_centers_limit"]["polygon"])

    c["eyes_grid"] = {}
    c["eyes_grid"]["spacing"] = eyes_spacing
    _xout = np.max(np.abs(c["eyes_centers_limit"]["xlim"]))
    _yout = np.max(np.abs(c["eyes_centers_limit"]["ylim"]))
    c["eyes_grid"][
        "radius_to_safely_include_all_potential_eyes"
    ] = 2 * np.hypot(_xout, _yout)

    c["eyes_grid"]["fn"] = int(
        np.ceil(
            c["eyes_grid"]["radius_to_safely_include_all_potential_eyes"]
            / c["eyes_grid"]["spacing"]
        )
    )
    _eyes_centers = optcad.geometry.grid.hexagonal.init_from_spacing(
        spacing=c["eyes_grid"]["spacing"],
        ref="",
        fN=c["eyes_grid"]["fn"],
    )
    alpha = np.deg2rad(30)
    sinA = np.sin(alpha)
    cosA = np.cos(alpha)
    for key in _eyes_centers:
        x, y, z = _eyes_centers[key]
        npos = np.array([cosA * x - sinA * y, sinA * x + cosA * y, 0.0])
        _eyes_centers[key] = npos

    c["eyes_grid"]["centers"] = optcad.polygon.get_vertices_inside(
        vertices=_eyes_centers,
        polygon=c["eyes_centers_limit"]["polygon"],
    )

    _eyes_voronoi_cells = (
        optcad.geometry.grid.hexagonal.init_voronoi_cells_from_centers(
            centers=c["eyes_grid"]["centers"],
            centers_spacing=c["eyes_grid"]["spacing"],
            rot=alpha,
        )
    )
    c["hull"] = {}
    c["hull"][
        "polygon"
    ] = optcad.geometry.grid.hexagonal.find_hull_of_voronoi_cells(
        voronoi_cells=_eyes_voronoi_cells,
        centers=c["eyes_grid"]["centers"],
        centers_spacing=c["eyes_grid"]["spacing"],
    )

    c["eyes"] = []
    for conid, hexid in enumerate(c["eyes_grid"]["centers"]):
        eye = {}
        eye["continuous_eye_id"] = conid
        eye["hexagonal_grid_id"] = hexid

        eye["pos"] = c["eyes_grid"]["centers"][hexid]

        eye["pos"][2] = optcad.geometry.sphere.surface_height(
            x=eye["pos"][0],
            y=eye["pos"][1],
            curvature_radius=c["curvature_radius"],
        )

        _axis, _angle = segmented_mirror.facet_rotation_axis_and_angle(
            facet_center=eye["pos"],
            target_point=eyes_facing_point,
            direction_incoming_light=np.array([0.0, 0.0, 1.0]),
        )
        if np.abs(_angle) > 0.0:
            eye["rot"] = {
                "repr": "axis_angle",
                "axis": _axis,
                "angle_deg": np.rad2deg(_angle),
            }
        else:
            eye["rot"] = {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]}

        c["eyes"].append(eye)

    return c


def _init_housing(
    screen_eyes_centers_limit_polygon,
    screen_hull_polygon,
    body_overhead_radius,
    curvature_radius,
    vs_hex_grid,
    min_hull_fn=36,
):
    assert body_overhead_radius >= 0.0
    assert min_hull_fn >= 3

    c = {}
    c["hull"] = {}
    c["hull"]["minkowski_regular_polygon_fn"] = max(
        [min_hull_fn, len(screen_hull_polygon)]
    )

    _circle = optcad.geometry.regular_polygon.make_vertices_xy(
        outer_radius=body_overhead_radius,
        fn=c["hull"]["minkowski_regular_polygon_fn"],
    )

    c["hull"]["polygon"] = optcad.minkowski.minkowski_hull_xy(
        poly1=screen_eyes_centers_limit_polygon,
        poly2=_circle,
        ref="body/outer_bound",
    )

    c["hull"]["mesh"] = optcad.primitives.template_curved_surface.init(
        outer_polygon=c["hull"]["polygon"],
        curvature_config={"curvature_radius": curvature_radius},
        curvature_height_function=optcad.geometry.sphere.surface_height,
        curvature_surface_normal_function=optcad.geometry.sphere.surface_normal,
        inner_polygon=screen_hull_polygon,
        vs_hex_grid=vs_hex_grid,
        fn_hex_grid=None,
        ref="housing",
    )
    return c


def _init_eye(
    eyes_spacing,
    lens_fill_factor,
    photosensor_fill_factor,
    lens_curvature_radius,
    photosensor_num_on_diagonal,
):
    c = {}
    c["geometry"] = optcad.primitives.light_field_eye.make_geometry(
        housing_outer_radius=0.05,
        housing_wall_width=2e-3,
        housing_height=0.01,
        lens_curvature_radius=lens_curvature_radius,
        lens_fn=3,
        photosensor_num_on_diagonal=photosensor_num_on_diagonal,
        photosensor_gap=0.0e-3,
        photosensor_plane_distance=0.12,
    )
    return c


def make_geometry(
    principal_aperture_focal_length,
    principal_aperture_diameter,
    screen_polygon,
    screen_curvature_radius,
    eyes_spacing,
    eyes_facing_point,
    eyes_photosensor_num_on_diagonal,
    eyes_photosensor_gap,
    eyes_lens_curvature_radius,
    eyes_lens_gap,
    eyes_lens_fn,
    eyes_housing_wall_width,
    body_overhead_radius,
):
    """
    principal_aperture_focal_length : float
        The expected focal-length of the primary optics (probably a mirror).
    principal_aperture_diameter : float
        The expected diameter of the primary optics.
    screen_polygon : OrderedDict[(str, [float,float,float])
        The outer perimeter of the screen. Only eyes which fit into it will
        be added to the camera.
    screen_curvature_radius : float
        The screen's radius of curvature.
    eyes_spacing : float
        Spacing of eyes in the x-y plane of the screen.
    eyes_facing_point : 3D np.array
        Eyes will face to this point. It effects the orientation of the eyes.
    eyes_lens_gap : float
        Gap in between lenses and housings of neighboring eyes.
    eyes_photosensor_num_on_diagonal:
        This many photosensor will be on the diagonal of an individual eye.
    eyes_photosensor_gap : float
        Gap in between neighboring photosensors inside the eye.
    eyes_housing_wall_width : float
        Width of the wall surrounding each eye.
    eyes_lens_curvature_radius : float
        Must match the lens' refractive index.
    body_overhead_radius : float
        Radius.
    """
    assert principal_aperture_focal_length > 0.0
    assert principal_aperture_diameter > 0.0
    assert len(screen_polygon) >= 3
    assert eyes_spacing > 0.0
    assert body_overhead_radius > 0.0
    HEXIN2OUT = 2.0 / np.sqrt(3.0)
    HEXOUT2IN = 1.0 / HEXIN2OUT

    c = {}
    c["principal_aperture"] = _init_principal_aperture(
        focal_length=principal_aperture_focal_length,
        inner_diameter=principal_aperture_diameter,
    )

    c["screen"] = _init_screen(
        polygon=screen_polygon,
        eyes_spacing=eyes_spacing,
        curvature_radius=screen_curvature_radius,
        eyes_facing_point=eyes_facing_point,
    )

    c["housing"] = _init_housing(
        screen_eyes_centers_limit_polygon=c["screen"]["eyes_centers_limit"][
            "polygon"
        ],
        screen_hull_polygon=c["screen"]["hull"]["polygon"],
        body_overhead_radius=body_overhead_radius,
        curvature_radius=screen_curvature_radius,
        vs_hex_grid=eyes_spacing / 2.0,
        min_hull_fn=36,
    )

    assert eyes_spacing > eyes_lens_gap
    eyes_housing_inner_radius = 0.5 * (eyes_spacing - eyes_lens_gap)
    eyes_housing_outer_radius = eyes_housing_inner_radius * HEXIN2OUT
    assert eyes_lens_curvature_radius >= eyes_housing_outer_radius

    eyes_fstop = c["principal_aperture"]["fstop"]
    eyes_photosensor_plane_distance = eyes_housing_inner_radius * eyes_fstop
    eyes_housing_height = 0.5 * eyes_photosensor_plane_distance

    c["eye"] = {}
    c["eye"]["geometry"] = optcad.primitives.light_field_eye.make_geometry(
        housing_outer_radius=eyes_housing_outer_radius,
        housing_wall_width=eyes_housing_wall_width,
        housing_height=eyes_housing_height,
        lens_curvature_radius=eyes_lens_curvature_radius,
        lens_fn=eyes_lens_fn,
        photosensor_num_on_diagonal=eyes_photosensor_num_on_diagonal,
        photosensor_gap=eyes_photosensor_gap,
        photosensor_plane_distance=eyes_photosensor_plane_distance,
    )

    c["eye"]["mesh"] = optcad.primitives.light_field_eye.init(
        geometry=c["eye"]["geometry"],
        ref="eye",
    )

    return c


def _assert_int(a):
    assert a % 1.0 == 0.0
    return int(a)


def add_to_frame(
    frame,
    light_field_camera_geometry,
    boundary_layer_eye_lens="eye/lens",
    boundary_layer_eye_housing_inner="eye/housing/inner",
    boundary_layer_eye_housing_outer="eye/housing/outer",
    boundary_layer_eye_photosesnsor="eye/photosensor",
    boundary_layer_housing="housing",
    continuous_eye_id_offset=1000 * 1000 * 1000,
    housing_id_offset=1000,
    ref="light_field_camera",
):
    """
    Parameters
    ----------
    frame : dict
        This light-field camera will be made a child of this frame. The frame
        is itself a child of scenery["geometry"]["relations"].
    light_field_camera_geometry : dict
        The light-field camera's geometry. A dict with all the relevant
        positions and orientations pre calculated.
    """
    continuous_eye_id_offset = _assert_int(continuous_eye_id_offset)
    housing_id_offset = _assert_int(housing_id_offset)

    c = light_field_camera_geometry
    objs = {}

    # eyes
    # ----
    eye_obj_key = posixpath.join(ref, "eye")
    objs[eye_obj_key] = optcad.export.reduce_mesh_to_obj(c["eye"]["mesh"])

    eye_mtl = {}
    for mtl in c["eye"]["mesh"]["materials"]:
        eye_mtl[mtl] = boundary_layer_eye_photosesnsor

    eye_mtl["eye/lens/top"] = boundary_layer_eye_lens
    eye_mtl["eye/lens/side"] = boundary_layer_eye_housing_outer
    eye_mtl["eye/lens/bot"] = boundary_layer_eye_lens

    eye_mtl["eye/housing/inner"] = boundary_layer_eye_housing_inner
    eye_mtl["eye/housing/outer"] = boundary_layer_eye_housing_outer
    eye_mtl["eye/housing/bottom"] = boundary_layer_eye_housing_outer
    eye_mtl["eye/housing/top"] = boundary_layer_eye_housing_outer

    for mtl in eye_mtl:
        eye_mtl[mtl] = posixpath.join(ref, eye_mtl[mtl])

    for i in range(len(c["screen"]["eyes"])):
        eye = c["screen"]["eyes"][i]
        child = {
            "id": continuous_eye_id_offset + eye["continuous_eye_id"],
            "pos": eye["pos"],
            "rot": eye["rot"],
            "obj": eye_obj_key,
            "mtl": eye_mtl,
        }
        frame["children"].append(child)

    # body
    # ----
    housing_obj_key = posixpath.join(ref, "housing")
    objs[housing_obj_key] = optcad.export.reduce_mesh_to_obj(
        c["housing"]["hull"]["mesh"]
    )

    housing_mtl = {"housing": boundary_layer_housing}
    for mtl in housing_mtl:
        housing_mtl[mtl] = posixpath.join(ref, housing_mtl[mtl])

    child = {
        "id": housing_id_offset + 0,
        "pos": [0, 0, 0],
        "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
        "obj": housing_obj_key,
        "mtl": housing_mtl,
    }
    frame["children"].append(child)

    return objs
