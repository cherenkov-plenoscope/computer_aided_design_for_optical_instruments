"""
Schmidt Camera
"""
import numpy as np
import json_utils
import merlict
import optic_object_wavefronts as optcad
import computer_aided_design_for_optical_instruments as cadoptins
import triangle_mesh_io
from optic_object_wavefronts import plot
import sebastians_matplotlib_addons as sebplt
import thin_lens


n_pmma_at_400nm = 1.505
n_air_2000m_asl = 1.000215821189739
mirror_reflectivity_at_400nm = 0.9
typical_distance_to_airshower = 99e3

config = {}
config["object_distance"] = typical_distance_to_airshower

# too low 0.94
# too big 1.05

config["fatchfactor"] = 0.96
config["refractive_index_corrector"] = n_pmma_at_400nm
config["corrector_mounting_radius"] = 600e-3
config["corrector_aperture_radius"] = 570e-3
config["corrector_width"] = 20e-3
config["corrector_stop_gap"] = 1e-3
config["focal_ratio"] = 2 / 3
config["focal_length"] = (2 * config["corrector_aperture_radius"]) * config[
    "focal_ratio"
]
config["relative_focus_delta"] = 1.25
config["mirror_aperture_radius"] = 750e-3
config["mirror_radius_of_curvature"] = 2.0 * config["focal_length"]
config["mirror_facet_inner_hex_radius"] = 120e-3
config["mirror_facet_body_width"] = 30e-3
config["mirror_gap_between_facets"] = 5e-3
config["screen_field_of_view_half_angle_rad"] = np.deg2rad(5.5)
config["screen_outer_radius"] = config["focal_length"] * np.tan(
    config["screen_field_of_view_half_angle_rad"]
)
config["screen_radius_of_curvature"] = (
    config["focal_length"] * config["fatchfactor"]
)
config["screen_distance"] = config["fatchfactor"] * config["focal_length"]

"""
thin_lens.compute_image_distance_for_object_distance(
    focal_length=config["focal_length"],
    object_distance=config["object_distance"],
)
"""
config["camera_body_radius"] = 2.0 * config["screen_outer_radius"]
config["camera_body_length"] = 2.0 * config["camera_body_radius"]


scenery = merlict.scenery.init()


scenery["materials"]["default_medium"] = "air"

# =========
# MATERIALS
# =========

# -----
# media
# -----
scenery["materials"]["media"]["air"] = {
    "refraction": [
        [200e-9, n_air_2000m_asl],
        [1.2e-6, n_air_2000m_asl],
    ],
    "absorbtion": [[200e-9, 1e9], [1.2e-6, 1e9]],
}
scenery["materials"]["media"]["pmma"] = {
    "refraction": [[200e-9, n_pmma_at_400nm], [1.2e-6, n_pmma_at_400nm]],
    "absorbtion": [[200e-9, 1e9], [1.2e-6, 1e9]],
}
scenery["materials"]["media"]["absorber"] = {
    "refraction": [[200e-9, 1], [1.2e-6, 1]],
    "absorbtion": [[200e-9, 0], [1.2e-6, 0]],
}

# --------
# surfaces
# --------
scenery["materials"]["surfaces"]["polished_aluminium_air"] = {
    "material": "Phong",
    "specular_reflection": [
        [200e-9, mirror_reflectivity_at_400nm],
        [1.2e-6, mirror_reflectivity_at_400nm],
    ],
    "diffuse_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "color": [255, 255, 255],
}
scenery["materials"]["surfaces"]["polished_pmma_air"] = {
    "material": "transparent",
    "specular_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "diffuse_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "color": [22, 91, 149],
}
scenery["materials"]["surfaces"]["dull_pmma_air"] = {
    "material": "Phong",
    "specular_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "diffuse_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "color": [222, 191, 149],
}
scenery["materials"]["surfaces"]["polished_silicon_air"] = {
    "material": "Phong",
    "specular_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "diffuse_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "color": [220, 180, 30],
}
scenery["materials"]["surfaces"]["absorber"] = {
    "material": "Phong",
    "specular_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "diffuse_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "color": [24, 24, 24],
}


# ---------------
# boundary layers
# ---------------
scenery["materials"]["boundary_layers"]["mirror_air"] = {
    "inner": {"medium": "air", "surface": "absorber"},
    "outer": {"medium": "air", "surface": "polished_aluminium_air"},
}
scenery["materials"]["boundary_layers"]["absorber_air"] = {
    "inner": {"medium": "absorber", "surface": "absorber"},
    "outer": {"medium": "air", "surface": "absorber"},
}
scenery["materials"]["boundary_layers"]["silicon_air"] = {
    "inner": {"medium": "air", "surface": "polished_silicon_air"},
    "outer": {"medium": "absorber", "surface": "absorber"},
}
scenery["materials"]["boundary_layers"]["pmma_air"] = {
    "inner": {"medium": "pmma", "surface": "polished_pmma_air"},
    "outer": {"medium": "air", "surface": "polished_pmma_air"},
}
scenery["materials"]["boundary_layers"]["dull_pmma_air"] = {
    "inner": {"medium": "pmma", "surface": "dull_pmma_air"},
    "outer": {"medium": "air", "surface": "dull_pmma_air"},
}

# ========
# GEOMETRY
# ========

# ------
# mirror
# ------
frame = {
    "id": 3000,
    "pos": [0, 0, -config["mirror_radius_of_curvature"]],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
    "children": [],
}
scenery["geometry"]["relations"]["children"].append(frame)

mirror_outer_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=config["mirror_aperture_radius"],
    fn=360,
    rot=np.pi / 6,
)
mirror_objs = cadoptins.segmented_mirror.add_to_frame(
    frame=frame,
    focal_length=config["focal_length"],
    aperture_outer_polygon=mirror_outer_polygon,
    aperture_inner_polygon=None,
    facet_inner_hex_radius=config["mirror_facet_inner_hex_radius"],
    gap_between_facets=config["mirror_gap_between_facets"],
    facet_body_width=config["mirror_facet_body_width"],
    davies_cotton_weight=0.0,
    parabola_weight=0.0,
    sphere_weight=1.0,
    facet_rotation="sphere",
    mean_distance_of_facet_centers_to_focal_point_is_focal_length=False,
    boundary_layer_facet_front="mirror_air",
    boundary_layer_facet_body="absorber_air",
    ref="mirror",
    facet_fn=6,
    facet_id_start=3100,
)
scenery["geometry"]["objects"].update(mirror_objs)

# ------
# screen
# ------
screen_mesh = optcad.primitives.spherical_cap_regular.init(
    outer_radius=config["screen_outer_radius"],
    curvature_radius=config["screen_radius_of_curvature"],
    inner_radius=None,
    fn_polygon=17,
    fn_hex_grid=10,
    ref="screen",
    rot=0.0,
)
screen_z = (
    -1.0 * config["mirror_radius_of_curvature"] + config["screen_distance"]
)
screen_obj = optcad.export.reduce_mesh_to_obj(screen_mesh)
screen_frame = {
    "id": 1,
    "pos": [0, 0, screen_z],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
    "obj": "screen",
    "mtl": {
        "screen": "silicon_air",
    },
}

scenery["geometry"]["relations"]["children"].append(screen_frame)
scenery["geometry"]["objects"]["screen"] = screen_obj

# ---------
# corrector
# ---------
corrector_outer_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=config["corrector_mounting_radius"],
    fn=90,
    rot=np.pi / 6,
    ref="outer_bound",
)
corrector_inner_polygon = None
"""
optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=(60e-3),
    fn=90,
    rot=np.pi / 6,
    ref="inner_bound",
)
"""

config[
    "schmidt_corrector_curvature_config"
] = optcad.geometry.schmidt_corrector.init_schmidt_corrector_curvature_config(
    corrector_aperture_radius=config["corrector_aperture_radius"],
    mirror_radius_of_curvature=-1.0 * config["mirror_radius_of_curvature"],
    refractive_index_corrector=config["refractive_index_corrector"],
    refractive_index_surrounding=n_air_2000m_asl,
    relative_focus_delta=config["relative_focus_delta"],
)


corrector_mesh = optcad.primitives.schmidt_corrector_plate.init(
    outer_polygon=corrector_outer_polygon,
    inner_polygon=corrector_inner_polygon,
    schmidt_corrector_curvature_config=config[
        "schmidt_corrector_curvature_config"
    ],
    offset=config["corrector_width"],
    fn_hex_grid=24,
    ref="corrector",
)

fig = sebplt.figure()
ax = sebplt.add_axes(fig, [0.1, 0.1, 0.8, 0.8])
ax.set_aspect("equal")
optcad.plot.ax_add_mesh_xy(
    ax=ax,
    mesh=corrector_mesh,
    vertex_color="red",
    vertex_marker="o",
    vertex_marker_size=0.3,
)
fig.savefig("corrector_mesh.jpg")

corrector_obj = optcad.export.reduce_mesh_to_obj(corrector_mesh)
with open("schmidt_corrector_plate.obj", "wt") as f:
    f.write(triangle_mesh_io.obj.dumps(corrector_obj))

fig = sebplt.figure()
ax = sebplt.add_axes(fig, [0.1, 0.1, 0.8, 0.8])
rr = np.linspace(0.0, 1.1 * config["corrector_aperture_radius"], 101)
hh = optcad.geometry.schmidt_corrector.surface_height(
    x=rr, y=0, **config["schmidt_corrector_curvature_config"]
)
ax.plot(
    rr * 1e3,
    hh * 1e3,
    "k-",
)
ax.plot(
    [
        config["corrector_aperture_radius"] * 1e3,
        config["corrector_aperture_radius"] * 1e3,
    ],
    [min(hh) * 1e3, max(hh) * 1e3],
    "k-",
    alpha=0.5,
)
ax.set_xlabel("radial / mm")
ax.set_ylabel("axial / mm")

ax.set_aspect(10)
fig.savefig("schmidt_corrector_plate_figure.jpg")


corrector_frame = {
    "id": 101,
    "pos": [0, 0, 0],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 180, 0]},
    "obj": "corrector",
    "mtl": {
        "corrector/bot": "pmma_air",
        "corrector/top": "pmma_air",
        "corrector/outer": "dull_pmma_air",
        "corrector/inner": "dull_pmma_air",
    },
}

scenery["geometry"]["relations"]["children"].append(corrector_frame)
scenery["geometry"]["objects"]["corrector"] = corrector_obj

# -------------
# aperture stop
# -------------
stop_outer_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=2.0 * config["corrector_aperture_radius"],
    fn=90,
    rot=0,
    ref="outer_bound",
)
stop_inner_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=config["corrector_aperture_radius"],
    fn=90,
    rot=0,
    ref="inner_bound",
)
stop_mesh = optcad.primitives.plane.init(
    outer_polygon=stop_outer_polygon,
    inner_polygon=stop_inner_polygon,
    fn_hex_grid=6,
    ref="stop",
)
stop_obj = optcad.export.reduce_mesh_to_obj(stop_mesh)
stop_frame = {
    "id": 102,
    "pos": [0, 0, config["corrector_width"] + config["corrector_stop_gap"]],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
    "obj": "stop",
    "mtl": {
        "stop": "absorber_air",
    },
}
scenery["geometry"]["relations"]["children"].append(stop_frame)
scenery["geometry"]["objects"]["stop"] = stop_obj


fig = sebplt.figure()
ax = sebplt.add_axes(fig, [0.1, 0.1, 0.8, 0.8])
ax.set_aspect("equal")
optcad.plot.ax_add_mesh_xy(
    ax=ax,
    mesh=stop_mesh,
    vertex_color="red",
    vertex_marker="o",
    vertex_marker_size=0.3,
)
fig.savefig("stop.jpg")

# -----------
# camera body
# -----------
camera_mesh = optcad.primitives.cylinder.init(
    outer_radius=config["camera_body_radius"],
    length=config["camera_body_length"],
    fn=24,
    rot=0.0,
    ref="camera",
)
camera_obj = optcad.export.reduce_mesh_to_obj(camera_mesh)
camera_frame = {
    "id": 103,
    "pos": [0, 0, screen_z + 0.1 * config["screen_outer_radius"]],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
    "obj": "camera",
    "mtl": {
        "camera/top": "absorber_air",
        "camera/outer": "absorber_air",
        "camera/bot": "absorber_air",
    },
}
scenery["geometry"]["relations"]["children"].append(camera_frame)
scenery["geometry"]["objects"]["camera"] = camera_obj


merlict.scenery.write_tar(sceneryPy=scenery, path="SchmidtOne.tar")

with open("SchmidtOne.json", "wt") as f:
    f.write(json_utils.dumps(config, indent=4))
