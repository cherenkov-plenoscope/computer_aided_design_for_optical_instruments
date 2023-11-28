"""
Create segmented mirrors from parameters.
"""
import numpy as np
import scipy
from scipy import spatial
import merlict
import optic_object_wavefronts as optcad
import optical_instruments_for_merlict as optmli
import triangle_mesh_io
from optic_object_wavefronts import plot
import sebastians_matplotlib_addons as sebplt


scenery = merlict.scenery.init(default_medium="vacuum")

frame = {
    "id": 0,
    "pos": [0, 0, 0],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
    "children": [],
}

scenery["geometry"]["relations"]["children"].append(frame)


scenery["materials"]["surfaces"][
    "perfect_mirror"
] = merlict.materials.surfaces.init("perfect_absorber/rgb_210_210_210")
scenery["materials"]["surfaces"][
    "perfect_absorber"
] = merlict.materials.surfaces.init("perfect_absorber/rgb_30_30_30")


aperture_outer_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=355e-3,
    fn=6,
    rot=np.pi / 6,
)

scenery = optmli.segmented_mirror.add_to_frame_in_scenery(
    frame=frame,
    scenery=scenery,
    focal_length=2000e-3,
    aperture_outer_polygon=aperture_outer_polygon,
    aperture_inner_polygon=None,
    facet_inner_hex_radius=54e-3,
    gap_between_facets=4e-3,
    davies_cotton_weight=0.0,
    parabola_weight=0.0,
    sphere_weight=1.0,
    mean_distance_of_facet_centers_to_focal_point_is_focal_length=False,
    facet_rotation="sphere",
    outer_medium="vacuum",
    inner_medium="vacuum",
    facet_surface_mirror="perfect_mirror",
    facet_surface_body="perfect_absorber",
    ref="a",
    facet_fn=3,
    facet_body_width=15e-3,
)

for objkey in scenery["geometry"]["objects"]:
    scenery["geometry"]["objects"][objkey] = optcad.export.reduce_mesh_to_obj(
        scenery["geometry"]["objects"][objkey]
    )

merlict.scenery.write_tar(sceneryPy=scenery, path="segmir.tar")


camera_geometry = optcad.primitives.light_field_eye.make_geometry(
    housing_outer_radius=0.05,
    housing_wall_width=2e-3,
    housing_height=0.01,
    lens_curvature_radius=0.125,
    lens_fn=3,
    photosensor_num_on_diagonal=9,
    photosensor_gap=0.0e-3,
    photosensor_plane_distance=0.12,
)
camera = optcad.primitives.light_field_eye.init(
    camera_geometry=camera_geometry,
    ref="cam",
)
camera_obj = optcad.export.reduce_mesh_to_obj(camera)

with open("eye.obj", "wt") as f:
    f.write(triangle_mesh_io.obj.dumps(camera_obj))


field_of_view_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=np.deg2rad(2.5 / 2 - 0.06667 / 2),
    fn=5,
)

lfs = optmli.light_field_camera.make_geometry(
    primary_optics_focal_length=106.5,
    primary_optics_diameter=71,
    camera_field_of_view_polygon=field_of_view_polygon,
    eye_field_of_view_full_angle=np.deg2rad(0.06667),
    eyes_point_towards_center_of_primary_optics=True,
    body_overhead=0.25,
)

fig = sebplt.figure(style={"rows": 1440, "cols": 2560, "fontsize": 1})
ax = sebplt.add_axes(fig=fig, span=[0.1, 0.1, 0.8, 0.8])
optcad.plot.ax_add_mesh_xy(ax=ax, mesh=lfs["body"]["mesh"])
optcad.plot.ax_add_polygon(
    ax=ax,
    polygon=lfs["screen"]["polygon"],
    closed=True,
    linestyle="-",
    color="r",
    alpha=0.3,
)
ax.set_aspect("equal")
fig.savefig("body.png")
sebplt.close(fig)


with open("body.obj", "wt") as f:
    mesh = lfs["body"]["mesh"]
    mesh = optcad.mesh.remove_unused_vertices_and_vertex_normals(mesh)
    f.write(triangle_mesh_io.obj.dumps(optcad.export.reduce_mesh_to_obj(mesh)))

print("Plot")
fig = sebplt.figure()
ax = sebplt.add_axes(fig=fig, span=[0.1, 0.1, 0.8, 0.8])

for eye in lfs["screen"]["eyes"]:
    ax.plot(
        eye["pos"][0],
        eye["pos"][1],
        ",r",
    )

optcad.plot.ax_add_polygon(
    ax=ax,
    polygon=lfs["screen"]["polygon"],
    closed=True,
    linestyle="-",
    color="r",
)

ax.set_aspect("equal")
fig.savefig("face_verts.jpg")
sebplt.close(fig)

lfc_frame = {
    "id": 133742,
    "pos": [0, 0, 5],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
    "children": [],
}
scenery["geometry"]["relations"]["children"].append(lfc_frame)

scenery = optmli.light_field_camera.add_to_frame_in_scenery(
    frame=lfc_frame,
    scenery=scenery,
    light_field_camera_geometry=lfs,
)
