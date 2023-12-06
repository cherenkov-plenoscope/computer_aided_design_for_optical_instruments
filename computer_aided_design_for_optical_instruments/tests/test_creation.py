"""
Create segmented mirrors from parameters.
"""
import numpy as np
import scipy
from scipy import spatial
import merlict
import optic_object_wavefronts as optcad
import computer_aided_design_for_optical_instruments as cadoptins
import triangle_mesh_io
from optic_object_wavefronts import plot
import sebastians_matplotlib_addons as sebplt


FOCAL_LENGTH = 106.5

scenery = merlict.scenery.init(default_medium="vacuum")

frame = {
    "id": 0,
    "pos": [0, 0, 0],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
    "children": [],
}

scenery["geometry"]["relations"]["children"].append(frame)


aperture_outer_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=0.5 * (FOCAL_LENGTH / 1.5) * (np.sqrt(3.0) / 2.0),
    fn=6,
    rot=np.pi / 6,
)

mirror_objs = cadoptins.segmented_mirror.add_to_frame(
    frame=frame,
    focal_length=FOCAL_LENGTH,
    aperture_outer_polygon=aperture_outer_polygon,
    aperture_inner_polygon=None,
    facet_inner_hex_radius=800e-3,
    gap_between_facets=20e-3,
    facet_body_width=25e-3,
    davies_cotton_weight=0.0,
    parabola_weight=1.0,
    sphere_weight=0.0,
    facet_rotation="sphere",
    mean_distance_of_facet_centers_to_focal_point_is_focal_length=False,
    boundary_layer_facet_front="facet/front",
    boundary_layer_facet_body="facet/body",
    ref="mirror",
    facet_fn=6,
)
scenery["geometry"]["objects"].update(mirror_objs)


scenery["materials"]["surfaces"][
    "perfect_mirror"
] = merlict.materials.surfaces.init("perfect_absorber/rgb_210_210_210")
scenery["materials"]["surfaces"][
    "perfect_absorber"
] = merlict.materials.surfaces.init("perfect_absorber/rgb_30_30_30")
scenery["materials"]["boundary_layers"]["mirror/facet/front"] = {
    "inner": {"medium": "vacuum", "surface": "perfect_absorber"},
    "outer": {"medium": "vacuum", "surface": "perfect_mirror"},
}
scenery["materials"]["boundary_layers"]["mirror/facet/body"] = {
    "inner": {"medium": "vacuum", "surface": "perfect_absorber"},
    "outer": {"medium": "vacuum", "surface": "perfect_absorber"},
}

merlict.scenery.write_tar(sceneryPy=scenery, path="segmir.tar")

screen_polygon = optcad.geometry.regular_polygon.make_vertices_xy(
    outer_radius=FOCAL_LENGTH * np.deg2rad(6.5 / 2 - 0.06667 / 2),
    fn=5,
)

lfc_geometry = cadoptins.light_field_camera.make_geometry(
    principal_aperture_focal_length=FOCAL_LENGTH,
    principal_aperture_diameter=71,
    screen_polygon=screen_polygon,
    screen_curvature_radius=-FOCAL_LENGTH,
    eyes_spacing=FOCAL_LENGTH * np.deg2rad(0.0667),
    eyes_facing_point=[0, 0, -FOCAL_LENGTH],
    eyes_photosensor_num_on_diagonal=5,
    eyes_photosensor_gap=1e-3,
    eyes_lens_curvature_radius=1.0,
    eyes_lens_gap=2e-3,
    eyes_lens_fn=7,
    eyes_housing_wall_width=2e-3,
    body_overhead_radius=0.33,
)

fig = sebplt.figure(style={"rows": 1440, "cols": 2560, "fontsize": 1})
ax = sebplt.add_axes(fig=fig, span=[0.1, 0.1, 0.8, 0.8])
optcad.plot.ax_add_mesh_xy(ax=ax, mesh=lfc_geometry["housing"]["hull"]["mesh"])
optcad.plot.ax_add_polygon(
    ax=ax,
    polygon=lfc_geometry["screen"]["hull"]["polygon"],
    closed=True,
    linestyle="-",
    color="r",
    alpha=0.3,
)
ax.set_aspect("equal")
fig.savefig("housing.png")
sebplt.close(fig)


with open("housing.obj", "wt") as f:
    mesh = lfc_geometry["housing"]["hull"]["mesh"]
    mesh = optcad.mesh.remove_unused_vertices_and_vertex_normals(mesh)
    f.write(triangle_mesh_io.obj.dumps(optcad.export.reduce_mesh_to_obj(mesh)))

print("Plot")
fig = sebplt.figure()
ax = sebplt.add_axes(fig=fig, span=[0.1, 0.1, 0.8, 0.8])

for eye in lfc_geometry["screen"]["eyes"]:
    ax.plot(
        eye["pos"][0],
        eye["pos"][1],
        ",r",
    )

optcad.plot.ax_add_polygon(
    ax=ax,
    polygon=lfc_geometry["screen"]["hull"]["polygon"],
    closed=True,
    linestyle="-",
    color="r",
)

ax.set_aspect("equal")
fig.savefig("face_verts.jpg")
sebplt.close(fig)


lfc_frame = {
    "id": 1337,
    "pos": [0, 0, FOCAL_LENGTH],
    "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
    "children": [],
}
scenery["geometry"]["relations"]["children"].append(lfc_frame)

lfc_objs = cadoptins.light_field_camera.add_to_frame(
    frame=lfc_frame,
    light_field_camera_geometry=lfc_geometry,
)
scenery["geometry"]["objects"].update(lfc_objs)


scenery["materials"]["surfaces"]["glass"] = {
    "material": "transparent",
    "specular_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "diffuse_reflection": [[200e-9, 0.0], [1.2e-6, 0.0]],
    "color": [200, 200, 200],
}
scenery["materials"]["media"]["schott_n_kzfs2"] = {
    "refraction": [
        [334.1e-9, 1.59259],
        [365.0e-9, 1.58382],
        [404.7e-9, 1.57580],
        [435.8e-9, 1.57114],
        [480.0e-9, 1.56612],
        [486.1e-9, 1.56553],
        [587.6e-9, 1.55836],
        [589.3e-9, 1.55827],
        [632.8e-9, 1.55617],
        [643.8e-9, 1.55570],
        [656.3e-9, 1.55519],
        [706.5e-9, 1.55337],
        [852.1e-9, 1.54944],
        [1014.0e-9, 1.54625],
        [1060.0e-9, 1.54546],
        [1529.6e-9, 1.53798],
        [1970.1e-9, 1.53011],
        [2325.4e-9, 1.52239],
    ],
    "absorbtion": [[200e-9, 1e99], [1.2e-6, 1e99]],
}

scenery["materials"]["boundary_layers"]["light_field_camera/eye/lens"] = {
    "inner": {"medium": "schott_n_kzfs2", "surface": "glass"},
    "outer": {"medium": "vacuum", "surface": "glass"},
}
scenery["materials"]["boundary_layers"][
    "light_field_camera/eye/housing/inner"
] = {
    "inner": {"medium": "vacuum", "surface": "perfect_mirror"},
    "outer": {"medium": "vacuum", "surface": "perfect_absorber"},
}
scenery["materials"]["boundary_layers"][
    "light_field_camera/eye/housing/outer"
] = {
    "inner": {"medium": "vacuum", "surface": "perfect_absorber"},
    "outer": {"medium": "vacuum", "surface": "perfect_absorber"},
}
scenery["materials"]["boundary_layers"][
    "light_field_camera/eye/photosensor"
] = {
    "inner": {"medium": "vacuum", "surface": "perfect_absorber"},
    "outer": {"medium": "vacuum", "surface": "perfect_absorber"},
}
scenery["materials"]["boundary_layers"]["light_field_camera/housing"] = {
    "inner": {"medium": "vacuum", "surface": "perfect_absorber"},
    "outer": {"medium": "vacuum", "surface": "perfect_absorber"},
}
merlict.scenery.write_tar(sceneryPy=scenery, path="plenoscope.tar")
