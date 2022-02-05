import numpy as np
import polyscope as ps

src = '../petiole/coords/'
bname = 'col-0_comp_0'

filename = src + bname + '_verts_inv.csv'
vertices = np.loadtxt(filename, delimiter=',')

filename = src + bname + '_triangles.csv'
triangles = np.loadtxt(filename, delimiter=',')

ps.set_ground_plane_mode("none")
ps.set_navigation_style("turntable")
ps.set_up_dir("z_up")

ps.init()

# alternately:
ps.register_surface_mesh("leaf", vertices, triangles, enabled=True,
                         color=(0.08,.5, 0.0), edge_color=((0.8, 0.8, 0.8)),
                         edge_width=0.7, smooth_shade=True,
                         material='clay', transparency=1)

ps.set_transparency_mode('none')
ps.set_screenshot_extension(".jpg");
ps.show()
