#%%
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
#%%

def resample(curve, dx, nmin=21):
    dl = np.linalg.norm(np.diff(curve, axis=0), axis=1)
    t = np.zeros(len(dl)+1)
    t[1:] = np.cumsum(dl)
    l = t[-1]
    t /= l
    nnew = max(int(np.ceil(l/dx))+1, nmin)
    tnew = np.linspace(0, 1, nnew)
    newcurve = np.zeros((nnew, 2))
    newcurve[:,0] = np.interp(tnew, t, curve[:,0])
    newcurve[:,1] = np.interp(tnew, t, curve[:,1])
    return newcurve

def add_arrow(line, positions=[0.1, 0.5, 0.9], direction='right', size=10):
    color = line.get_color()
    xdata = line.get_xdata()
    ydata = line.get_ydata()
    
    for position in positions:
        start_ind = int(position * len(xdata))
        if direction == 'right':
            end_ind = start_ind + 1
        else:
            end_ind = start_ind - 1

        line.axes.annotate('',
            xytext=(xdata[start_ind], ydata[start_ind]),
            xy=((xdata[start_ind]+xdata[end_ind])/2, (ydata[start_ind]+ydata[end_ind])/2),
            arrowprops=dict(arrowstyle="-|>", color=color),
            size=size
        )

import gmsh
def make_mesh(curves, sizes, loops, resample_curves=True, visualize=True):
    print("Entering mesh generation")
    print("loops are made of the following curves:")
    for i, loop in enumerate(loops):
        print(f"    loop {i}: {loop}")
    print("Make sure to check the orientration of the loops in the visualization!!")

    if resample_curves:
        curves = [resample(c, s/2) for c, s in zip(curves, sizes)]
    else:
        curves = [np.array(c) for c in curves]
    all_curves = np.vstack(curves)
    all_points, iindex = np.unique(all_curves, axis=0, return_inverse=True)
    
    # visualize boundary curves for validation
    if visualize:
        for i, loop in enumerate(loops): 
            colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
            color = colors[i%len(colors)]
            for icurve in loop:
                curve = curves[icurve]
                line = plt.plot(curve[:,0], curve[:,1], color=color)[0]
                plt.plot(curve[[0,-1],0], curve[[0,-1],1], ".", color=color)
                pts = np.linspace(0, 1, max(min(30, len(curve)//12), 3))[1:-1]
                add_arrow(line, pts)
                t = plt.text(curve[2*len(curve)//5,0], curve[2*len(curve)//5,1], str(icurve), ha='center', va='center', color=color)
                t.set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0, boxstyle='round'))
        plt.tight_layout()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()


    idxcurves = []
    start = 0
    for curve in curves:
        idxcurves.append(iindex[start:start+len(curve)]+1)
        start += len(curve)

    gmsh.initialize()
    for i, (x, y) in enumerate(all_points):
        gmsh.model.geo.addPoint(x, y, 0, tag=i+1)
    
    for i, line in enumerate(idxcurves):
        gmsh.model.geo.add_polyline(line, i+1)
    
    looptags = []
    for i, l in enumerate(loops):
        looptags.append(gmsh.model.geo.addCurveLoop(np.array(l)+1, i+1))

    gmsh.model.geo.addPlaneSurface(looptags, 1)
    gmsh.model.geo.synchronize()

    for i, line in enumerate(idxcurves):
        gmsh.model.mesh.set_size_at_parametric_points(1, i+1, [0, 1],  [sizes[i], sizes[i]])


    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)

    # set size callback
    # def cb(dim, tag, x, y, z, lc):
    #     return ((y < 0.4) + 1.0) * lc
    # cb = lambda dim, tag, x, y, z, lc: lc*(1+ x > 0.5) * 100*1e-3
    # gmsh.model.mesh.setSizeCallback(cb)
    
    # triangles
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    
    # quads (pah hyper robuste)
    # gmsh.option.setNumber("Mesh.SaveAll", 1)
    # gmsh.option.setNumber("Mesh.RecombineAll", 1)
    # gmsh.option.setNumber("Mesh.Algorithm", 11)
    # gmsh.option.setNumber("Mesh.SmoothRatio", 21.5)
    # gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
    # gmsh.model.geo.mesh.setRecombine(2, 1, 45)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write("mesh.msh")
    
    # visualize
    gmsh.fltk.run()
    gmsh.finalize()


curves = np.load('curves.npy', allow_pickle=True)
if curves.ndim == 1:
    curves = list(curves)
else:
    curves = [curves[0].astype(float)]

# Split the first curve (contour) into 3 parts (base0, countour0, base1)
contour = curves[0]
base0 = contour[0:100+1]
countour0 = contour[100:-100+1]
base1 = contour[-100:]

# Replace the first curve with the 3 parts
curves = [base0[::-1], base1[::-1], countour0[::-1]] + curves[1:]

# Mesh size for each curve
sizes = [20e-3]*len(curves)

# Loops : each loop is a list of curve indices. The first loop is the contour, made of the first 3 curves. The other loops are made of a single curve.
loops = [[i] for i in range(len(curves))]
loops = [[0, 1, 2]] + loops[3:]

make_mesh(curves, sizes, loops)
# %%
