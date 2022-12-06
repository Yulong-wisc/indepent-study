import pybullet as p
#import pybullet_data as pd
import os
import wavefront as wf
import util
import optimization as opt
import numpy as np
import argparse, sys
import vtk
import meshplex
import trimesh
import math
parser = argparse.ArgumentParser(description=__doc__, formatter_class=
        argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--input_obj", "-i", type=str, dest="input_obj", default="mesh.obj", help="File name of the .obj input")
parser.add_argument("--output_csv", "-o", type=str, dest="output_csv", default="result.csv", help="File name of the .csv output (after optimization)")
parser.add_argument("--approx_csv", "-a", type=str, dest="approx_csv", default="quick_approx.csv", help="File name of a quick approximation of the decomposition (before optimization)")
parser.add_argument("--convex_obj", "-c", type=str, dest="convex_obj", default="parts.obj", help="File name of the intermediate .obj convex shapes")
parser.add_argument("--convex_log", type=str, dest="convex_log", default="log.txt", help="File name of the intermediate convex decomposition logs")
parser.add_argument("--interpolate_ratio", "-r", type=float, dest="ratio", default=0.5,
			help="Must be in [0,1]; set it close to 1 to make the final output spheres large, close to 0 to make them small")
parser.add_argument("--voxel_resolution", "-v", type=int, dest="voxel_resolution", default=50000,
			help="The resolution at which the convex decomposition takes place; larger number usually leads to more final spheres")
parser.add_argument("--max_hull_verts", type=int, dest="max_hull_verts", default=64, help="The max number of vertices of one convex hull")
parser.add_argument("--budget", "-b", type=int, dest="budget", default=5, help="The target number of spheres we want")
args = parser.parse_args(sys.argv[1:])

original_meshes = wf.load_obj(args.input_obj)
assert len(original_meshes) == 1, "This script handles OBJ with one mesh group only."
mesh = original_meshes[0]

p.connect(p.DIRECT)

############## IF OTHER CONVEX DECOMPOSITION PARAM NEEDED TO CHANGE: PLEASE DIRECTLY CHANGE THEM HERE #######################
p.vhacd(args.input_obj, args.convex_obj, args.convex_log, concavity=0.0025, alpha=0.04,
        resolution=args.voxel_resolution, maxNumVerticesPerCH=args.max_hull_verts)
#############################################################################################################################
parts = wf.load_obj(args.convex_obj)#, triangulate=True)

'''if (len(parts) == 1):
    parts = []
    parts.append(wf.load_obj(args.convex_obj)[0])

    planes = [[0,0,100],[0,0,80],[0,0,40],[0,0,0],[0,0,-40],[0,0,-80],[0,0,-100]]
    normals = [[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]]
    for i in range(len(planes)):
        plane = np.array(planes[i])
        normal = np.array(normals[i])
        cut_result = []
        for part in parts:
            halfA = util.sliceMesh(part, normal, plane)
            halfB = util.sliceMesh(part, -normal, plane)
            if len(halfA.vertices) > 0:
                cut_result.append(halfA)
            if len(halfB.vertices) > 0:
                cut_result.append(halfB)
        parts = cut_result'''

xyzr = np.zeros((len(parts), 4))
part_id = 0
print(xyzr)
#1. the ratio of length and width
box_length =[]
box_width = []
box_ratio=[]
box_diference=[]
volume_percentage = []


# 3.the percentage that the part takes refer to the whole object

whole_mesh = trimesh.load('parts.obj')
whole_volume=whole_mesh.volume
#part_volume=whole_mesh[0].volume

#print(whole_volume)
#print(part_volume)
#max_convex_num=input();
#while(len(parts)<=max_convex_num):

for part in parts:

    bounding_box = util.bbox(part.vertices)

    #1. get the length, width of the box and their ratio
    box_length.append( math.sqrt( (bounding_box[1,0]-bounding_box[0,0])**2 + (bounding_box[1,1]-bounding_box[0,1])**2 + (bounding_box[1,2]-bounding_box[0,2])**2)) #the difference of the z axis value
    a = [bounding_box[1,1]-bounding_box[0,1],bounding_box[1,0]-bounding_box[0,0]] #the difference of x axis and y axis
    box_width.append(min(a))
    box_ratio.append(box_length[part_id]/box_width[part_id])


    big_sphere = util.box2ball(bounding_box)

    mesh_center = util.coord_avg(part.vertices)

    small_sphere = util.inscribedSphereViaPointMesh(mesh_center, part)

    decomp_sphere = util.interpolateSphere(big_sphere, small_sphere, args.ratio)
    xyzr[part_id,:3] = decomp_sphere.center
    xyzr[part_id,-1] = decomp_sphere.radius

    # 2. the value of modification (suspected)
    box_diference.append((box_length[part_id] - box_width[part_id]) / decomp_sphere.radius)

    # 3.the percentage that the part takes refer to the whole object

    '''verts = np.array(part.vertices)
    faces = np.empty((len(part.polygons), 3), dtype=int)
    mesh = meshplex.MeshTri(np.array(verts), np.array(faces))
    
    V = mesh.volume
    
    print("Volume =", V)'''

    '''verts = np.array(part.vertices)
    faces = np.empty((len(part.polygons), 3), dtype=int)
    
    mesh = meshplex.MeshTri(verts, faces)
    
    V = np.sum(mesh.cell_volumes)
    
    part_volume = trimesh.visual'''

    #p_mesh = trimesh.load(args.convex_obj[0][0])
    #p_volume = p.volume

    filename = "my_temp_part" + str(part_id) + ".obj"

    wf.save_obj(part, filename)
    part_mesh = trimesh.load(filename)
    part_volume = part_mesh.volume
    #print(part_volume)
    volume_percentage.append(part_volume/whole_volume)
    part_id += 1
    pass





weight1 = 1/3
weight2 = 1/3
weight3 = 2 / 3

box_judge = []
for i in range(len(parts)):
    box_judge.append(weight1 * box_ratio[i] + weight2 * box_diference[i] + weight3 * volume_percentage[i])
# print(box_judge[i])
pass
##outputting critical convex
output_convex_id = box_judge.index(max(box_judge))
print("%d",output_convex_id)
out_put_convex = parts[output_convex_id+1]
print("%d",out_put_convex)






##cutting the critical convex into two parts
'''planes = [[-0.91073949, 10.05606754, 2.68031705]]
normals = [[-3.54744241, -2.447003, -0.12203631]]
cut_result = []
halfA = util.sliceMesh(out_put_convex, normals, planes)
halfB = util.sliceMesh(out_put_convex, -normals, planes)
if len(halfA.vertices) > 0:
    cut_result.append(halfA)
if len(halfB.vertices) > 0:
    cut_result.append(halfB)'''

#parts.remove(out_put_convex)
print("%d",len(parts))
#parts.append(cut_result)
#xyzr = np.zeros((len(parts), 4))
print("%d",len(parts))





##try
part_id=0;
for part in parts:

    bounding_box = util.bbox(part.vertices)
    big_sphere = util.box2ball(bounding_box)

    mesh_center = util.coord_avg(part.vertices)

    small_sphere = util.inscribedSphereViaPointMesh(mesh_center, part)

    decomp_sphere = util.interpolateSphere(big_sphere, small_sphere, args.ratio)
    xyzr[part_id, :3] = decomp_sphere.center
    xyzr[part_id, -1] = decomp_sphere.radius
    part_id += 1





#parts.append(input)


"""Max = box_ratio[0]
for i in range(len(box_ratio)):
    if box_ratio[i] > Maxnums:

        Maxnums = box_ratio[i]
        index_part = i"""

# index_part is the part that gives the max value of length/width
#print(Maxnums,index_part)

#2. the value of modification (suspected)

'''for part in parts:
    box_diference.append((box_length[part]-box_width[part])/decomp_sphere.radius)
    part_id += 1'''


'''#3. the percentage that the part takes refer to the whole object
part_id = 0
part_portion = []

for part in parts:
    part_portion.append(/)
    part_id +=1'''


#assign the weight to the three criteria
'''weight1= 2/3
weight2= 1/6
weight3= 1/6

box_judge= []
for i in range(len(parts)):
    box_judge.append(weight1*box_ratio[i]+weight2*box_diference[i]+weight3*volume_percentage[i])
    #print(box_judge[i])
    pass

print('the most critical mesh is convex %d' %(box_judge.index(max(box_judge))+1))'''




'''print("%d",len(args.convex_obj)) #parts.obj
print("%d",len(args.approx_csv)) #quickapproximation.cvs
print("%d",len(args.output_csv)) #result.cvs'''


np.savetxt(args.approx_csv, xyzr, header = "x,y,z,r", delimiter=",")

# The ball number which this vertex got assigned to.
assign_list = util.findClosestSphere(mesh.vertices, xyzr)
opt_spheres = opt.optimizeAsgdSpheresFromVert(mesh.vertices, xyzr, assign_list)
np.savetxt(args.output_csv, opt_spheres, header = "x,y,z,r", delimiter=",")





