import pybullet as p
# import pybullet_data as pd
import os
import wavefront as wf
import util
import testforGDS
import Newsurf
import optimization as opt
import numpy as np
import argparse, sys
import vtk
import meshplex
import trimesh
import math

parser = argparse.ArgumentParser(description=__doc__, formatter_class=
argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--input_obj", "-i", type=str, dest="input_obj", default="mesh.obj",
                    help="File name of the .obj input")
parser.add_argument("--output_csv", "-o", type=str, dest="output_csv", default="result.csv",
                    help="File name of the .csv output (after optimization)")
parser.add_argument("--approx_csv", "-a", type=str, dest="approx_csv", default="quick_approx.csv",
                    help="File name of a quick approximation of the decomposition (before optimization)")
parser.add_argument("--convex_obj", "-c", type=str, dest="convex_obj", default="parts.obj",
                    help="File name of the intermediate .obj convex shapes")
parser.add_argument("--convex_log", type=str, dest="convex_log", default="log.txt",
                    help="File name of the intermediate convex decomposition logs")
parser.add_argument("--interpolate_ratio", "-r", type=float, dest="ratio", default=0.5,
                    help="Must be in [0,1]; set it close to 1 to make the final output spheres large, close to 0 to make them small")
parser.add_argument("--voxel_resolution", "-v", type=int, dest="voxel_resolution", default=50000,
                    help="The resolution at which the convex decomposition takes place; larger number usually leads to more final spheres")
parser.add_argument("--max_hull_verts", type=int, dest="max_hull_verts", default=64,
                    help="The max number of vertices of one convex hull")
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
parts = wf.load_obj(args.convex_obj)  # , triangulate=True)

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



#print(xyzr)
# 1. the ratio of length and width
box_length = []
box_width = []
box_ratio = []
box_diference = []
volume_percentage = []

# 3.the percentage that the part takes refer to the whole object

whole_mesh = trimesh.load('parts.obj')
whole_volume = whole_mesh.volume
# part_volume=whole_mesh[0].volume

# print(whole_volume)
# print(part_volume)

loop_time=0
max_num_loops=int(input())
print(max_num_loops)
while(loop_time<max_num_loops):
    part_id = 0
    for part in parts:

        bounding_box = util.bbox(part.vertices)

        # 1. get the length, width of the box and their ratio
        box_length.append(math.sqrt(
            (bounding_box[1, 0] - bounding_box[0, 0]) ** 2 + (bounding_box[1, 1] - bounding_box[0, 1]) ** 2 + (
                        bounding_box[1, 2] - bounding_box[0, 2]) ** 2))  # the difference of the z axis value
        a = [bounding_box[1, 1] - bounding_box[0, 1],
             bounding_box[1, 0] - bounding_box[0, 0]]  # the difference of x axis and y axis
        box_width.append(min(a))
        box_ratio.append(box_length[part_id] / box_width[part_id])

        big_sphere = util.box2ball(bounding_box)

        mesh_center = util.coord_avg(part.vertices)

        small_sphere = util.inscribedSphereViaPointMesh(mesh_center, part)

        decomp_sphere = util.interpolateSphere(big_sphere, small_sphere, args.ratio)
        #xyzr[part_id, :3] = decomp_sphere.center
       # xyzr[part_id, -1] = decomp_sphere.radius

        # 2. the value of modification (suspected)
        box_diference.append((box_length[part_id] - box_width[part_id]) / decomp_sphere.radius)

        # 3.the percentage that the part takes refer to the whole object

        filename = "my_temp_part" + str(part_id) + ".obj"

        wf.save_obj(part, filename)
        part_mesh = trimesh.load(filename)
        part_volume = part_mesh.volume
        # print(part_volume)
        volume_percentage.append(part_volume / whole_volume)
        part_id += 1



    #getting new parts arr
    weight1 = 1 / 3
    weight2 = 1 / 3
    weight3 = 1 / 3

    box_judge = []
    for i in range(len(parts)):
        box_judge.append(weight1 * box_ratio[i] + weight2 * box_diference[i] + weight3 * volume_percentage[i])
    # print(box_judge[i])
    pass
    ##outputting critical convex
    output_convex_id = box_judge.index(max(box_judge))
    print("%d", output_convex_id)
    out_put_convex = parts[output_convex_id]
    print("%d", len(parts))

    #print("%d", out_put_convex)

    '''
    input a convex into testforGDS
    '''
    objFilePath = "my_temp_part" + str( output_convex_id) + ".obj"
    with open(objFilePath) as file:
        points = []
        face = []
        while 1:
            line = file.readline()
            if not line:
                break
            strs = line.split(" ")
            if strs[0] == "v":
                points.append((float(strs[1]), float(strs[2]), float(strs[3])))
            if strs[0] == "f":
                face.append((int(strs[1]), int(strs[2]), int(strs[3])))
    # I get points set and surfaces set
    points = np.array(points)
    face = np.array(face)
    # get the centerpoint of outersphere, and the two most distal points
    centerpoint, radiusoutersphere, point1, point2 = testforGDS.outersphere(points)
    centerpoint = np.array(centerpoint)
    point1 = np.array(point1)
    point2 = np.array(point2)

 #   print('pointa:', point1, 'pointb:', point2)

    Volume = 0
   # print('radiusoutersphere:', radiusoutersphere)
    # get the volume of object
    for i in range(len(face)):
        delta = testforGDS.volume(points[face[i][0] - 1], points[face[i][1] - 1], points[face[i][2] - 1], centerpoint)
        Volume += delta
   # print('Volume:', Volume)
    # interplate the point on each side for the preparation of computing variance
    a = []
    for i in range(len(face)):
        set = testforGDS.surface(points[face[i][0] - 1], points[face[i][1] - 1], points[face[i][2] - 1])
        a = a + set
    # split the object with the same volume, and compare the variance to get the best one.
    # spherepoint1,spherepoint2,Radius,radius,Vari=split(Volume,point1,point2,set)
    spherepoint1, spherepoint2, Radius, radius, Vari = testforGDS.split(Volume, point1, point2, points)
    print('point1:', spherepoint1, 'point2:', spherepoint2)
    print('Radius:', Radius, 'radius:', radius, 'Vari:', Vari)

    R = Radius
    r = radius
    a = spherepoint1
    b = spherepoint2
    d = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5

    if d > R + r:
        a = np.array(a)
        b = np.array(b)
        origin = (a - b) * (r / d) + b
        normal = a - b
    if d <= R + r:
        origin, normal = Newsurf.Newsurf(R, r, d, a, b)
    normal = np.array(normal)
    normal=normal/((normal[0]**2+normal[1]**2+normal[2]**2)**0.5)

    print('origin:', origin, 'normal:', normal)



    #####inputting the origin and normal to planes and normal
    #####
    planes = origin
    normals = normal
    cut_result = []
    halfA = util.sliceMesh(out_put_convex, normals, planes)
    halfB = util.sliceMesh(out_put_convex, -(normals), planes)
    if len(halfA.vertices) > 0:
        cut_result.append(halfA)
    if len(halfB.vertices) > 0:
        cut_result.append(halfB)

    #print(parts)

    parts.remove(out_put_convex)
    #print(parts)

    '''print("%d \n", len(parts))
    filename = "my_temp_part_halfA" + ".obj"

    wf.save_obj(halfA, filename)
    filename = "my_temp_part_halfB" + ".obj"

    wf.save_obj(halfB, filename)'''



    parts.append(halfA)
    #print(parts)

    parts.append(halfB)
    #print(parts)

    #filename = "my_temp_part_parts" + ".obj"

    #wf.save_obj(parts[18], filename)

    print("%d", len(parts))

    #getting new xyzr arr
    xyzr = np.zeros((len(parts), 4))

    part_id=0
    for part in parts:
        bounding_box = util.bbox(part.vertices)
        big_sphere = util.box2ball(bounding_box)

        mesh_center = util.coord_avg(part.vertices)

        small_sphere = util.inscribedSphereViaPointMesh(mesh_center, part)

        decomp_sphere = util.interpolateSphere(big_sphere, small_sphere, args.ratio)
        xyzr[part_id, :3] = decomp_sphere.center
        xyzr[part_id, -1] = decomp_sphere.radius
        part_id += 1

    loop_time+=1















##optimization
np.savetxt(args.approx_csv, xyzr, header="x,y,z,r", delimiter=",")

# The ball number which this vertex got assigned to.
assign_list = util.findClosestSphere(mesh.vertices, xyzr)
opt_spheres = opt.optimizeAsgdSpheresFromVert(mesh.vertices, xyzr, assign_list)
np.savetxt(args.output_csv, opt_spheres, header="x,y,z,r", delimiter=",")





