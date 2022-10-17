import vtk

filename = "bear.obj"
reader = vtk.vtkOBJReader()
reader.SetFileName(filename)
reader.Update()

polydata = reader.GetOutput()

# normals = polydata.GetPointData().GetNormals()
# print(type(normals))
# print(normals.GetNumberOfTuples())

# for i in range(normals.GetNumberOfTuples()):
#     print(normals.GetTuple3(i))


#print(polydata)

triangleFilter = vtk.vtkTriangleFilter()
triangleFilter.SetInputData(reader.GetOutput())
triangleFilter.Update()

polygonProperties = vtk.vtkMassProperties()
polygonProperties.SetInputData(triangleFilter.GetOutput())
polygonProperties.Update()

vol = polygonProperties.GetVolume()
print(vol)