#!/usr/bin/env python
# coding: utf-8

# # Importing necessary libraries:

# **********************************************

# In[1]:


from Mesh_Preprocess import *


# In[2]:


from Source_Term import *


# In[3]:


from LS_Solver import *


# In[4]:


from Post_Process import *


# In[5]:


from Advec_Diff_Solver import *


# In[6]:


from Analytical_Solution_Contour_Plot import *


# In[7]:


from Analytical_Sol import *


# In[8]:


from Error_Estimate import *


# In[9]:


from Metric_Calculate import *


# In[10]:


from Error_Evaluate import *


# ## Reading the contents of mesh file

# In[11]:


Mesh_File = "mesh.msh"


# In[12]:


Refine_mesh = "mesh_refine.msh"


# In[13]:


subprocess.run(f"cp {Mesh_File} ../../gmsh-4.12.1-Linux64/bin/{Mesh_File}",shell =  True)


# In[14]:


subprocess.run(f"cd ../../gmsh-4.12.1-Linux64/bin; ./gmsh {Mesh_File} -2 -refine -o {Refine_mesh} -save_all",shell =  True)


# In[15]:


subprocess.run(f"cp ../../gmsh-4.12.1-Linux64/bin/{Refine_mesh} {Refine_mesh}",shell =  True)


# In[16]:


Mesh_File = Refine_mesh


# In[17]:


# Mesh file format
Mesh_Ascii_Format, Mesh_File_Type, Mesh_Data_Size_Native = Mesh_Format(Mesh_File)
print(f"Mesh File Format: ")
print(f"Mesh Ascii Format: {Mesh_Ascii_Format}")
print(f"Mesh File Type: {Mesh_File_Type}")
print(f"Mesh Data Size_Native: {Mesh_Data_Size_Native}")


# In[18]:


# Mesh Elements
Mesh_Elements_Data = Mesh_Elements(Mesh_File)


# In[19]:


# This will count the number of triangles in the mesh
Num_Triangles = np.count_nonzero(Mesh_Elements_Data[:,-1] == 2)


# In[20]:


Element_Node_Connectivity = Element_Node_Connectivity_Calculate(Num_Triangles,Mesh_Elements_Data)


# In[21]:


Edge_Node_Connectivity,Boundary_Edges = Edge_Node_Connectivity_Calculate(Element_Node_Connectivity,Mesh_Elements_Data)


# In[22]:


Element_Edge_Connectivity = Element_Edge_Connectivity_Calculate(Num_Triangles,Element_Node_Connectivity,Edge_Node_Connectivity)


# In[23]:


Node_Coordinates,Point_Nodes,Curve_Nodes,Surface_Nodes = Mesh_Nodes(Mesh_File)


# In[24]:


Num_Nodes = Node_Coordinates.shape[0]


# In[25]:


Element_Element_Connectivity = Element_Element_Connectivity_Calculate_fast(Num_Triangles,Num_Nodes,Element_Node_Connectivity)
# Element_Element_Connectivity = Element_Element_Connectivity_Calculate(Num_Triangles,Element_Node_Connectivity)


# ### We need to renumber the elements

# In[26]:


Element_Element_Connectivity_new,Element_Edge_Connectivity_new,Element_Node_Connectivity_new,Edge_Node_Connectivity_new = Renumbering(Element_Element_Connectivity,Element_Edge_Connectivity,Element_Node_Connectivity,Edge_Node_Connectivity)


# In[27]:


Face_Centroid = Face_Centroid_Calculate(Edge_Node_Connectivity_new,Node_Coordinates)


# In[28]:


# Number of vertices in each element
Num_vertices = 3


# In[29]:


# Dimension of the problem
dim = 2


# In[30]:


Element_Node_Connectivity[0]
Element_Node_Coordinates_l = np.zeros((3,2))
Element_Node_Coordinates_l[0,:] = Node_Coordinates[int(Element_Node_Connectivity[0,1])][0:2]
Element_Node_Coordinates_l[1,:] = Node_Coordinates[int(Element_Node_Connectivity[0,2])][0:2]
Element_Node_Coordinates_l[2,:] = Node_Coordinates[int(Element_Node_Connectivity[0,3])][0:2]


# In[31]:


Anticlock_vertices = Anti_Clock_Triangle_vertices(Element_Node_Coordinates_l)


# In[32]:


# Check_Element_all_prop(Element_Node_Connectivity_new,Node_Coordinates,Edge_Node_Connectivity_new)


# In[33]:


Edge_Element_Connectivity = Edge_Element_Connectivity_Calculate(Edge_Node_Connectivity_new,Element_Element_Connectivity_new,Element_Edge_Connectivity_new,Boundary_Edges)


# In[34]:


Diffusion_mesh_data,Element_cen = Diffusion_mesh_data_Calculate(Num_Triangles,Element_Node_Connectivity_new,Element_Edge_Connectivity_new,Node_Coordinates,Edge_Node_Connectivity_new,Edge_Element_Connectivity)


# In[35]:


Element_Area = Element_Area_Calculate(Num_Triangles,Element_Node_Connectivity_new,Node_Coordinates)


# In[36]:


Num_Edges = Edge_Node_Connectivity_new.shape[0]


# In[37]:


Edge_Len = Edge_Len_Calculate(Num_Edges,Edge_Node_Connectivity_new,Node_Coordinates)


# ****************************************************

# # Physical Variables:

# *******************************************

# In[38]:


rho = 1
Gamma_phi = 0.4


# In[39]:


u = 10
v = 10
V = np.array([u,v])


# *************************************************

# # Calculating the source term:

# ******************************************

# In[40]:


Element_Mass = Elements_Mass_Calculate(rho,Element_Area)


# In[41]:


Element_Source_diff = Source_Cal_diff_Elements(Num_Triangles,Element_cen,Gamma_phi,u,v)


# In[42]:


Source_Term_diff = Source_Term_Elements(rho,Element_Area,Element_Source_diff)


# In[43]:


Element_Source_advec = Source_Cal_advec_Elements(Num_Triangles,Element_cen,rho,Gamma_phi,u,v)


# In[44]:


Source_Term_advec = Source_Term_Elements(rho,Element_Area,Element_Source_advec)


# ******************************************

# # Advection Diffusion solver:

# *****************************************************

# In[45]:


A,b = Advec_Diff_Solver_Steady_State_A_b(0,rho,V,Gamma_phi,1e-10,Element_Element_Connectivity_new,Source_Term_diff,Source_Term_advec,Element_cen,Element_Edge_Connectivity_new,Edge_Len,Diffusion_mesh_data,Boundary_Edges)


# In[47]:


def A_b_Element_Cen_Info_Refine():
    return A,b,Element_cen


# In[ ]:




