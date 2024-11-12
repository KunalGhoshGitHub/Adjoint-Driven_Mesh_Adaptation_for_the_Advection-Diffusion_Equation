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


# In[11]:


from Solver_Refine_System import *


# In[12]:


from Residual_Coarse_Projection import *


# ## Reading the contents of mesh file

# In[13]:


Mesh_File = "mesh.msh"


# In[14]:


# Mesh file format
Mesh_Ascii_Format, Mesh_File_Type, Mesh_Data_Size_Native = Mesh_Format(Mesh_File)
print(f"Mesh File Format: ")
print(f"Mesh Ascii Format: {Mesh_Ascii_Format}")
print(f"Mesh File Type: {Mesh_File_Type}")
print(f"Mesh Data Size_Native: {Mesh_Data_Size_Native}")


# In[15]:


# Mesh Elements
Mesh_Elements_Data = Mesh_Elements(Mesh_File)


# In[16]:


# This will count the number of triangles in the mesh
Num_Triangles = np.count_nonzero(Mesh_Elements_Data[:,-1] == 2)


# In[17]:


Element_Node_Connectivity = Element_Node_Connectivity_Calculate(Num_Triangles,Mesh_Elements_Data)


# In[18]:


Edge_Node_Connectivity,Boundary_Edges = Edge_Node_Connectivity_Calculate(Element_Node_Connectivity,Mesh_Elements_Data)


# In[19]:


Element_Edge_Connectivity = Element_Edge_Connectivity_Calculate(Num_Triangles,Element_Node_Connectivity,Edge_Node_Connectivity)


# In[20]:


Node_Coordinates,Point_Nodes,Curve_Nodes,Surface_Nodes = Mesh_Nodes(Mesh_File)


# In[21]:


Num_Nodes = Node_Coordinates.shape[0]


# In[22]:


Element_Element_Connectivity = Element_Element_Connectivity_Calculate_fast(Num_Triangles,Num_Nodes,Element_Node_Connectivity)
# Element_Element_Connectivity = Element_Element_Connectivity_Calculate(Num_Triangles,Element_Node_Connectivity)


# ### We need to renumber the elements

# In[23]:


Element_Element_Connectivity_new,Element_Edge_Connectivity_new,Element_Node_Connectivity_new,Edge_Node_Connectivity_new = Renumbering(Element_Element_Connectivity,Element_Edge_Connectivity,Element_Node_Connectivity,Edge_Node_Connectivity)


# In[24]:


Face_Centroid = Face_Centroid_Calculate(Edge_Node_Connectivity_new,Node_Coordinates)


# In[25]:


# Number of vertices in each element
Num_vertices = 3


# In[26]:


# Dimension of the problem
dim = 2


# In[27]:


Element_Node_Connectivity[0]
Element_Node_Coordinates_l = np.zeros((3,2))
Element_Node_Coordinates_l[0,:] = Node_Coordinates[int(Element_Node_Connectivity[0,1])][0:2]
Element_Node_Coordinates_l[1,:] = Node_Coordinates[int(Element_Node_Connectivity[0,2])][0:2]
Element_Node_Coordinates_l[2,:] = Node_Coordinates[int(Element_Node_Connectivity[0,3])][0:2]


# In[28]:


Anticlock_vertices = Anti_Clock_Triangle_vertices(Element_Node_Coordinates_l)


# In[29]:


# Check_Element_all_prop(Element_Node_Connectivity_new,Node_Coordinates,Edge_Node_Connectivity_new)


# In[30]:


Edge_Element_Connectivity = Edge_Element_Connectivity_Calculate(Edge_Node_Connectivity_new,Element_Element_Connectivity_new,Element_Edge_Connectivity_new,Boundary_Edges)


# In[31]:


Diffusion_mesh_data,Element_cen = Diffusion_mesh_data_Calculate(Num_Triangles,Element_Node_Connectivity_new,Element_Edge_Connectivity_new,Node_Coordinates,Edge_Node_Connectivity_new,Edge_Element_Connectivity)


# In[32]:


Element_Area = Element_Area_Calculate(Num_Triangles,Element_Node_Connectivity_new,Node_Coordinates)


# In[33]:


Num_Edges = Edge_Node_Connectivity_new.shape[0]


# In[34]:


Edge_Len = Edge_Len_Calculate(Num_Edges,Edge_Node_Connectivity_new,Node_Coordinates)


# ****************************************************

# # Physical Variables:

# *******************************************

# In[35]:


rho = 1
Gamma_phi = 0.4


# In[36]:


u = 10
v = 10
V = np.array([u,v])


# *************************************************

# # Calculating the source term:

# ******************************************

# In[37]:


Element_Mass = Elements_Mass_Calculate(rho,Element_Area)


# In[38]:


Element_Source_diff = Source_Cal_diff_Elements(Num_Triangles,Element_cen,Gamma_phi,u,v)


# In[39]:


Source_Term_diff = Source_Term_Elements(rho,Element_Area,Element_Source_diff)


# In[40]:


Element_Source_advec = Source_Cal_advec_Elements(Num_Triangles,Element_cen,rho,Gamma_phi,u,v)


# In[41]:


Source_Term_advec = Source_Term_Elements(rho,Element_Area,Element_Source_advec)


# ******************************************

# # Advection Diffusion solver:

# *****************************************************

# In[42]:


A,RHS,phi_0 = Advec_Diff_Solver_Steady_State(0,rho,V,Gamma_phi,1e-10,Element_Element_Connectivity_new,Source_Term_diff,Source_Term_advec,Element_cen,Element_Edge_Connectivity_new,Edge_Len,Diffusion_mesh_data,Boundary_Edges)


# In[43]:


A_psi,RHS_psi,psi = Advec_Diff_Adjoint_Solver_Steady_State(0,rho,V,Gamma_phi,1e-10,Element_Element_Connectivity_new,Source_Term_diff,Source_Term_advec,Element_cen,Element_Edge_Connectivity_new,Edge_Len,Diffusion_mesh_data,Boundary_Edges,phi_0)


# In[44]:


A_u,RHS_u,phi_0_u = Advec_Diff_Solver_Steady_State(1,rho,V,Gamma_phi,1e-10,Element_Element_Connectivity_new,Source_Term_diff,Source_Term_advec,Element_cen,Element_Edge_Connectivity_new,Edge_Len,Diffusion_mesh_data,Boundary_Edges)


# *******************************************

# In[45]:


np.linalg.cond(A)


# In[46]:


A_r,b_r,Element_cen_r = A_b_Element_Cen_Info_Refine()


# In[47]:


Residual_Coarse_Projection = Residual_Coarse_Projection_Create(Element_Node_Connectivity_new,Node_Coordinates,Element_cen_r,Num_Triangles,A_r,b_r,phi_0)


# Coarse_Refine_Map = np.zeros((Num_Triangles,5),dtype = np.int32)

# for i in range(Num_Triangles):
#     Element = Element_Node_Connectivity_new[i,0]
#     Nodes = Element_Node_Connectivity_new[i,1:]
#     Node_Coordinates_temp = np.zeros((3,2))
#     for j in range(3):
#         Node_Coordinates_temp[j,:] = Node_Coordinates[int(Nodes[j]),:-1]
#     temp = Tri_area(Node_Coordinates_temp)
# 
#     Coarse_Refine_Map[i,0] = Element
#     ctr = 1
#     for j in range(Element_cen_r.shape[0]):
#         temp_r = 0
#         Element_r = Element_cen_r[j,0]
#         Element_r_cen = Element_cen_r[j,1:]
#         for k in range(3):
#             Node_Coordinates_temp_1 = Node_Coordinates_temp.copy()
#             Node_Coordinates_temp_1[k,:] = Element_r_cen
#             temp_r = temp_r + Tri_area(Node_Coordinates_temp_1)
#             
#         if np.abs(temp_r - temp) < 1e-9:
#             Coarse_Refine_Map[i,ctr] = Element_r
#             # print(f"i: {i}, j: {j}")
#             ctr = ctr +1 

# Coarse_Refine_Solution_Projection = np.zeros(Element_cen_r.shape[0])

# for i in range(Num_Triangles):
#     for j in range(1,5):
#         Element_r = Coarse_Refine_Map[i,j]
#         Coarse_Refine_Solution_Projection[Element_r] = phi_0[i]

# Residual_r = (A_r@Coarse_Refine_Solution_Projection) - b_r

# Residual_Coarse_Projection = np.zeros(Num_Triangles)

# for i in range(Num_Triangles):
#     temp  = 0
#     for j in range(1,5):
#         Element_r = Coarse_Refine_Map[i,j]
#         temp = temp + np.abs(Residual_r[Element_r])
#     Residual_Coarse_Projection[i] = temp

# Residual_Coarse_Projection

# np.abs(psi*(A@phi_0 - RHS))

# # Analytical Solution

# ***********************************************

# In[48]:


Anal_sol_over_Area, Num_sol_over_Area = Sol_Over_Area(u,v,Gamma_phi,phi_0,Element_Area)


# In[49]:


Element_cen_phi_Actual = Analytical_Solution(Element_cen,V,Gamma_phi)


# In[50]:


Element_Vertex_Avg_sol = Element_Vertex_Avg_Anal_Sol(Num_Triangles,Element_Node_Connectivity_new,Node_Coordinates,V,Gamma_phi)


# *****************************************************

# # Analytical Solution Contour Plot

# ********************************

# In[51]:


pts = 250
vmax = np.round(max(Element_cen_phi_Actual),1)
vmin = np.round(min(Element_cen_phi_Actual),1)
cmap = "turbo"
Analytical_Solution_Contour_Plot(pts,u,v,Gamma_phi,vmin,vmax,cmap,"Analytical_Contour.eps")


# ***************************************

# # Error

# **************************************

# In[52]:


# Calculates the exact error using the solution at the centroid of the element
# error = Error_Estimate_CDAS(Element_cen_phi_Actual,phi_0,Element_Area)

# Calculates the exact error using the average of the solution at the vertices of the element
# error = Error_Estimate_CDEVAAS(Element_Vertex_Avg_sol,phi_0,Element_Area)

# Calculates the error estimate using: (phi_0 - phi_0_u)*(|Gradient(phi_0)|)*(Element_Area)
error,mod_grad = Error_Estimate_EMGA(phi_0,phi_0_u,Element_Element_Connectivity_new,Element_cen,Element_Area)


# In[ ]:





# In[53]:


psi_shifted = np.abs(psi -min(psi) + 1)
error = np.abs(psi_shifted*(Residual_Coarse_Projection))
# error = np.abs(psi -min(psi) + 1)


# In[54]:


-min(psi) + 1


# In[55]:


# Post_Process_Without_Grid(mod_grad,r"$|\nabla\Phi|$",Element_Node_Connectivity_new,Node_Coordinates,Img_file)


# In[56]:


Error_Discrete = abs(error)


# In[57]:


# Post_Process(Error_Discrete,"Error Estimate",Element_Node_Connectivity_new,Node_Coordinates,Img_file)


# In[58]:


# Post_Process_Without_Grid(Error_Discrete,"Error Estimate",Element_Node_Connectivity_new,Node_Coordinates,Img_file)


# In[59]:


# Post_Process((Error_Discrete/Error_Discrete.mean())*Element_Area[:,1],"Mean Ratio True Error",Element_Node_Connectivity_new,Node_Coordinates,Img_file)


# In[60]:


Post_Process(phi_0,r"",Element_Node_Connectivity_new,Node_Coordinates,"phi_0.eps")


# In[61]:


x_min = 0
x_max = 1
y_min = 0
y_max = 1
x_points = 1000
y_points = 1000
Smooth_Contour(x_min,x_max,y_min,y_max,x_points,y_points,Element_cen,phi_0,"linear","turbo",vmin,vmax,"Smooth_Contour_Phi_0.eps")


# In[62]:


Post_Process_Without_Grid(phi_0,r"",Element_Node_Connectivity_new,Node_Coordinates,"phi_0_without_grid.eps")


# In[63]:


Post_Process_Without_Grid(psi,r"",Element_Node_Connectivity_new,Node_Coordinates,"psi.eps")
Post_Process_Without_Grid(psi_shifted,r"",Element_Node_Connectivity_new,Node_Coordinates,"psi_shifted.eps")
Post_Process_Without_Grid(error,r"",Element_Node_Connectivity_new,Node_Coordinates,"error.eps")
Post_Process_Without_Grid(np.abs(Residual_Coarse_Projection),r"",Element_Node_Connectivity_new,Node_Coordinates,"residual.eps")


# In[64]:


# Post_Process(Element_cen_phi_Actual,r"$\Phi: Actual$",Element_Node_Connectivity_new,Node_Coordinates,Img_file)


# *******************************

# # Metric

# *********************************************

# In[65]:


Target_Dof = 1024


# In[66]:


p = 0
D = 2
q = 1
d = Calculate_Optimal_Mesh_Density(Error_Discrete,Element_Area,Target_Dof,p,D,q)


# In[67]:


metric_term_a_element = d


# ## We need to scale the metric

# In[68]:


metric_scale = Metric_Scale_Calculate(Target_Dof,Num_Triangles,metric_term_a_element,Element_Area)


# In[69]:


metric_term_a_element = metric_term_a_element*metric_scale


# In[70]:


Node_metric = Node_Wise_Metric(metric_term_a_element,Node_Coordinates,Element_Node_Connectivity_new,Element_Area)


# In[71]:


adj_metric_file = "adj_metric_file.mtr"


# In[72]:


metric_file_writer(adj_metric_file,Node_Coordinates,Node_metric)


# In[73]:


Adapted_mesh_file_name = "bamg_adapted.mesh"


# In[74]:


Bamg_mesh_filename = "Bamg.mesh"


# In[75]:


subprocess.run(f"bamg -b {Bamg_mesh_filename} -M {adj_metric_file} -v 3 -o {Adapted_mesh_file_name} -nbv 200000",shell = True)


# In[76]:


Adapted_mesh_file_writer(Adapted_mesh_file_name)


# In[77]:


subprocess.run(f"cp {Adapted_mesh_file_name} ../../gmsh-4.12.1-Linux64/bin/{Adapted_mesh_file_name}",shell =  True)


# In[78]:


subprocess.run(f"cd ../../gmsh-4.12.1-Linux64/bin; ./gmsh {Adapted_mesh_file_name} -2 -o {Mesh_File} -save_all",shell =  True)


# In[79]:


subprocess.run(f"cp ../../gmsh-4.12.1-Linux64/bin/{Mesh_File} {Mesh_File}",shell =  True)


# In[80]:


subprocess.run(f"cp {Adapted_mesh_file_name} {Bamg_mesh_filename}",shell =  True)


# In[81]:


# Plot_Edge_Number_Cell_Number(Node_Coordinates,Element_Node_Connectivity,Face_Centroid,Element_cen)


# In[82]:


Error_Discrete = np.abs(Element_Vertex_Avg_sol[:,1] - phi_0)


# In[83]:


Cen_Error_Data_Writer(Error_Discrete,Element_Area,Anal_sol_over_Area,Num_sol_over_Area,Num_Triangles,u,v,Gamma_phi)


# In[84]:


Vertex_Based_Error_Data_Writer(Element_Vertex_Avg_sol,phi_0,Num_Triangles,Element_Area,u,v,Gamma_phi)


# In[85]:


# Refinement by splitting
# ./gmsh mesh.msh -2 -refine -o out.msh -save_all

