#!/usr/bin/env python
# coding: utf-8

# # Importing necessary libraries:

# *************************

# In[4]:


from Mesh_Preprocess import *


# In[5]:


import numpy as np


# In[7]:


from numba import jit


# **********************

# In[8]:


@jit(nopython=True)
def Coasre_Refine_Map_Create(Num_Triangles,Element_Node_Connectivity_new,Node_Coordinates,Element_cen_r):
    Coarse_Refine_Map = np.zeros((Num_Triangles,5),dtype = np.int32)
    
    for i in range(Num_Triangles):
        Element = Element_Node_Connectivity_new[i,0]
        Nodes = Element_Node_Connectivity_new[i,1:]
        Node_Coordinates_temp = np.zeros((3,2))
        for j in range(3):
            Node_Coordinates_temp[j,:] = Node_Coordinates[int(Nodes[j]),:-1]
        temp = Tri_area(Node_Coordinates_temp)
    
        Coarse_Refine_Map[i,0] = Element
        ctr = 1
        for j in range(Element_cen_r.shape[0]):
            temp_r = 0
            Element_r = Element_cen_r[j,0]
            Element_r_cen = Element_cen_r[j,1:]
            for k in range(3):
                Node_Coordinates_temp_1 = Node_Coordinates_temp.copy()
                Node_Coordinates_temp_1[k,:] = Element_r_cen
                temp_r = temp_r + Tri_area(Node_Coordinates_temp_1)
                
            if np.abs(temp_r - temp) < 1e-9:
                Coarse_Refine_Map[i,ctr] = Element_r
                # print(f"i: {i}, j: {j}")
                ctr = ctr +1 
    return Coarse_Refine_Map


# In[9]:


@jit(nopython=True)
def Residual_Coarse_Projection_Create(Element_Node_Connectivity_new,Node_Coordinates,Element_cen_r,Num_Triangles,A_r,b_r,phi_0):
    Coarse_Refine_Map = Coasre_Refine_Map_Create(Num_Triangles,Element_Node_Connectivity_new,Node_Coordinates,Element_cen_r)
    Coarse_Refine_Solution_Projection = np.zeros(Element_cen_r.shape[0])
    for i in range(Num_Triangles):
        for j in range(1,5):
            Element_r = Coarse_Refine_Map[i,j]
            Coarse_Refine_Solution_Projection[Element_r] = phi_0[i]
    Residual_r = (A_r@Coarse_Refine_Solution_Projection) - b_r
    Residual_Coarse_Projection = np.zeros(Num_Triangles)
    for i in range(Num_Triangles):
        temp  = 0
        for j in range(1,5):
            Element_r = Coarse_Refine_Map[i,j]
            temp = temp + np.abs(Residual_r[Element_r])
        Residual_Coarse_Projection[i] = temp

    return Residual_Coarse_Projection


# In[ ]:




