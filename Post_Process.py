#!/usr/bin/env python
# coding: utf-8

# # Importing necessary libraries:

# ******************************************************

# In[1]:


import numpy as np


# In[2]:


import matplotlib.pyplot as plt


# In[3]:


import matplotlib


# In[4]:


from matplotlib.colors import Normalize


# In[5]:


from scipy.interpolate import griddata


# # Post Processing:

# **************************************************

# In[6]:


def Post_Process(phi,string,Element_Node_Connectivity_new,Node_Coordinates,Img_file):
    """
    Input:
    phi: Scalar
    string: String for the title
    Element_Node_Connectivity_new:
    Node_Coordinates:
    
    Output:
    A contour plot showing the value of phi in different elements
    """
    
    temp = np.zeros((4,2))
    x = Node_Coordinates[:,0]
    y = Node_Coordinates[:,1]
    # cmap = "YlGn"
    # cmap = "viridis"
    # cmap = "rainbow"
    # cmap = "jet"
    cmap = "turbo"

    cmap = matplotlib.colormaps.get_cmap(cmap)
    
    norm = Normalize(vmin=min(phi), vmax=max(phi))  # Normalize the scalar values

    fig, ax = plt.subplots()
    
    for i in range(Element_Node_Connectivity_new.shape[0]):
        Nodes = Element_Node_Connectivity_new[i,1:]
        Nodes = np.array(Nodes,dtype = int)
        temp[:3,0] = x[Nodes]
        temp[:3,1] = y[Nodes]
        temp[-1,0] = x[Nodes[0]]
        temp[-1,1] = y[Nodes[0]]
        X = temp[:,0]
        Y = temp[:,1]
        ax.plot(X,Y,"k",linewidth = 0.15)
        # color = plt.cm.get_cmap(cmap)((phi[i] - min(phi)) / (max(phi) - min(phi)))
        color = cmap(norm(phi[i]))
        
        ax.fill(X,Y,color=color)
        
        # Debug Script:
        # plt.plot(x[Nodes[0]],y[Nodes[0]],"-b*")
        
    scalar_map = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=min(phi), vmax=max(phi)))
    # scalar_map = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=0, vmax=0.7))
    scalar_map.set_array(phi)
    plt.colorbar(scalar_map,ax = ax)
    # cbar = plt.colorbar(scalar_map)
    # cbar.set_label(r'$\Phi$: '+string)
    plt.title(string)
    plt.axis("scaled")
    # plt.arrow(-0.1,-0.1,0.1,0,color = "k",width=0.0035)
    # plt.arrow(-0.1,-0.1,0,0.1,color = "k",width=0.0035)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    # Adding arrows outside the contour plot
    plt.annotate('', xy=(-0.1, 0.1), xytext=(-0.1,-0.11),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),xycoords='axes fraction', textcoords='axes fraction')
    
    plt.annotate('', xy=(0.1, -0.1), xytext=(-0.11,-0.1),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),xycoords='axes fraction', textcoords='axes fraction')
    # plt.xlim(-0.2, 1.2)
    # plt.ylim(-0.2, 1.2)

    # plt.text(-0.1, 0.11, 'y', ha='center', va='center', rotation=90, fontsize=12, transform=plt.gca().transAxes)
    plt.text(-0.1, 0.12, 'y', ha='center', va='center', fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.11, -0.1, 'x', ha='center', va='center', fontsize=12, transform=plt.gca().transAxes)

    plt.savefig(Img_file)
    # plt.show()
    plt.close()


# In[7]:


def Post_Process_Without_Grid(phi,string,Element_Node_Connectivity_new,Node_Coordinates,Img_file):
    """
    Input:
    phi: Scalar
    string: String for the title
    Element_Node_Connectivity_new:
    Node_Coordinates:
    
    Output:
    A contour plot showing the value of phi in different elements
    """
    
    temp = np.zeros((4,2))
    x = Node_Coordinates[:,0]
    y = Node_Coordinates[:,1]
    # cmap = "YlGn"
    # cmap = "viridis"
    # cmap = "rainbow"
    # cmap = "jet"
    cmap = "turbo"

    cmap = matplotlib.colormaps.get_cmap(cmap)
    
    norm = Normalize(vmin=min(phi), vmax=max(phi))  # Normalize the scalar values

    fig, ax = plt.subplots()
    
    for i in range(Element_Node_Connectivity_new.shape[0]):
        Nodes = Element_Node_Connectivity_new[i,1:]
        Nodes = np.array(Nodes,dtype = int)
        temp[:3,0] = x[Nodes]
        temp[:3,1] = y[Nodes]
        temp[-1,0] = x[Nodes[0]]
        temp[-1,1] = y[Nodes[0]]
        X = temp[:,0]
        Y = temp[:,1]
        # ax.plot(X,Y,"k")
        # color = plt.cm.get_cmap(cmap)((phi[i] - min(phi)) / (max(phi) - min(phi)))
        color = cmap(norm(phi[i]))
        
        ax.fill(X,Y,color=color)
        
        # Debug Script:
        # plt.plot(x[Nodes[0]],y[Nodes[0]],"-b*")
        
    scalar_map = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=min(phi), vmax=max(phi)))
    # scalar_map = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=0, vmax=0.7))
    
    scalar_map.set_array(phi)
    plt.colorbar(scalar_map,ax = ax)
    # cbar = plt.colorbar(scalar_map)
    # cbar.set_label(r'$\Phi$: '+string)
    plt.title(string)
    plt.axis("scaled")
    # Add arrows to indicate axis direction
    # plt.arrow(-0.1,-0.1,0.1,0,color = "k",width=0.0035)
    # plt.arrow(-0.1,-0.1,0,0.1,color = "k",width=0.0035)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    # Adding arrows outside the contour plot
    plt.annotate('', xy=(-0.1, 0.1), xytext=(-0.1,-0.11),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),xycoords='axes fraction', textcoords='axes fraction')
    
    plt.annotate('', xy=(0.1, -0.1), xytext=(-0.11,-0.1),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),xycoords='axes fraction', textcoords='axes fraction')
    # plt.xlim(-0.2, 1.2)
    # plt.ylim(-0.2, 1.2)

    # plt.text(-0.1, 0.11, 'y', ha='center', va='center', rotation=90, fontsize=12, transform=plt.gca().transAxes)
    plt.text(-0.1, 0.12, 'y', ha='center', va='center', fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.11, -0.1, 'x', ha='center', va='center', fontsize=12, transform=plt.gca().transAxes)
    
    plt.savefig(Img_file)
    # plt.show()
    plt.close()


# In[8]:


def Smooth_Contour(x_min,x_max,y_min,y_max,x_points,y_points,Element_cen,phi_0,interpolation_scheme,cmap,vmin,vmax,Img_file):
    
    # Create grid points
    grid_x, grid_y = np.meshgrid(np.linspace(x_min, x_max, x_points),np.linspace(y_min, y_max, y_points))

    # Getting the element centroid
    x_Element_cen = Element_cen[:,1]
    y_Element_cen = Element_cen[:,2]

    # Hard coding the corners
    x_Element_cen = np.concatenate((x_Element_cen,np.array([x_min,x_min,x_max,x_max])))
    y_Element_cen = np.concatenate((y_Element_cen,np.array([y_min,y_max,y_min,y_max])))
    phi_0_smooth = np.concatenate((phi_0,np.array([0,0,0,0])))
    grid_z = griddata((x_Element_cen, y_Element_cen), phi_0_smooth, (grid_x, grid_y), method=interpolation_scheme)
    plt.contourf(grid_x, grid_y, grid_z, levels=100, cmap=cmap)
    # plt.colorbar()
    
    scalar_map = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=vmin, vmax=vmax))
    scalar_map.set_array(grid_z)
    ax = plt.gca()
    plt.colorbar(scalar_map,ax = ax)
    plt.axis("scaled")
    
    # plt.scatter(Element_cen[:,1], Element_cen[:,2], c=phi_0, edgecolors='black')
    # plt.arrow(-0.1,-0.1,0.1,0,color = "k",width=0.0035)
    # plt.arrow(-0.1,-0.1,0,0.1,color = "k",width=0.0035)

    # Adding arrows outside the contour plot
    plt.annotate('', xy=(-0.1, 0.1), xytext=(-0.1,-0.11),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),xycoords='axes fraction', textcoords='axes fraction')
    
    plt.annotate('', xy=(0.1, -0.1), xytext=(-0.11,-0.1),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),xycoords='axes fraction', textcoords='axes fraction')
    # plt.xlim(-0.2, 1.2)
    # plt.ylim(-0.2, 1.2)

    # plt.text(-0.1, 0.11, 'y', ha='center', va='center', rotation=90, fontsize=12, transform=plt.gca().transAxes)
    plt.text(-0.1, 0.12, 'y', ha='center', va='center', fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.11, -0.1, 'x', ha='center', va='center', fontsize=12, transform=plt.gca().transAxes)
    
    plt.savefig(Img_file)
    # plt.show()
    plt.close()

