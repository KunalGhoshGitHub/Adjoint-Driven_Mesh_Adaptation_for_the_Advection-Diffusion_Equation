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


# ****************************************

# # Analytical Solution Contour Plot

# *************************************************

# In[5]:


def Analytical_Solution_Contour_Plot(pts,u,v,Gamma_phi,vmin,vmax,cmap,Img_File):
    x = np.linspace(0,1,pts)
    y = np.linspace(0,1,pts)
    xv, yv = np.meshgrid(x, y, indexing='ij')
    term_x = (xv + (((np.exp(u*xv/Gamma_phi) - 1))/(1-(np.exp(u/Gamma_phi)))))
    term_y = (yv + (((np.exp(v*yv/Gamma_phi) - 1))/(1-(np.exp(v/Gamma_phi)))))
    Ans = term_x*term_y
    norm = Normalize(vmin=vmin, vmax=vmax)
    plt.contourf(x,y,Ans,100, cmap = cmap, norm=norm)
    cmap = cmap
    scalar_map = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=vmin, vmax=vmax))
    scalar_map.set_array(Ans)
    ax = plt.gca()
    plt.colorbar(scalar_map,ax = ax)
    plt.axis("scaled")

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
    
    plt.savefig(Img_File)
    # plt.show()
    plt.close()

