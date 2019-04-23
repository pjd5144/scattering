# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:27:32 2019

@author: pjdudenas
"""
import fabio
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':22})

        
class reduction:
    
    def __init__(self,name):
        self.name = name
        
    def geometry(self,SDD,center_x, center_y, xpixels = 1475, ypixels = 1679, 
                 pixel_size = .172,wavelength = 1.23894):
        self.SDD = SDD
        self.center_x = center_x
        self.center_y = center_y
        self.xpixels = xpixels
        self.ypixels = ypixels
        self.pixel_size = pixel_size
        self.qpx = 2*np.pi*(pixel_size*10e6)/(self.SDD*10e6*wavelength/10)
        self.qp = (np.arange(1,xpixels+1)-center_x)*self.qpx
        self.qz = -(np.arange(1,ypixels+1)-center_y)*self.qpx
        self.xran = np.arange(self.center_x,self.xpixels)
        self.yran = np.arange(0,self.center_y+1)
        
    def load(self):
        self.data = np.array(fabio.open(self.name).data)
        self.data[self.data < 1] = 1

    def raw_plot(self,size=(12,9)):
        fig, ax = plt.subplots(figsize=size)
        plt.imshow(np.log(self.data),cmap='jet')
        plt.colorbar()
        plt.show()
    
    def plot(self,size=(9,9),show_cbar='False'):
        fig, ax = plt.subplots(figsize=size)
        plt.imshow(np.log(self.data),cmap='jet',extent=[self.qp[0],self.qp[-1],self.qz[-1],self.qz[0]],
                   aspect='auto')
        plt.xlim([-2,2])
        plt.ylim([0,4])
        plt.yticks(np.arange(1,5))
        plt.xticks(np.arange(-2,3))
#        tick_loc, tick_label = plt.yticks()
#        ax.set_yticklabels(map(str,(np.abs(tick_loc))))
        plt.xlabel(r'$q_p$ $[nm^{-1}]$')
        plt.ylabel(r'$q_z$ $[nm^{-1}]$')
        if show_cbar == 'True':
            plt.colorbar()
        plt.show()
    
    def qp_linecut(self,ypixel1=1200,ypixel2=1215):
        I = np.mean(self.data[ypixel1:ypixel2,self.xran],axis=0)
        qp = self.qp[self.xran]
#        print(self.xran.shape)
        return qp, I

    def qz_linecut(self,xpixel1=650,xpixel2=660):
        I = np.mean(self.data[self.yran,xpixel1:xpixel2],axis=1)
        qz = self.qz[self.yran]
        return qz, I


if __name__ == '__main__':
    

    image1 = reduction('PAAGNP1_A0p160_2m.edf')
    image1.load()
    data = image1.data
    image1.raw_plot((12,10))
    image1.geometry(1990,622,1320)
#    print(image1.qz)
    image1.plot((12,10))
    qz, I = image1.qz_linecut()
    plt.figure()
    plt.loglog(qz,I)
    plt.show()
