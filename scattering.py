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
                 pixel_size = .172,wavelength = .123894):
        self.SDD = SDD
        self.center_x = center_x
        self.center_y = center_y
        self.xpixels = xpixels
        self.ypixels = ypixels
        self.pixel_size = pixel_size
        self.qpx = 2*np.pi*(pixel_size*10e6)/(self.SDD*10e6*wavelength)
        self.qp = (np.arange(1,xpixels+1)-center_x)*self.qpx
        self.qz = -(np.arange(1,ypixels+1)-center_y)*self.qpx
        self.xran = np.arange(self.center_x,self.xpixels)
        self.yran = np.arange(0,self.center_y+1)
        self.wavelength = wavelength
        
    def load(self):
        self.data = np.array(fabio.open(self.name).data,dtype=float)
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
    
    def SRFconvert(self,alphai):
        self.alphai = alphai*np.pi/180
        ypos = (self.center_y - np.arange(1,self.ypixels+1))*self.pixel_size
        ypos = ypos.reshape(len(ypos),1)
        xpos = (np.arange(1,self.xpixels+1)-self.center_x)*self.pixel_size
        xpos = xpos.reshape((1,len(xpos)))
        gamma = np.arctan(xpos/self.SDD)
        delta = np.arctan(ypos/np.sqrt(self.SDD**2+xpos**2))
        qx = 2*np.pi/self.wavelength*(-np.cos(delta)*np.sin(gamma));
        qy = 2*np.pi/self.wavelength*(np.cos(self.alphai)*(np.cos(delta)*np.cos(gamma)-1)+np.sin(self.alphai)*np.sin(delta));
        self.qz = 2*np.pi/self.wavelength*(np.cos(self.alphai)*np.sin(delta) + np.sin(self.alphai)*(1 - np.cos(delta)*np.cos(gamma)));
        self.q = np.sqrt(qx**2+qy**2+self.qz**2);
        self.qr = np.sqrt(qx**2+qy**2);
        self.chi = np.arccos(self.qz/self.q);
        self.chi[:,0:self.center_x] *= -1
        self.qr[:,0:self.center_x] *= -1
        placeholder = self.data
        placeholder[:,self.center_x-1] = np.nan
        self.dataN = placeholder
        self.q[self.center_y:-1,:] *= -1
        self.chi[np.isnan(self.chi)] = 0
    
    def remesh_qrqz(self):
        self.qrmin = np.min(self.qr)
        self.qrmax = np.max(self.qr)
        self.qzmax = np.max(self.qz)
        self.qzmin = np.min(self.qz)
        qr1d = self.qr.ravel()
        qz1d = self.qz.ravel()
        D = self.data.ravel()
        
        # Remesh with no pixel-splitting or interpolation
        num_bins_r = int((self.qrmax-self.qrmin)/self.qpx) # need to calculate dchi if want to use q v chi plot
        num_bins_z = int((self.qzmax-self.qzmin)/self.qpx)
        bins=[num_bins_r,num_bins_z]
        remesh_data, rbin, zbin = np.histogram2d(qr1d, qz1d, bins=bins, normed=False, weights=D)
        num_per_bin, rbin, zbin = np.histogram2d(qr1d,qz1d,bins=bins, normed=False,weights=None)
        with np.errstate(divide='ignore',invalid='ignore'):
            remesh_data =  remesh_data/num_per_bin
            remesh_data[np.isnan(remesh_data)] = 1
        self.remesh_data = np.rot90(remesh_data)
        
    def remesh_chiq(self):
        self.chimin = np.min(self.chi)
        self.chimax = np.max(self.chi)
        self.qmax = np.max(self.q)
        self.qmin = np.min(self.q)
        chi1d = self.chi.ravel()
        q1d = self.q.ravel()
        D = self.data.ravel()
        
        # Remesh with no pixel-splitting or interpolation
        num_bins_r = int((self.chimax-self.chimin)/np.arctan(self.qpx))/2.5 # need to calculate dchi if want to use q v chi plot
        num_bins_z = int((self.qmax-self.qmin)/self.qpx)
        bins=[num_bins_r,num_bins_z]
        remesh_data, rbin, zbin = np.histogram2d(chi1d, q1d, bins=bins, normed=False, weights=D)
        num_per_bin, rbin, zbin = np.histogram2d(chi1d,q1d,bins=bins, normed=False,weights=None)
        with np.errstate(divide='ignore',invalid='ignore'):
            remesh_data =  remesh_data/num_per_bin
            remesh_data[np.isnan(remesh_data)] = 1
        self.remesh_data = np.rot90(remesh_data)

    def meshplot(self,plottype='qrqz',cmap='jet',raster='True'):
        Zm = np.ma.masked_invalid(self.dataN)
        if plottype == 'qrqz':
            cs = plt.pcolormesh(self.qr,self.qz,np.log(Zm),cmap=cmap,antialiased='True',rasterized=raster)
        elif plottype == 'chiq':
            cs = plt.pcolormesh(self.chi,self.q,np.log(Zm),cmap=cmap,antialiased='True',rasterized=raster)
        elif plottype == 'chiqz':
            cs = plt.pcolormesh(self.chi,self.qz,np.log(Zm),cmap=cmap,antialiased='True',rasterized=raster)
        else:
            print('Incorrect plot type')
        return cs
    
    def implot(self,plottype='qrqz',cmap='jet'):
        if plottype == 'qrqz':
            cs = plt.imshow(np.log(self.remesh_data),interpolation='nearest',cmap=cmap,extent=(self.qrmin, self.qrmax, self.qzmin, self.qzmax))
        elif plottype == 'chiq':
            cs = plt.imshow(np.log(self.remesh_data),interpolation='nearest',cmap=cmap,extent=(self.chimin, self.chimax, self.qmin, self.qmax))
        else:
            print('Incorrect plot type')
        return cs
    
    
    
if __name__ == '__main__':
    
    # import packages
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
#    import time
    
    # update parameters
    mpl.rcParams.update({'font.sans-serif':'Arial','font.family':'sans-serif'})
    mpl.rcParams['figure.dpi'] = 150
    sns.set_context(context='paper')
    mpl.rcParams['svg.fonttype'] = 'none'
    mapvalue = mpl.cm.get_cmap('jet')
    fcolor = mapvalue(0)
    
    # import and process data
    image1 = reduction('PAAGNP1_A0p160_2m.edf')
    image1.load()
    data = image1.data
#    image1.raw_plot((12,10))
    image1.geometry(300,622,1320)
    image1.SRFconvert(0.16)
    
    # plotting
    fig, ax = plt.subplots(figsize=(6,4.5))
    ax.set_facecolor(fcolor)
    cs = image1.meshplot('qrqz')
#    plt.xlim((-1.5,1.5))
#    plt.ylim((0, 30))
    plt.colorbar(cs)
    plt.xlabel(r'Azimuthal Angle - $\chi$ [rad]')
    plt.ylabel(r'q$_z$ [nm$^{-1}$]')
#    time1 = time.time()
    fig.savefig("testimage.png",dpi=600)
#    time2 = time.time()
#    print(time2-time1)
    plt.show()

