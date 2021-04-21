'''
plotwise.py
Maryam Hami 
7/30/18

PLot WISE catalogue sources on VLA radio continuum image
'''

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style 
plt.style.use(astropy_mpl_style)
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import astropy.units as u
from matplotlib.patches import Ellipse


def read_fits(fitsfile):
    '''
    reads the fits file, returns the data
    '''
    hdus=fits.open(fitsfile) #remember that the fits files are in the fits directory! 
    hdu=hdus[0] #selecting the first and only image in the file 
    data=hdu.data
    header=hdu.header
    #need to fix CRVAL1, because it is negative and causing problems everywhere else:
    if header['CRVAL1'] < 0:
        header['CRVAL1']=header['CRVAL1']+360
    return data,header


def read_csv(csvfile): 
    '''
    read the csv file, return the data
    '''
    tab=np.genfromtxt(csvfile,delimiter=',',usecols=(0,2,3,4,5),skip_header=1, names=None, autostrip=True, dtype=[('GName','U16'),('Catalog','U1'),('GLong','f8'),('GLat','f8'),('Size','f8')])
    return tab


def plot(data,header,tab,plotfile):
    '''
    plot the image
    '''
    wcs=WCS(header)
    #fits file in celestial coords (world sky coords)
    wcs_celest = wcs.sub(['celestial'])
    #wise catalogue lat and long is in galactic coords
    wise_coords = SkyCoord(l=tab['GLong']*u.degree, b=tab['GLat']*u.degree, frame='galactic')
    #transform wise coords into celestial coords 
    wise_celest = wise_coords.transform_to('fk5')
    #in order to plot in right coords, need to convert pixels to world 
    #xpos,ypos = wcs_celest.wcs_world2pix(wise_celest.ra.deg,wise_celest.dec.deg,1)
    wise_RA = wise_celest.ra
    wise_DEC = wise_celest.dec
    corners = wcs_celest.calc_footprint() #find corners of image, to limit the regions of wise data
    min_RA = np.min(corners[:,0])
    max_RA = np.max(corners[:,0])
    RA_range = max_RA - min_RA
    min_Dec = np.min(corners[:,1])
    max_Dec = np.max(corners[:,1])
    Dec_range = max_Dec - min_Dec
    good = (min_RA < wise_RA.deg)&(wise_RA.deg < max_RA)&(min_Dec < wise_DEC.deg)&(wise_DEC.deg < max_Dec)
    #Plotting both fits file and wise pts 
    fig=plt.figure(figsize=(7,7))
    ax=plt.subplot(projection=wcs_celest) #world coord system, to use RA & dec in plot
    cax=ax.imshow(data[0,0],origin='lower', cmap=plt.cm.gray_r) #add a color map! 
    cbar=fig.colorbar(cax) #adding color bar
    cbar.set_label('flux density (Jy/Beam)') #labeling the color bar 
    #make CASA region file
    with open('wisesouth.rgn','w') as f:
        f2=open('wise_notes_og.txt','w')
        f.write('#CRTFv0 CASA Region Text Format version 0\n')  
        #plotting WISE regions through for loops to make the ellipses of diff sizes
        for dat,RA,DEC in zip(tab[good],wise_RA[good],wise_DEC[good]):
            xpos,ypos = wcs_celest.wcs_world2pix(RA.deg,DEC.deg,1)
            size = dat['Size']*2. #should also be divided by pixsize, but my pixsize is 1 pix/arcsec 
            if 'K'==dat['Catalog']: 
                color='red'
            elif 'C'==dat['Catalog']:
                color='blue'
            elif 'G'==dat['Catalog']:
                color='green'
            elif 'Q'==dat['Catalog']:
                color='yellow' 
            ell = Ellipse((xpos,ypos),size,size,
                    color=color, fill=False, linestyle='dashed',zorder=100)
            ax.add_patch(ell)  
            #add a line for each object in wise catalogue, so i can load the regions into CASA  
            sign='+'
            if DEC.deg < 0:
                sign='-'
            line='ellipse [[{0:02}:{1:02}:{2:08.5f}, {3}{4:03}.{5:02}.{6:07.4f}],'.format(int(RA.hms.h),int(RA.hms.m),RA.hms.s,sign,abs(int(DEC.dms.d)),abs(int(DEC.dms.m)),abs(DEC.dms.s))
            line=line+' [{0:.4f}arcsec, {1:.4f}arcsec],'.format(size/2,size/2)
            line=line+' 90.00000000deg] coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color={0}, font="DejaVu Sans", fontsize=11, fontstyle=normal, usetex=false\n'.format(color)
            #line=line+", label='{0}'\n".format(dat['GName'])
            f.write(line)
            #add the f2 so that I get GName, size, catalogue, notes
            line='{0},{1:7.2f},{2},\n'.format(dat['GName'],dat['Size'],dat['Catalog'])
            f2.write(line)
        f2.close()
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    fig.savefig(plotfile)

def main(fitsfile,csvfile,plotfile):
    data,header=read_fits(fitsfile)
    tab=read_csv(csvfile)
    plot(data,header,tab,plotfile)


if __name__=='__main__':
    main('southbar.linmos.sault.fits','wise_hii_V2.1_hrds.csv','mosaic_south.pdf')


