#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr 17 11:29:14 2019
function:
@author: leizhao
"""
import os,imageio
import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import time
import zlconversions as zl
import sys
import pandas as pd
try:
    import cPickle as pickle
except ImportError:
    import pickle
import glob

def get_doppio_url(dtime):
    '''dtime ids gmt time'''
    date=dtime.strftime('%Y-%m-%d')
    url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/runs/History_RUN_2018-11-12T00:00:00Z'
    return url.replace('2018-11-12',date)


def plot(lons,lats,slons,slats,temp,depth,time_str,path_save,dpi=80):
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_axes([0.01,0.05,0.98,0.87])
    # create polar stereographic Basemap instance.
    m = Basemap(projection='stere',lon_0=-70.25,lat_0=40.5,lat_ts=0,llcrnrlat=37,urcrnrlat=44,\
                llcrnrlon=-75.5,urcrnrlon=-65,rsphere=6371200.,resolution='i',area_thresh=100)
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    if len(slons)!=0:
        x1,y1=m(slons,slats)
        ax.plot(x1,y1,'ro',markersize=10)
    # draw parallels.
    parallels = np.arange(0.,90,2.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=20)
    # draw meridians
    meridians = np.arange(180.,360.,2.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20)
    x, y = m(lons, lats) # compute map proj coordinates.
    # draw filled contours.
    
    
    dept_clevs=[50,150,300,1000,2000]
    dept_cs=m.contour(x,y,depth,dept_clevs,colors='black')
    plt.clabel(dept_cs, inline = True, fontsize =15,fmt="%1.0f")
    
    
    clevs=np.arange(34,57.5,0.5)  #for all year:np.arange(34,84,1) or np.arange(34,68,1)
    cs = m.contourf(x,y,temp,clevs,cmap=plt.get_cmap('rainbow'))
    # add colorbar.
    cbar = m.colorbar(cs,location='right',pad="2%",size="5%")
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label('Fahrenheit',fontsize=25)
    # add title
    plt.title('DOPPIO MODEL BOTTOM TEMPERATURE '+time_str,fontsize=30)
    if not os.path.exists(path_save):
        os.makedirs(path_save)
    plt.savefig(os.path.join(path_save,time_str.replace(' ','t')+'.png'),dpi=dpi)
def mean_temp(temps):
    mean_temp=temps[0,0]
    for i in range(1,24):
        mean_temp+=temps[i,0]
    return mean_temp/24.0
        
def make_images(dpath,path,dt=datetime(2019,5,1,0,0,0),interval=31):
    '''dpath: the path of dictionary, use to store telemetered data
        path: use to store images
        dt: start time
        interval: how many days we need make 
    '''
    with open(dpath,'rb') as fp:
         telemetered_dict=pickle.load(fp)
    for j in range(interval):
        dtime=dt+timedelta(days=j)
        print(dtime)
        url=get_doppio_url(dtime)
        while True:
            if zl.isConnected(address=url):
                break
            print('check the website is well or internet is connected?')
            time.sleep(5)
        skip=0
        while True: 
            try:
                nc = NetCDFFile(url)
                lons=nc.variables['lon_rho'][:]
                lats=nc.variables['lat_rho'][:]
                temps=nc.variables['temp']
                depth=nc.variables['h'][:]
                break
            except RuntimeError:
                print(str(url)+': need reread')
            except OSError:
                if zl.isConnected(address=url):
                    print(str(url)+': file not exit.')
                    skip=1
                    break
            except KeyboardInterrupt:
                sys.exit()
        if skip==1:
            continue
                
        m_temp=mean_temp(temps)
        ntime=dtime
        time_str=ntime.strftime('%Y-%m-%d')
        temp=m_temp*1.8+32
        Year=str(ntime.year)
        Month=str(ntime.month)
        Day=str(ntime.day)
        slons,slats=[],[]
        try:
            slons,slats=[],[]
            for i in telemetered_dict[Year][Month][Day].index:
                slons.append(telemetered_dict[Year][Month][Day]['lon'].iloc[i])
                slats.append(telemetered_dict[Year][Month][Day]['lat'].iloc[i])
        except:
            slons,slats=[],[]
        plot(lons,lats,slons,slats,temp,depth,time_str,path)
        
def read_telemetry(path):
    """read the telemetered data and fix a standard format, the return the standard data"""
    tele_df=pd.read_csv(path,sep='\s+',names=['vessel_n','esn','month','day','Hours','minates','fracyrday',\
                                          'lon','lat','dum1','dum2','depth','rangedepth','timerange','temp','stdtemp','year'])
    if len(tele_df)<6000:
        print('check the telemetered website!')
        sys.exit()
        
    return tele_df


def seperate(filepathsave,ptelemetered='https://www.nefsc.noaa.gov/drifter/emolt.dat'):
    '''create a dictionary use to store the data from telemetered, index series is year, month, day and hour
    ptelemetered: the path of telemetered
    '''
    dfdict={}
    df=read_telemetry(ptelemetered)
    for i in df.index:
        if df['depth'][i]<2.0:
            continue
        if df['minates'].iloc[i]<=30:
            Ctime=datetime.strptime(str(df['year'].iloc[i])+'-'+str(df['month'].iloc[i])+'-'+str(df['day'].iloc[i])+' '+\
                                         str(df['Hours'].iloc[i])+':'+str(df['minates'].iloc[i])+':'+'00','%Y-%m-%d %H:%M:%S')
        else:
            Ctime=datetime.strptime(str(df['year'].iloc[i])+'-'+str(df['month'].iloc[i])+'-'+str(df['day'].iloc[i])+' '+\
                                         str(df['Hours'].iloc[i])+':'+str(df['minates'].iloc[i])+':'+'00','%Y-%m-%d %H:%M:%S')+timedelta(seconds=1800)
        Year=str(Ctime.year)
        Month=str(Ctime.month)
        Day=str(Ctime.day)
        if not Year in dfdict:
            dfdict[Year]={}
        if not Month in dfdict[Year]:
            dfdict[Year][Month]={}
        if not Day in dfdict[Year][Month]:
            dfdict[Year][Month][Day]={}

        if len(dfdict[Year][Month][Day])!=0:
            dfdict[Year][Month][Day]=dfdict[Year][Month][Day].append(pd.DataFrame(data=[[df['lat'].iloc[i],df['lon'].iloc[i],df['temp'].iloc[i]]],columns=['lat','lon','temp']).iloc[0])
            dfdict[Year][Month][Day].index=range(len(dfdict[Year][Month][Day]))
        else:
            dfdict[Year][Month][Day]=pd.DataFrame(data=[[df['lat'].iloc[i],df['lon'].iloc[i],df['temp'].iloc[i]]],columns=['lat','lon','temp'])
    with open(filepathsave,'wb') as fp:
        pickle.dump(dfdict,fp,protocol=pickle.HIGHEST_PROTOCOL)


def make_gif(gif_name,png_dir,start_time=False,end_time=False,frame_length = 0.2,end_pause = 4 ):
    '''use images to make the gif
    frame_length: seconds between frames
    end_pause: seconds to stay on last frame
    the format of start_time and end time is string, for example: %Y-%m-%d(YYYY-MM-DD)'''
    
    if not os.path.exists(os.path.dirname(gif_name)):
        os.makedirs(os.path.dirname(gif_name))
    allfile_list = glob.glob(os.path.join(png_dir,'*.png')) # Get all the pngs in the current directory
    print(allfile_list)
    file_list=[]
    if start_time:    
        for file in allfile_list:
            if start_time<=os.path.basename(file).split('.')[0]<=end_time:
                file_list.append(file)
    else:
        file_list=allfile_list
    list.sort(file_list, key=lambda x: x.split('/')[-1].split('t')[0]) # Sort the images by time, this may need to be tweaked for your use case
    images=[]
    # loop through files, join them to image array, and write to GIF called 'wind_turbine_dist.gif'
    for ii in range(0,len(file_list)):       
        file_path = os.path.join(png_dir, file_list[ii])
        if ii==len(file_list)-1:
            for jj in range(0,int(end_pause/frame_length)):
                images.append(imageio.imread(file_path))
        else:
            images.append(imageio.imread(file_path))
    # the duration is the time spent on each image (1/duration is frame rate)
    imageio.mimsave(gif_name, images,'GIF',duration=frame_length)
def main():  
    #hardcodes
    realpath=os.path.dirname(os.path.abspath(__file__))
    dpath=realpath[::-1].replace('py'[::-1],'result/Doppio'[::-1],1)[::-1]  # the directory of the result
    if not os.path.exists(dpath):
        os.makedirs(dpath)
    dictionary=os.path.join(dpath,'dictionary_emolt.p')
    gif_path=os.path.join(dpath,'gif')
    map_save=os.path.join(dpath,'map')
    gif_name =os.path.join(gif_path,'2019-08Doppio.gif')
    
    
    #############################
    #run functions
    seperate(filepathsave=dictionary)
#    make_images(dpath=dictionary,path=map_save)
#    make_gif(gif_name,map_save)
    make_images(dpath=dictionary,path=map_save,dt=datetime(2019,8,1,0,0,0))
    make_gif(gif_name,map_save,start_time='2019-08-01',end_time='2019-08-15')
    
if __name__ == "__main__":
    main()

