from mpl_toolkits.basemap import Basemap
from numpy import *
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------------------------------------------------------------------
#- input ----------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------

#- directory where error files are located
path_errorfile='C:/errorfiles/'

#- list of events for which errorfile is to be read
#event_index_list=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,1002,1005,1009,1010,1012,1014,1017,1018,1019,1020,1021,1027,1029,1030,1037,1039,1048,1054,1062,1065,1066,1068,1069,1071,1072,1073,1074,1075,1076,2001,2004,2006,2011,2012,2014,2017,2020,2021,2022,2024,2025,2032,2035,2036,5008,5011,5022,5023,5024,5026,5028,5033,5035,5046,5047]
event_index_list=[1,2]

#- plot ray paths?
plot_ray_paths=True

#- parameters of Lambert Conical Projection
L_width=10000000
L_height=8000000
L_res='c'
L_lat1=45.0
L_lat2=55.0
L_lat0=60.0
L_lon0=20.0
lat_spacing=10.0
lon_spacing=10.0

#--------------------------------------------------------------------------------------------------------------------------------------
#- define errorfile class ------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------

class errorfile:

    def __init__(self):

        #- event index
        self.idx=0
        #- source latitude and longitude
        self.src_lat=0.0
        self.src_lon=0.0
        #- number of stations
        self.n_rec=0
        #- list of station names
        self.station_name=[]
        #- arrays of station latitudes and longitudes
        self.rec_lat=array([],dtype=float)
        self.rec_lon=array([],dtype=float)
        #- array of weights
        self.weight=array([],dtype=float)
        #- array of unweighted misfits
        self.misfit_unweighted=array([],dtype=float)

#--------------------------------------------------------------------------------------------------------------------------------------
#- march through the list of errorfiles -----------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------

N=len(event_index_list)
event_list=[]

for n in range(N):

    #- make class

    event=errorfile()
    event.idx=int(event_index_list[n])

    #- open file

    fn=path_errorfile+str(event_index_list[n])+'/errorfile'
    fid=open(fn,'r')

    #- read event coordinates

    line=fid.readline().strip().split(':')[1]
    event.src_lat=float(line.strip().split(' ')[0])
    event.src_lon=float(line.strip().split(' ')[1])

    #- march through all the stations

    eof=False
    k=0
    fid.readline()

    while eof==False:

        event.n_rec=event.n_rec+1
        #- station name
        event.station_name.append(fid.readline().strip())
        #- station coordinates
        line=fid.readline().strip().split(':')[1]
        lat=float(line.strip().split(' ')[0])
        lon=float(line.strip().split(' ')[1])
        event.rec_lat=hstack((event.rec_lat,[lat]))
        event.rec_lon=hstack((event.rec_lon,[lon]))
        #- misfit and weights
        line=fid.readline()
        while line[0:5]=='taper':
            line=fid.readline()

        line=line.strip().split(':')
        misfit=float(line[1].strip().split(' ')[0])
        event.misfit_unweighted=hstack((event.misfit_unweighted,[misfit]))
        total_weight=float(line[1].strip().split(' ')[1])+float(line[1].strip().split(' ')[2])
        event.weight=hstack((event.weight,[total_weight]))

         #- next receiver index
        idx=fid.readline()

        #- check end of file
        if idx=='':
            eof=True

    #- close file

    fid.close()

    #- append to event_list

    event_list.append(event)

#--------------------------------------------------------------------------------------------------------------------------------------
#- plot source, receiver, ray path map ---------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------

#- generate a map

plt.subplot(121)
m = Basemap(width=L_width,height=L_height,rsphere=(6378137.00,6356752.3142),resolution=L_res,projection='lcc',lat_1=L_lat1,lat_2=L_lat2,lat_0=L_lat0,lon_0=L_lon0)

#- plot receiver positions and great-circle paths

for n in range(N):

    slat=event_list[n].src_lat
    slon=event_list[n].src_lon
    if slon>180:
        slon=slon-360

    for k in range(event_list[n].n_rec):
        rlat=float(event_list[n].rec_lat[k])
        rlon=float(event_list[n].rec_lon[k])
        if rlon>180:
            rlon=rlon-360

        if plot_ray_paths==True:
            m.drawgreatcircle(slon,slat,rlon,rlat,linewidth=0.5,color='k')

        x,y=m(rlon,rlat)
        m.plot(x,y,'bo',markersize=5.0)

#- plot source positions

for n in range(N):

    slat=event_list[n].src_lat
    slon=event_list[n].src_lon
    if slon>180:
        slon=slon-360

    x,y=m(slon,slat)
    if (event_list[n].idx>=1000) and (event_list[n].idx<2000):
        m.plot(x,y,'y*',markersize=15.0)
    else:
        m.plot(x,y,'r*',markersize=15.0)


m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color=[0.9,0.9,0.9],lake_color=[1.0,1.0,1.0])
m.drawparallels(arange(-90,90,lat_spacing),labels=[1,0,0,1])
m.drawmeridians(arange(-180,180,lon_spacing),labels=[1,0,0,1])
m.drawmapboundary(fill_color=[1.0,1.0,1.0])
plt.title('source-receiver geometry')

#--------------------------------------------------------------------------------------------------------------------------------------
#- plot misfit histograms ----------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------

#- assemble weighted and unweighted misfits into vectors

misfit_unweighted=array([],dtype=float)
misfit_weighted=array([],dtype=float)

for n in range(N):

        misfit_unweighted=hstack((misfit_unweighted,event_list[n].misfit_unweighted))
        misfit_weighted=hstack((misfit_weighted,event_list[n].misfit_unweighted*event_list[n].weight))

plt.subplot(122)
plt.hist(misfit_weighted,50,facecolor='blue',alpha=0.20)
plt.hist(misfit_unweighted,50,facecolor='black',alpha=0.75)
plt.title('misfit unweighted (black) and weighted (blue)')

#--------------------------------------------------------------------------------------------------------------------------------------
#- print statistics --------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------

print '- statistics -------------------------------------------------------'
print '-----------------------------------------------------------------------'
print 'number of measurements: {}'.format(len(misfit_weighted))
print '-----------------------------------------------------------------------'
print 'minimum misfit unweighted: {}'.format(min(misfit_unweighted))
print 'maximum misfit unweighted: {}'.format(max(misfit_unweighted))
print 'mean misfit unweighted: {}'.format(mean(misfit_unweighted))
print '-----------------------------------------------------------------------'
print 'minimum misfit weighted: {}'.format(min(misfit_weighted))
print 'maximum misfit weighted: {}'.format(max(misfit_weighted))
print 'mean misfit weighted: {}'.format(mean(misfit_weighted))
print '-----------------------------------------------------------------------'

plt.show()
