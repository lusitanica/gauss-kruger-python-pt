
from math import *

############################################################################
# Meridian Arc
############################################################################
def arcmer(a,equad,lat1,lat2):
    b=a*sqrt(1-equad)
    n=(a-b)/(a+b)
    a0=1.+((n**2)/4.)+((n**4)/64.)
    a2=(3./2.)*(n-((n**3)/8.))
    a4=(15./16.)*((n**2)-((n**4)/4.))
    a6=(35./48.)*(n**3)

    s1=a/(1+n)*(a0*lat1-a2*sin(2.*lat1)+a4*sin(4.*lat1)-a6*sin(6.*lat1))
    s2=a/(1+n)*(a0*lat2-a2*sin(2.*lat2)+a4*sin(4.*lat2)-a6*sin(6.*lat2))
    return s2-s1
#############################################################################
# Gauss-Kruger Projection
#############################################################################
def geogauss(lat,lon,a,equad,lat0,lon0):

    lat0=radians(lat0)
    lon0=radians(lon0)

    lat=radians(lat)
    lon=radians(lon)

    lon=lon-lon0

    N=a/sqrt(1-equad*(sin(lat))**2)
    RO=a*(1-equad)/((1-equad*(sin(lat)**2))**(3./2.))

    k1=(N/RO)+(4.*(N**2)/(RO**2))-((tan(lat))**2)

    k2=(N/RO)-((tan(lat))**2)

    k3=N/RO*(14.-58.*((tan(lat)))**2)+40.*((tan(lat))**2)+((tan(lat))**4)-9.

    x=lon*N*cos(lat)+(lon**3)/6.*N*((cos(lat))**3)*k2+(lon**5)/120.*N*((cos(lat))**5)*k3

    y=arcmer(a,equad,lat0,lat)+(lon**2)/2.*N*sin(lat)*cos(lat)+((lon**4)/24.)*N*sin(lat)*((cos(lat))**3)*k1

    return x,y
#############################################################################
# Irish Transverse Mercator
#############################################################################
def itm(lat,lon):
    # GRS-80
    a = 6378137.
    equad =0.00669437999

    # Natural Origin 
    lat0=53.5
    lon0=-8.
    
    coords = geogauss(lat,lon,a,equad,lat0,lon0)
    
    k0=0.999820
    x = (k0*coords[0])+600000.
    y = (k0*coords[1])+750000.
    return x,y
#############################################################################
# Test values 
#############################################################################    

lat=53.5
lon=-8.

xy = itm(lat,lon)
print lat,lon
print "x= %.16f" %xy[0]
print "y= %.16f" %xy[1]

