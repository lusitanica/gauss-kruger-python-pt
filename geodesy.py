from math import *

# Elipsoide pre-definidos
def elipsoides(elip):
    
# Hayford(Internacional)
    if elip == 1:
        a = 6378388.
        equad = 0.006722670022
# WGS84
    elif elip == 2:
        a = 6378137.
        equad = 0.00669437999013
# Bessel
    elif elip == 3:
        a = 6377397.155
        equad = 0.0066743723
        
    return a,equad
############################################################################
# Comprimento de arco de meridiano
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
# Projeccao directa de gauss
def geogauss(lat,lon,a,equad):

# Ponto Central da projeccao
    lat0=radians((39.+(40./60.)))
    lon0=radians(-(8.+(7./60.)+(54.862/3600.)))

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

###############################################################################
#Projeccao inversa de gauss
def gaussgeo(m,p,a,equad):

    lat0=radians((39.+(40./60.)))
    lon0=radians(-(8.+(7./60.)+(54.862/3600.)))

    sigma1=p

    fil=lat0+sigma1/(a*(1-equad))

    deltafi=1

    while deltafi > 0.0000000001:

        sigma2=arcmer(a,equad,lat0,fil)

        RO=a*(1-equad)/((1-equad*(sin(fil)**2))**(3./2.))

        deltafi=(sigma1-sigma2)/RO

        fil=fil+deltafi 


    N=a/sqrt(1-equad*(sin(fil))**2)

    RO=a*(1-equad)/((1-equad*(sin(fil)**2))**(3./2.))

    t=tan(fil)

    psi=N/RO

    lat=fil-(t/RO)*((m**2)/(2.*N))+(t/RO)*((m**4)/(24.*(N**3)))*(-4.*(psi**2)-9.*psi*(1.-t**2)+12.*(t**2))-(t/RO)*(m**6/(720.*(N**5)))*(8.*(psi**4)*(11.-24.*(t**2))-12.*(psi**3)*(21.-71.*(t**2))+15.*(psi**2)*(15.-98.*(t**2)+15.*(t**4))+180.*psi*(5.*(t**2)-3.*(t**4))-360.*(t**4))+(t/RO)*((m**8)/(40320.*(N**7)))*(1385.+3633.*(t**2)+4095.*(t**4)+1575.*(t**6))

    lon=(m/(N))-((m**3)/(6.*(N**3)))*(psi+2.*(t**2))+((m**5)/(120.*(N**5)))*(-4.*(psi**3)*(1.-6.*(t**2))+(psi**2)*(9.-68.*(t**2))+72.*psi*(t**2)+24.*(t**4))-((m**7)/(5040.*(N**7)))*(61.+662.*(t**2)+1320.*(t**4)+720.*(t**6))
	
    lon=lon0+lon/cos(fil)
	
    lat=degrees(lat)
    lon=degrees(lon)
	     
    return lat,lon

###############################################################################
# Codigo para testar as funcoes
a=elipsoides(1)[0]
equad=elipsoides(1)[1]
print a,equad
lat=38.01132704
lon=-7.87167688
print lat,lon
m=geogauss(lat,lon,a,equad)[0]
p=geogauss(lat,lon,a,equad)[1]
print "M= %.3f" %m
print "P= %.3f" %p
m=22854.099
p=-183736.383
lat=gaussgeo(m,p,a,equad) [0]
lon=gaussgeo(m,p,a,equad) [1]
print "Latitude= %.8f" %lat
print "Longitude= %.8f" %lon
