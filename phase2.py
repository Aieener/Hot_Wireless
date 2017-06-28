# phase2.py
# CIS 553
# data process phase 2
# data visualization phase 3
# Author: Yuding Ai
# date: Feb 10
from mpl_toolkits.basemap import Basemap
import math
import gmplot
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
import matplotlib.ticker as mtick
from matplotlib.ticker import LinearLocator, FormatStrFormatter


#=============================================================================
# Original data 
#=============================================================================
# since I realized that the origin point of each access point I record at phase
#1 is not actually the accurate one, so I decided to redo the origin location
#of each access point, the BSSID is the same as in my phase 1 though

# the original data of lattitude and longitude for the library WAP 84:d4:7e:dc:2a:40
liboriginsignal = -41
liborigin=[39.95233362,-75.19346528]

liblat=[39.95233362,39.95230982,39.95228379,39.95227624,39.95232304,39.95236796,
        39.95237908,39.9523782,39.95237557,39.9523449,39.95231382,39.95228106,
        39.95228917,39.95231623,39.95232877,39.9522439,39.95221351,39.95216466,
        39.95214442,39.95216179,39.9523527,39.95239043,39.95242305,39.95242487,
        39.95240718,39.95240224,39.95232515,39.95229446,39.95233105]

liblon=[-75.19346528,-75.19354176,-75.19354936,-75.19355259,-75.19359157,
        -75.19357628,-75.19354717,-75.19348767,-75.19339785,-75.1933452,
        -75.19332773,-75.19327585,-75.19323223,-75.19326823,-75.19337275,
        -75.19339994,-75.1934048,-75.19332293,-75.1932546,-75.19316639,
        -75.1934005,-75.19345689,-75.19352098,-75.19343724,-75.19331876,
        -75.19327910,-75.19335194,-75.19333902,-75.19334085]

libsignal=[-41,-55,-61,-61,-67,-68,-63,-58,-60,-70,-75,-76,-76,-70,-70,
        -82,-83,-78,-80,-85,-69,-62,-52,-61,-67,-69,-72,-66,-77]

# the original data of lattitude and longitude (in degree) for the Starbucks WAP 94:b4:0f:53:07:00
staroriginsignal = -43 
starorigin = [39.95322759,-75.19217689]
starlat = [39.95300328,39.95314201,39.95308746,39.95322759,39.95331679,39.95341871,
        39.95341624,39.95336173,39.95333451,39.95329762,39.95323437,39.95316914,
        39.95305603,39.9530076,39.95296533,39.9529472,39.95293145,39.95293909,
        39.95294534,39.95302457,39.95307636,39.95329533,39.95334517,39.953197,
        39.95318632,39.95323969,39.95331454,39.95331931,39.95334841,39.95335327,
        39.95333676,39.9531888,39.95318109,39.95315349,39.95332706,39.95334124,
        39.9533159,39.95329751,39.95325922]
starlon = [-75.19218521,-75.19215558,-75.19214409,-75.19217689,-75.19219967,
        -75.19223796,-75.19221528,-75.19219009,-75.19220592,-75.19221843,
        -75.19222063,-75.19222198,-75.19222518,-75.19222164,-75.19222223,
        -75.19219909,-75.19221608,-75.19220119,-75.19199002,-75.19194761,
        -75.1918776,-75.19169873,-75.19191655,-75.191923,-75.19193303,
        -75.19191041,-75.19192654,-75.19193617,-75.19195773,-75.19196058,
        -75.19197971,-75.19204021,-75.19197923,-75.19183212,-75.19149968,
        -75.1914769,-75.19142933,-75.19145057,-75.19149104]
starsignal = [-67,-54,-54,-43,-68,-77,-69,-64,-67,-65,-55,-54,-56,-69,-68,
        -70,-68,-70,-70,-72,-68,-68,-73,-59,-59,-61,-65,-68,-67,-67,-71,-59,
        -70,-71,-73,-77,-75,-76,-76]


# the original data of lattitude and longitude (in degree) for the united by blue  access point ac:86:74:57:d6:82
unitedoriginsignal = -41  
unitedorigin = [39.95303359,-75.1931862]
unitedlat = [39.9531262,39.95306805,39.95303359,39.95300903,39.95308204,39.95307472,
        39.95305298,39.95304666,39.95309792,39.95311104,39.95314332,39.95315069,
        39.95318367,39.95316839,39.9529689,39.95295116,39.95298173,39.95299407,
        39.95297741,39.9530528,39.95305761,39.95311303,39.95304962,39.95300972,
        39.95293554,39.95287939]
unitedlon = [-75.19316856,-75.19322581,-75.1931862,-75.19329689,-75.19333511,
        -75.19350965,-75.19355542,-75.19358466,-75.19371891,-75.19376055,
        -75.19378227,-75.19381858,-75.1938578,-75.19387127,-75.19327191,
        -75.1932809,-75.19337015,-75.19344377,-75.19354724,-75.19362915,
        -75.19366132,-75.19377652,-75.19363119,-75.19348948,-75.19321491,
        -75.19300501]
unitedsignal = [-60,-61,-41,-59,-62,-64,-77,-69,-79,-75,-81,-82,-82,-82,
        -74,-70,-66,-64,-76,-67,-76,-74,-71,-71,-69,-79]


# the original data of lattitude and longitude (in degree) for the bluemercury  access point e2:55:7d:46:ba:70
blueoriginsignal = -42 
blueorigin = [39.95348435,-75.1954269]
bluelat = [39.95344475,39.95346047,39.95351389,39.9535116,39.95351465,39.95348435,
        39.95350455,39.95353418,39.95356256,39.95353174,39.95344449,39.95345135,
        39.95342285,39.95340621,39.95338485,39.9533588,39.95336648,39.95336993,
        39.95327412,39.95328151,39.95328667,39.95329931,39.95333245,39.95343911,
        39.95343569,39.95322091,39.95330729,39.95339435,39.95330516]
bluelon = [-75.1953241,-75.19534486,-75.19536677,-75.19539087,-75.19540162,
        -75.1954269,-75.19552856,-75.19554469,-75.1956038,-75.19586533,
        -75.19582356,-75.19575759,-75.19569923,-75.19565746,-75.19551169,
        -75.19546526,-75.19540467,-75.19533335,-75.1955647,-75.19562271,
        -75.19567539,-75.19576064,-75.19582484,-75.19562141,-75.19557638,
        -75.19541343,-75.19556592,-75.19585788,-75.19589166]
bluesignal = [-67,-68,-62,-62,-56,-42,-61,-61,-68,-69,-68,-69,-57,-69,-60,-60,
        -72,-72,-65,-65,-76,-71,-70,-69,-69,-66,-66,-70,-71]

# the original data of lattitude and longitude (in degree) for the Penne-bar access point e8:40:40:80:52:e8
penneoriginsignal = -43
penneorigin = [39.95339425,-75.1960868] 
pennelat = [39.95337588,39.95337808,39.95334292,39.95331375,39.95329611,39.953289,
        39.95329387,39.95333236,39.95327273,39.95336876,39.95344523,39.95342261,
        39.9534374,39.9534569,39.95350044,39.95345781,39.9534526,39.95336277,
        39.95336332,39.95340967,39.95338347,39.95337122,39.95338304,39.95338875,
        39.95339425,39.95339517,39.95340651,39.95343159,39.95349377]
pennelon = [-75.19609743,-75.19618277,-75.19610598,-75.19602633,-75.19600751,
        -75.19598095,-75.1959173,-75.19589952,-75.19608164,-75.19570731,-75.19592838,
        -75.19593443,-75.19599858,-75.19608773,-75.19651072,-75.19636381,-75.19627185,
        -75.19591907,-75.19594091,-75.19603491,-75.19618058,-75.19617366,-75.19612399,
        -75.19603628,-75.1960868,-75.19612735,-75.19616441,-75.19620826,-75.19634645]
pennesignal = [-70,-77,-72,-80,-72,-76,-70,-70,-80,-70,-64,-82,-71,-72,-82,-78,-69,
        -65,-64,-64,-71,-70,-71,-63,-43,-66,-66,-76,-79]
#=============================================================================
# helper functions
#=============================================================================

def discalc(lat,lon,origin):
    """Calculate the distance (unit = meter) of all locations from the origin of one WAP"""
    # calculation following online sourse:http://www.movable-type.co.uk/scripts/latlong.html
    # using the Harversine formula
    R = 6373  #unit = km the radius of earth, reference: http://andrew.hedges.name/experiments/haversine/
    # create a list to store all the distances
    Distance = []
    for idx in range(len(lat)):
        #convert degree to radian using numpy
        deltalat = np.radians(lat[idx]-origin[0])
        deltalon = np.radians(lon[idx]-origin[1])
        rlat1 = np.radians(origin[0])
        rlat2 = np.radians(lat[idx])
        a = (math.sin(deltalat/2))**2 + math.cos(rlat1)*math.cos(rlat2)*(math.sin(deltalon/2))**2
        c = 2*math.atan2(math.sqrt(a),math.sqrt(1-a))
        d = R* c * 1000
        Distance.append(d)

    return Distance

#-----------------------------------------------------------------------------
#unit converter
def dbmtomw(signal):
    """Convert dBm to mW"""
    return 10.0**(signal/10.0) 

def mwtodbm(power):
    """Convert mW to dBm"""
    return 10.0*math.log10(power/1.0)
#-----------------------------------------------------------------------------
#calculate Rpower (the "ideal power") in free space
def Rpower(dis,originsignal):
    """Calculate the ideal power for each access point and store
    them into a list"""
    rpower = []
    for r in dis:
        if r == 0:
            rpower.append(originsignal)
        else:
            result = dbmtomw(float(originsignal)) * (1.0)/(1.0*r**2)
            result = mwtodbm(float(result))
            rpower.append(result)
    return rpower
#-----------------------------------------------------------------------------


#=============================================================================
#Plot locations direct onto google map as for better visualization
#=============================================================================

#reference:https://github.com/vgm64/gmplot
#-----------------------------------------------------------------------------
# gmap = gmplot.GoogleMapPlotter(39.95211, -75.19368, 17)

# # gmap.plot(starlat, starlon, 'cornflowerblue', edge_width=10)
# gmap.scatter(liblat, liblon, '#3B0B39', size=5, marker=False)
# gmap.scatter(liblat, liblon, 'r', marker=True)
# gmap.draw("lib.html")

def gdraw(lat,lon,origin,name):
    """Draw all the measured location as scatters (based on 
    their lat and long) onto google map and output the 
    result as a html file"""
    gmap = gmplot.GoogleMapPlotter(origin[0], origin[1], 17)
    gmap.scatter(lat, lon, '#3B0B39', size=2, marker=False)
    # gmap.scatter(lat, lon, 'r', marker=True)
    gmap.draw(name)

#=============================================================================
#Plot Heatmap
#=============================================================================

def plotheatmap(lat,lon,origin,signal,name):
    """Plot the heap map for each WAP"""
    fig, ax = plt.subplots(1)
    maxs =1.000* max(signal)
    # color = [str(maxs/(1.00*item)) for item in signal]
    color = [(1.00*item)/(100.0) for item in signal]
    plt.scatter(lon,lat,s = 400,marker = '8', c = color,alpha=0.5)

    ax.scatter(origin[1],origin[0],marker = 'd',color = 'g', s = 60,alpha = 1,label = 'Wifi Access Point')
    ax.legend(scatterpoints = 1)

    # plt.gray()
    cbar = plt.colorbar()
    plt.ylim(min(lat)-1E-4,max(lat)+1E-4)
    plt.xlim(min(lon)-1E-4,max(lon)+1E-4)
    for i, txt in enumerate(signal):
        ax.annotate(txt,(lon[i],lat[i]),fontsize = 10)
    xlabel = r'Longitude'
    ylabel = r'Latitude'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    filename = name + '.png'
    titlename = name + ' Wifi heatmap'
    plt.title(titlename)
    fig.savefig(filename,dpi=300)

#=============================================================================
#Plot Power VS Distance
#=============================================================================
def plotpower(dis,Rpower,signal,name):
    """Plot Power VS Distance as to compare the measured average bulk propagation loss
    with the ideal propagation loss for a free space"""

    fig2 = plt.figure()
    plt.plot(dis,Rpower,'^',color='r',alpha = 0.5, label = r'$Ideal\ signal\ \backsim$ 1/r^2')
    plt.plot(dis,signal,'o',color='g', alpha = 0.5, label = r'$Measured\ signal$')
    ylabel = 'Signal Power; unit = dBm'
    plt.ylabel(ylabel)
    xlabel = 'distance from origin; unit = meters'
    plt.xlabel(xlabel)
    titlename = name + "'s Wifi Power VS distance from WAP origin"
    filename = name +"power.png"
    plt.legend()
    plt.title(titlename)
    fig2.savefig(filename, dpi=300, bbox_inches='tight')

#=============================================================================
#Plot Delta Power VS Distance
#=============================================================================
def deltaplotpower(dis,Rpower,signal,originsignal,name):
    """Plot Power VS Distance as to compare the measured average bulk propagation loss
    with the ideal propagation loss for a free space"""
    for i in range(len(dis)):
        Rpower[i] = originsignal - Rpower[i]
        signal[i] = originsignal - signal[i]

    fig2 = plt.figure()
    plt.plot(dis,Rpower,'^',color='r',alpha = 0.5, label = r'$Ideal\ signal\ \backsim$ 1/r^2')
    plt.plot(dis,signal,'o',color='g', alpha = 0.5, label = r'$Measured\ signal$')
    ylabel = 'Propagation Loss; unit = dBm'
    plt.ylabel(ylabel)
    xlabel = 'distance from origin; unit = meters'
    plt.xlabel(xlabel)
    titlename = name + "'s Porpapgation loss VS distance (r)"
    filename = name +"deltapower.png"
    plt.legend(loc = 4)
    plt.title(titlename)
    fig2.savefig(filename, dpi=300, bbox_inches='tight')

#=============================================================================
# calculate k value the average propagation loss
#=============================================================================
# Assume noise power to be C  = -100dBm
def avepropaloss(signal,dis,Rpower,name):
    # calculation follows Piazza post @64
    C = -100
    K = []
    R = []
    ronly = []
    konly = []
    for i in range(len(signal)):
        snr = 10*math.log10(dbmtomw(signal[i])/dbmtomw(C)) #unit = dB
        rsnr = 10*math.log10(dbmtomw(Rpower[i])/dbmtomw(C)) #unit = dB
        if dis[i] !=0:
            r = rsnr*dis[i]**2
            k = snr*dis[i]**2
            lst = [k,dis[i],r]
            K.append(lst)
            ronly.append(r)
            konly.append(k)
    ave = sum(konly)/float(len(konly))
    rave = sum(ronly)/float(len(ronly))
    
    """out put the average propagation loss"""
    outfile = open(name,'w')
    outfile.write("Average propagation loss is "+ str(ave) +" dB*m^2\n")
    outfile.write("Average propagation loss according to 1/r^2 is " + str(rave) +" dB*m^2\n")
    outfile.write("For each location, the propagation loss k and the 1/r^2 prediction k_r is shown below along with the distance to the origin WAP\n")

    for data in K:
        content = "k = "+ str(round(data[0])) +" dB*m^2;\t"+ "k_r = "+str(round(data[2])) +" dB*m^2;\t"+"dis = "+ str(data[1])+ " m\n"
        outfile.write(content)
    outfile.close()

    return K

#=============================================================================
# write the outputdata into a text file
#=============================================================================
def outputdata(lat,lon,signal,name,bssid,dis,rpower):
    """write the output data into a text file"""
    outfile = open(name,'w')
    d = []
    for i in range(len(lat)):
        t = (lat[i],lon[i],signal[i],dis[i],rpower[i])
        d.append(t)
    outfile.write("BSSID " + "  \t\t" + "Latitude" + " \t" + "longitude" + " \t" + "Power" + " \t" +"Dis"+" \t"+"Rpower" + "\n")

    for data in d: 
        content = bssid + "\t" + str(data[0])+ "\t"+ str(data[1]) +"\t" + str(data[2]) +"\t" + str(round(data[3],2)) + "\t"+str(round(data[4],2))+"\n"
        outfile.write(content)

    outfile.close()



#-----------------------------------------------------------------------------
def main():

    # location visualization
    # locate each location onto googlemap one by one as html file
    gdraw(liblat,liblon,liborigin,"library.html")
    gdraw(starlat,starlon,starorigin,"starbucks.html")
    gdraw(unitedlat,unitedlon,unitedorigin,"united.html")
    gdraw(bluelat,bluelon,blueorigin,"bluemercury.html")
    gdraw(pennelat,pennelon,penneorigin,"penne.html")

    # locate all locations onto googlemap as one html file 
    gmap = gmplot.GoogleMapPlotter(liborigin[0], liborigin[1], 17)
    gmap.scatter(liblat, liblon, '#3B0B39', size=2, marker=False)
    gmap.scatter(starlat,starlon,'g',size = 2 ,marker = False)
    gmap.scatter(bluelat,bluelon,'b',size = 2 ,marker = False)
    gmap.scatter(unitedlat,unitedlon,'r',size = 2 ,marker = False)
    gmap.scatter(pennelat,pennelon,'m',size = 2 ,marker = False)
    gmap.draw("all.html")

    # plot the heatmap
    plotheatmap(liblat,liblon,liborigin,libsignal,"library")
    plotheatmap(starlat,starlon,starorigin,starsignal,"starbucks")
    plotheatmap(unitedlat,unitedlon,unitedorigin,unitedsignal,"unitedblue")
    plotheatmap(bluelat,bluelon,blueorigin,bluesignal,"bluemercury")
    plotheatmap(pennelat,pennelon,penneorigin,pennesignal,"penne")

    # calculate the distance:
    libdis = discalc(liblat,liblon,liborigin)
    stardis = discalc(starlat,starlon,starorigin)
    unitedis = discalc(unitedlat,unitedlon,unitedorigin)
    bluedis = discalc(bluelat,bluelon,blueorigin)
    pennedis = discalc(pennelat,pennelon,penneorigin)
    
    # # calculate Rpower
    Rlib = Rpower(libdis,liboriginsignal) 
    Rstar = Rpower(stardis,staroriginsignal) 
    Runited = Rpower(unitedis,unitedoriginsignal) 
    Rblue = Rpower(bluedis,blueoriginsignal) 
    Rpenne = Rpower(pennedis,penneoriginsignal) 

    #write the data for phase 2 in to a text file
    outputdata(liblat,liblon,libsignal,"library.txt","84:d4:7e:dc:2a:40",libdis,Rlib)
    outputdata(starlat,starlon,starsignal,"starbucks.txt","94:b4:0f:53:07:00",stardis,Rstar)
    outputdata(unitedlat,unitedlon,unitedsignal,"unitedblue.txt","ac:86:74:57:d6:82",unitedis,Runited)
    outputdata(bluelat,bluelon,bluesignal,"bluemercury.txt","e2:55:7d:46:ba:70",bluedis,Rblue)
    outputdata(pennelat,pennelon,pennesignal,"penne.txt","e8:40:40:80:52:e8",pennedis,Rpenne)

    #Plot distance vs signal strength
    plotpower(libdis,Rlib,libsignal,"Library")
    plotpower(stardis,Rstar,starsignal,"Starbucks")
    plotpower(unitedis,Runited,unitedsignal,"Unitedblue")
    plotpower(bluedis,Rblue,bluesignal,"Bluemercury")
    plotpower(pennedis,Rpenne,pennesignal,"Penne")

    #plot the distance vs propagation loss
    deltaplotpower(libdis,Rlib,libsignal,liboriginsignal,"Library")
    deltaplotpower(stardis,Rstar,starsignal,staroriginsignal,"Starbucks")
    deltaplotpower(unitedis,Runited,unitedsignal,unitedoriginsignal,"Unitedblue")
    deltaplotpower(bluedis,Rblue,bluesignal,blueoriginsignal,"Bluemercury")
    deltaplotpower(pennedis,Rpenne,pennesignal,penneoriginsignal,"Penne")

    #calculate k
    avepropaloss(libsignal,libdis,Rlib,'libaveloss.txt')
    avepropaloss(starsignal,stardis,Rstar,'staraveloss.txt')
    avepropaloss(unitedsignal,unitedis,Runited,'unitedaveloss.txt')
    avepropaloss(bluesignal,bluedis,Rblue,'blueaveloss.txt')
    avepropaloss(pennesignal,pennedis,Rpenne,'penneaveloss.txt')


main()


