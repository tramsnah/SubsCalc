import os
import numpy as np
import pandas as pd
import shapefile # This version uses PyShp rather than geopandas
from shapely.geometry import Polygon, MultiPolygon, Point

if __name__ == "__main__": # ugly
    import vanopstal as VO
else:
    from . import vanopstal as VO

###########################################
# Read shapefile into polygon
def read_shapefile(shapefolder, shapefl):
    # Allow for multiple shapes
    if (isinstance(shapefl, str)):
        shapefl = [shapefl]

    # collate the shapes
    shape=[]
    for shapef in shapefl:
        # There may be spaces...
        shapef = shapef.strip()

        # Check extension
        filename, file_extension = os.path.splitext(shapef)
        
        # Check if path ends in slash
        if (shapefolder[-1] != "\\"):
            shapefolder += "\\"
        
        # Compound filename
        fname = shapefolder + shapef
        print("    Loading:", fname)
        
        # Load into list of Polygon
        with shapefile.Reader(fname) as sf:
            shps = sf.shapes()
            for shp in shps:
                pts=shp.parts
                pts.append(len(shp.points)+1)
                # This is ugly
                for i in range(len(pts)-1):
                    #print(pts[i],pts[i+1])
                    p = Polygon(shp.points[pts[i]:pts[i+1]])
                    shape.append(p)
                    
    # And convert to multipolygon
    shape=MultiPolygon(shape)
    
    print()
    
    return shape

############################################################
# Calculate subsidence
def calc_subsidence(shape, dxy, depth, rb, tComp):
    ###########################################
    # Grid resolution
    dxy=200 # m
    
    # Effective volume of grid cell
    tCompV = tComp*dxy*dxy

    ###########################################
    # Bounding box
    (xmin,ymin,xmax,ymax)=shape.bounds
    #print(xmin, xmax)
    #print(ymin, ymax)

    ###########################################
    # Grid points inside the polygon(s)
    x=[]
    y=[]
    nx=int((xmax-xmin)/dxy+1)
    ny=int((ymax-ymin)/dxy+1)
    print("nx,ny=", nx,ny)
    for ii in range(nx):
        for jj in range(ny):
            rx=xmin+ii*dxy
            ry=ymin+jj*dxy
            # circular reservoir
            inside=shape.contains(Point(rx,ry))
            if (inside):
                x.append(rx)
                y.append(ry)
    x=np.array(x)
    y=np.array(y)

    ###########################################
    # Grid for calculation
    # Leave room for the slope of the side of the bowl
    nxy_xtra=int(depth*3/dxy)+1
    xming=xmin-nxy_xtra*dxy
    yming=ymin-nxy_xtra*dxy
    nxg=nx+2*nxy_xtra
    nyg=ny+2*nxy_xtra
    xpm=np.zeros(nxg*nyg)
    ypm=np.zeros(nxg*nyg)
    kk=0
    for ii in range(nxg):
        for jj in range(nyg):
            xpm[kk]=xming+ii*dxy
            ypm[kk]=yming+jj*dxy
            kk+=1

    ###########################################
    # Do the deed
    print("starting the calculation...")
    import time
    start_time = time.time()
    subsidence=VO.vanopstal(x,y,depth,depth*rb,xpm,ypm,tCompV)
    subsidence *= 1000 # m-->mm
    dt = (time.time() - start_time)
    print("CPU time:", dt)
    print("max. subs.", np.amax(subsidence))
    print("compact   ", tComp*1000)
    print("---------------")
    print("max/comp = ",max(subsidence)/(tComp*1000)*100,"%" )

    # Convert 1D back to 2D
    X2=np.reshape(xpm, (nxg,nyg))
    Y2=np.reshape(ypm, (nxg,nyg))
    S2=np.reshape(subsidence, (nxg,nyg))
    
    print()

    return X2, Y2, S2

###########################################
# Plot
# ESPG code 23031 is ED50/UTM-31-NS
# S2 and c_interval is in mm
def plot_subs_countours(X2, Y2, S2, shape, epsg=23031, c_interval=10, first_level=None):
    import matplotlib
    import matplotlib.pyplot as plt

    # Init plot
    fig, ax = plt.subplots()

    # Add field shape if provided
    if (shape is not None):
        for s in shape.geoms:    
            xs, ys = s.exterior.xy
            ax.plot(xs,ys,color="black",linewidth=0.3)
        
    # Determine contour levels & calculate them
    if (first_level is None):
        first_level = c_interval
    max_subs = np.amax(S2)
    n_levels = int((max_subs-first_level)/c_interval)+1
    levels = range(first_level, first_level + c_interval*n_levels, c_interval)
    print("Calculating contours...")
    CS = ax.contour(X2, Y2, S2, levels)

    # Plot the plot
    ax.clabel(CS, CS.levels, inline=True)

    # Base map
    import contextily as cx
    from pyproj import CRS
    crs = CRS.from_user_input(epsg)
    cx.add_basemap(ax, crs=crs);

    print()

    return ax
    
if __name__ == "__main__":
    import shapefile # This version uses PyShp rather than geopandas
    from shapely.geometry import Polygon, MultiPolygon, Point

    ###########################################
    # Read shapefile into polygon
    shape = read_shapefile("test_data", "test.shp")

    ###########################################
    # Parameters
    depth=1886.74 # m depth
    rb=1.2 # fraction
    Cm=1.41e-5 # 1/bar
    thick=24.3 # m
    Pini=206 # bar
    Paban=15 # bar
    ntg = 0.7

    # Net thickness
    thick *= ntg

    # Pressure drop to abandonment pressure
    dP=Pini-Paban
            
    # put into compound thickness multiplier
    tComp=thick*dP*Cm;

    # Grid resolution
    dxy=200 # m

    # Do the calc
    X2, Y2, S2 = calc_subsidence(shape, dxy, depth, rb, tComp)

    ###########################################
    # Plot
    import matplotlib.pyplot as plt

    ax = plot_subs_countours(X2, Y2, S2, shape)

    plt.show()