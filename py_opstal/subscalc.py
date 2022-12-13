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
def read_shapefile(shapefolder, shapefl):
    '''
    Read shapefile into shapely MultiPolygon
    Allow for multiple shapes (i.e.shapefl can be a list of shapefile names)
    '''
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
def calc_subsidence(shape, dxy, depth, rb, tComp, snap_to_grid = True):
    '''
    Calculate subsidence from compaction 'tComp'
    '''
    # Effective volume of grid cell
    tCompV = tComp*dxy*dxy

    ###########################################
    # Bounding box
    (xmin,ymin,xmax,ymax)=shape.bounds
    #print(xmin, xmax)
    #print(ymin, ymax)

    # Snap to grid (makes combining easier)
    if (snap_to_grid):
        xmin=dxy*int(xmin/dxy)
        ymin=dxy*int(ymin/dxy)
        xmax=dxy*int(xmax/dxy+1)
        ymax=dxy*int(ymax/dxy+1)
    
    ###########################################
    # Grid points inside the polygon(s)
    x=[]
    y=[]
    nx=int((xmax-xmin)/dxy+1)
    ny=int((ymax-ymin)/dxy+1)
    #print("    nx,ny=", nx,ny)
    for ii in range(nx):
        for jj in range(ny):
            rx=xmin+ii*dxy
            ry=ymin+jj*dxy
            # Check if point is inside
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
    print("    Starting the calculation...")
    import time
    start_time = time.time()
    subsidence=VO.vanopstal(x,y,depth,depth*rb,xpm,ypm,tCompV)
    subsidence *= 1000 # m-->mm
    dt = (time.time() - start_time)
    print("        CPU time:", dt)
    print("        max. subs.", np.amax(subsidence))
    print("        compact   ", tComp*1000)
    print("        ---------------")
    print("        Volume = ", np.sum(subsidence)/1000*dxy*dxy, "m3")
    print("        max/comp = ",np.amax(subsidence)/(tComp*1000)*100,"%" )

    
    # Convert 1D back to 2D
    X2=np.reshape(xpm, (nxg,nyg))
    Y2=np.reshape(ypm, (nxg,nyg))
    S2=np.reshape(subsidence, (nxg,nyg))
    
    print()
    
    return X2, Y2, S2

def plot_subs_countours(X2, Y2, S2, shapes, epsg=23031, c_interval=10, first_level=None):
    '''
    Plot countours of S2, assuming grid is in X2,Y2
    (all2D numpy arrays).
    A background map is added (OSM).
    ESPG code 23031 is ED50/UTM-31-NS.
    S2, first_level and c_interval are in mm.
    if first_level is not provided, it is set equal to c_interval.
    
    returns ax
    '''
    import matplotlib
    import matplotlib.pyplot as plt
    
    if (not isinstance(shapes, list)):
        shapes=[shapes]
        
    # Init plot
    fig, ax = plt.subplots()

    # Add field shape if provided
    for shape in shapes:
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
    
    # Smaller fonts work better
    ax.tick_params(axis='both', which='major', labelsize=5)

    # We're done
    print()

    # To be able to add stuff to the plot, if desired
    return ax

#################################################################    
def get_merged_grid(xs, ys):
    '''
    Determine the grid to use to merge a the grids in lists xs and ys.
    Assumes all grids snapped to integer multiples of dxy

    returns X2big, Y2big, di, dj
    
    (di[i] and dj[i] are integer shifts needed to map array #i to the 
    merged grid.
    '''
    # Of arrays provided
    nshape = len(xs)
    assert(len(ys)==nshape)
    
    # Assume dx, dy are the fixed & same for all grids
    dx = xs[0][1,0] - xs[0][0,0]
    dy = ys[0][0,1] - ys[0][0,0]
    
    xmin = 1e37
    xmax = -1e37
    ymin = 1e37
    ymax = -1e37
    for i in range(nshape):
        ydims = xs[i].shape
        xmin = min(xmin, xs[i][0,0],xs[i][ydims[0]-1, ydims[1]-1])
        xmax = max(xmax, xs[i][0,0],xs[i][ydims[0]-1, ydims[1]-1])
        ymin = min(ymin, ys[i][0,0],ys[i][ydims[0]-1, ydims[1]-1])
        ymax = max(ymax, ys[i][0,0],ys[i][ydims[0]-1, ydims[1]-1])

    # Init the X/Y arrays to the merged size
    nxbig = int((xmax-xmin)/dx + 1)
    nybig = int((ymax-ymin)/dy + 1)
    X2big=np.zeros([nxbig,nybig])
    Y2big=np.zeros([nxbig,nybig])
    for i in range(nxbig):
        for j in range(nybig):
            X2big[i,j]=xmin+i*dx # clumsy
            Y2big[i,j]=ymin+j*dy # clumsy

    # Determine shifts from small grid to large one
    di = []
    dj = []
    for i in range(nshape):
        # This should be an integer to begin with (all grids snapped to
        # multiples of dx, dy), but be wary of roundoff
        di.append(int((xs[i][0,0] - xmin)/dx + 0.5))
        dj.append(int((ys[i][0,0] - ymin)/dy + 0.5))
    
    # Return results
    return X2big, Y2big, di, dj

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

