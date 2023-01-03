'''
Supporting code for handling subsidence cases. These are stored in an excel
sheet (database). The case contains geomechanical parameters, shapes,
profiles.

Example:
    Case	    Parameter	Value
    test_case_1	depth	    2985
    test_case_1	rb	        1.05
    test_case_1	cm	        7.00E-06
    test_case_1	thickness	212
    test_case_1	ntg	        0.5
    test_case_1	pini	    344
    test_case_1	paban	    89
    test_case_1	tau	        3
    test_case_1	shape	    test, test_2, test_3, test_4,test_5, test_6
    test_case_1	shapefolder	.\test_data\
    test_case_1	comment	    Test case for merged time dependence. Assume ROSLL is undepleted
    test_case_1	profile_h	A HIST 21, D6 HIST 21, C Hist 21, D8 HIST 21,,
    test_case_1	profile_fc	A FC 21 2P, D6 FC 21 2P, C FC 21 2P, D8 FC 21 2P, D9 FC 21 2P, B FC 21 2P
    test_case_1	IGIP	    9.3, 1.01, 1.4, 2.6,0.9,0.9

In the xls file, this is expected to be stored in the tab 'Parameters'.
'''
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# UGLY! One import from extrernal source, one if this file is run as-is (for testing)
try: # ugly, but we want it part of the library as well as stand-alone test
    from . import vanopstal as VO
    from . import profile_utils as PU
    from . import subscalc as SC
except ImportError:
    import vanopstal as VO
    import profile_utils as PU
    import subscalc as SC

def _get_float_parameter(df, key, defval=None, n=None):
    '''
    Get float parameter from database. Returns a list.
    If one parameter is supplied for all n cases,
    n identical values are in the list.
    '''
    try:
        val = df[df["Parameter"]==key]["Value"].values[0]     # m depth
        if (isinstance(val, str)):
            val = val.split(',')  # Can have multiple, if comma-separated
            val = [float(g) for g in val] # Convert to float
    except IndexError:
        val = defval

    if (not isinstance(val, list)):
        val=[val]
        
    # Check # of parameters is OK
    if (n is not None):
        if not (len(val)==1 or len(val)==n):
            print("# of data items provided for ", key,  " is ", 
                    len(val), ", but should be 1 or ", n,"\n")
            assert(0)
        if (len(val) == 1):
            val = val*n
        
    return val

def _get_string_parameter(df, key, defval=None, n=None):
    '''
    Get string parameter from database. Returns a list.
    If one parameter is supplied for all n cases,
    n identical values are in the list.
    '''
    try:
        val = df[df["Parameter"]==key]["Value"].values[0]     # m depth
        val = val.split(',')  # Can have multiple, if comma-separated
        val = [v.strip() for v in val] 
    except IndexError:
        val = defval

    if (not isinstance(val, list)):
        val=[val]

    # Check # of parameters is OK
    if (n is not None):
        if not (len(val)==1 or len(val)==n):
            print("# of data items provided for ", key,  " is ", 
                    len(val), ", but should be 1 or ", n,"\n")
            assert(0)
        if (len(val) == 1):
            val = val*n

    return val

###########################################
def handle_subs_case(subspar_xls, case, dxy=200, cmult=1.0, par_sheet_name="Parameters", 
                    prof_sheet_name="Profiles", gas=None):
    '''
    handle case 'case' from the database 'subspar_xls"
    returns X2big, Y2big, S2big, df_out, df_subs, df_prod_all, lshapes
    
    If conversion from Gp to volume needs to take place a 'gas' object needs to be supplied, e.g.
        my_gas = PU.dry_gas(110, 0.6, 0.0, 0.01, 0.01)
    '''
    # Read Parameter DB
    print("Loading case database ", subspar_xls)
    df = pd.read_excel(subspar_xls, sheet_name=par_sheet_name)
    print()

    # Get specific case data
    df = df[df["Case"]==case]
    print("Case '"+case+"' parameters:")
    print(df)
    print()
    if (len(df)==0):
        print("No data found for case '"+case+"'\n")
        assert(0)

    # Get list of shapes. This determines the # of subcases we need to run
    shapef = _get_string_parameter(df, "shape")
    nshape = len(shapef)

    # Base parameters
    shapefolder = _get_string_parameter(df, "shapefolder",".", n=nshape)
    depth = _get_float_parameter(df, "depth", n=nshape)           # m depth
    rb = _get_float_parameter(df, "rb", 1, n=nshape)              # fraction
    Cm = _get_float_parameter(df, "cm", n=nshape)                 # 1/bar
    thick = _get_float_parameter(df, "thickness", n=nshape)       # m
    ntg = _get_float_parameter(df, "ntg", 0, n=nshape)            # fraction
    Pini = _get_float_parameter(df, "pini", n=nshape)             # bar
    Paban = _get_float_parameter(df, "paban", 0, n=nshape)        # bar
    
    # Relative path? Assume it is to my path
    def get_relative_path(fname):
        if (isinstance(fname, list)):
            ofnames = []
            for lfname in fname:
                ofnames.append(get_relative_path(lfname))
            return ofnames
            
        if (not os.path.isabs(fname)):
            me = sys.argv[0]
            path = os.path.dirname (me)  
            fname = path + "\\"+ fname
            
        return fname
    shapefolder = get_relative_path(shapefolder)

    # Profile data is optional. Consists of a historic part and a forecast part
    profileh= _get_string_parameter(df, "profile_h", n=nshape)
    profilefc= _get_string_parameter(df, "profile_fc", n=nshape)
    # IGIP (for Gp-->pressure conversion)
    IGIP = _get_float_parameter(df, "IGIP", 0, n=nshape)          # bar
    tau = _get_float_parameter(df, "tau", 5, n=nshape)            # y

    # Log data along the way
    lnames=[]
    lshapes=[]
    xs=[]
    ys=[]
    ss=[]
    ts=[]
    deepest=[]
    ddata=[]

    # Loop over the shapes
    for ishape in range(nshape):
        # Read the shapefile
        shape = SC.read_shapefile(shapefolder[ishape], shapef[ishape])

        # Do the calc for unit compaction
        X2, Y2, S2 = SC.calc_subsidence(shape, dxy, depth[ishape], rb[ishape], 1)
        
        # Log the results
        d={}
        d["shape"]=[shapef[ishape]]
        d["area"]=[shape.area] # m2
        d["depth"]=[depth[ishape]]
        d["rb"]=[rb[ishape]]
        d["cm"]=[Cm[ishape]]
        d["thickness"]=[thick[ishape]]
        d["ntg"]=[ntg[ishape]]
        d["Pini"]=[Pini[ishape]]
        d["SubsFrac"]=[np.amax(S2)/1000] # Subsidence as fraction of compaction

        # If no profile supplied, no time dependence. Third dimension = 1
        if (profileh is None):
            xdims = list(S2.shape)
            xdims.append(1)
            S2t = np.zeros(xdims)
            S2t[:,:,1] = S2 * (Pini[ishape] - Paban[ishape])* \
                        thick[ishape]*Cm[ishape]*ntg[ishape]*cmult
            
            # Log the results specific to this branch
            d["Paban"]=[Paban[ishape]]
        else:
            ###########################################
            # Read Profile
            df = pd.read_excel(subspar_xls, sheet_name=prof_sheet_name)

            # Then extract the profile
            df_prod = PU.get_merged_profile(df, profileh[ishape], 
                                profilefc[ishape], extend_year=2060)

            # Define our gas (XXX TE DOEN)
            d["gas"]=gas.name

            # Convert to pressure
            df_prod, l_p_aban, l_IGIP = \
                PU.convert_gp_to_pressure(df_prod, Pini[ishape], Paban[ishape], 
                                            IGIP[ishape], gas = gas)

            # Apply time lag
            df_prod = PU.apply_time_lag(df_prod, tau[ishape])

            # Calculate time dependent compaction
            df_prod["Comp"] = (Pini[ishape] - df_prod["Pressure_D"])*\
                        thick[ishape]*Cm[ishape]*ntg[ishape]*cmult

            # Cull to first month of the year
            df_prod_y = df_prod.iloc[::12, :].copy()

            # Multiply by compaction (time dependent), adding a timedimension 
            # to the 2D subsidence array
            nt = len(df_prod_y)
            xdims = list(S2.shape)
            xdims.append(nt)
            S2t = np.zeros(xdims)
            for i in range(nt):
                S2t[:,:,i] = S2 * df_prod_y["Comp"].values[i]
                
            # Log data
            lnames.append(shapef[ishape])
            lshapes.append(shape)
            xs.append(X2)
            ys.append(Y2)
            ss.append(S2t)
            ts.append(df_prod_y.copy()) 
            
            # Log the deepest point
            ind = np.argmax(S2t[:,:,-1], axis=None)
            ind = np.unravel_index(ind, S2t[:,:,-1].shape)
            deepest.append(ind)
            
            # Log the results specific to this branch
            d["Paban"]=[l_p_aban]
            d["Gp"]=[df_prod["Gp"].values[-1]]
            d["profileh"]=[profileh[ishape]]
            d["profilefc"]=[profilefc[ishape]]
            d["IGIP"]=[l_IGIP]
            d["tau"]=[tau[ishape]]

        # Log the generic results
        d["max_subsidence"]=[np.amax(S2t[:,:,-1])] # mm
        d["volume"]=[np.sum(S2t[:,:,-1])/1000*dxy*dxy] # m3
        ddata.append(pd.DataFrame(d))
        
        print()
        
    #####################################################################
    print("    Adding bowls")
    # Merge the X,Y grids
    X2big, Y2big, di, dj = SC.get_merged_grid(xs,ys)
    nxbig = X2big.shape[0]
    nybig = Y2big.shape[1]

    # Merge the time grid
    # Determine the time list to use. 
    # End time and intervals are the same for all, so we only need to 
    # look at # of values (or, equivalently, start time)
    ntbig=0
    tmax=None
    for i in range(nshape):
        ydims = ss[i].shape
        if (ydims[2]>ntbig):
            ntbig = ydims[2]
            tmax = ts[i]["Date"].values

    # S2 will also have a time dimension
    S2big=np.zeros([nxbig,nybig,ntbig]) 
    
    #####################################################################
    # Merge the profiles. 
    #df_prod_all = pd.DataFrame({"Date": tmax})
    for i in range(nshape):
        label = lnames[i]
        label = label.replace(".shp","")
        # Label the individual profiles
        ts[i]["Label"] = label
        ## Make sure the columns have different names
        #dft = ts[i][["Date","Gp","Pressure","Comp"]]
        #cols = list(dft.columns)
        #for j,coln in enumerate(cols):
        #    if (coln != "Date"):
        #        coln += "_"+label
        #        cols[j] = coln
        #dft.columns = cols
        #df_prod_all = pd.merge(df_prod_all, dft, how="outer", on="Date")
    df_prod_all = pd.concat(ts)
    
    # Add them together
    for i in range(nshape):
        # This should be an integer to begin with (all grids snapped to
        # multiples of dxy), but be wary of roundoff
        ydims = ss[i].shape
        dk = ntbig - ydims[2]
        S2big[di[i]:di[i]+ydims[0],dj[i]:dj[i]+ydims[1],dk:] += ss[i]
        
        # Also shift the indices of the deepest point. Log in list (for convienience), as 
        # well as in dataframe
        (i0,j0) = deepest[i]
        deepest[i] = (i0+di[i],j0+dj[i])
        ddata[i]["Deepest"]="s" # Init to non-float, so we can store tuple
        ddata[i].at[0,"Deepest"]=deepest[i]
        
    print()
    
    ###########################################
    # Plot the final state
    ax = SC.plot_subs_countours(X2big, Y2big, S2big[:,:,-1], lshapes, title=case+" (final)")
    plt.show()
    
    # Get time dependence, if provided
    if (profileh is not None):
        # Also deepest point of total bowl
        if (nshape>1):
            ind = np.argmax(S2big[:,:,-1], axis=None)
            ind = np.unravel_index(ind, S2big[:,:,-1].shape)
            deepest.append(ind)
            lnames.append(case)
            
            d={}
            d["shape"]=["Total"]
            d["max_subsidence"]=[np.amax(S2big[:,:,-1])] # mm
            d["volume"]=[np.sum(S2big[:,:,-1])/1000*dxy*dxy] # m3
            ddata.append(pd.DataFrame(d))
            ddata[-1]["Deepest"]="s" # Init to non-float, so we can store tuple
            ddata[-1].at[0,"Deepest"]=(ind[0],ind[1])

        # Get time series at all the deepest points
        df_subs = pd.DataFrame({"Date": tmax})
        for i in range(len(deepest)):
            label=lnames[i]
            label = label.replace(".shp","")
            ind=deepest[i]
            cst = S2big[ind[0],ind[1],:]
            cname="Subs_"+label
            df_subs[cname]= cst
            ax = df_subs.plot(x="Date",y=cname)
            ax.set_title(label)
            ax.invert_yaxis()
            ax.set_ylabel('Subsidence [mm]')
            plt.show()
            
        # Calculate the volumes
        nt = len(tmax)
        v1 = 0 
        df_subs["Volume"] = 0
        df_subs["dV/dt"] = 0
        for i in range(1,nt):
            v0 = v1
            v1 = np.sum(S2big[:,:,i])/1000*dxy*dxy # m3
            df_subs["Volume"].values[i] = v1
            df_subs["dV/dt"].values[i] = (v1-v0) # m3/y, since table has been resampled to yearly
        
        ax = df_subs.plot(x="Date",y="dV/dt")
        ax = df_subs.plot(x="Date",y="Volume", secondary_y=True, ax=ax)
        plt.show()

    # Collate data for all the shapes processed into 1 dataframe
    df_out=pd.concat(ddata, ignore_index=True)
    
    return X2big, Y2big, S2big, df_out, df_subs, df_prod_all, lshapes

###################################################################
# Test code
if __name__ == "__main__":
    # Relative path? Assume it is to my path
    def get_relative_path(fname):
        if (not os.path.isabs(fname)):
            me = sys.argv[0]
            path = os.path.dirname (me)  
            fname = path + "\\"+ fname
        return fname

    ############################################################
    # Master data
    subspar_xls = get_relative_path(r".\test_data\test_data.xlsx")
    dxy = 200 # m
    
    my_gas = PU.dry_gas("Test", 110, 0.6, 0.0, 0.01, 0.01)

    ############################################################
    # test case
    case="test_case_1"

    # Do the deed
    X2big3, Y2big3, S2big3, df_out, df_subs, df_prod_all, lshapes3 = \
            handle_subs_case(subspar_xls, case, dxy=dxy, gas=my_gas)
    print(df_out)

    # Also plot current
    t_idx_now = 29
    ax = SC.plot_subs_countours(X2big3, Y2big3, S2big3[:,:,t_idx_now], lshapes3)
    ts_now = (df_prod_all["Date"].values[t_idx_now])
    ts_now = np.datetime_as_string(ts_now, unit='D')
    ax.set_title("Current ("+ts_now+"), test case")
    plt.show()

    # Write key outputs to excel
    xls_out = case.replace(" ","_").replace("/","_").replace("\\","_") + "_out.xlsx"
    print("Exporting key data", xls_out)
    with pd.ExcelWriter(xls_out) as writer:  
        df_out.to_excel(writer, sheet_name="Parameters")
        df_subs.to_excel(writer, sheet_name="Deepest_Subs")
        df_prod_all.to_excel(writer, sheet_name="Profiles")
    print()
