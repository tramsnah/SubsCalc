import numpy as np
import pandas as pd
from datetime import datetime

if __name__ == "__main__": # ugly
    import pvtcorrelation as PVT
else:
    from . import pvtcorrelation as PVT

class dry_gas:
    def __init__(self, temp, sg, x_h2s, x_co2, x_n2):
        # degF to degC
        self._temp = temp*9/5+32
        self._sg = sg
        self._x_h2s = x_h2s
        self._x_co2 = x_co2
        self._x_n2 = x_n2

    def zfactor(self, pressure):
        # bar to psi
        pressure *=14.5038
        
        # calculate gas pseudoproperties 
        P_pc, T_pc, P_pr, T_pr = PVT.gas_pseudoprops(self._temp, pressure, self._sg, self._x_h2s, self._x_co2)
        
        # calculate gas z-factor 
        pseudo_rho, z_factor = PVT.gas_zfactor(T_pr, P_pr)
        return z_factor
        
    def p_from_pz(self, pz):
        # Use successive substitution
        p = pz
        p1 = 1e30
        while (abs(p-p1)>0.01):
            p1 = p
            z = self.zfactor(p)
            p = pz*z

        return p
        

def get_merged_profile(df, c_hist, c_fc, extend_year=2050, force_january_start=True):
    '''
    Get profiles with keys 'c_hist' and 'c-fc', merge them,
    resample them to monthly, and extend them to 31-12 of the
    year extend_year, if needed.
    Dataframe structure assumed is
        index Profile Date Gp
    '''
    # Get the two profiles from the overarching df
    print("    Loading profiles ", "'"+c_hist+"', '" +c_fc+"'")
    #print(df)
    df_hist=None
    if (c_hist is not None) and (c_hist != ""):
        df_hist=df[df["Profile"]==c_hist]
        if (len(df_hist) == 0): 
            df_hist = None
    df_fc=None
    if (c_fc is not None) and (c_fc != ""):    
        df_fc=df[df["Profile"]==c_fc]
        if (len(df_fc) == 0): 
            df_fc = None

    # Merge them. if needed\
    if (df_hist is not None) and (df_fc is not None):
        # Need to find the final history date 
        # that matches a date in the FC
        # They need to match *exactly* (XXX: be less sensitive)
        i = 0
        while(True):
            i += 1
            final_date = df_hist["Date"].values[-i]
            final_gp = df_hist["Gp"].values[-i]
            df_fc_suture = df_fc[df_fc["Date"] == final_date]
            if (len(df_fc_suture)>0): break
        indx_fc_suture = df_fc_suture.index.values[0]
        gp_fc_suture=df_fc_suture["Gp"].values[0]

        # Small adaptation may be needed to merge them seamlessly
        dgp = final_gp-gp_fc_suture

        # Then merge
        df_prod = pd.concat([df_hist.iloc[:-i,:], df_fc.iloc[indx_fc_suture:,:]])
    elif (df_fc is not None):
        print("        No history data found ", "'"+c_hist+"'")
        df_prod = df_fc.copy()
    elif (df_hist is not None):
        print("        No forecast data found ", "'"+c_fc+"'")
        df_prod = df_hist.copy()
    else:
        print("        No profile data found ", "'"+c_hist+"', '" +c_fc+"'")
        assert(0)

    # Extend it to desired date, if needed
    t_end = pd.to_datetime('31-dec-'+str(int(extend_year)))
    if (t_end > df_prod['Date'].values[-1]):
        #print("extend to", t_end)
        final_gp_fc = df_prod['Gp'].values[-1]
        df_extend = pd.DataFrame({"Profile": ["Extension"], "Date": [t_end], "Gp": [final_gp_fc]})
        df_prod = pd.concat([df_prod, df_extend])

    # Prepend zeros, if needed, so it starts in January
    if (force_january_start):
        # Get first date
        d0 = df_prod["Date"].values[0]
        # Get month number of first date
        m0 = np.datetime64(d0,'M').astype(int) % 12
        ts=[]
        gps=[]
        while (m0 != 0):    
            # Shift to first day of month
            d1 = np.datetime64(d0, 'M')
            d1 = np.datetime64(d1, 'D')
            # Then subtract one more day (convention is EOM)
            d1 -= np.timedelta64(1,'D')
            
            # Check the month again
            d0 = d1
            m0 = np.datetime64(d0,'M').astype(int) % 12
            
            # Append, so we get a list of months to be prefixed
            ts.append(d0)
        # Prepand the times found
        if (len(ts)>0):
            gp0 = df_prod["Gp"].values[0]
            ts.sort()
            df_extend = pd.DataFrame({"Profile": ["Prefix"]*len(ts), "Date": ts, "Gp": gp0*len(ts)})
            df_prod = pd.concat([df_extend, df_prod])

    # Set date as index, then resample to monthly
    df_prod.set_index("Date", inplace=True)
    df_prod_m = df_prod.resample('M').interpolate()
    df_prod_m.reset_index(inplace=True)
    
    # Also fill profile column
    df_prod_m['Profile'] = df_prod_m['Profile'].fillna(method='bfill')
    
    return df_prod_m

def apply_time_lag(df, tau, value_key="Pressure", extend_year=2050):
    '''
    Apply 1st order delay to column 'value_key', as a simple way to implement 
    delayed subsidence.
    Value column can be either pressure or Gp (as desired).
    Delay is in years.
    Dataframe structure assumed is
        index Profile Date 'value_key'
    '''
    # Add empty column
    df[value_key+"_D"] = float(0)
    ovals = df[value_key+"_D"].values
    
    # For now do it the ugly way, with a loop
    n = len(df)
    v1 = df[value_key].values[0]
    t1 = df['Date'].values[0]
    w1 = float(v1)
    dt1 = np.timedelta64(1, 'D')*365.25 # Leap year cycle should be 365.2425, but year 2000 was a leap year
    ovals[0] = w1
    for i in range(1,n):
        w0 = w1
        v0 = v1
        t0 = t1
        v1 = df[value_key].values[i]
        t1 = df['Date'].values[i]
        dt = (t1-t0)/dt1
        w1 = w0 + (v1 - w0)*dt/tau
        ovals[i] = w1
 
    return df
    
def convert_gp_to_pressure(df, p_ini, p_aban, IGIP, gas = None):
    '''
    Convert Gp to field average pressure.
    Assume ideal gas if gas=None
    Dataframe structure assumed is
        index Profile Date Fp
    
    Paban is used (and IGIP calculated), unless None is provided, in which case IGIP is used
    and Paban calculated.
    
    returns a modified gp, p_aban, IGIP
    
    Note IGIP is *apparent* *connected* IGIP, not static IGIP
    '''
    
    gp_final = df["Gp"].values[-1]
    
    if (gas is None):
        # Back calculate IGIP, if p_aban provided
        if (p_aban is not None and p_aban>0):
            IGIP = gp_final/(1-p_aban/p_ini)

        df["Pressure"] = p_ini*(1 - df["Gp"]/IGIP)
    else:
        # For now, ugly
        z_ini = gas.zfactor(p_ini)
        pz_ini = p_ini/z_ini

        # Back calculate IGIP, if p_aban provided
        if (p_aban is not None and p_aban>0):
            z_aban = gas.zfactor(p_aban)
            pz_aban = p_aban/z_aban
            IGIP = gp_final/(1-pz_aban/pz_ini)
        
        # Linear relationship between Gp and pz
        df["Pressure"] = pz_ini*(1 - df["Gp"]/IGIP)
        
        # Then convert to pressure
        n = len(df)
        p_vals = df["Pressure"].values
        for i in range(n):
            pz = p_vals[i]
            p = gas.p_from_pz(pz)
            p_vals[i] = p
            
    p_aban=df["Pressure"].values[-1]
    
    return df, p_aban, IGIP

if __name__ == "__main__":
    # Read Profile
    df = pd.read_excel(r".\test_data\test_data.xlsx", sheet_name="Profiles")

    # Get the relevant ones
    c_fc = "Test FC 21 2P"
    c_hist = "Test HIST 21"

    # Then extract the profile
    df_prod = get_merged_profile(df, c_hist, c_fc)

    # Define our gas
    L15_gas = dry_gas(110, 0.6, 0.0, 0.01, 0.01)

    # Convert to pressure
    df_prod = convert_gp_to_pressure(df_prod, 340, 10, gas = L15_gas)

    # Apply time lag
    df_prod = apply_time_lag(df_prod, 5)

    # Show off
    print(df_prod)

    import matplotlib
    import matplotlib.pyplot as plt
    
    ax = df_prod.plot(x="Date", y="Pressure")
    ax = df_prod.plot(x="Date", y="Pressure_D", ax=ax)
    plt.show()
