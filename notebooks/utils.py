from scipy.constants import pi
import numpy as np

def get_gs(df, Tl = 'Tliq', jez=False):
    df['K_W_Tc'] = (df.Tc - df.Tg) / df[Tl] # best one in the paper
    df['K_W_Tx'] = (df.Tx - df.Tg) / df[Tl]
    df['gamma_Tc'] = df.Tc / (df.Tg+df[Tl])
    df['H_prime_Tx'] = (df.Tx - df.Tg) / df.Tg
    df['K_H_Tc'] = (df.Tc - df.Tg) / (df[Tl] - df.Tc) # replaced Tmelt with Tliq
    df['H_prime_Tc'] = (df.Tc - df.Tg) / df.Tg
    df['K_H_Tx'] = (df.Tx - df.Tg) / (df[Tl] - df.Tx) # replaced Tmelt with Tliq
    df['deltaT_rg'] = (df.Tx - df.Tg) / (df[Tl] - df.Tg)
    df['K_cr'] = (df[Tl] - df.Tx) / (df[Tl] - df.Tg)
    if jez:
        df['Jezica'] = (df.ViscosityAtTl) - 2 * np.log10(df[Tl])
    return df
    
def get_eta_tl(df, Tl = 'Tliq'):
    return df['log10 (η∞)'] + (12-df['log10 (η∞)'])*(df.T12/df[Tl])*np.exp((df.m/(12-df['log10 (η∞)'])-1)*(df.T12/df[Tl] - 1))

def get_gfa(df, logXs = -2, logNs = 3, g=pi,  Tl = 'Tliq', **kw):
    Umax = 10 ** df.log_Umax
    
    tn = (10**logXs / (g * 10**logNs * Umax**2))**(1 / 2)

    df['GFA'] = -np.log10((df[Tl] - df.T_Umax) / tn)
    return df
