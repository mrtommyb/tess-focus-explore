import numpy as np



def calctmagv(v, j, k, ve = None, je = None, ke = None):
    ''' calculate tess magnitudes from V and JHK
        Guillermo Torres private email 5/27/2016
    '''

    def sccalctmagv(v, j, k, ve, je, ke):
        if v is None or j is None or k is None:
            return (None, None, None)

        if ve is None:
            ve = 0.0
        if je is None:
            je = 0.0
        if ke is None:
            ke = 0.0

        vk = v - k
        vke = np.sqrt(ve**2 + ke**2)
        t = j + 0.0293 + 0.38523 * vk - 0.01862 * vk**2 + 0.00152 * vk**3
        te = np.sqrt(je**2 + 0.021**2 +
                   ((0.38523 - 2 * 0.01862 * vk + 3 * 0.00152 * vk**2) * vke)**2)
        ts = 'tor_vjk'
        return (round(t, 3), round(te, 3), ts)
    
    
    def npcalctmagv(v, j, k, ve, je, ke):
        if ve is None:
            ve = np.zeros(len(j))
        if je is None:
            je = np.zeros(len(j))
        if ke is None:
            ke = np.zeros(len(k))

        bd = np.isnan(ve)
        ve[bd] = 0.0
        bd = np.isnan(je)
        je[bd] = 0.0
        bd = np.isnan(ke)
        ke[bd] = 0.0
        
        t  = np.zeros(len(j))
        te = np.zeros(len(j))
        t.fill(np.nan)
        te.fill(np.nan)
        ts = len(j) * [None]
        ts = np.array(ts)
        
        vk = np.zeros(len(j))
        vke = np.zeros(len(j))

        idx = np.where(~(np.isnan(v + j + k)))
        nidx = np.shape(idx)[1]
        if nidx <= 0:
            return (t, te, ts)
        
        vk[idx] = v[idx]  - k[idx]
        vke[idx] = np.sqrt(ve[idx]**2 + ke[idx]**2)
        t[idx]  = j[idx] + 0.0293 + 0.38523 * vk[idx] - 0.01862 * vk[idx]**2 + \
                   0.00152 * vk[idx]**3
        te[idx] = np.sqrt(je[idx]**2 + 0.021**2 +
                           ((0.38523 - 2 * 0.01862 * vk[idx] + 
                             3 * 0.00152 * vk[idx]**2) * vke[idx])**2)
        ts[idx] = 'tor_vjk'
        return (np.round(t, 3), np.round(te, 3), ts)

    
    if type(j) == np.ndarray or type(k) == np.ndarray:
        return npcalctmagv(v, j, k, ve, je, ke)
    return sccalctmagv(v, j, k, ve, je, ke)



def calctmagbphot(b, j, k, be = None, je = None, ke = None):
    ''' calculate tess magnitudes from photographic B and JHK
        Guillermo Torres private email 5/27/2016
    '''
    
    def sccalctmagbphot(b, j, k, be, je, ke):
        if b is None or j is None or k is None:
            return (None, None, None)
        
        if be is None:
            be = 0.0
        if je is None:
            je = 0.0
        if ke is None:
            ke = 0.0

        bk = b - k
        bke = sqrt(be**2 + ke**2)
        t  = j + 0.0381 + 0.31926 * bk - 0.01780 * bk**2 + 0.00178 * bk**3
        te = sqrt(je**2 + 0.020**2 +
                   ((0.31926 - 2 * 0.01780 * bk + 3 * 0.00178 * bk**2) * bke)**2)
        ts = 'tor_bpjk'
        return (round(t, 3), round(te, 3), ts)
    
    
    def npcalctmagbphot(b, j, k, be, je, ke):
        if be is None:
            be = np.zeros(len(j))
        if je is None:
            je = np.zeros(len(j))
        if ke is None:
            ke = np.zeros(len(k))

        bd = np.isnan(be)
        be[bd] = 0.0
        bd = np.isnan(je)
        je[bd] = 0.0
        bd = np.isnan(ke)
        ke[bd] = 0.0
        
        t  = np.zeros(len(j))
        te = np.zeros(len(j))
        t.fill(np.nan)
        te.fill(np.nan)
        ts = len(j) * [None]
        ts = np.array(ts)
        
        bk = np.zeros(len(j))
        bke = np.zeros(len(j))

        idx = np.where(~(np.isnan(b + j + k)))
        bk[idx] = b[idx]  - k[idx]
        bke[idx] = np.sqrt(be[idx]**2 + ke[idx]**2)
        t[idx]  = j[idx] + 0.0381 + 0.31926 * bk[idx] - 0.01780 * bk[idx]**2 + \
                   0.00178 * bk[idx]**3
        te[idx] = np.sqrt(je[idx]**2 + 0.020**2 +
                           ((0.31926 - 2 * 0.01780 * bk[idx] + 
                             3 * 0.00178 * bk[idx]**2) * bke[idx])**2)
        ts = 'tor_bpjk'
        return (np.round(t, 3), np.round(te, 3), ts)

    
    if type(j) == np.ndarray or type(k) == np.ndarray:
        return npcalctmagbphot(b, j, k, be, je, ke)
    return sccalctmagbphot(b, j, k, be, je, ke)



def calctmagb(b, j, k, be = None, je = None, ke = None):
    ''' calculate tess magnitudes from Johnson B and JHK
        Guillermo Torres private email 5/27/2016
    '''
    
    def sccalctmagb(b, j, k, be, je, ke):
        if b is None or j is None or k is None:
            return (None, None, None)
        
        if be is None:
            be = 0.0
        if je is None:
            je = 0.0
        if ke is None:
            ke = 0.0

        bk = b - k
        bke = sqrt(be**2 + ke**2)
        t  = j + 0.0407 +  0.29688 * bk - 0.02313 * bk**2 + 0.00226 * bk**3
        te = sqrt(je**2 + 0.031**2 +
                   ((0.29688 - 2 * 0.02313 * bk + 3 * 0.00226 * bk**2) * bke)**2)
        ts = 'tor_bjk'
        return (round(t, 3), round(te, 3), ts)
    
    
    def npcalctmagb(b, j, k, be, je, ke):
        if be is None:
            be = np.zeros(len(b))
        if je is None:
            je = np.zeros(len(j))
        if ke is None:
            ke = np.zeros(len(k))

        bd = np.isnan(be)
        be[bd] = 0.0
        bd = np.isnan(je)
        je[bd] = 0.0
        bd = np.isnan(ke)
        ke[bd] = 0.0
        
        t  = np.zeros(len(j))
        te = np.zeros(len(j))
        t.fill(np.nan)
        te.fill(np.nan)
        ts = len(j) * [None]
        ts = np.array(ts)
        
        bk = np.zeros(len(j))
        bke = np.zeros(len(j))

        idx = np.where(~(np.isnan(b + j + k)))
        bk[idx] = b[idx]  - k[idx]
        bke[idx] = np.sqrt(be[idx]**2 + ke[idx]**2)
        t[idx]  = j[idx] + 0.0407 + 0.29688 * bk[idx] - 0.02313 * bk[idx]**2 + \
                   0.00226 * bk[idx]**3
        te[idx] = sqrt(je[idx]**2 + 0.031**2 +
                        ((0.29688 - 2 * 0.02313 * bk[idx] + 
                          3 * 0.00226 * bk[idx]**2) * bke[idx])**2)
        ts = 'tor_bjk'
        return (np.round(t, 3), np.round(te, 3), ts)

    
    if type(b) == np.ndarray or type(j) == np.ndarray or type(k) == np.ndarray:
        return npcalctmagb(b, j, k, be, je, ke)
    return sccalctmagb(b, j, k, be, je, ke)



def calctmagjhk(j, k, je = None, ke = None):
    ''' calculate tess magnitudes from JHK
        Guillermo Torres private email 5/27/2016
    '''
    
    def sccalctmagjhk(j, k, je, ke):
        if j is None or k is None:
            return (None, None, None)
        
        if je is None:
            je = 0.0
        if ke is None:
            ke = 0.0
        t = None
        te = None
        ts = None
        
        jk = j - k
        jke = sqrt(je**2 + ke**2)
        
        if jk <= 0.70:
            t = j + 0.0563 + 1.89115 * jk - 1.74299 * jk**2 + 1.22163 * jk**3
            te = sqrt(je**2 + 0.008**2 +
                      ((1.89115 - 2 * 1.74299 * jk + 3 * 1.22163 * jk**2) * jke)**2)
        elif jk > 0.7 and jk <= 1.0:
            t = j + 147.811 - 545.64 * jk + 668.453 * jk**2 - 269.372 * jk**3
            te = sqrt(je**2 + 0.17**2 +
                      ((-545.64 + 2 * 668.453 * jk - 3 * 269.372 * jk**2) * jke)**2)
        elif jk > 1.0:
            t = j + 1.75
            te = 1.0
        ts = 'tor_jhk'
        
        return (round(t, 3), round(te, 3), ts)


    def npcalctmagjhk(j, k, je, ke):
        if je is None:
            je = np.zeros(len(j))
        if ke is None:
            ke = np.zeros(len(k))

        bd = np.isnan(je)
        je[bd] = 0.0
        bd = np.isnan(ke)
        ke[bd] = 0.0
        
        t  = np.zeros(len(j))
        te = np.zeros(len(j))
        t.fill(np.nan)
        te.fill(np.nan)
        ts = len(j) * [None]
        ts = np.array(ts)
        
        jk = np.zeros(len(j))
        jke = np.zeros(len(j))

        jk = j - k
        mask = ~np.isnan(jk)
        gd  = np.where(mask)
        mask[gd] = (jk[gd] <= 0.7)
        gd  = np.where(mask)
        ngd = np.shape(gd)[1]
        if ngd > 0:
            jk[gd] = j[gd] - k[gd]
            jke[gd] = np.sqrt(je[gd]**2 + ke[gd]**2)
            t[gd]  = j[gd] + 0.0563 + 1.89115 * jk[gd] - 1.74299 * jk[gd]**2 + \
                      1.22163 * jk[gd]**3
            te[gd] = np.sqrt(je[gd]**2 + 0.008**2 +
                              ((1.89115 - 2 * 1.74299 * jk[gd] + 
                                3 * 1.22163 * jk[gd]**2) * jke[gd])**2)
            ts[gd] = 'tor_jhk'
        
        mask = ~np.isnan(jk)
        gd  = np.where(mask)
        mask[gd] = (jk[gd] > 0.7 & jk <= 1.0)
        gd  = np.where(mask)
        ngd = np.shape(gd)[1]
        if ngd > 0:
            jk[gd] = j[gd] - k[gd]
            jke[gd] = np.sqrt(je[gd]**2 + ke[gd]**2)
            t[gd]  = j[gd] + 147.811 - 545.64 * jk[gd] + 668.453 * jk[gd]**2 - \
                      269.372 * jk[gd]**3
            te[gd] = np.sqrt(je[gd]**2 + 0.17**2 +
                              ((-545.64 + 2 * 668.453 * jk[gd] - \
                                3 * 269.372 * jk[gd]**2) * jke[gd])**2)
            ts[gd] = 'tor_jhk'

        mask = ~np.isnan(jk)
        gd  = np.where(mask)
        mask[gd] = (jk[gd] > 1.0)
        gd  = np.where(mask)
        ngd = np.shape(gd)[1]
        if ngd > 0:
            t[gd] = j[gd] + 1.75
            te[gd] = 1.0        
            ts[gd] = 'tor_jhk'

        return (np.round(t, 3), np.round(te, 3), ts)
        

    if type(j) == np.ndarray or type(k) == np.ndarray:
        return npcalctmagjhk(j, k, je, ke)
    return sccalctmagjhk(j, k, je, ke)



def calctmagvjh(v, j, h, ve = None, je = None, he = None):
    ''' calculate tess magnitudes from V and JHK
        Guillermo Torres private email 5/27/2016
    '''
    
    def sccalctmagvjh(v, j, h, ve, je, he):
        if v is None or j is None or h is None:
            return (None, None, None)
        
        if ve is None or np.isnan(ve):
            ve = 0.0
        if je is None or np.isnan(je):
            je = 0.0
        if he is None or np.isnan(he):
            he = 0.0

        jh = j - h
        jhe = sqrt(je**2 + he**2)
        t  = v - 0.1140 - 1.96827 * jh + 0.75955 * jh**2 - 0.28408 * jh**3
        te = sqrt(ve**2 + 0.063**2 + 
                   ((-1.96827 + 2 * 0.75955 * jh - 3 * 0.28408 * jh**2) * jhe)**2)
        ts = 'tor_vjh'
        
        if (v - t) > 1.5:
            t = v - 1.5
            te = 1.5
            ts = 'voffset'
        
        return (round(t, 3), round(te, 3), ts)
    
    
    def npcalctmagvjh(v, j, h, ve, je, he):
        if ve is None:
            ve = np.zeros(len(v))
        if je is None:
            je = np.zeros(len(j))
        if he is None:
            he = np.zeros(len(h))

        bd = np.isnan(ve)
        ve[bd] = 0.0
        bd = np.isnan(je)
        je[bd] = 0.0
        bd = np.isnan(he)
        he[bd] = 0.0
        
        t  = np.zeros(len(j))
        te = np.zeros(len(j))
        t.fill(np.nan)
        te.fill(np.nan)
        ts = len(j) * [None]
        ts = np.array(ts)
        
        jh = np.zeros(len(j))
        jhe = np.zeros(len(j))

        idx = np.where(~(np.isnan(v + j + h)))
        nidx = np.shape(idx)[1]
        if nidx <= 0:
            return (t, te, ts)
        
        jh[idx] = j[idx]  - h[idx]
        jhe[idx] = np.sqrt(ve[idx]**2 + he[idx]**2)
        t[idx]  = v[idx] - 0.1140 + - 1.96827 * jh[idx] + 0.75955 * jh[idx]**2 - \
                   0.28408 * jh[idx]**3
        te[idx] = np.sqrt(ve[idx]**2 + 0.063**2 +
                           ((-1.96827 + 2 * 0.75955 * jh[idx] - 
                             3 * 0.28408 * jh[idx]**2) * jhe[idx])**2)
        ts[idx] = 'tor_vjh'
        
        bd = np.where(((~np.isnan(v)) & (np.isnan(t))) | 
                      ((v - t) > 1.5))
        nbd = np.shape(bd)[1]
        if nbd > 0:
            t[bd] = v[bd] - 1.5
            te[bd] = 1.5
            ts[bd] = 'voffset'
        
        return (np.round(t, 3), np.round(te, 3), ts)

    
    if type(j) == np.ndarray or type(h) == np.ndarray:
        return npcalctmagvjh(v, j, h, ve, je, he)
    return sccalctmagvjh(v, j, h, ve, je, he)



def calctmagjh(j, h, je = None, he = None):
    ''' calculate tess magnitudes from V and JHK
        Guillermo Torres private email 5/27/2016
    '''
    
    def sccalctmagjh(j, h, je, he):
        if j is None or np.isnan(j) or h is None or np.isnan(h):
            return (None, None, None)
        
        if je is None or np.isnan(je):
            je = 0.0
        if he is None or np.isnan(he):
            he = 0.0

        jh = j - h
        jhe = sqrt(je**2 + he**2)
        t  = j + 0.1561 + 1.93384 * jh - 1.49220 * jh**2 - 0.99995 * jh**3
        te = sqrt(je**2 + 0.040**2 +
                   ((1.93384 - 2 * 1.49220 * jh - 3 * 0.99995 * jh**2) * jhe)**2)
        ts = 'tor_jh'
        
        if (t - j) > 0.5 + 0.8:    # offset + sigma for t - j
            t = j + 0.4
            te = 0.5
            ts = 'joffset'
        
        return (round(t, 3), round(te, 3), ts)
    
    
    def npcalctmagjh(j, h, je, he):
        if je is None:
            je = np.zeros(len(j))
        if he is None:
            he = np.zeros(len(h))

        bd = np.isnan(je)
        je[bd] = 0.0
        bd = np.isnan(he)
        he[bd] = 0.0
        
        t  = np.zeros(len(j))
        te = np.zeros(len(j))
        t.fill(np.nan)
        te.fill(np.nan)
        ts = len(j) * [None]
        ts = np.array(ts)
        
        jh = np.zeros(len(j))
        jhe = np.zeros(len(j))

        idx = np.where(~(np.isnan(j + h)))
        nidx = np.shape(idx)[1]
        if nidx <= 0:
            return (t, te, ts)
        
        jh[idx] = j[idx]  - h[idx]
        jhe[idx] = np.sqrt(je[idx]**2 + he[idx]**2)
        t[idx]  = j[idx] + 0.1561 + 1.93384 * jh[idx] - 1.49220 * jh[idx]**2 - \
                   0.99995 * jh[idx]**3
        te[idx] = np.sqrt(je[idx]**2 + 0.040**2 +
                           ((1.93384 - 2 * 1.49220 * jh[idx] - \
                             3 * 0.99995 * jh[idx]**2) * jhe[idx])**2)
        ts[idx] = 'tor_jh'
        
#         bd = np.where(((~np.isnan(j)) & (np.isnan(h))) | 
#                       ((t - j) > 0.4 + 0.5))
#         nbd = np.shape(bd)[1]
#         if nbd > 0:
#             t[bd] = j[bd] + 0.4
#             te[bd] = 0.5
#             ts[bd] = 'joffset'
        
        return (np.round(t, 3), np.round(te, 3), ts)

    
    if type(j) == np.ndarray or type(h) == np.ndarray:
        return npcalctmagjh(j, h, je, he)
    return sccalctmagjh(j, h, je, he)


def tessmag_err(tessmag):
    '''Tessmag precision in ppm
        total photon error from the source using the Ricker 2015 paper, Buoma 2016 private comm.
    '''

    if (tessmag is None):
        return(0.0)
    else:
        # 
        F = 4.73508403525e-05
        E = -0.0022308015894
        D = 0.0395908321369
        C = -0.285041632435
        B = 0.850021465753
        A = exp(3.29685004771)
        pht_er_src = A * np.exp(B * tessmag + C * tessmag ** 2 + 
                                D * tessmag ** 3 + E * tessmag ** 4 + 
                                F * tessmag ** 5)
        
        if (tessmag < 5):
            pht_er_src = 61.75
        return(pht_er_src)
    
def contam(x0, y0, xb, yb, s, dx0 = None, dy0 = None):
    ''' computes the flux of a star sitting at x0, y0 falling into a box
        sized xb, yb (radius) around the origin (0, 0)
        The target star sits at (0, 0), the contaminant at (x0, y0)
        (NOTE: we get a - sign from integrating the Gaussian which is why
               we have +- x0 instead of (+- xb) in the argument of errf)
        We compute how much of the normalized flux in a box xb, yb around (x0, y0) 
        falls into the box xb, yb around (0, 0)
        x0, y0 - point around which to compute the contamination (origin of flux)
        xb, yb - width of the box around (0, 0) in x and y direction
        s      - sigma of gaussian over which we integrate (int(gauss) = erf)
        dx0, dy0 - error for x0, y0
    '''
    # use internal functions for single values and numpy arrays, providing
    # the same call signature no matter if called with scalars or arrays
    
    def sccontam(x0, y0, xb, yb, s, dx0 = None, dy0 = None):
        sq2 = sqrt(2)
        contx = erf((xb + x0) / (sq2 * s)) + erf((xb - x0) / (sq2 * s))
        conty = erf((yb + y0) / (sq2 * s)) + erf((yb - y0) / (sq2 * s))
        cont = 0.25 * contx * conty
        
        dcont = None
        if dx0 is not None or dy0 is not None:
            dcont = 0.0
            tmp = 0.25 * sqrt(2.0 / pi) / s
            dnm = 2.0 * s * s
            if dx0 is not None:
                dcont += (tmp * (exp(-(xb + x0)**2 / dnm) - exp(-(xb - x0)**2 / dnm)) * 
                          conty * dx0)**2
            if dy0 is not None:
                dcont += (tmp * (exp(-(yb + y0)**2 / dnm) - exp(-(yb - y0)**2 / dnm)) * 
                          contx * dy0)**2
            dcont = sqrt(dcont)
    
        return (cont, dcont)

    def npcontam(x0, y0, xb, yb, s, dx0 = None, dy0 = None):
        sq2 = sqrt(2)
        contx = npspec.erf((xb + x0) / (sq2 * s)) + npspec.erf((xb - x0) / (sq2 * s))
        conty = npspec.erf((yb + y0) / (sq2 * s)) + npspec.erf((yb - y0) / (sq2 * s))
        cont = 0.25 * contx * conty
    
        dcont = None
        if dx0 is not None or dy0 is not None:
            dcont = 0.0
            tmp = 0.25 * sqrt(2.0 / pi) / s
            dnm = 2.0 * s * s
            if dx0 is not None:
                dcont += (tmp * (np.exp(-(xb + x0)**2 / dnm) - np.exp(-(xb - x0)**2 / dnm)) * 
                          conty * dx0)**2
            if dy0 is not None:
                dcont += (tmp * (np.exp(-(yb + y0)**2 / dnm) - np.exp(-(yb - y0)**2 / dnm)) * 
                          contx * dy0)**2
            dcont = sqrt(dcont)
        
        return (cont, dcont)

    if type(x0) == np.ndarray:
        return npcontam(x0, y0, xb, yb, s, dx0, dy0)
    return sccontam(x0, y0, xb, yb, s, dx0 = None, dy0 = None)

def tmag2flux(tmag, tmagerr = None):
    ''' we do not have a 0 magnitude flux for tessmags, so we assume it to be
        the same as for Cousins I mags. It does not really play a role, 
        because we use flux ratios in the end 
    '''
    return iflux(tmag, tmagerr)

def iflux(imag, imagerr = None):
    ''' compute flux and error for given I magnitude and error
    '''
    # use internal functions for single values and numpy arrays, providing
    # the same call signature no matter if called with scalars or arrays
    def sciflux(imag, imagerr = None):
        if (imag is None):
            return (None, None)
        flux = 2635.0 * pow(10.0, -imag / 2.5)
        ferr = None
        if (imagerr is not None):
            ferr = abs(-log(10) / 2.5 * flux * imagerr)
        return (flux, ferr) 

    def npiflux(imag, imagerr = None):
        if (imag is None):
            return (None, None)
        nmag = len(imag)
        flux = 2635.0 * np.power(10.0, -imag / 2.5)
        ferr = np.zeros(nmag)
        if (imagerr is not None):
            ferr = abs(-log(10) / 2.5 * np.multiply(flux, imagerr))
        return (flux, ferr) 

    if type(imag) == np.ndarray:
        return npiflux(imag, imagerr)
    return sciflux(imag, imagerr)

def aperture(tmag, tmage = None):
    # get tmag dependant quadratic aperture, 
    # polynomial from email by Jon Jenkins 2017/02/28
    
    def scaperture(tmag, tmage):
        if tmage is None:
            tmage = 0.25
        
        npix = None
        npixe = None

        if tmag <= 4.0:
            npix = 274.2898 - 77.7918 * 4.0 + 7.7410 * 4.0**2 - 0.2592 * 4.0**3
            npixe = abs((77.7918 + 2 * 7.7410 * 4.0 - 3 * 0.2592 * 4.0**2) * tmage)
        else:
            npix = (274.2898 - 77.7918 * tmag + 
                     7.7410 * tmag**2 - 0.2592 * tmag**3)
            npixe = abs(77.7918 + 2 * 7.7410 * tmag - 3 * 0.2592 * tmag**2) * tmage
            if npix < 1:
                npix = 1
                npixe = 0.05
        
        # square of the same area - this will underestimate the contaminating flux
        aper1 = sqrt(npix)
        apere1 = abs(0.5 / sqrt(npix) * npixe)
        
        # surrounding square - area by factor 4 / pi bigger than circle
        # this will overestimate the contaminating flux
        aper2 = 2 * sqrt(npix / pi)
        apere2 = 2 * sqrt(1.0 / (npix * pi)) * npixe
        
        # average over surrounding square and square of the same area
        # this should give a more realistic estimate
        aper  = (aper1 + aper2) / 2
        apere = (apere1 + apere2) / 2
        
        return (aper, apere)
    
    
    def npaperture(tmag, tmage):
        npix = np.zeros(len(tmag))
        npixe = np.copy(npix)
        
        if tmage is None:
            tmage = np.copy(npix)
        else:
            bd = np.where(np.isnan(tmage))
            tmage[bd] = 0.25
            
        bd  = np.where(tmag <= 4.0)
        nbd = np.shape(bd)[1]
        if nbd > 0:
            npix[bd] = 274.2898 - 77.7918 * 4.0 + 7.7410 * 4.0**2 - 0.2592 * 4.0**3
            npixe[bd] = np.abs((77.7918 + 2 * 7.7410 * 4.0 - 3 * 0.2592 * 4.0**2) *
                               tmage[bd])
            
        gd  = np.where(tmag > 4.0)
        ngd = np.shape(gd)[1]
        if ngd > 0:
            npix[gd] = (274.2898 - 77.7918 * tmag[gd] + 
                         7.7410 * tmag[gd]**2 - 0.2592 * tmag[gd]**3)
            npixe[gd] = np.abs(77.7918 + 2 * 7.7410 * tmag[gd] - 3 * 0.2592 * tmag[gd]**2) \
                        * tmage
        bd = np.where(npix < 1.0)
        npix[bd] = 1.0
        bd = np.where(npixe < 0.25)
        npixe[bd] = 0.25
        
        aper = np.sqrt(npix)
        apere = np.abs(0.5 / aper * npixe)
        return (aper, apere)

    if type(tmag) == np.ndarray:
        return npaperture(tmag, tmage)
    return scaperture(tmag, tmage)

def haversine(ra1, dec1, ra2, dec2, 
              dra1 = None, ddec1 = None, dra2 = None, ddec2 = None):
    ''' compute distance between ra1, dec1 and ra2, dec2 using haversine 
        approximation, dra1, ddec2, dra2, ddec2 errors for point 1 and 2
        return - distance, error for distance
    '''
    # use internal functions for single values and numpy arrays, providing
    # the same call signature no matter if called with scalars or arrays
    def schaversine(ra1, dec1, ra2, dec2, 
                    dra1 = None, ddec1 = None, dra2 = None, ddec2 = None):
        ras = ra1  * deg2rad
        raf = ra2  * deg2rad
        decs = dec1 * deg2rad
        decf = dec2 * deg2rad
    
        dra = raf - ras
        ddec = decf - decs
        delta = 2 * asin(sqrt(sin(ddec / 2)**2 + 
                              cos(decs) * cos(decf) * sin(dra / 2)**2))
        delta = delta * rad2deg
    
        ddelta = None
        if (dra1 is not None or ddec1 is not None or 
            dra2 is not None or ddec2 is not None):
            tmp  = sin(0.5 * (dec2 - dec1))**2 * cos(ra1) * cos(ra2) + \
                    sin(0.5 * (ra2 - ra1))**2
            denom = 4 * (1 - tmp) * tmp
            ddelta = 0.0
            if dra1 is not None:
                ddelta += (cos(dec1 - dec2) * cos(ra2) * sin(ra1) - 
                           cos(ra1) * sin(ra2))**2 / denom * (dra1 * deg2rad)**2
            if dra2 is not None:
                ddelta += (sin(ra1 - ra2) + 
                           2 * cos(ra1) * sin(ra2) * sin(0.5 * (dec2 - dec1))**2)**2 / \
                             denom * (dra2 * deg2rad)**2
            if ddec1 is not None:
                ddelta += (cos(ra1) * cos(ra2) * sin(dec1 - dec2))**2 / \
                          denom * (ddec1 * deg2rad)**2
            if ddec2 is not None:
                ddelta += (cos(ra1) * cos(ra2) * sin(dec1 - dec2))**2 / \
                          denom * (ddec2 * deg2rad)**2
            ddelta = sqrt(ddelta) * rad2deg
            
        return (delta, ddelta)


    def nphaversine(ra1, dec1, ra2, dec2,
                    dra1 = None, ddec1 = None, dra2 = None, ddec2 = None):
        ras = ra1 * deg2rad
        raf = ra2 * deg2rad
        decs = dec1 * deg2rad
        decf = dec2 * deg2rad
    
        dra = raf - ras
        ddec = decf - decs
        tmp1 = np.sin(ddec / 2)**2
        tmp2 = np.cos(decs) * np.cos(decf) * np.sin(dra / 2)**2
        delta = 2.0 * np.arcsin(np.sqrt(tmp1 + tmp2))
        delta = delta * rad2deg
        
        ddelta = None
        if (dra1 is not None or ddec1 is not None or 
            dra2 is not None or ddec2 is not None):
            tmp  = sin(0.5 * (dec2 - dec1))**2 * cos(ra1) * cos(ra2) + \
                    sin(0.5 * (ra2 - ra1))**2
            denom = 4 * (1 - tmp) * tmp
            ddelta = np.array(len(ra1) * [0.0])
            if dra1 is not None:
                ddelta += (np.cos(dec1 - dec2) * np.cos(ra2) * np.sin(ra1) - 
                           np.cos(ra1) * np.sin(ra2))**2 / denom * (dra1 * deg2rad)**2
            if dra2 is not None:
                ddelta += (np.sin(ra1 - ra2) + 
                           2 * np.cos(ra1) * np.sin(ra2) * 
                           np.sin(0.5 * (dec2 - dec1))**2)**2 / \
                             denom * (dra2 * deg2rad)**2
            if ddec1 is not None:
                ddelta += (np.cos(ra1) * np.cos(ra2) * np.sin(dec1 - dec2))**2 / \
                          denom * (ddec1 * deg2rad)**2
            if ddec2 is not None:
                ddelta += (np.cos(ra1) * np.cos(ra2) * np.sin(dec1 - dec2))**2 / \
                          denom * (ddec2 * deg2rad)**2
            ddelta = np.sqrt(ddelta) * rad2deg
        
        return (delta, ddelta)

    if type(ra1) == np.ndarray or type(dec1) == np.ndarray or \
       type(ra2) == np.ndarray or type(dec2) == np.ndarray:
        return nphaversine(ra1, dec1, ra2, dec2, dra1, ddec1, dra2, ddec2)
    return schaversine(ra1, dec1, ra2, dec2, dra1, ddec1, dra2, ddec2)

