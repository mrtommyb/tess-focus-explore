'''


Description

Compute flux contamination for every ctl star using
    Pixel scale:  20.25 arcsec/pixel
    PSF FWHM: 2 pixels
    Aperture size: 4 pixel diameter, mag varying aperture
    Search radius for contaminants: 10 pixel radius 

'''

# original values from email by Josh Pepper from 06/03/2016
tess_scale = 20.25 / 3600.0        # arcsec / pixel --> deg / pixel
# tess_fwhm  = 2.0  * tess_scale     # 2 pixel
tess_fwhm  = 1.88  * tess_scale    # new fwhm from camera image (sigma = 0.80)
tess_aper  = 4.0  * tess_scale     # 4 pixel aperture
tess_srad  = 10.0 * tess_scale     # 10 pixel search radius 
tess_sigma = tess_fwhm / (2.0 * sqrt(2.0 * log(2)))



if __name__ == '__main__':
    watch = Stopwatch()
    watch.start()
    
    parser = get_stdopts()
    
    # add local options
    parser.add_option('--tmin', dest='tmin', type='float', default='1',
                      help='tessmag min (default: 1')
    parser.add_option('--tmax', dest='tmax', type='float', default='-1',
                      help='tessmag min (default: -1')
    
    options = get_config(parser)

    options.lf = Logfile(options.lfname)
    lf = options.lf


        
    nblocks = int((options.stop - options.start) / options.blen) + 1
    lf.write(str(nblocks - 1) + ' blocks to process')

    rmydb  = DBReflect(options.catdburl)
    rconn  = rmydb.getConnection()
    sconn  = rmydb.getNewConnection()
    
    wmydb  = DBReflect(options.catdburl)
    wconn  = wmydb.getConnection()

    # no tesmag limits
    selcmd  = sql.text('select ticlinkpk, ra, dec, tessmag, tessmag_e '
                       'from ticlink_ctl '
                       'where ticlinkpk >= :start and ticlinkpk < :stop')

    starcmd = sql.text('select pk, ra, dec, tessmag, tessmag_e, '
                              '3600 * q3c_dist(ra, dec, :ctlra, :ctldec) as rdist '
                       'from ticlink '
                       'where q3c_radial_query(ra, dec, :ctlra, :ctldec, :srad) and '
                             'pk != :ticlinkpk '
                       'order by rdist asc')
    
    galcmd  = sql.text('select pk, ra, dec, tmag as tessmag, e_tmag as tessmag_e, '
                              '3600 * q3c_dist(ra, dec, :ctlra, :ctldec) as rdist '
                       'from ticnasa_ext '
                       'where q3c_radial_query(ra, dec, :ctlra, :ctldec, :srad) '
                       'order by rdist')
    
    updcmd  = sql.text('update ticvals set numcont = :nrcont, '
                                          'contratio = :flxratio, '
                                          'contratio_e = :flxratioe '
                       'where ticlinkpk = :ticlinkpk ')
    
    nrstars = 0
    blockno = 0
    nrzeros = 1
    nrcont  = 0
    contcnt = 0
    updates = []
#     parms   = {'tmin' : options.tmin, 'tmax' : options.tmax}
#     for row in rmydb.windowed(selcmd, options, single = True):
    for row in rmydb.fetchall(selcmd, single = True):
        nrstars += 1
        tmage = row['tessmag_e']
        
        if tmage is None or tmage == np.nan:
            tmage = 0.25
            
        (aper, apere) = aperture(row['tessmag'], tmage)
        aper  = aper * tess_scale
        apere = apere * tess_scale
        if options.debug > 1:
            print 'target star'
            print 'keys  : ', row.keys()
            print 'values: ', row
            print 'apert.: ', aper / tess_scale, ' [px]'
            
        try:
            (tflux, tfluxe) = tmag2flux(row['tessmag'], tmage)
            if options.debug > 1:
                print 'tflux, tfluxe = ', tflux, tfluxe, ' (target flux and error)'
        except:
            e, v = sys.exc_info()[:2]
            lf.write()
            lf.write('while computing tflux for ticlinkpk ' + str(row['ticlinkpk']))
            lf.write('ERROR: ' + str(e))
            lf.write('VALUE: ' + str(v))
            raise

        if np.isnan(tflux):
            lf.write('nan flux for pk ' + str(row['ticlinkpk']))
            continue
        
        # contamination from stellar sources
        pars   = {'ctlra' : row['ra'], 'ctldec' : row['dec'], 'srad' : tess_srad,
                  'ticlinkpk' : row['ticlinkpk']}
        crows  = sconn.execute(starcmd, pars).fetchall()
        nrcont = len(crows)
        if nrcont == 0:
            lf.write('INVESTIGATE: no contaminants for ticlinkpk ' + 
                     str(row['ticlinkpk']))
            continue 
        if options.debug > 1:
            print ''
            print nrcont, ' stellar contaminants in order of distance [asec]'
            print 'keys :', crows[0].keys()
            for ii in range(nrcont):
                print ii, ' : ', crows[ii]
        cdict = rmydb.getDictOfCols(crows)
        bd    = np.where(np.isnan(cdict['tessmag_e']))
        nbd   = np.shape(bd)[1]
        if nbd > 0:
            cdict['tessmag_e'][bd] = 0.25 
        (cflux, cfluxe) = tmag2flux(cdict['tessmag'], cdict['tessmag_e'])
        if options.debug > 1:
            print ''
            print 'stellar contaminant fluxes and errors'
            print cflux
            print cfluxe
        
        # compile information for contamfluxes
#         (dist, diste) = haversine(row.ra, row.dec, cdict['ra'], cdict['dec'])
        (dra, drae)   = haversine(row.ra, row.dec, cdict['ra'], row.dec)
        (ddec, ddece) = haversine(row.ra, row.dec, row.ra, cdict['dec'])
        (tmp, dtmp)   = contam(dra, ddec, aper / 2.0, aper / 2.0, tess_sigma)
        if options.debug > 1:
            print ''
            print 'geometric fractions of contaminant fluxes in target aperture'
            print tmp
        cflx = tmp * cflux
        cfle = abs(tmp * cfluxe)
        if options.debug > 1:
            print ''
            print 'contaminating fluxes and errors (fraction * contaminant fluxes)'
            print cflx
            print cfle

        # sum it up
        sumcflx = np.sum(cflx)
        sumcfle = np.sum(cfle)
        if np.isnan(sumcflx):
            msg = 'sumcflx is nan for ticlinkpk %d' % (row['ticlinkpk'],) 
            lf.write(msg)
        if options.debug > 1:
            print ''
            print 'sum of stellar contaminating flux and errors'
            print sumcflx, sumcfle

        # contamination from galactic sources
        del pars['ticlinkpk']
        crows  = sconn.execute(galcmd, pars).fetchall()
        nrgals = len(crows)
        if nrgals != 0:
            nrcont += nrgals
            cdict  = rmydb.getDictOfCols(crows)
            (cflux, cfluxe) = tmag2flux(cdict['tessmag'], cdict['tessmag_e'])
            bd    = np.where(np.isnan(cdict['tessmag_e']))
            nbd   = np.shape(bd)[1]
            if nbd > 0:
                cdict['tessmag_e'][bd] = 0.25 
                
            if options.debug > 1:
                print ''
                print '-------------------------------------------------------'
                print ''
                print nrgals, ' galactic contaminants in order of distance [asec]'
                print 'keys :', crows[0].keys()
                for ii in range(nrgals):
                    print ii, ' : ', crows[ii]
                print ''
                print 'galactic contaminant fluxes and errors'
                print cflux
                print cfluxe
            
            # compile information for contamfluxes
#             (dist, diste) = haversine(row.ra, row.dec, cdict['ra'], cdict['dec'])
            (dra, drae)   = haversine(row.ra, row.dec, cdict['ra'], row.dec)
            (ddec, ddece) = haversine(row.ra, row.dec, row.ra, cdict['dec'])
            (tmp, dtmp)   = contam(dra, ddec, aper / 2.0, aper / 2.0, tess_sigma)
            if options.debug > 1:
                print ''
                print 'geometric fractions of contaminant fluxes in target aperture'
                print tmp
            cflx = tmp * cflux
            cfle = abs(tmp * cfluxe)
            if options.debug > 1:
                print ''
                print 'contaminating fluxes and errors (fraction * contaminant fluxes)'
                print cflx
                print cfle
    
            # sum it up
            if options.debug > 1:
                print ''
                print 'sum of galactic contaminating flux and error'
                print np.sum(cflx), np.sum(cfle)
            sumcflx = sumcflx + np.sum(cflx)
            sumcfle = sumcfle + np.sum(cfle)
            if np.isnan(sumcflx):
                msg = 'sumcflx is nan after gal for ticlinkpk %d' % (row['ticlinkpk'],) 
                lf.write(msg)
                
        if options.debug > 1:
            print ''
            print 'candidate flux and error:'
            print tflux, tfluxe
            print 'total contaminating flux and error (stellar and galactic)'
            print sumcflx, sumcfle
        flxratio  = sumcflx / tflux
        flxratioe = sqrt((sumcfle / tflux)**2 + (sumcflx * tfluxe / tflux**2)**2)
        if np.isnan(flxratio):
            msg = 'flxratio is nan for ticlinkpk %d' % (row['ticlinkpk'],) 
            lf.write(msg)
        if options.debug > 1:
            print ''
            print 'contamination ratio = contaminating flux / target flux with error'
            print flxratio, flxratioe
            
        if np.isinf(flxratio):
            lf.write('INVESTIGATE: flxratio is inf for ticlinkpk ' + str(row['ticlinkpk']))

        updict  = {'nrcont' : nrcont, 'flxratio' : flxratio, 
                   'flxratioe' : flxratioe, 'ticlinkpk' : row['ticlinkpk']}
        updates.append(updict)
        
        nrupd = len(updates)
        if nrupd >= options.blen:
            blockno += 1
            contcnt += nrupd
            wmydb.singleTrans(updcmd, updates, lf, commit = True)
            ltime = watch.lap()
            msg = '%d, %d: %d stars in %f s, %f stars/s' % \
                  (blockno, row['ticlinkpk'], nrupd, ltime, nrupd / ltime)
            lf.write(msg)
            updates = []

    # check for remaining updates
    nrupd = len(updates)
    if nrupd > 0:
        blockno += 1
        contcnt  += nrupd
        wmydb.singleTrans(updcmd, updates, lf, commit = True)
        ltime = watch.lap()
        msg = '%d: %d stars in %f s, %f stars/s' % \
              (blockno, nrupd, ltime, nrupd / ltime)
        lf.write(msg)
        updates = []
        
    rmydb.close()
    wmydb.close()
    
    ltime = watch.stop()
    msg = '%d stars in %f s, %f stars/s' % \
          (nrstars, ltime, contcnt / ltime)
    lf.write(msg)
    
    lf.write('done')
