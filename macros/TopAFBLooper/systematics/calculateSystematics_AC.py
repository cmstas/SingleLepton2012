#!/usr/bin/env python

import sys,getopt,urllib2,json,math,fnmatch,numpy
from optparse import OptionParser
from numpy import *


#def calculateVariation(refNominal, systematic, binname):
#    """
#    Calculates the variation due to a particular systematic/subtype
#    """
#
#    vartype = systematic['vartype']
#    usenotes = systematic['usenotes']
#
#    if 'newnom' in usenotes: nominal = systematic['nominal'][binname]
#    else: nominal = refNominal
#
#    #Calculate the contribution from this systematic.
#    #Calculation method depends on the vartype flag, determined in parseSystematics.py
#
#    # Max of the 'up' and 'down' variations (potentially with respect to a special nominal value)
#    if vartype == 'updown' or vartype == 'updown_newnom':
#        upvar   = systematic['up'][binname][0]   - nominal[0]
#        downvar = systematic['down'][binname][0] - nominal[0]
#        varsign = sign(systematic['up'][binname][0] - systematic['down'][binname][0])
#        if abs(upvar) >= abs(downvar): maxvar  = abs(upvar)*varsign
#        elif abs(upvar) < abs(downvar): maxvar = abs(downvar)*varsign
#
#    # Half the difference between the 'up' and the 'down' variation
#    elif vartype == 'half':
#        maxvar = ( systematic['up'][binname][0] - systematic['down'][binname][0] ) /2
#
#    # Simple difference from nominal
#    elif vartype == 'diffnom':
#        maxvar = systematic['diffnom'][binname][0] - nominal[0]
#
#    # Anything else we don't understand
#    else:
#        print ""
#        print "I don't know what to do with this variation type:", vartype
#        #print systematic, subtype
#        sys.exit(1)
#
#    if 'divideby' in usenotes:
#        idx = usenotes.index('divideby')
#        factor = float( usenotes[idx+1] )
#        maxvar /= factor
#
#    return maxvar





def calculateVariation_combine(systematics, plotname, systematicname, subtypename, binname, allowmaxstat, combine, binlist, binlistcombine):
    """
    Calculates the variation due to a particular systematic/subtype
    """

    systematic = systematics[plotname][systematicname][subtypename]
    refNominal = systematics[plotname]['Nominal']['default']

    if(combine==1):
        systematiccombine = systematics['lepPlusCosTheta'][systematicname][subtypename]
        refNominalcombine = systematics['lepPlusCosTheta']['Nominal']['default']
    if(combine==-1):
        systematiccombine = systematics['lepMinusCosTheta'][systematicname][subtypename]
        refNominalcombine = systematics['lepMinusCosTheta']['Nominal']['default']

    binnamecombine = binname
    if 'bin' in binname and 'tt' not in binname and combine != 0: binnamecombine = binlist[ binlistcombine.index(binname) ]

    vartype = systematic['vartype']
    usenotes = systematic['usenotes']

    if 'newnom' in usenotes:
        nominal = systematic['nominal']
        if(combine != 0): nominalcombine = systematiccombine['nominal']
    else:
        nominal = refNominal['nominal']
        if(combine != 0): nominalcombine = refNominalcombine['nominal']

    nominalinclusive = nominal['Unfolded']
    if(combine != 0): nominalinclusivecombine = nominalcombine['Unfolded']

    #Calculate the contribution from this systematic.
    #Calculation method depends on the vartype flag, determined in parseSystematics.py

    # Max of the 'up' and 'down' variations (potentially with respect to a special nominal value)
    if vartype == 'updown' or vartype == 'updown_newnom':
        useupvariation = False
        upvarinclusive   = systematic['up']['Unfolded'][0]   - nominalinclusive[0]
        downvarinclusive = systematic['down']['Unfolded'][0] - nominalinclusive[0]
        halfdiffinclusive = (systematic['up']['Unfolded'][0] - systematic['down']['Unfolded'][0])/2.0
        if(combine != 0):
            upvarinclusive   += combine*(systematiccombine['up']['Unfolded'][0]   - nominalinclusivecombine[0])
            downvarinclusive += combine*(systematiccombine['down']['Unfolded'][0] - nominalinclusivecombine[0])
            halfdiffinclusive += combine*(systematiccombine['up']['Unfolded'][0] - systematiccombine['down']['Unfolded'][0])/2.0
        if abs(upvarinclusive) >= abs(downvarinclusive): useupvariation = True
        halfdiffSF = 1.0
        usehalfdiffwithSF = False

        #use both up and down variations where possible [except when both move in same direction] to get a more accurate measure of the individual bin variations
        #and abs(upvarinclusive/downvarinclusive )<100 and abs(downvarinclusive/upvarinclusive )<100 
        if( upvarinclusive!=0 and downvarinclusive!=0 and upvarinclusive/downvarinclusive < 0): 
            if useupvariation: halfdiffSF = abs(upvarinclusive/halfdiffinclusive)
            else: halfdiffSF = abs(downvarinclusive/halfdiffinclusive)
            usehalfdiffwithSF = True
            #print halfdiffSF

        if usehalfdiffwithSF:
            maxvar = ( systematic['up'][binname][0] - systematic['down'][binname][0] ) /2
            if(combine != 0):
                if 'bin' in binname and 'tt' not in binname:
                    maxvar  += ( systematiccombine['up'][binnamecombine][0] - systematiccombine['down'][binnamecombine][0] ) /2
                else:
                    maxvar  += combine*(systematiccombine['up'][binnamecombine][0] - systematiccombine['down'][binnamecombine][0] )/2
            maxvar *= halfdiffSF

        else:
            upvar   = systematic['up'][binname][0]   - nominal[binname][0]
            downvar = systematic['down'][binname][0] - nominal[binname][0]
            if(combine != 0):
                if 'bin' in binname and 'tt' not in binname:
                    upvar   += systematiccombine['up'][binnamecombine][0]   - nominalcombine[binnamecombine][0]
                    downvar += systematiccombine['down'][binnamecombine][0] - nominalcombine[binnamecombine][0]
                else:
                    upvar   += combine*(systematiccombine['up'][binnamecombine][0]   - nominalcombine[binnamecombine][0])
                    downvar += combine*(systematiccombine['down'][binnamecombine][0] - nominalcombine[binnamecombine][0])
            if useupvariation: maxvar  = upvar
            else: maxvar = -downvar

    # Half the difference between the 'up' and the 'down' variation
    elif vartype == 'half':
        maxvar = ( systematic['up'][binname][0] - systematic['down'][binname][0] ) /2
        if(combine != 0):
            if 'bin' in binname and 'tt' not in binname:
                maxvar  += ( systematiccombine['up'][binnamecombine][0] - systematiccombine['down'][binnamecombine][0] ) /2
            else:
                maxvar  += combine*(systematiccombine['up'][binnamecombine][0] - systematiccombine['down'][binnamecombine][0] )/2


    # Simple difference from nominal
    elif vartype == 'diffnom':
        maxvar = systematic['diffnom'][binname][0] - nominal[binname][0]
        if(combine != 0):
            if 'bin' in binname and 'tt' not in binname:
                maxvar  += systematiccombine['diffnom'][binnamecombine][0]   - nominalcombine[binnamecombine][0]
            else:
                maxvar  += combine*(systematiccombine['diffnom'][binnamecombine][0]   - nominalcombine[binnamecombine][0])



    # Anything else we don't understand
    else:
        print ""
        print "I don't know what to do with this variation type:", vartype
        #print systematic, subtype
        sys.exit(1)

    # Find the max stat error on all the directions, then take the max of that and the systematic error
    # Only do this for the inclusive asymmetry (because this is the best criterion for the statistical significance of the uncertainty), and multiply the covariance matrix by corresponding factor if necessary
    if allowmaxstat and 'maxstat' in usenotes and binname == 'Unfolded':
        varsyst = abs(maxvar)
        directionlist = systematic.keys()
        directionlist.remove('vartype')
        directionlist.remove('usenotes')
        staterrlist = [ systematic[i][binname][1] for i in directionlist ]
        if(combine != 0): staterrlistcombine = [ systematiccombine[i][binnamecombine][1] for i in directionlist ]
        varstat = 0
        #print vartype
        if vartype == 'half': varstat = sqrt(systematic['up'][binname][1]*systematic['up'][binname][1] + systematic['down'][binname][1]*systematic['down'][binname][1])/2.
        elif vartype == 'diffnom': varstat = systematic['diffnom'][binname][1] # nominal MC stat uncertainty is accounted for separately so don't add it here too
        else: varstat = max( staterrlist ) / sqrt(2.)  #this is an approximate upper limit on the stat uncertainty
        if(combine != 0):
            if vartype == 'half': varstat = sqrt(varstat*varstat + (systematiccombine['up'][binnamecombine][1]*systematiccombine['up'][binnamecombine][1] + systematiccombine['down'][binnamecombine][1]*systematiccombine['down'][binnamecombine][1])/4. )
            elif vartype == 'diffnom': varstat = sqrt(varstat*varstat + systematiccombine['diffnom'][binnamecombine][1]*systematiccombine['diffnom'][binnamecombine][1] ) # nominal MC stat uncertainty is accounted for separately so don't add it here too
            else: varstat = sqrt(varstat*varstat + max( staterrlistcombine )*max( staterrlistcombine )/2. )  #this is an approximate upper limit on the stat uncertainty
        #print  max( staterrlist )
        #print varstat
        #print varsyst
        if(maxvar>=0): maxvar = max(varsyst, varstat)
        else: maxvar = -max(varsyst, varstat)
        if systematicname == 'Scale' and (plotname == 'rapiditydiffMarco' or plotname == 'lepChargeAsym'):
            #the observed scaledown charge asymmetries seem to have unphysical fluctuations, so limit them to the MC stat uncertainty
            if(maxvar>=0): maxvar = varstat
            else: maxvar = -varstat
        #print maxvar


    if 'divideby' in usenotes:
        idx = usenotes.index('divideby')
        factor = float( usenotes[idx+1] )
        maxvar /= factor

    if(combine != 0):
        maxvar /= 2.;

    return maxvar








def calculateVariation(refNominal, systematic, binname, allowmaxstat = 0):
    """
    Calculates the variation due to a particular systematic/subtype
    """

    vartype = systematic['vartype']
    usenotes = systematic['usenotes']

    if 'newnom' in usenotes: nominal = systematic['nominal'][binname]
    else: nominal = refNominal['nominal'][binname]

    if 'newnom' in usenotes: nominalinclusive = systematic['nominal']['Unfolded']
    else: nominalinclusive = refNominal['nominal']['Unfolded']

    #Calculate the contribution from this systematic.
    #Calculation method depends on the vartype flag, determined in parseSystematics.py

    # Max of the 'up' and 'down' variations (potentially with respect to a special nominal value)
    if vartype == 'updown' or vartype == 'updown_newnom':
        useupvariation = False
        upvarinclusive   = systematic['up']['Unfolded'][0]   - nominalinclusive[0]
        downvarinclusive = systematic['down']['Unfolded'][0] - nominalinclusive[0]
        if abs(upvarinclusive) >= abs(downvarinclusive): useupvariation = True
        upvar   = systematic['up'][binname][0]   - nominal[0]
        downvar = systematic['down'][binname][0] - nominal[0]
        if useupvariation: maxvar  = upvar
        else: maxvar = -downvar

    # Half the difference between the 'up' and the 'down' variation
    elif vartype == 'half':
        maxvar = ( systematic['up'][binname][0] - systematic['down'][binname][0] ) /2

    # Simple difference from nominal
    elif vartype == 'diffnom':
        maxvar = systematic['diffnom'][binname][0] - nominal[0]

    # Anything else we don't understand
    else:
        print ""
        print "I don't know what to do with this variation type:", vartype
        #print systematic, subtype
        sys.exit(1)

    # Find the max stat error on all the directions, then take the max of that and the systematic error
    # Only do this for the inclusive asymmetry (because this is the best criterion for the statistical significance of the uncertainty), and multiply the covariance matrix by corresponding factor if necessary
    if allowmaxstat and 'maxstat' in usenotes and binname == 'Unfolded':
        varsyst = abs(maxvar)
        directionlist = systematic.keys()
        directionlist.remove('vartype')
        directionlist.remove('usenotes')
        staterrlist = [ systematic[i][binname][1] for i in directionlist ]
        varstat = 0
        #print vartype
        if vartype == 'half': varstat = sqrt(systematic['up'][binname][1]*systematic['up'][binname][1] + systematic['down'][binname][1]*systematic['down'][binname][1])/2.
        elif vartype == 'diffnom': varstat = systematic['diffnom'][binname][1] # nominal MC stat uncertainty is accounted for separately so don't add it here too
        else: varstat = max( staterrlist ) / sqrt(2.)  #this is an approximate upper limit on the stat uncertainty
        #print  max( staterrlist )
        #print varstat
        #print varsyst
        if(maxvar>=0): maxvar = max(varsyst, varstat)
        else: maxvar = -max(varsyst, varstat)
        #print maxvar


    if 'divideby' in usenotes:
        idx = usenotes.index('divideby')
        factor = float( usenotes[idx+1] )
        maxvar /= factor

    return maxvar





def GetCorrectedAfb(covarianceM, nbins, n):
    """
    WARNING: DOESN'T WORK FOR NON-UNIFORM BIN WIDTH
    """

    # Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    # Setup Alpha Vector
    alpha = zeros( [nbins] )
    for i in range(nbins):
        if i < nbins/2: alpha[i] = -1
        else: alpha[i] = 1

    # Components of the error calculation
    sum_n = 0.
    sum_alpha_n = 0.
    for i in range(nbins):
        sum_n += n[i]
        sum_alpha_n += alpha[i] * n[i]

    dfdn = zeros( [nbins] )
    for i in range(nbins):
        dfdn[i] = ( alpha[i] * sum_n - sum_alpha_n ) / pow(sum_n,2)

    # Error Calculation
    afberr = 0.
    for i in range(nbins):
        for j in range(nbins):
            afberr += covarianceM[i,j] * dfdn[i] * dfdn[j]

    afberr = sqrt(afberr)

    # Calculate Afb
    afb = sum_alpha_n / sum_n

    print "afb: %2.6f , afberr: %2.6f " % (afb, afberr )  
    return afberr



def GetCorrectedAfb_integratewidth_V(covarianceM, nbins, n, binwidth):

    # Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    # Setup Alpha Vector
    alpha = zeros( [nbins] )
    for i in range(nbins):
        if i < nbins/2: alpha[i] = -1
        else: alpha[i] = 1

    # Components of the error calculation
    sum_n = 0.
    sum_alpha_n = 0.
    for i in range(nbins):
        sum_n += n[i] * binwidth[i]
        sum_alpha_n += alpha[i] * n[i] * binwidth[i]

    dfdn = zeros( [nbins] )
    for i in range(nbins):
        dfdn[i] = ( alpha[i] * sum_n - sum_alpha_n ) / pow(sum_n,2)

    # Error Calculation
    afberr = 0.
    for i in range(nbins):
        for j in range(nbins):
            afberr += covarianceM[i,j] * dfdn[i] * dfdn[j] * binwidth[i] * binwidth[j]

    afberr = sqrt(afberr)

    # Calculate Afb
    afb = sum_alpha_n / sum_n

    #print "afb: %2.6f , afberr: %2.6f " % (afb, afberr )  
    return (afb,afberr)



def GetNormalisedCovarianceMatrix(covarianceM, nbins, n, binwidth):

    # Components of the error calculation
    sum_n = 0.
    for i in range(nbins):
        if binwidth[i]==0: binwidth[i]=1 #set dummy bin widths for the 2D distributions where the bin widths are not initialised because the bin contents are absolute numbers of events
        sum_n += n[i] * binwidth[i]

    #print sum_n
    if abs(sum_n-1.)>2e-6: print "WARNING: input distribution not normalised to unit area: %2.6g" % (sum_n-1.)
    sum_n = 1. #since sum_n should always = 1, set it exactly to minimise rounding errors

    newbinerr = zeros( [nbins] )
    newcov = zeros( [nbins,nbins] )
    #sum_n_minus_n = zeros( [nbins] )

    dfdn = zeros( [nbins,nbins] )
    for k in range(nbins):
        #sum_n_minus_n[k] = sum_n - n[k]
        for j in range(nbins):
            if k==j: dfdn[k,j] = (sum_n-n[k]*binwidth[k]) / pow(sum_n,2)
            else: dfdn[k,j] = (-n[k]*binwidth[k]) / pow(sum_n,2)

    # Error Calculation
    for k in range(nbins):
        for i in range(nbins):
            for j in range(nbins):
                newbinerr[k] += covarianceM[i,j] * dfdn[k,i] * dfdn[k,j] * binwidth[i] * binwidth[j]
                for l in range(nbins):
                    newcov[k][l] += covarianceM[i,j] * dfdn[k,i] * dfdn[l,j] * binwidth[i] * binwidth[j]
        newbinerr[k] = sqrt(newbinerr[k]/binwidth[k]/binwidth[k])

    for k in range(nbins):
        for l in range(nbins):
            newcov[k][l] /= (binwidth[k] * binwidth[l])

    #print newbinerr
    #print newcov
    binsumerr = GetBinSumUncertainty(newcov, nbins, n, binwidth)
    if (abs(binsumerr))>1e-15: print "WARNING: output of GetNormalisedCovarianceMatrix doesn't give zero uncertainty on the sum of all bins %2.3g" % (binsumerr)

    return newcov



def GetBinSumUncertainty(covarianceM, nbins, n, binwidth):

    #check output of GetNormalisedCovarianceMatrix really does give zero uncertainty on the sum of all bins
    for i in range(nbins):
        if binwidth[i]==0: binwidth[i]=1 #set dummy bin widths for the 2D distributions where the bin widths are not initialised because the bin contents are absolute numbers of events

    dfdn = zeros( [nbins] )
    for k in range(nbins): dfdn[k] = 1

    binsumerr = 0.

    # Error Calculation
    for i in range(nbins):
        for j in range(nbins):
            binsumerr += covarianceM[i,j] * dfdn[i] * dfdn[j] * binwidth[i] * binwidth[j]
    #print sqrt(abs(binsumerr))

    return binsumerr




def GetCorrectedAfb2D(covarianceM, nbins, nbins2D, n):

    # Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    numbinsx = nbins/nbins2D
    numbinsy = nbins2D

    # Setup Alpha Vector
    alpha = zeros( [numbinsx] )
    for i in range(numbinsx):
        if i < numbinsx/2: alpha[i] = -1
        else: alpha[i] = 1

    # Components of the error calculation
    sum_n = zeros( [numbinsy+1] )
    sum_alpha_n = zeros( [numbinsy+1] )
    for i in range(nbins):
        i_2di = i % numbinsx
        i_2dj = i / numbinsx
        sum_n[i_2dj] += n[i]
        sum_alpha_n[i_2dj] += alpha[i_2di] * n[i]
        sum_n[numbinsy] += n[i]
        sum_alpha_n[numbinsy] += alpha[i_2di] * n[i]

    dfdn = zeros( [nbins,numbinsy+1] )
    for i in range(nbins):
        i_2di = i % numbinsx
        i_2dj = i / numbinsx
        dfdn[i][i_2dj] = ( alpha[i_2di] * sum_n[i_2dj] - sum_alpha_n[i_2dj] ) / pow(sum_n[i_2dj],2)
        dfdn[i][numbinsy] = ( alpha[i_2di] * sum_n[numbinsy] - sum_alpha_n[numbinsy] ) / pow(sum_n[numbinsy],2)

    # Error Calculation
    afb = zeros( [numbinsy+1] )
    afberr = zeros( [numbinsy+1] )
    afbcov = zeros( [numbinsy,numbinsy] )
    for i in range(nbins):
        for j in range(nbins):
            i_2di = i % numbinsx
            i_2dj = i / numbinsx
            j_2di = j % numbinsx
            j_2dj = j / numbinsx
            afberr[numbinsy] += covarianceM[i,j] * dfdn[i][numbinsy] * dfdn[j][numbinsy]
            if(i_2dj==j_2dj): afberr[i_2dj] += covarianceM[i,j] * dfdn[i][i_2dj] * dfdn[j][i_2dj]
            afbcov[i_2dj][j_2dj] += covarianceM[i,j] * dfdn[i][i_2dj] * dfdn[j][j_2dj]

    # Calculate Afb
    for i in range(numbinsy+1):
        if afberr[i]<0: print "index: %i err %1.6g %1.6g %1.6g %1.6g" % (i, afberr[0], afberr[1], afberr[2], afberr[3])
        afberr[i] = sqrt(afberr[i])
        afb[i] = sum_alpha_n[i] / sum_n[i]

    #print afberr
    return (afb,afberr,afbcov)




def Get1DProjection(covarianceM, nbins, nbins2D, n):

    # Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    numbinsx = nbins/nbins2D
    numbinsy = nbins2D

    # Components of the error calculation
    sum_n = zeros( [numbinsy+1,numbinsx] )
    sum_n_err = zeros( [numbinsy+1,numbinsx] )
    sum_n_cov = zeros( [numbinsy+1,numbinsx,numbinsx] )
    for i in range(nbins):
        i_2di = i % numbinsx
        i_2dj = i / numbinsx
        sum_n[i_2dj][i_2di] += n[i]
        sum_n[numbinsy][i_2di] += n[i]

    # Error Calculation
    for i in range(nbins):
        for j in range(nbins):
            i_2di = i % numbinsx
            i_2dj = i / numbinsx
            j_2di = j % numbinsx
            j_2dj = j / numbinsx
            if(i_2di==j_2di): sum_n_err[numbinsy][i_2di] += covarianceM[i,j] 
            if(i_2di==j_2di and i_2dj==j_2dj): sum_n_err[i_2dj][i_2di] += covarianceM[i,j] 
            sum_n_cov[numbinsy][i_2di][j_2di] += covarianceM[i,j] 
            if(i_2dj==j_2dj): sum_n_cov[i_2dj][i_2di][j_2di] += covarianceM[i,j] 

    sum_n_err = sqrt(sum_n_err)
    # Calculate Afb to check covariance matrix calculation
    #for i in range(numbinsy+1):
    #    afberr = GetCorrectedAfb(sum_n_cov[i], numbinsx, sum_n[i])
    #    print afberr

    return (sum_n,sum_n_err,sum_n_cov)


def printCplusplusarray(covarianceM):
    numrows = len(covarianceM)
    numcols = len(covarianceM[0])
    print "float covmatrix[] = {"
    for row in range(numrows):
        for col in range(numcols):
            if col==0: print "{ ",
            if col!=numcols-1: print " %2.8g , " % (covarianceM[row,col]),
            else: print " %2.8g " % (covarianceM[row,col]),
            if col==numcols-1:
                if row!=numrows-1: print "},"
                else: print "}"
    print "};"


def increasevariance_integratewidth_V(covarianceM, nbins, n, binwidth, factor):

    #create new array for covarianceM0. covarianceM0 = covarianceM doesn't work because it binds the new name to the same object.
    covarianceM0 = zeros( [nbins] )
    for i in range(nbins):
        covarianceM0[i] = covarianceM[i,i]

    (afb,afberr0) = GetCorrectedAfb_integratewidth_V(covarianceM, nbins, n, binwidth)
    afberr = afberr0
    nit = 0
    while afberr < afberr0*factor:
        nit += 1
        for i in range(nbins):
            covarianceM[i,i] = covarianceM0[i] + 0.00000001*nit*n[i]*n[i]
            #covarianceM[i,i] = covarianceM0[i]*pow(1.0001,nit)
            #covarianceM[i,i] *= 1.0001
        (afb,afberr) = GetCorrectedAfb_integratewidth_V(covarianceM, nbins, n, binwidth)
    #if nit>0: print " %i iterations, %2.2f%% uncertainty added to covariance matrix" % (nit,100.*sqrt(0.00000001*nit))


def increasevariance2D(covarianceM, nbins, nbins2D, n, factor):

    #create new array for covarianceM0. covarianceM0 = covarianceM doesn't work because it binds the new name to the same object.
    covarianceM0 = zeros( [nbins] )
    for i in range(nbins):
        covarianceM0[i] = covarianceM[i,i]

    (afb,afberr0,afbcov) = GetCorrectedAfb2D(covarianceM, nbins, nbins2D, n)
    (afb,afberr,afbcov) = GetCorrectedAfb2D(covarianceM, nbins, nbins2D, n)
    nit = 0
    while afberr[nbins2D] < afberr0[nbins2D]*factor:
        nit += 1
        for i in range(nbins):
            covarianceM[i,i] = covarianceM0[i] + 0.0000001*nit*n[i]*n[i]
            #covarianceM[i,i] = covarianceM0[i]*pow(1.001,nit)
            #covarianceM[i,i] *= 1.0001
        (afb,afberr,afbcov) = GetCorrectedAfb2D(covarianceM, nbins, nbins2D, n)
    #if nit>0: print " %i iterations, %2.2f%% uncertainty added to covariance matrix" % (nit,100.*sqrt(0.0000001*nit))




def main():

    usage  = "Usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="verbose output")
    parser.add_option("-f", "--file", action="store", type="string", default=None, dest="systematicsjsonfilename", help="Name of systematics file in json format")
    parser.add_option("-n", "--nomatrix", action="store_true", default=False, dest="nomatrix", help="Suppress printout of covariance and correlation matrices")
    parser.add_option("-t", "--type", action="store", type="string", default=None, dest="unfoldingtype", help="None gives 1D, options are mtt, ttpt, ttrapidity2")
    parser.add_option("-c", "--combine", action="store_true", default=False, dest="allowcombine", help="Allow combination of A_P+ and A_P-")

    (opts, args) = parser.parse_args()

    numpy.set_printoptions(linewidth=800)
    numpy.set_printoptions(threshold=10000)
    numpy.set_printoptions(precision=6)
    
    if ( opts.systematicsjsonfilename == None ) :
        print ""
        print "Please specify systematics file in json format!"
        print ""
        parser.print_help()
        sys.exit(2)

    verbose = opts.verbose
    systematicsjsonfilename = opts.systematicsjsonfilename
    nomatrix = opts.nomatrix
    allowcombine = opts.allowcombine

    try:
        systematicsjsonfile = open(systematicsjsonfilename)
    except:
        print ""
        print "Systematics file in json format",systematicsjsonfilename," cannot be opened"
        print ""
        sys.exit(1)
        
    try:
        systematics = json.load(open(systematicsjsonfilename))
    except:
        print ""
        print "Systematics in systematics file",systematicsjsonfilename,"are not in json format"
        print ""
        sys.exit(1)
        
    #results = {}
    print ""

    systorderlist = ['JESMET','JER','LES','Btag','LepSelEff','PUshape','Background','Mass','Scale','PDF','hadronization','Regularisation','pt_reweight']

    #if allowcombine: asymorderlist = ['rapiditydiffMarco','lepChargeAsym','lepAzimAsym2','lepCosOpeningAngle','topSpinCorr','lepCosTheta','lepPlusCosTheta']
    #else: asymorderlist = ['rapiditydiffMarco','lepChargeAsym','lepAzimAsym2','lepCosOpeningAngle','topSpinCorr','lepPlusCosTheta','lepMinusCosTheta']

    if allowcombine: asymorderlist = ['rapiditydiffMarco','lepChargeAsym']
    else: asymorderlist = ['rapiditydiffMarco','lepChargeAsym']

    #if allowcombine: asymorderlist = ['lepChargeAsym']
    #else: asymorderlist = ['lepChargeAsym']

    #if allowcombine: asymorderlist = ['lepAzimAsym2','lepCosOpeningAngle','topSpinCorr','lepCosTheta','lepPlusCosTheta']
    #else: asymorderlist = ['lepAzimAsym2','lepCosOpeningAngle','topSpinCorr','lepPlusCosTheta','lepMinusCosTheta']

    #if allowcombine: asymorderlist = ['lepAzimAsym2','lepCosOpeningAngle','lepCosTheta','lepPlusCosTheta','topSpinCorr']
    #else: asymorderlist = ['lepAzimAsym2','lepCosOpeningAngle','lepPlusCosTheta','lepMinusCosTheta','topSpinCorr']

    afberrs = {}
    afberrs2D = {}
    afberr_MCstat = {}
    afberr_MCstat2D = {}
    afberrdatastat = {}
    afberrsyst = {}
    afberrtotal = {}
    afbs = {}
    plotKeyword = {}
    plotKeywordAsym = {}
    plotLatex = {}
    plotLatexAsym = {}
    plotLatexUnit = {}
    binlabels = {}
    printhepdata = False


#    for plot in sorted(systematics.keys()):
#        for systematic in sorted(systematics[plot].keys()):
#            afberrs[plot,systematic] = 0
#            afberrs2D[plot,systematic] = zeros( [3+1] )


    #set up strings for hepdata
    if opts.unfoldingtype == 'mtt':
        secvarLatex = '$M_\mathrm{t\\bar{t}}$'
        secvarLatexUnit = '$M_\mathrm{t\\bar{t}}$ (GeV)'
        secvarKeyword = 'MTTBAR'
        secvarbinlabels = [0., 430., 530., 9999.]
    if opts.unfoldingtype == 'ttpt':
        secvarLatex = '$p_\mathrm{T}^\mathrm{t\\bar{t}}$'
        secvarLatexUnit = '$p_\mathrm{T}^\mathrm{t\\bar{t}}$ (GeV)'
        secvarKeyword = 'PTTTBAR'
        secvarbinlabels = [0., 41., 92., 9999.]
    if opts.unfoldingtype == 'ttrapidity2':
        secvarLatex = '$|y_\mathrm{t\\bar{t}}|$'
        secvarLatexUnit = '$|y_\mathrm{t\\bar{t}}|$'
        secvarKeyword = 'ABSYTTBAR'
        secvarbinlabels = [0., 0.34, 0.75, 9999.]
        

    for plot in sorted(systematics.keys()):
        if plot == 'rapiditydiffMarco': plotKeyword[plot] = 'DELTAABSYTOP'
        elif plot == 'lepChargeAsym': plotKeyword[plot] = 'DELTAABSETALEP'
        elif plot == 'lepAzimAsym2': plotKeyword[plot] = 'DELTAPHI'
        elif plot == 'lepCosOpeningAngle': plotKeyword[plot] = 'COSOPENINGANGLE'
        elif plot == 'topSpinCorr': plotKeyword[plot] = 'C1C2'
        elif plot == 'lepCosTheta': plotKeyword[plot] = 'COSTHETASTAR'
        elif plot == 'lepPlusCosTheta':
            if allowcombine: plotKeyword[plot] = 'CPVCOSTHETASTAR'
            else:  plotKeyword[plot] = 'COSTHETASTARPLUS'
        elif plot == 'lepMinusCosTheta':
            if allowcombine: plotKeyword[plot] = 'CPCCOSTHETASTAR'
            else:  plotKeyword[plot] = 'COSTHETASTARMINUS'
        else: plotKeyword[plot] = ''

        if plot == 'rapiditydiffMarco': plotKeywordAsym[plot] = 'CHARGEASYM'
        elif plot == 'lepChargeAsym': plotKeywordAsym[plot] = 'LEPCHARGEASYM'
        elif plot == 'lepAzimAsym2': plotKeywordAsym[plot] = 'DELTAPHIASYM'
        elif plot == 'lepCosOpeningAngle': plotKeywordAsym[plot] = 'OPENINGANGLEASYM'
        elif plot == 'topSpinCorr': plotKeywordAsym[plot] = 'C1C2ASYM'
        elif plot == 'lepCosTheta': plotKeywordAsym[plot] = 'POLARIZATION'
        elif plot == 'lepPlusCosTheta':
            if allowcombine: plotKeywordAsym[plot] = 'CPVPOLARIZATION'
            else:  plotKeywordAsym[plot] = 'POLARIZATIONPLUS'
        elif plot == 'lepMinusCosTheta':
            if allowcombine: plotKeywordAsym[plot] = 'CPCPOLARIZATION'
            else:  plotKeywordAsym[plot] = 'POLARIZATIONMINUS'
        else: plotKeywordAsym[plot] = ''

        if plot == 'rapiditydiffMarco': plotLatex[plot] = '$\Delta|y_\mathrm{t}|$'
        elif plot == 'lepChargeAsym': plotLatex[plot] = '$\Delta|\eta_{\ell}|$'
        elif plot == 'lepAzimAsym2': plotLatex[plot] = '$\left|\Delta \phi_{\ell^+\ell^-}\\right|$'
        elif plot == 'lepCosOpeningAngle': plotLatex[plot] = '$\cos\\varphi$'
        elif plot == 'topSpinCorr': plotLatex[plot] = '$\cos\\theta^{\star}_{\ell^+} \cos\\theta^{\star}_{\ell^-}$'
        elif plot == 'lepCosTheta': plotLatex[plot] = '$\cos\\theta^{\star}_\ell$'
        elif plot == 'lepPlusCosTheta':
            if allowcombine: plotLatex[plot] = '$\cos\\theta^{\star}_{\ell,\mathrm{CPV}}$'
            else: plotLatex[plot] = '$\cos\\theta^{\star}_{\ell+}$'
        elif plot == 'lepMinusCosTheta':
            if allowcombine: plotLatex[plot] = '$\cos\\theta^{\star}_{\ell,\mathrm{CPC}}$'
            else: plotLatex[plot] = '$\cos\\theta^{\star}_{\ell-}$'
        else: plotLatex[plot] = ''

        if plot == 'rapiditydiffMarco': plotLatexUnit[plot] = '$\Delta|y_\mathrm{t}|$'
        elif plot == 'lepChargeAsym': plotLatexUnit[plot] = '$\Delta|\eta_{\ell}|$'
        elif plot == 'lepAzimAsym2': plotLatexUnit[plot] = '$\left|\Delta \phi_{\ell^+\ell^-}\\right|$ ($\pi$ radians)'
        elif plot == 'lepCosOpeningAngle': plotLatexUnit[plot] = '$\cos\\varphi$'
        elif plot == 'topSpinCorr': plotLatexUnit[plot] = '$\cos\\theta^{\star}_{\ell^+} \cos\\theta^{\star}_{\ell^-}$'
        elif plot == 'lepCosTheta': plotLatexUnit[plot] = '$\cos\\theta^{\star}_\ell$'
        elif plot == 'lepPlusCosTheta':
            if allowcombine: plotLatexUnit[plot] = '$\cos\\theta^{\star}_{\ell,\mathrm{CPV}}$'
            else: plotLatexUnit[plot] = '$\cos\\theta^{\star}_{\ell+}$'
        elif plot == 'lepMinusCosTheta':
            if allowcombine: plotLatexUnit[plot] = '$\cos\\theta^{\star}_{\ell,\mathrm{CPC}}$'
            else: plotLatexUnit[plot] = '$\cos\\theta^{\star}_{\ell-}$'
        else: plotLatexUnit[plot] = ''

        if plot == 'rapiditydiffMarco': plotLatexAsym[plot] = '$A_{\mathrm{C}}$'
        elif plot == 'lepChargeAsym': plotLatexAsym[plot] = '$A^{\mathrm{lep}}_{\mathrm{C}}$'
        elif plot == 'lepAzimAsym2': plotLatexAsym[plot] = '$A_{\Delta\phi}$'
        elif plot == 'lepCosOpeningAngle': plotLatexAsym[plot] = '$A_{\cos\\varphi}$'
        elif plot == 'topSpinCorr': plotLatexAsym[plot] = '$A_{c1c2}$'
        elif plot == 'lepCosTheta': plotLatexAsym[plot] = '$A_{P}$'
        elif plot == 'lepPlusCosTheta':
            if allowcombine: plotLatexAsym[plot] = '$A_{P}^{\\rm{CPV}}$'
            else: plotLatexAsym[plot] = '$A_{P+}$'
        elif plot == 'lepMinusCosTheta':
            if allowcombine: plotLatexAsym[plot] = '$A_{P}^{\\rm{CPC}}$'
            else: plotLatexAsym[plot] = '$A_{P-}$'
        else: plotLatexAsym[plot] = ''


        if plot == 'lepChargeAsym': binlabels[plot] =  [ -2., -68./60., -48./60., -32./60., -20./60., -8./60., 0., 8./60., 20./60., 32./60., 48./60., 68./60., 2.]
        elif plot == 'lepAzimAsym2': binlabels[plot] = [0., 5./60., 10./60., 15./60., 20./60., 25./60., 30./60., 35./60., 40./60., 45./60., 50./60., 55./60., 1.]
        elif plot == 'lepAzimAsym': binlabels[plot] = [-1., -50./60., -40./60., -30./60., -20./60., -10./60.,  0., 10./60., 20./60., 30./60., 40./60., 50./60., 1.]
        elif plot == 'topCosTheta': binlabels[plot] = [-1., -2./3., -1./3., 0., 1./3., 2./3., 1.]
        elif plot == 'pseudorapiditydiff': binlabels[plot] =  [ -2., -1.0, -28./60., 0., 28./60., 1.0, 2.]
        elif plot == 'rapiditydiff': binlabels[plot] =  [ -2., -1.0, -20./60., 0., 20./60., 1.0, 2.]
        elif plot == 'rapiditydiffMarco': binlabels[plot] =  [ -2., -44./60., -20./60., 0., 20./60., 44./60., 2.]
        elif plot == 'lepCosTheta': binlabels[plot] = [-1., -2./3., -1./3., 0., 1./3., 2./3., 1.]
        elif plot == 'lepPlusCosTheta': binlabels[plot] = [-1., -2./3., -1./3., 0., 1./3., 2./3., 1.]
        elif plot == 'lepMinusCosTheta': binlabels[plot] = [-1., -2./3., -1./3., 0., 1./3., 2./3., 1.]
        elif plot == 'topSpinCorr': binlabels[plot] = [-1., -0.4, -10./60., 0., 10./60., 0.4, 1.]
        elif plot == 'lepCosOpeningAngle': binlabels[plot] = [-1., -2./3., -1./3., 0., 1./3., 2./3., 1.]




    #Loop over the different asymmetry variables...
    #for plot in sorted(systematics.keys()):
    for plot in asymorderlist:

        if plot == 'lepCosThetaCPV': continue

        #if plot not in results.keys(): results[plot] = {}

        if plot == '': printhepdata = False
        #elif plot == 'rapiditydiffMarco': printhepdata = True
        #elif plot == 'lepChargeAsym': printhepdata = True
        elif plot == 'lepAzimAsym2': printhepdata = True
        elif plot == 'lepCosOpeningAngle': printhepdata = True
        elif plot == 'topSpinCorr': printhepdata = True
        elif plot == 'lepCosTheta': printhepdata = True
        elif plot == 'lepPlusCosTheta': printhepdata = True
        #elif plot == 'lepMinusCosTheta': printhepdata = True
        else: printhepdata = False

        #if not printhepdata: continue


        print 'Variable:', plot

        if 'Nominal' not in systematics[plot].keys(): 
            print 'Could not find Nominal values'
            sys.exit(1)

        combine = 0
        if(allowcombine):
            if plot == 'lepMinusCosTheta':
                combine = 1 #CP+ (add lepMinusCosTheta)
                combinevar = 'lepPlusCosTheta'
            if plot == 'lepPlusCosTheta':
                combine = -1 #CP- (add reflected lepMinusCosTheta)
                combinevar = 'lepMinusCosTheta'

        #Count how many bins we're dealing with
        typelist = list( systematics[plot]['Nominal']['default']['nominal'].keys() )
        binnames = fnmatch.filter(typelist, 'bin*')
        #binnames.sort()
        nbins = len(binnames)

        binlist = []
        binlistcombine = []
        if 'x' in binnames[1]:  #2D binning
            for i in range(1,3+1):
                for j in range(1, (nbins/3)+1):
                    binlist.append('bin' + str(i) + 'x' + str(j))
                    if(combine == -1): binlistcombine.append('bin' + str(i) + 'x' + str((nbins/3)+1 - j))
                    if(combine == 1): binlistcombine.append('bin' + str(i) + 'x' + str(j))
        else:
            binlist = ['bin'+str(i) for i in range(1, nbins+1)]   #1D binning
            if(combine == -1): binlistcombine = ['bin'+str(nbins+1-i) for i in range(1, nbins+1)]   #1D binning
            if(combine == 1): binlistcombine = ['bin'+str(i) for i in range(1, nbins+1)]   #1D binning

        if sorted(binlist) != sorted(binnames):
            print 'Problem: these lists of bin names don\'t match:'
            print sorted(binnames)
            print sorted(binlist)
            sys.exit(1)


        binlist2D = []
        if opts.unfoldingtype: binlist2D = [str(opts.unfoldingtype)+"bin1",str(opts.unfoldingtype)+"bin2",str(opts.unfoldingtype)+"bin3"]
        nbins2D = len(binlist2D)

        binlisttrue = []
        if opts.unfoldingtype: binlisttrue = ["True_Top_from_acceptance_denominator_bin0","True_Top_from_acceptance_denominator_bin1","True_Top_from_acceptance_denominator_bin2","True_Top_from_acceptance_denominator"]
        else: binlisttrue = ["True_Top_from_acceptance_denominator"]
        nbinstrue = len(binlisttrue)


        #Now read bin widths we're dealing with
        binwidthnames = fnmatch.filter(typelist, 'bwidth*')
        #binwidthnames.sort()
        binwidthlist = []
        if nbins2D==0: binwidthlist = ['bwidth'+str(i) for i in range(1, nbins+1)] 

        #print binwidthlist

        if sorted(binwidthlist) != sorted(binwidthnames) and nbins2D==0:
            print 'Problem: these lists of binwidth names don\'t match:'
            print sorted(binwidthnames)
            print sorted(binwidthlist)
            sys.exit(1)

        covnames = fnmatch.filter(typelist, 'cov*')
        #print covnames
        #print sorted(covnames)



        sumsq_total = 0
        sumsq_total_2D = zeros( [nbins2D] )
        denominator = zeros( [nbinstrue] )
        denominator_err = zeros( [nbinstrue] )
        covar_total = zeros( [nbins,nbins] )
        covar_total_incDataStat = zeros( [nbins,nbins] )
        covMCStat = zeros( [nbins,nbins] )
        covDataStat = zeros( [nbins,nbins] )
        #errMCStat = zeros( [nbins] )
        #errDataStat = zeros( [nbins] )
        corr_total  = zeros( [nbins,nbins] )
        #binsyst_preMCStat = zeros( [nbins] )

        bin_nominals_i = zeros( [nbins] )
        binwidth = zeros( [nbins] )

        #get the MC truth values
        binindex = 0
        for i in binlisttrue:
            denominator[binindex] = systematics[plot]['Nominal']['default']['nominal'][i][0]
            denominator_err[binindex] = systematics[plot]['Nominal']['default']['nominal'][i][1]
            if(combine == 1):
                denominator[binindex] += systematics[combinevar]['Nominal']['default']['nominal'][i][0]
                denominator[binindex] /= 2.
                denominator_err[binindex] = sqrt(denominator_err[binindex]*denominator_err[binindex] + systematics[combinevar]['Nominal']['default']['nominal'][i][1]*systematics[combinevar]['Nominal']['default']['nominal'][i][1])
                denominator_err[binindex] /= 2.
            if(combine == -1):
                denominator[binindex] -= systematics[combinevar]['Nominal']['default']['nominal'][i][0]
                denominator[binindex] /= 2.
                denominator_err[binindex] = sqrt(denominator_err[binindex]*denominator_err[binindex] + systematics[combinevar]['Nominal']['default']['nominal'][i][1]*systematics[combinevar]['Nominal']['default']['nominal'][i][1])
                denominator_err[binindex] /= 2.
            binindex += 1


        #Get the nominal values for each bin, and overall
        binindex = 0
        for i in binlist:
            bin_nominals_i[binindex] = systematics[plot]['Nominal']['default']['nominal'][i][0]
            if(combine != 0):
                bin_nominals_i[binindex] += systematics[combinevar]['Nominal']['default']['nominal'][binlist[ binlistcombine.index(i) ]][0]
                bin_nominals_i[binindex] /= 2.
            binindex += 1

        binindex = 0
        for i in binwidthlist:
            binwidth[binindex] = systematics[plot]['Nominal']['default']['nominal'][i][0]
            binindex += 1

        binindex = 0
        for i in sorted(covnames):
            rowindex = binindex/nbins
            colindex = binindex%nbins
            covMCStat[rowindex,colindex] = systematics[plot]['Nominal']['default']['nominal'][i][0]
            covDataStat[rowindex,colindex] = systematics[plot]['NominalDataStat']['default']['nominal'][i][0]
            if(combine == 1):
                covMCStat[rowindex,colindex] += systematics[combinevar]['Nominal']['default']['nominal'][i][0]
                covDataStat[rowindex,colindex] += systematics[combinevar]['NominalDataStat']['default']['nominal'][i][0]
            if(combine == -1):
                covMCStat[nbins-rowindex-1,nbins-colindex-1] += systematics[combinevar]['Nominal']['default']['nominal'][i][0]
                covDataStat[nbins-rowindex-1,nbins-colindex-1] += systematics[combinevar]['NominalDataStat']['default']['nominal'][i][0]
            #if rowindex==colindex:
                #errMCStat[rowindex] = sqrt(covMCStat[rowindex,rowindex])
                #errDataStat[rowindex] = sqrt(covDataStat[rowindex,rowindex])
            binindex += 1

        if(combine != 0):
                covMCStat /= 4.
                covDataStat /= 4.
                #errMCStat /= 2.
                #errDataStat /= 2.


        #Loop over the different systematics...
        for systematic in sorted(systematics[plot].keys()):
            if systematic == 'Nominal': continue
            if systematic == 'NominalDataStat': continue
            if systematic == 'RegularisationFullKIT': continue
            if systematic == 'name': continue

            sumsq_syst = 0
            sumsq_syst_2D = zeros( [nbins2D] )
            covar_syst = zeros( [nbins,nbins] )

            maxstat_factor = 1

            #Loop over the subtypes within each systematic...
            #print 'test0'
            #print  systematics[plot][systematic].keys()
            #for subtype in systematics[plot][systematic].keys(): print systematics[plot][systematic][subtype]['usenotes']
            if 'extrapolate' in (systematics[plot][systematic][subtype]['usenotes'][0] for subtype in systematics[plot][systematic].keys() ):
                #print 'extrapolating'
                #now do the extrapolation
                subtypenum =0
                weightedaveragegradient = 0.
                bin_weightedaveragegradient = zeros( [nbins] )
                bin_variations_0 = {}
                #sort the subtypes into numerical order before doing the extrapolation
                for subtype in sorted(systematics[plot][systematic].keys()):
                    #use only the first 6 points for the extrapolation
                    if subtypenum > 5: continue
                    sumsq_subtype = 0
                    sumsq_subtype_2D = zeros( [nbins2D] )
                    covar_subtype = zeros( [nbins,nbins] )
                    bin_variations = {}

                    #Calculate the variation in the overall asymmetry
                    var_overall = calculateVariation_combine( systematics, plot, systematic, subtype, 'Unfolded' , 0, combine, binlist, binlistcombine)
                    sumsq_subtype = var_overall*var_overall
                    if subtypenum == 0:
                        var_overall_maxstat = calculateVariation_combine( systematics, plot, systematic, subtype, 'Unfolded' , 1, combine, binlist, binlistcombine)
                        maxstat_factor = abs(var_overall_maxstat/var_overall)
                        if(maxstat_factor<1): print "something went wrong with maxstat_factor"

                    #Calculate the variation in the 2D asymmetries
                    for binindex in range(nbins2D):
                        var_overall_2D = calculateVariation_combine( systematics, plot, systematic, subtype, sorted(binlist2D)[binindex] , 0, combine, binlist, binlistcombine)
                        if subtypenum == 0: sumsq_subtype_2D[binindex] = var_overall_2D*var_overall_2D

                    #Calculate the variation in individual bins
                    for i in binlist: bin_variations[i] = calculateVariation_combine( systematics, plot, systematic, subtype, i , 0, combine, binlist, binlistcombine)
                    if subtypenum == 0: bin_variations_0 = bin_variations

                    #Calculate the covariance matrix using the individual bin variations:
                    #covar_ij = var(bin i) * var(bin j)
                    if subtypenum == 0:
                        for row in range(nbins):
                            for col in range(nbins):
                                covar_subtype[row, col] = bin_variations[binlist[row]] * bin_variations[binlist[col]]
                    else:
                        for i in range(nbins):
                            bin_weightedaveragegradient[i] += ( sqrt(covar_syst[i, i]) -  sqrt( bin_variations[binlist[i]] * bin_variations[binlist[i]] ) )/subtypenum


                    #Add to the running total of the variance and covariance
                    if subtypenum == 0:
                        sumsq_syst = sumsq_subtype
                        sumsq_syst_2D = sumsq_subtype_2D
                        covar_syst = covar_subtype
                    else:
                        weightedaveragegradient += (sqrt(sumsq_syst) - sqrt(sumsq_subtype))/subtypenum
                        #print (sqrt(sumsq_syst) - sqrt(sumsq_subtype))/subtypenum
                        #print 15.*weightedaveragegradient/subtypenum
                    subtypenum += 1
                weightedaveragegradient/=subtypenum
                bin_weightedaveragegradient/=subtypenum
                #print 15.*weightedaveragegradient

                extrapfactor = 1 + 15.*weightedaveragegradient / sqrt(sumsq_syst)
                #take maximum from extrapolation or maxstat_factor due to the MC statistical uncertainty (guaranteed to be >=1)
                print "          %s extrapolation factor: %2.3f , maxstat_factor: %2.3f " % (systematic,extrapfactor,maxstat_factor)
                #if maxstat_factor>1. and (plot == 'rapiditydiffMarco' or plot == 'lepChargeAsym'):
                #    maxstat_factor = 1.
                #    print "limiting maxstat_factor factor to 1.0 for charge asymmetry"
                if (plot == 'rapiditydiffMarco' or plot == 'lepChargeAsym'):
                    extrapfactor = maxstat_factor
                    print "setting extrapolation factor to maxstat_factor for charge asymmetry"
                if(extrapfactor<maxstat_factor):  extrapfactor = maxstat_factor
                if(extrapfactor>1000.00):
                    extrapfactor = 1000.00  #don't extrapolate more than a factor of 1000.00 to avoid producing highly (anti)correlated covariance matrix
                    print "limiting uncertainty scaling factor to 1000.00"
                #print "extrapolation factor: %2.3f , extrapolation amount: %2.3f " % (extrapfactor,100.*15.*weightedaveragegradient )
                sumsq_syst *= extrapfactor*extrapfactor
                sumsq_syst_2D *= extrapfactor*extrapfactor
                #instead of extrapolating for every bin, use same SF as for asymmetry (automatically ensures the covariance matrix is consistent with the uncertainty on the asymmetry, so no need to recalculate it)
                #covar_syst *= extrapfactor*extrapfactor
                #(afbtemp,afberrtemp,afbcovtemp) = GetCorrectedAfb2D(covar_syst, nbins, nbins2D, bin_nominals_i)
                #print "before increase factor %6.2f: %6.6g" % ( extrapfactor, afberrtemp[nbins2D] )
                if extrapfactor > 1.:
                    if nbins2D==0: increasevariance_integratewidth_V(covar_syst, nbins, bin_nominals_i, binwidth, extrapfactor)
                    else: increasevariance2D(covar_syst, nbins, nbins2D, bin_nominals_i, extrapfactor)
                else: covar_syst *= extrapfactor*extrapfactor
                #(afbtemp,afberrtemp,afbcovtemp) = GetCorrectedAfb2D(covar_syst, nbins, nbins2D, bin_nominals_i)
                #print "after increase factor %6.2f: %6.6g" % ( extrapfactor, afberrtemp[nbins2D] )

#The code below extrapolates the uncertainty on each bin individually, calculates the resulting new covariance matrix, then recalculates the inclusive uncertainty. Currently only working for variables with uniform bin width.
#Results are typically the same or less conservative than the simple method above, so use that instead.
#                bin_extrapfactor = zeros( [nbins] )
#                for i in range(nbins):
#                    bin_extrapfactor[i] = 1 + 15.*bin_weightedaveragegradient[i] / sqrt(covar_syst[i, i])
#                    print "extrapolation factor bin %1.1i: %2.3f , extrapolation amount: %2.3f " % (i,bin_extrapfactor[i],100.*15.*bin_weightedaveragegradient[i]  )     
#                    if bin_extrapfactor[i] < 1: bin_extrapfactor[i] = 1
#                    if bin_extrapfactor[i] > 3: bin_extrapfactor[i] = 3
#
#                for row in range(nbins):
#                    for col in range(nbins):
#                        covar_syst[row, col] = bin_variations_0[binlist[row]] * bin_variations_0[binlist[col]] * bin_extrapfactor[row] * bin_extrapfactor[col]
#
#                afberr = GetCorrectedAfb(covar_syst, nbins, bin_nominals_i)


            else:
                for subtype in systematics[plot][systematic].keys():
                    sumsq_subtype = 0
                    sumsq_subtype_2D = zeros( [nbins2D] )
                    covar_subtype = zeros( [nbins,nbins] )
                    bin_variations = {}

                    #Calculate the variation in the overall asymmetry
                    var_overall = calculateVariation_combine( systematics, plot, systematic, subtype, 'Unfolded' , 0, combine, binlist, binlistcombine)
                    sumsq_subtype = var_overall*var_overall
                    var_overall_maxstat = calculateVariation_combine( systematics, plot, systematic, subtype, 'Unfolded' , 1, combine, binlist, binlistcombine)
                    if var_overall == 0.0: maxstat_factor = 1
                    else: maxstat_factor = abs(var_overall_maxstat/var_overall)
                    if(maxstat_factor<1): print "something went wrong with maxstat_factor"
                    if(systematic == 'hadronization'):
                        if(plot=='lepAzimAsym' or plot=='lepAzimAsym2' or plot=='lepChargeAsym'):
                            maxstat_factor = 1
                            #print "setting maxstat_factor to 1 for purely leptonic variable because hadronization systematic only affects S matrix"
                        elif (maxstat_factor > 1):
                            maxstat_factor = maxstat_factor * 0.9 #hack to remove component due to acceptance correction
                            if (maxstat_factor<1): maxstat_factor = 1



                    #Calculate the variation in the 2D asymmetries
                    for binindex in range(nbins2D):
                        var_overall_2D = calculateVariation_combine( systematics, plot, systematic, subtype, sorted(binlist2D)[binindex] , 0, combine, binlist, binlistcombine)
                        sumsq_subtype_2D[binindex] = var_overall_2D*var_overall_2D


                    #Calculate the variation in individual bins
                    for i in binlist: bin_variations[i] = calculateVariation_combine( systematics, plot, systematic, subtype, i , 0, combine, binlist, binlistcombine)

                    #Calculate the covariance matrix using the individual bin variations:
                    #covar_ij = var(bin i) * var(bin j)
                    for row in range(nbins):
                        for col in range(nbins):
                            covar_subtype[row, col] = bin_variations[binlist[row]] * bin_variations[binlist[col]]

                    if(maxstat_factor>1): print "          %s maxstat_factor: %2.3f " % (systematic,maxstat_factor)
                    if(maxstat_factor>1000.00):
                        maxstat_factor = 1000.00  #don't extrapolate more than a factor of 1000.00 to avoid producing highly (anti)correlated covariance matrix
                        print "limiting uncertainty scaling factor to 1000.00"

                    #if nbins2D==0:
                        #(afb,afberr) = GetCorrectedAfb_integratewidth_V(covar_subtype, nbins, bin_nominals_i, binwidth)
                        #print "%15s subtype: %2.6f" % (subtype, afberr)
                    #else:
                        #(afb,afberr,afbcov) = GetCorrectedAfb2D(covar_subtype, nbins, nbins2D, bin_nominals_i)   #the 2D bin values represent normalised number of events and don't need to be multiplied by the bin widths
                        #print "%15s subtype: %2.6f" % (subtype, afberr[nbins2D])


                    #double-check GetNormalisedCovarianceMatrix has no effect on total matrix
                    #print " "
                    #print GetBinSumUncertainty(covar_subtype, nbins, bin_nominals_i, binwidth)
                    #covar_subtype_test = GetNormalisedCovarianceMatrix(covar_subtype, nbins, bin_nominals_i, binwidth)
                    #print GetBinSumUncertainty(covar_subtype_test, nbins, bin_nominals_i, binwidth)
                    #for i in range (nbins):
                    #    for j in range (nbins):
                    #        if abs(covar_subtype_test[i,j] - covar_subtype[i,j])/sqrt(covar_subtype[i,i]*covar_subtype[j,j])>1e-1: print "%15s systematic, %15s subtype , %15s vartype, WARNING: if this is not a rounding error, something went wrong with GetNormalisedCovarianceMatrix: %2.3g  %2.3g" % (systematic,subtype,systematics[plot][systematic][subtype]['vartype'],covar_subtype_test[i,j] , covar_subtype[i,j])



                    #Add to the running total of the variance and covariance
                    sumsq_syst += sumsq_subtype*maxstat_factor*maxstat_factor
                    sumsq_syst_2D += sumsq_subtype_2D*maxstat_factor*maxstat_factor
                    if maxstat_factor > 1.:
                        if nbins2D==0: increasevariance_integratewidth_V(covar_subtype, nbins, bin_nominals_i, binwidth, maxstat_factor)
                        else: increasevariance2D(covar_subtype, nbins, nbins2D, bin_nominals_i, maxstat_factor)
                    covar_syst += covar_subtype


                    #if nbins2D==0:
                        #(afb,afberr) = GetCorrectedAfb_integratewidth_V(covar_subtype, nbins, bin_nominals_i, binwidth)
                        #print "%15s syubtype2: %2.6f" % (subtype, afberr)
                    #else:
                        #(afb,afberr,afbcov) = GetCorrectedAfb2D(covar_subtype, nbins, nbins2D, bin_nominals_i)   #the 2D bin values represent normalised number of events and don't need to be multiplied by the bin widths
                        #print "%15s syubtype2: %2.6f" % (subtype, afberr[nbins2D])


            #the slight differences between afberr and sqrt(sumsq_syst) and sqrt(sumsq_syst_2D), particularly for the 2D asymmetries, are because GetCorrectedAfb(2D) uses the nominal bins whereas sumsq_syst(_2D) effectively uses the average of the up and down variations, which are not necessarily centred at nominal
            if nbins2D==0:
                (afb,afberrs[plot,systematic]) = GetCorrectedAfb_integratewidth_V(covar_syst, nbins, bin_nominals_i, binwidth)
                print "%15s systematic: %2.6f" % (systematic, afberrs[plot,systematic])
            else:
                (afb,afberrs2D[plot,systematic],afbcov) = GetCorrectedAfb2D(covar_syst, nbins, nbins2D, bin_nominals_i)   #the 2D bin values represent normalised number of events and don't need to be multiplied by the bin widths
                print "%15s systematic: %2.6f" % (systematic, afberrs2D[plot,systematic][nbins2D])
            for binindex in range(nbins2D): print "%25s %5s: %2.6f" % (systematic, sorted(binlist2D)[binindex], afberrs2D[plot,systematic][binindex])
            #end loop over subtypes
            if systematic != 'RegularisationFullKIT':
                sumsq_total += sumsq_syst
                sumsq_total_2D += sumsq_syst_2D
                covar_total += covar_syst

        #end loop over systematics

        #check afberr from covariance matrix is consistent with sqrt(sumsq_total) and sqrt(sumsq_total_2D) - again, the slight differences are because we use the nominal bins with the covariance matrix        
        if nbins2D==0:
            (afb,afberr_preMCStat) = GetCorrectedAfb_integratewidth_V(covar_total, nbins, bin_nominals_i, binwidth)
            if abs(afberr_preMCStat-math.sqrt(sumsq_total))>0.000002: print "WARNING: inconsistent covariance matrix. DeltaAfberr = %2.4g" % ( afberr_preMCStat-math.sqrt(sumsq_total) )
            #print "%20s uncertainty: %2.6f" % ("Sum systs", afberr_preMCStat)
            (afb,afberr_MCstat[plot]) = GetCorrectedAfb_integratewidth_V(covMCStat, nbins, bin_nominals_i, binwidth)
            print "%20s uncertainty: %2.6f" % ("MC stat", afberr_MCstat[plot])
            #(afb,afberr) = GetCorrectedAfb_integratewidth_V(covDataStat, nbins, bin_nominals_i, binwidth)
            #print "%20s uncertainty: %2.6f" % ("Data stat", afberr)
        else:
            (afb,afberr_preMCStat,afbcov) = GetCorrectedAfb2D(covar_total, nbins, nbins2D, bin_nominals_i)
            if abs(afberr_preMCStat[nbins2D]-math.sqrt(sumsq_total))>0.000002: print "WARNING: inconsistent covariance matrix. DeltaAfberr = %2.4g" % ( afberr_preMCStat[nbins2D]-math.sqrt(sumsq_total) )
            #print "%20s uncertainty: %2.6f" % ("Sum systs", afberr_preMCStat[nbins2D])
            #for binindex in range(nbins2D): print "%25s %5s: %2.6f" % ("Sum systs", sorted(binlist2D)[binindex], afberr_preMCStat[binindex])
            (afb,afberr_MCstat2D[plot],afbcov) = GetCorrectedAfb2D(covMCStat, nbins, nbins2D, bin_nominals_i)
            print "%20s uncertainty: %2.6f" % ("MC stat", afberr_MCstat2D[plot][nbins2D])
            for binindex in range(nbins2D): print "%25s %5s: %2.6f" % ("MC stat", sorted(binlist2D)[binindex], afberr_MCstat2D[plot][binindex])
            #(afb,afberr,afbcov) = GetCorrectedAfb2D(covDataStat, nbins, nbins2D, bin_nominals_i)
            #print "%20s uncertainty: %2.6f" % ("Data stat", afberr[nbins2D])
            #for binindex in range(nbins2D): print "%25s %5s: %2.6f" % ("Data stat", sorted(binlist2D)[binindex], afberr[binindex])


        #save bin systs for unfolded plots (unfolding macro automatically adds MC stat uncertainty so don't add it here)
        #for row in range(nbins): binsyst_preMCStat[row] = sqrt(covar_total[row,row])

#        print ""
#        #print covar_total
#        print binsyst_preMCStat
#        if nbins2D==0:
#            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(covar_total, nbins, bin_nominals_i, binwidth)
#            print "systerrbefore: %s = %2.6f +/- %2.6f " % (plot, afb, afberrdatastat[plot])
#        else:
#            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(covar_total, nbins, nbins2D, bin_nominals_i)
#            print "systerrbefore: %s = %2.6f +/- %2.6f " % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D])
#        (newbinerr,newcov) = GetNormalisedCovarianceMatrix(covar_total, nbins, bin_nominals_i, binwidth)
#        if nbins2D==0:
#            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(newcov, nbins, bin_nominals_i, binwidth)
#            print "systerrafter: %s = %2.6f +/- %2.6f " % (plot, afb, afberrdatastat[plot])
#        else:
#            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(newcov, nbins, nbins2D, bin_nominals_i)
#            print "systerrafter: %s = %2.6f +/- %2.6f " % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D])
#        print "repeat and confirm there is no change this time:"
#        if nbins2D==0:
#            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(newcov, nbins, bin_nominals_i, binwidth)
#            print "systerrbefore: %s = %2.6f +/- %2.6f " % (plot, afb, afberrdatastat[plot])
#        else:
#            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(newcov, nbins, nbins2D, bin_nominals_i)
#            print "systerrbefore: %s = %2.6f +/- %2.6f " % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D])
#        (newbinerr,newcov) = GetNormalisedCovarianceMatrix(newcov, nbins, bin_nominals_i, binwidth)
#        if nbins2D==0:
#            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(newcov, nbins, bin_nominals_i, binwidth)
#            print "systerrafter: %s = %2.6f +/- %2.6f " % (plot, afb, afberrdatastat[plot])
#        else:
#            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(newcov, nbins, nbins2D, bin_nominals_i)
#            print "systerrafter: %s = %2.6f +/- %2.6f " % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D])
#        print ""
#        #print covMCStat
#        print errMCStat
#        if nbins2D==0:
#            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(covMCStat, nbins, bin_nominals_i, binwidth)
#            print "MCstaterrbefore: %s = %2.6f +/- %2.6f " % (plot, afb, afberrdatastat[plot])
#        else:
#            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(covMCStat, nbins, nbins2D, bin_nominals_i)
#            print "MCstaterrbefore: %s = %2.6f +/- %2.6f " % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D])
#        (newbinerr,newcov) = GetNormalisedCovarianceMatrix(covMCStat, nbins, bin_nominals_i, binwidth)
#        if nbins2D==0:
#            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(newcov, nbins, bin_nominals_i, binwidth)
#            print "MCstaterrafter: %s = %2.6f +/- %2.6f " % (plot, afb, afberrdatastat[plot])
#        else:
#            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(newcov, nbins, nbins2D, bin_nominals_i)
#            print "MCstaterrafter: %s = %2.6f +/- %2.6f " % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D])
#        print ""
#        #print covDataStat
#        print errDataStat
#        if nbins2D==0:
#            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(covDataStat, nbins, bin_nominals_i, binwidth)
#            print "Datastaterrbefore: %s = %2.6f +/- %2.6f " % (plot, afb, afberrdatastat[plot])
#        else:
#            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(covDataStat, nbins, nbins2D, bin_nominals_i)
#            print "Datastaterrbefore: %s = %2.6f +/- %2.6f " % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D])
#        (newbinerr,newcov) = GetNormalisedCovarianceMatrix(covDataStat, nbins, bin_nominals_i, binwidth)
#        if nbins2D==0:
#            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(newcov, nbins, bin_nominals_i, binwidth)
#            print "Datastaterrafter: %s = %2.6f +/- %2.6f " % (plot, afb, afberrdatastat[plot])
#        else:
#            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(newcov, nbins, nbins2D, bin_nominals_i)
#            print "Datastaterrafter: %s = %2.6f +/- %2.6f " % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D])



        #add stat uncertainties from TUnfold
        covar_total+=covMCStat

        #correct covariance matrices so they correspond to normalised differential cross section (i.e. error on sum of bins = 0)
        covar_total = GetNormalisedCovarianceMatrix(covar_total, nbins, bin_nominals_i, binwidth)
        covDataStat = GetNormalisedCovarianceMatrix(covDataStat, nbins, bin_nominals_i, binwidth)

        covar_total_incDataStat = covar_total + covDataStat

        binsumerr = GetBinSumUncertainty(covar_total_incDataStat, nbins, bin_nominals_i, binwidth)
        if (abs(binsumerr))>1e-15: print "WARNING: sum of all bins doesn't have zero uncertainty %2.3g" % (binsumerr)

        #double-check GetNormalisedCovarianceMatrix has no effect on total matrix
        covar_total_incDataStat_test = GetNormalisedCovarianceMatrix(covar_total_incDataStat, nbins, bin_nominals_i, binwidth)
        for i in range (nbins):
            for j in range (nbins):
                if abs(covar_total_incDataStat_test[i,j] - covar_total_incDataStat[i,j])/sqrt(covar_total_incDataStat[i,i]*covar_total_incDataStat[j,j])>1e-7: print "WARNING: if this is not a rounding error, something went wrong with GetNormalisedCovarianceMatrix: %2.3g" % ((covar_total_incDataStat_test[i,j] - covar_total_incDataStat[i,j])/sqrt(covar_total_incDataStat[i,i]*covar_total_incDataStat[j,j]))

        if nbins2D==0:
            (afb,afberrdatastat[plot]) = GetCorrectedAfb_integratewidth_V(covDataStat, nbins, bin_nominals_i, binwidth)
            (afb,afberrtotal[plot]) = GetCorrectedAfb_integratewidth_V(covar_total_incDataStat, nbins, bin_nominals_i, binwidth)
            (afbs[plot],afberrsyst[plot]) = GetCorrectedAfb_integratewidth_V(covar_total, nbins, bin_nominals_i, binwidth)
            print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)  (%2.6f total)" % (plot, afb, afberrdatastat[plot], afberrsyst[plot], afberrtotal[plot])
            if plot=='lepAzimAsym2' or plot=='lepCosOpeningAngle':
                print "1D %s numbers for Jeremy:" % plot
                print bin_nominals_i
                for row in range(nbins): print " %2.6f " % (sqrt(covDataStat[row,row])),
                print "(stat)"
                for row in range(nbins): print " %2.6f " % (sqrt(covar_total[row,row])),
                print "(syst)"
                for row in range(nbins): print " %2.6f " % (sqrt(covar_total_incDataStat[row,row])),
                print "(total)"
                printCplusplusarray(covar_total_incDataStat)

                #double-check printed covariance matrix
                #(afb_temp,afberr_temp) = GetCorrectedAfb_integratewidth_V(covar_total_incDataStat, nbins, bin_nominals_i, binwidth)
                #binsumunc_temp = GetBinSumUncertainty(covar_total_incDataStat, nbins, bin_nominals_i, binwidth)
                #print (afb_temp,afberr_temp,binsumunc_temp)
                #normcov_temp = GetNormalisedCovarianceMatrix(covar_total_incDataStat, nbins, bin_nominals_i, binwidth)
                #(afb_temp,afberr_temp) = GetCorrectedAfb_integratewidth_V(normcov_temp, nbins, bin_nominals_i, binwidth)
                #binsumunc_temp = GetBinSumUncertainty(normcov_temp, nbins, bin_nominals_i, binwidth)
                #print (afb_temp,afberr_temp,binsumunc_temp)

        else:
            (afb,afberrdatastat[plot],afbcovdatastat) = GetCorrectedAfb2D(covDataStat, nbins, nbins2D, bin_nominals_i)
            (afb,afberrtotal[plot],afbcovtotalincDataStat) = GetCorrectedAfb2D(covar_total_incDataStat, nbins, nbins2D, bin_nominals_i)
            (afbs[plot],afberrsyst[plot],afbcovsyst) = GetCorrectedAfb2D(covar_total, nbins, nbins2D, bin_nominals_i)   #the 2D bin values represent normalised number of events and don't need to be multiplied by the bin widths
            print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)  (%2.6f total)" % (plot, afb[nbins2D], afberrdatastat[plot][nbins2D], afberrsyst[plot][nbins2D], afberrtotal[plot][nbins2D])
            (sum_n,sum_n_errdatastat,sum_n_covdatastat) = Get1DProjection(covDataStat, nbins, nbins2D, bin_nominals_i)
            (sum_n,sum_n_errtotal,sum_n_covtotalincDataStat) = Get1DProjection(covar_total_incDataStat, nbins, nbins2D, bin_nominals_i)
            (sum_n,sum_n_errsyst,sum_n_covsyst) = Get1DProjection(covar_total, nbins, nbins2D, bin_nominals_i)

            if plot=='lepAzimAsym2' or plot=='lepCosOpeningAngle':
                print "2D %s numbers for Jeremy:" % plot
                nbins_temp = nbins/nbins2D
                tempbinwidthcorrection = 1.
                if plot=='lepAzimAsym2': tempbinwidthcorrection = nbins_temp/math.pi  #divide deltaPhi projection by bin width (like we do for the pure 1D distribution) to produce numbers for Jeremy
                if plot=='lepCosOpeningAngle': tempbinwidthcorrection = nbins_temp/2.  #divide opening angle projection by bin width (like we do for the pure 1D distribution) to produce numbers for Jeremy
                #print tempbinwidthcorrection
                sum_n*=tempbinwidthcorrection
                sum_n_errdatastat*=tempbinwidthcorrection
                sum_n_covdatastat*=tempbinwidthcorrection*tempbinwidthcorrection
                sum_n_errtotal*=tempbinwidthcorrection
                sum_n_covtotalincDataStat*=tempbinwidthcorrection*tempbinwidthcorrection
                sum_n_errsyst*=tempbinwidthcorrection
                sum_n_covsyst*=tempbinwidthcorrection*tempbinwidthcorrection

                print sum_n[nbins2D]
                print sum_n_errdatastat[nbins2D]
                print sum_n_errsyst[nbins2D]
                print sum_n_errtotal[nbins2D]
                printCplusplusarray(sum_n_covtotalincDataStat[nbins2D])

                #double-check printed covariance matrix
                #binwidth_temp = zeros( [nbins_temp] )
                #for row in range(nbins_temp): binwidth_temp[row] = 1./tempbinwidthcorrection
                #(afb_temp,afberr_temp) = GetCorrectedAfb_integratewidth_V(sum_n_covtotalincDataStat[nbins2D], nbins_temp, sum_n[nbins2D], binwidth_temp)
                #binsumunc_temp = GetBinSumUncertainty(sum_n_covtotalincDataStat[nbins2D], nbins_temp, sum_n[nbins2D], binwidth_temp)
                #print (afb_temp,afberr_temp,binsumunc_temp)
                #normcov_temp = GetNormalisedCovarianceMatrix(sum_n_covtotalincDataStat[nbins2D], nbins_temp, sum_n[nbins2D], binwidth_temp)
                #(afb_temp,afberr_temp) = GetCorrectedAfb_integratewidth_V(normcov_temp, nbins_temp, sum_n[nbins2D], binwidth_temp)
                #binsumunc_temp = GetBinSumUncertainty(normcov_temp, nbins_temp, sum_n[nbins2D], binwidth_temp)
                #print (afb_temp,afberr_temp,binsumunc_temp)

        #print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)" % (plot, systematics[plot]['Nominal']['default']['nominal']['Unfolded'][0], systematics[plot]['Nominal']['default']['nominal']['Unfolded'][1], math.sqrt(sumsq_total))

        binindex = 0
        for i in sorted(binlist2D):
            #print "%s %s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)" % (plot, i, systematics[plot]['Nominal']['default']['nominal'][i][0], systematics[plot]['Nominal']['default']['nominal'][i][1], math.sqrt(sumsq_total_2D[binindex]))
            print "%s %s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)  (%2.6f total)" % (plot, i, afb[binindex], afberrdatastat[plot][binindex], afberrsyst[plot][binindex], afberrtotal[plot][binindex])
            binindex += 1
        print ""

        for i in range (nbinstrue):
            print "%s %s = %2.6f +/- %2.6f (stat)" % (plot, binlisttrue[i], denominator[i], denominator_err[i])
        print ""



        #Now calculate correlation matrix
        for row in range(nbins):
            for col in range(nbins):
                corr_total[row,col] = covar_total[row,col] / math.sqrt( covar_total[row,row] * covar_total[col,col] )

        afbcorsyst = zeros( [nbins2D,nbins2D] )
        afbcordatastat = zeros( [nbins2D,nbins2D] )
        afbcortotalincDataStat = zeros( [nbins2D,nbins2D] )

        for row in range(nbins2D):
            for col in range(nbins2D):
                afbcorsyst[row,col] = afbcovsyst[row,col] / math.sqrt( afbcovsyst[row,row] * afbcovsyst[col,col] )
                afbcordatastat[row,col] = afbcovdatastat[row,col] / math.sqrt( afbcovdatastat[row,row] * afbcovdatastat[col,col] )
                afbcortotalincDataStat[row,col] = afbcovtotalincDataStat[row,col] / math.sqrt( afbcovtotalincDataStat[row,row] * afbcovtotalincDataStat[col,col] )

        if nomatrix == False:
            print "%s syst covariance matrix:" % plot
            print binlist
            print covar_total
            print ""
            print "%s stat covariance matrix:" % plot
            print binlist
            print covDataStat
            print ""
            #print "%s bin systematics:" % plot
            #print binlist
            #print binsyst_preMCStat
            #print ""
            if nbins2D>0:
                print "%s total covariance matrix:" % plot
                print binlist2D
                print afbcovtotalincDataStat
                print ""
                print "%s total correlation matrix:" % plot
                print binlist2D
                print afbcortotalincDataStat
                print ""
                #print "%s stat correlation matrix:" % plot
                #print binlist2D
                #print afbcordatastat
                #print ""
                #print "%s syst correlation matrix:" % plot
                #print binlist2D
                #print afbcorsyst
                #print ""
        print "code fragment to paste in AfbFinalUnfold.h:"
        #if nbins2D==0: 
        #    for row in range(nbins): print "syst_corr[%i] = %2.6f;" % (row, binsyst_preMCStat[row])
        #else:
        #    for binindex in range(nbins2D): print "syst_corr[%i] = %2.6f;" % (binindex, afberr_preMCStat[binindex] )
        #print ""

        
        if nbins2D==0: 
            if combine!=0 and plot == "lepPlusCosTheta":
                print " special hard-coded bin values for A_P^CPV: "
                for row in range(nbins): print "hardcodeddata[%i] = %2.6f;" % (row, bin_nominals_i[row])
            for row in range(nbins): print "stat_corr[%i] = %2.6f;" % (row, sqrt(covDataStat[row,row]))
            for row in range(nbins): print "syst_corr[%i] = %2.6f;" % (row, sqrt(covar_total[row,row]))
        else:
            if combine!=0 and plot == "lepPlusCosTheta":
                print " special hard-coded bin values for A_P^CPV: "
                for binindex in range(nbins2D): print "hardcodeddata[%i] = %2.6f;" % (binindex, afbs[plot][binindex] )
            for binindex in range(nbins2D): print "stat_corr[%i] = %2.6f;" % (binindex, afberrdatastat[plot][binindex] )
            for binindex in range(nbins2D): print "syst_corr[%i] = %2.6f;" % (binindex, afberrsyst[plot][binindex] )
        print ""



        #hepdata printouts

        if nbins2D==0: 
            #bin values
            print "*dataset:"
            print "*location: Figure 3"
            print "*dscomment: Values of the %i bins of the normalized differential cross section as a function of %s." % (nbins, plotLatex[plot])
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s" % (plotKeyword[plot])
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader: Value"
            print "*xheader: Bin of %s" % (plotLatexUnit[plot])
            print "*data: x : y"
            for col in range(nbins):
                print "** %2.4f TO %2.4f; %2.4f;" % (binlabels[plot][col],binlabels[plot][col+1], bin_nominals_i[col])
            print "*dataend:\n"


            #stat covariance matrix
            print "*dataset:"
            print "*location: Figure 3"
            print "*dscomment: Statistical covariance matrix for the %i bins of the normalized differential cross section as a function of %s." % (nbins, plotLatex[plot])
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s" % (plotKeyword[plot])
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader:",
            for col in range(nbins):
                if col<nbins-1: print "$%i$ :" % (col+1),
                else: print "$%i$" % (col+1)
            print "*xheader: Bin of %s" % (plotLatex[plot])
            print "*data: x",
            for col in range(nbins): print ": y",
            print ""
            for row in range(nbins):
                print "** $%i$;" % (row+1),
                for col in range(nbins): print " %2.4g;" % (covDataStat[row,col]),
                print ""
            print "*dataend:\n"

            #syst covariance matrix
            print "*dataset:"
            print "*location: Figure 3"
            print "*dscomment: Systematic covariance matrix for the %i bins of the normalized differential cross section as a function of %s." % (nbins, plotLatex[plot])
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s" % (plotKeyword[plot])
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader:",
            for col in range(nbins):
                if col<nbins-1: print "$%i$ :" % (col+1),
                else: print "$%i$" % (col+1)
            print "*xheader: Bin of %s" % (plotLatex[plot])
            print "*data: x",
            for col in range(nbins): print ": y",
            print ""
            for row in range(nbins):
                print "** $%i$;" % (row+1),
                for col in range(nbins): print " %2.4g;" % (covar_total[row,col]),
                print ""
            print "*dataend:\n"

        else:

            #bin values
            print "*dataset:"
            print "*location: Figure 4"
            print "*dscomment: Values of %s in the %i bins of %s, and inclusively (bottom row). The value 9999 is used as a placeholder for infinity." % (plotLatexAsym[plot], nbins2D, secvarLatex)
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s %s" % (plotKeywordAsym[plot], secvarKeyword)
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader: %s" % (plotLatexAsym[plot])
            print "*xheader: Bin of %s" % (secvarLatexUnit)
            print "*data: x : y"
            for col in range(nbins2D):
                if opts.unfoldingtype == 'ttrapidity2': print "** %2.2f TO %2.2f; %2.4f +- %2.4f (DSYS=%2.4f);" % (secvarbinlabels[col],secvarbinlabels[col+1], afbs[plot][col], afberrdatastat[plot][col], afberrsyst[plot][col])
                else: print "** %2.0f TO %2.0f; %2.4f +- %2.4f (DSYS=%2.4f);" % (secvarbinlabels[col],secvarbinlabels[col+1], afbs[plot][col], afberrdatastat[plot][col], afberrsyst[plot][col])
            if opts.unfoldingtype == 'ttrapidity2': print "** %2.2f TO %2.2f; %2.4f +- %2.4f (DSYS=%2.4f);" % (0.,9999., afbs[plot][nbins2D], afberrdatastat[plot][nbins2D], afberrsyst[plot][nbins2D])
            else: print "** %2.0f TO %2.0f; %2.4f +- %2.4f (DSYS=%2.4f);" % (0.,9999., afbs[plot][nbins2D], afberrdatastat[plot][nbins2D], afberrsyst[plot][nbins2D])
            print "*dataend:\n"


            #stat correlation matrix
            print "*dataset:"
            print "*location: Figure 4"
            print "*dscomment: Statistical correlation matrix for %s in the %i bins of %s." % (plotLatexAsym[plot], nbins2D, secvarLatex)
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s %s" % (plotKeywordAsym[plot], secvarKeyword)
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader:",
            for col in range(nbins2D):
                if col<nbins2D-1: print "$%i$ :" % (col+1),
                else: print "$%i$" % (col+1)
            print "*xheader: Bin of %s" % (secvarLatex)
            print "*data: x",
            for col in range(nbins2D): print ": y",
            print ""
            for row in range(nbins2D):
                print "** $%i$;" % (row+1),
                for col in range(nbins2D): print " %2.4g;" % (afbcordatastat[row,col]),
                print ""
            print "*dataend:\n"

            #syst correlation matrix
            print "*dataset:"
            print "*location: Figure 4"
            print "*dscomment: Systematic correlation matrix for %s in the %i bins of %s." % (plotLatexAsym[plot], nbins2D, secvarLatex)
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s %s" % (plotKeywordAsym[plot], secvarKeyword)
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader:",
            for col in range(nbins2D):
                if col<nbins2D-1: print "$%i$ :" % (col+1),
                else: print "$%i$" % (col+1)
            print "*xheader: Bin of %s" % (secvarLatex)
            print "*data: x",
            for col in range(nbins2D): print ": y",
            print ""
            for row in range(nbins2D):
                print "** $%i$;" % (row+1),
                for col in range(nbins2D): print " %2.4g;" % (afbcorsyst[row,col]),
                print ""
            print "*dataend:\n"




            #2D bin values
            print "*dataset:"
            print "*location: Figure 4"
            print "*dscomment: Fraction of events in each of the $%i\\times%i$ bins of the normalized differential cross section as a function of %s and %s. The value 9999 is used as a placeholder for infinity." % (nbins/nbins2D, nbins2D, plotLatex[plot], secvarLatex)
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s %s" % (plotKeyword[plot], secvarKeyword)
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader: Fraction of events"
            print "*xheader: Bin of %s : Bin of %s" % (plotLatexUnit[plot], secvarLatexUnit)
            print "*data: x : x : y"
            for col in range(nbins):
                if opts.unfoldingtype == 'ttrapidity2': print "** %2.4f TO %2.4f; %2.2f TO %2.2f; %2.4f;" % (binlabels[plot][col%(nbins/nbins2D)],binlabels[plot][1+(col)%(nbins/nbins2D)], secvarbinlabels[col/(nbins/nbins2D)], secvarbinlabels[1+col/(nbins/nbins2D)], bin_nominals_i[col])
                else: print "** %2.4f TO %2.4f; %2.0f TO %2.0f; %2.4f;" % (binlabels[plot][col%(nbins/nbins2D)],binlabels[plot][1+(col)%(nbins/nbins2D)], secvarbinlabels[col/(nbins/nbins2D)], secvarbinlabels[1+col/(nbins/nbins2D)], bin_nominals_i[col])
            print "*dataend:\n"


            #stat covariance matrix
            print "*dataset:"
            print "*location: Figure 4"
            print "*dscomment: Statistical covariance matrix for the $%i\\times%i$ bins of the normalized differential cross section as a function of %s and %s." % (nbins/nbins2D, nbins2D, plotLatex[plot], secvarLatex)
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s %s" % (plotKeyword[plot], secvarKeyword)
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader:",
            for col in range(nbins):
                if col<nbins-1: print "$%i\\times%i$ :" % (1+(col)%(nbins/nbins2D),1+(col)/(nbins/nbins2D)),
                else: print "$%i\\times%i$" % (1+(col)%(nbins/nbins2D),1+(col)/(nbins/nbins2D))
            print "*xheader: Bin of %s and %s" % (plotLatex[plot], secvarLatex)
            print "*data: x",
            for col in range(nbins): print ": y",
            print ""
            for row in range(nbins):
                print "** $%i\\times%i$;" % (1+(row)%(nbins/nbins2D),1+(row)/(nbins/nbins2D)),
                for col in range(nbins): print " %2.4g;" % (covDataStat[row,col]),
                print ""
            print "*dataend:\n"

            #syst covariance matrix
            print "*dataset:"
            print "*location: Figure 4"
            print "*dscomment: Systematic covariance matrix for the $%i\\times%i$ bins of the normalized differential cross section as a function of %s and %s." % (nbins/nbins2D, nbins2D, plotLatex[plot], secvarLatex)
            print "*reackey: P P --> TOP TOPBAR X"
            #print "*obskey: %s %s" % (plotKeyword[plot], secvarKeyword)
            print "*qual: RE : P P --> TOP TOPBAR X"
            print "*qual: SQRT(S) IN GEV : 8000.0"
            print "*yheader:",
            for col in range(nbins):
                if col<nbins-1: print "$%i\\times%i$ :" % (1+(col)%(nbins/nbins2D),1+(col)/(nbins/nbins2D)),
                else: print "$%i\\times%i$" % (1+(col)%(nbins/nbins2D),1+(col)/(nbins/nbins2D))
            print "*xheader: Bin of %s and %s" % (plotLatex[plot], secvarLatex)
            print "*data: x",
            for col in range(nbins): print ": y",
            print ""
            for row in range(nbins):
                print "** $%i\\times%i$;" % (1+(row)%(nbins/nbins2D),1+(row)/(nbins/nbins2D)),
                for col in range(nbins): print " %2.4g;" % (covar_total[row,col]),
                print ""
            print "*dataend:\n"            


    #end loop over asymmetries

    if nbins2D==0:
        #table of asymmetries
        print "*dataset:"
        print "*location: Table 5"
        print "*dscomment: Inclusive values of the asymmetry variables."
        print "*reackey: P P --> TOP TOPBAR X"
        #print "*obskey: "
        print "*qual: RE : P P --> TOP TOPBAR X"
        print "*qual: SQRT(S) IN GEV : 8000.0"
        print "*yheader: Data (Unfolded)"
        print "*xheader: Asymmetry variable"
        print "*data: x : y"
        for plot in asymorderlist:
            print "** %s; %2.4f +- %2.4f (DSYS=%2.4f);" % (plotLatexAsym[plot], afbs[plot], afberrdatastat[plot], afberrsyst[plot])
        print "*dataend:\n"



    #print table of systematics
    #print "\\begin{table}[!h]"
    #print "\\begin{center}"
    #print "\\begin{tabular}{c||r|r|r|r|r}"
    print sorted(systematics['lepAzimAsym2'].keys())

    if nbins2D==0:
        #for systematic in sorted(systematics['lepAzimAsym2'].keys()):
        for systematic in systorderlist:
            if systematic == 'Nominal': continue
            if systematic == 'NominalDataStat': continue
            if systematic == 'RegularisationFullKIT': continue
            if systematic == 'name': continue
            print systematic,
            #print " & ",
            for plot in asymorderlist:
                print " & $%2.3f$ " % afberrs[plot,systematic],
            print " \\\\ "
        print "MCstat",
        for plot in asymorderlist:
            print " & $%2.3f$ " % afberr_MCstat[plot],
        print " \\\\ "
        print "totalsyst",
        for plot in asymorderlist:
            print " & $%2.3f$ " % afberrsyst[plot],
        print " \\\\ "
        print " "
        print "(total uncertainty)",
        for plot in asymorderlist:
            print " & $%2.4f$ " % afberrtotal[plot],
        print " \\\\ "
        print " "
    else:
        for binindex in range(nbins2D+1):
            if binindex<nbins2D: print binlist2D[binindex]
            else: print "inclusive"
            #for systematic in sorted(systematics['lepAzimAsym2'].keys()):
            for systematic in systorderlist:
                if systematic == 'Nominal': continue
                if systematic == 'NominalDataStat': continue
                if systematic == 'RegularisationFullKIT': continue
                if systematic == 'name': continue
                print systematic,
                #print " & ",
                #for plot in sorted(systematics.keys()):
                for plot in asymorderlist:
                    print " & $%2.3f$ " % afberrs2D[plot,systematic][binindex],
                print " \\\\ "
            print "MCstat",
            for plot in asymorderlist:
                print " & $%2.3f$ " % afberr_MCstat2D[plot][binindex],
            print " \\\\ "
            print "totalsyst",
            for plot in asymorderlist:
                print " & $%2.3f$ " % afberrsyst[plot][binindex],
            print " \\\\ "
            print " "
            print "(total uncertainty)",
            for plot in asymorderlist:
                print " & $%2.4f$ " % afberrtotal[plot][binindex],
            print " \\\\ "
            print " "

    if nbins2D==0:
        for plot in asymorderlist:
            print "%s & $%2.3f \\pm %2.3f \\pm %2.3f$ " % (plotLatexAsym[plot], afbs[plot], afberrdatastat[plot], afberrsyst[plot])
    else:
        for binindex in range(nbins2D+1):
            if binindex<nbins2D: print binlist2D[binindex]
            else: print "inclusive"
            for plot in asymorderlist:
                print "%s & $%2.3f \\pm %2.3f \\pm %2.3f$ " % (plotLatexAsym[plot], afbs[plot][binindex], afberrdatastat[plot][binindex], afberrsyst[plot][binindex])

        print " "
        print " All bins: "
        for plot in asymorderlist:
            print plotLatexAsym[plot],
            for binindex in range(nbins2D): print " & $%2.3f \\pm %2.3f \\pm %2.3f$ " % (afbs[plot][binindex], afberrdatastat[plot][binindex], afberrsyst[plot][binindex]),
            print "\\\\"
 


if __name__ == '__main__':
    main()
