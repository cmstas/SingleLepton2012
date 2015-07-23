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
        if(combine != 0):
            upvarinclusive   += combine*(systematiccombine['up']['Unfolded'][0]   - nominalinclusivecombine[0])
            downvarinclusive += combine*(systematiccombine['down']['Unfolded'][0] - nominalinclusivecombine[0])
        if abs(upvarinclusive) >= abs(downvarinclusive): useupvariation = True
        upvar   = systematic['up'][binname][0]   - nominal[binname][0]
        downvar = systematic['down'][binname][0] - nominal[binname][0]
        if(combine != 0):
            if 'bin' in binname:
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
            if 'bin' in binname:
                maxvar  += ( systematiccombine['up'][binnamecombine][0] - systematiccombine['down'][binnamecombine][0] ) /2
            else:
                maxvar  += combine*(systematiccombine['up'][binnamecombine][0] - systematiccombine['down'][binnamecombine][0] )/2


    # Simple difference from nominal
    elif vartype == 'diffnom':
        maxvar = systematic['diffnom'][binname][0] - nominal[binname][0]
        if(combine != 0):
            if 'bin' in binname:
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
            covarianceM[i,i] = covarianceM0[i] + 0.0000001*nit*n[i]*n[i]
            #covarianceM[i,i] = covarianceM0[i]*pow(1.0001,nit)
            #covarianceM[i,i] *= 1.0001
        (afb,afberr) = GetCorrectedAfb_integratewidth_V(covarianceM, nbins, n, binwidth)


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
            covarianceM[i,i] = covarianceM0[i] + 0.000001*nit*n[i]*n[i]
            #covarianceM[i,i] = covarianceM0[i]*pow(1.001,nit)
            #covarianceM[i,i] *= 1.0001
        (afb,afberr,afbcov) = GetCorrectedAfb2D(covarianceM, nbins, nbins2D, n)




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

    #Loop over the different asymmetry variables...
    for plot in sorted(systematics.keys()):

        #if plot not in results.keys(): results[plot] = {}

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
        corr_total  = zeros( [nbins,nbins] )
        binsyst_preMCStat = zeros( [nbins] )

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
            binindex += 1

        if(combine != 0):
                covMCStat /= 4.
                covDataStat /= 4.


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
                (afb,afberr) = GetCorrectedAfb_integratewidth_V(covar_syst, nbins, bin_nominals_i, binwidth)
                print "%15s systematic: %2.6f" % (systematic, afberr)
            else:
                (afb,afberr,afbcov) = GetCorrectedAfb2D(covar_syst, nbins, nbins2D, bin_nominals_i)   #the 2D bin values represent normalised number of events and don't need to be multiplied by the bin widths
                print "%15s systematic: %2.6f" % (systematic, afberr[nbins2D])
            for binindex in range(nbins2D): print "%25s %5s: %2.6f" % (systematic, sorted(binlist2D)[binindex], afberr[binindex])
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
            (afb,afberr) = GetCorrectedAfb_integratewidth_V(covMCStat, nbins, bin_nominals_i, binwidth)
            print "%20s uncertainty: %2.6f" % ("MC stat", afberr)
            #(afb,afberr) = GetCorrectedAfb_integratewidth_V(covDataStat, nbins, bin_nominals_i, binwidth)
            #print "%20s uncertainty: %2.6f" % ("Data stat", afberr)
        else:
            (afb,afberr_preMCStat,afbcov) = GetCorrectedAfb2D(covar_total, nbins, nbins2D, bin_nominals_i)
            if abs(afberr_preMCStat[nbins2D]-math.sqrt(sumsq_total))>0.000002: print "WARNING: inconsistent covariance matrix. DeltaAfberr = %2.4g" % ( afberr_preMCStat[nbins2D]-math.sqrt(sumsq_total) )
            #print "%20s uncertainty: %2.6f" % ("Sum systs", afberr_preMCStat[nbins2D])
            #for binindex in range(nbins2D): print "%25s %5s: %2.6f" % ("Sum systs", sorted(binlist2D)[binindex], afberr_preMCStat[binindex])
            (afb,afberr,afbcov) = GetCorrectedAfb2D(covMCStat, nbins, nbins2D, bin_nominals_i)
            print "%20s uncertainty: %2.6f" % ("MC stat", afberr[nbins2D])
            for binindex in range(nbins2D): print "%25s %5s: %2.6f" % ("MC stat", sorted(binlist2D)[binindex], afberr[binindex])
            #(afb,afberr,afbcov) = GetCorrectedAfb2D(covDataStat, nbins, nbins2D, bin_nominals_i)
            #print "%20s uncertainty: %2.6f" % ("Data stat", afberr[nbins2D])
            #for binindex in range(nbins2D): print "%25s %5s: %2.6f" % ("Data stat", sorted(binlist2D)[binindex], afberr[binindex])


        #save bin systs for unfolded plots (unfolding macro automatically adds MC stat uncertainty so don't add it here)
        for row in range(nbins): binsyst_preMCStat[row] = sqrt(covar_total[row,row])
        #add stat uncertainties from TUnfold
        covar_total+=covMCStat
        covar_total_incDataStat = covar_total + covDataStat


        if nbins2D==0:
            (afb,afberrdatastat) = GetCorrectedAfb_integratewidth_V(covDataStat, nbins, bin_nominals_i, binwidth)
            (afb,afberrtotal) = GetCorrectedAfb_integratewidth_V(covar_total_incDataStat, nbins, bin_nominals_i, binwidth)
            (afb,afberrsyst) = GetCorrectedAfb_integratewidth_V(covar_total, nbins, bin_nominals_i, binwidth)
            print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)  (%2.6f total)" % (plot, afb, afberrdatastat, afberrsyst, afberrtotal)
        else:
            (afb,afberrdatastat,afbcovdatastat) = GetCorrectedAfb2D(covDataStat, nbins, nbins2D, bin_nominals_i)
            (afb,afberrtotal,afbcovtotalincDataStat) = GetCorrectedAfb2D(covar_total_incDataStat, nbins, nbins2D, bin_nominals_i)
            (afb,afberrsyst,afbcov) = GetCorrectedAfb2D(covar_total, nbins, nbins2D, bin_nominals_i)   #the 2D bin values represent normalised number of events and don't need to be multiplied by the bin widths
            print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)  (%2.6f total)" % (plot, afb[nbins2D], afberrdatastat[nbins2D], afberrsyst[nbins2D], afberrtotal[nbins2D])

        #print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)" % (plot, systematics[plot]['Nominal']['default']['nominal']['Unfolded'][0], systematics[plot]['Nominal']['default']['nominal']['Unfolded'][1], math.sqrt(sumsq_total))

        binindex = 0
        for i in sorted(binlist2D):
            #print "%s %s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)" % (plot, i, systematics[plot]['Nominal']['default']['nominal'][i][0], systematics[plot]['Nominal']['default']['nominal'][i][1], math.sqrt(sumsq_total_2D[binindex]))
            print "%s %s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)  (%2.6f total)" % (plot, i, afb[binindex], afberrdatastat[binindex], afberrsyst[binindex], afberrtotal[binindex])
            binindex += 1
        print ""

        for i in range (nbinstrue):
            print "%s %s = %2.6f +/- %2.6f (stat)" % (plot, binlisttrue[i], denominator[i], denominator_err[i])
        print ""



        #Now calculate correlation matrix
        for row in range(nbins):
            for col in range(nbins):
                corr_total[row,col] = covar_total[row,col] / math.sqrt( covar_total[row,row] * covar_total[col,col] )

        afbcor = zeros( [nbins2D,nbins2D] )
        afbcordatastat = zeros( [nbins2D,nbins2D] )
        afbcortotalincDataStat = zeros( [nbins2D,nbins2D] )

        for row in range(nbins2D):
            for col in range(nbins2D):
                afbcor[row,col] = afbcov[row,col] / math.sqrt( afbcov[row,row] * afbcov[col,col] )
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
                #print afbcor
                #print ""
        print "code fragment to paste in AfbFinalUnfold.h:"
        if nbins2D==0: 
            for row in range(nbins): print "syst_corr[%i] = %2.6f;" % (row, binsyst_preMCStat[row])
        else:
            for binindex in range(nbins2D): print "syst_corr[%i] = %2.6f;" % (binindex, afberr_preMCStat[binindex] )
        print ""
    #end loop over asymmetries
    

if __name__ == '__main__':
    main()
