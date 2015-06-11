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



def calculateVariation(refNominal, systematic, binname, allowmaxstat = 0):
    """
    Calculates the variation due to a particular systematic/subtype
    """

    numpy.set_printoptions(linewidth=500)
    numpy.set_printoptions(precision=6)

    vartype = systematic['vartype']
    usenotes = systematic['usenotes']

    if 'newnom' in usenotes: nominal = systematic['nominal'][binname]
    else: nominal = refNominal

    #Calculate the contribution from this systematic.
    #Calculation method depends on the vartype flag, determined in parseSystematics.py

    # Max of the 'up' and 'down' variations (potentially with respect to a special nominal value)
    if vartype == 'updown' or vartype == 'updown_newnom':
        upvar   = systematic['up'][binname][0]   - nominal[0]
        downvar = systematic['down'][binname][0] - nominal[0]
        varsign = sign(systematic['up'][binname][0] - systematic['down'][binname][0])
        if abs(upvar) >= abs(downvar): maxvar  = abs(upvar)*varsign
        elif abs(upvar) < abs(downvar): maxvar = abs(downvar)*varsign

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
    sum_n = 0.;
    sum_alpha_n = 0.;
    for i in range(nbins):
        sum_n += n[i]
        sum_alpha_n += alpha[i] * n[i]

    dfdn = zeros( [nbins] )
    for i in range(nbins):
        dfdn[i] = ( alpha[i] * sum_n - sum_alpha_n ) / pow(sum_n,2)

    # Error Calculation
    afberr = 0.;
    for i in range(nbins):
        for j in range(nbins):
            afberr += covarianceM[i,j] * dfdn[i] * dfdn[j];

    afberr = sqrt(afberr);

    # Calculate Afb
    afb = sum_alpha_n / sum_n;

    print "afb: %2.6f , afberr: %2.6f " % (afb, afberr )  
    return afberr






def main():

    usage  = "Usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="verbose output")
    parser.add_option("-f", "--file", action="store", type="string", default=None, dest="systematicsjsonfilename", help="Name of systematics file in json format")
    parser.add_option("-n", "--nomatrix", action="store_true", default=False, dest="nomatrix", help="Suppress printout of covariance and correlation matrices")

    (opts, args) = parser.parse_args()
    
    if ( opts.systematicsjsonfilename == None ) :
        print ""
        print "Please specify systematics file in json format!"
        print ""
        parser.print_help()
        sys.exit(2)

    verbose = opts.verbose
    systematicsjsonfilename = opts.systematicsjsonfilename
    nomatrix = opts.nomatrix

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

        #Count how many bins we're dealing with
        typelist = list( systematics[plot]['Nominal']['default']['nominal'].keys() )
        binnames = fnmatch.filter(typelist, 'bin*')
        binnames.sort()
        nbins = len(binnames)

        if 'x' in binnames[1]:  #2D binning
            binlist = []
            for i in range(1,3+1):
                for j in range(1, (nbins/3)+1):
                    binlist.append('bin' + str(i) + 'x' + str(j))
        else:
            binlist = ['bin'+str(i) for i in range(1, nbins+1)]   #1D binning

        if sorted(binlist) != sorted(binnames):
            print 'Problem: these lists of bin names don\'t match:'
            print sorted(binnames)
            print sorted(binlist)
            sys.exit(1)

        sumsq_total = 0
        covar_total = zeros( [nbins,nbins] )
        corr_total  = zeros( [nbins,nbins] )
        binsyst_total = zeros( [nbins] )
        bin_nominals = {}
        bin_nominals_i = zeros( [nbins] )

        #Get the nominal values for each bin, and overall
        for i in binlist: bin_nominals[i] = systematics[plot]['Nominal']['default']['nominal'][i]
        (nominal_unfolded, stat_unfolded) = systematics[plot]['Nominal']['default']['nominal']['Unfolded'] 
        #print nominal_unfolded

        for i in range(nbins): 
            #print bin_nominals[binlist[i]][0]
            bin_nominals_i[i] = bin_nominals[binlist[i]][0]

        #Loop over the different systematics...
        for systematic in sorted(systematics[plot].keys()):
            if systematic == 'Nominal': continue
            if systematic == 'name': continue

            sumsq_syst = 0
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
                weightedaveragegradient = 0.;
                bin_weightedaveragegradient = zeros( [nbins] )
                bin_variations_0 = {}
                #sort the subtypes into numerical order before doing the extrapolation
                for subtype in sorted(systematics[plot][systematic].keys()):
                    #use only the first 6 points for the extrapolation
                    if subtypenum > 5: continue
                    sumsq_subtype = 0
                    covar_subtype = zeros( [nbins,nbins] )
                    bin_variations = {}

                    #Calculate the variation in the overall asymmetry
                    var_overall = calculateVariation( [nominal_unfolded, stat_unfolded], systematics[plot][systematic][subtype], 'Unfolded' )
                    sumsq_subtype = var_overall*var_overall
                    if subtypenum == 0:
                        var_overall_maxstat = calculateVariation( [nominal_unfolded, stat_unfolded], systematics[plot][systematic][subtype], 'Unfolded' , 1)
                        maxstat_factor = abs(var_overall_maxstat/var_overall)
                        if(maxstat_factor<1): print "something went wrong with maxstat_factor"


                    #Calculate the variation in individual bins
                    for i in binlist: bin_variations[i] = calculateVariation( bin_nominals[i], systematics[plot][systematic][subtype], i )
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
                        covar_syst = covar_subtype
                    else:
                        weightedaveragegradient += (sqrt(sumsq_syst) - sqrt(sumsq_subtype))/subtypenum
                        #print (sqrt(sumsq_syst) - sqrt(sumsq_subtype))/subtypenum
                        #print 15.*weightedaveragegradient/subtypenum
                    subtypenum += 1
                weightedaveragegradient/=subtypenum
                bin_weightedaveragegradient/=subtypenum
                #print 15.*weightedaveragegradient

                #afberr = GetCorrectedAfb(covar_syst, nbins, bin_nominals_i)

                extrapfactor = 1 + 15.*weightedaveragegradient / sqrt(sumsq_syst)
                #take maximum from extrapolation or maxstat_factor due to the MC statistical uncertainty (guaranteed to be >=1)
                print "          %s extrapolation factor: %2.3f , maxstat_factor: %2.3f " % (systematic,extrapfactor,maxstat_factor)
                if(extrapfactor<maxstat_factor):  extrapfactor = maxstat_factor
                if(extrapfactor>3.0):
                    extrapfactor = 3.0  #don't extrapolate more than a factor of 3 to avoid producing highly (anti)correlated covariance matrix
                    print "limiting uncertainty scaling factor to 3.0 to maintain sensible covariance matrix"
                #print "extrapolation factor: %2.3f , extrapolation amount: %2.3f " % (extrapfactor,100.*15.*weightedaveragegradient )
                sumsq_syst *= extrapfactor*extrapfactor
                #instead of extrapolating for every bin, use same SF as for asymmetry (automatically ensures the covariance matrix is consistent with the uncertainty on the asymmetry, so no need to recalculate it)
                covar_syst *= extrapfactor*extrapfactor

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
                    covar_subtype = zeros( [nbins,nbins] )
                    bin_variations = {}

                    #Calculate the variation in the overall asymmetry
                    var_overall = calculateVariation( [nominal_unfolded, stat_unfolded], systematics[plot][systematic][subtype], 'Unfolded' )
                    sumsq_subtype = var_overall*var_overall
                    var_overall_maxstat = calculateVariation( [nominal_unfolded, stat_unfolded], systematics[plot][systematic][subtype], 'Unfolded' , 1)
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

                    #Calculate the variation in individual bins
                    for i in binlist: bin_variations[i] = calculateVariation( bin_nominals[i], systematics[plot][systematic][subtype], i )

                    #Calculate the covariance matrix using the individual bin variations:
                    #covar_ij = var(bin i) * var(bin j)
                    for row in range(nbins):
                        for col in range(nbins):
                            covar_subtype[row, col] = bin_variations[binlist[row]] * bin_variations[binlist[col]]

                    if(maxstat_factor>1): print "          %s maxstat_factor: %2.3f " % (systematic,maxstat_factor)
                    if(maxstat_factor>3.0):
                        maxstat_factor = 3.0  #don't extrapolate more than a factor of 3 to avoid producing highly (anti)correlated covariance matrix
                        print "limiting uncertainty scaling factor to 3.0 to maintain sensible covariance matrix"
                    #Add to the running total of the variance and covariance
                    sumsq_syst += sumsq_subtype*maxstat_factor*maxstat_factor
                    covar_syst += covar_subtype*maxstat_factor*maxstat_factor

            #end loop over subtypes
            print "%15s systematic: %2.6f" % (systematic, math.sqrt(sumsq_syst))
            sumsq_total += sumsq_syst
            covar_total += covar_syst

        #end loop over systematics
        #Now calculate correlation matrix
        for row in range(nbins):
            for col in range(nbins):
                corr_total[row,col] = covar_total[row,col] / math.sqrt( covar_total[row,row] * covar_total[col,col] )

        for row in range(nbins): binsyst_total[row] = sqrt(covar_total[row,row])

        print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)" % (plot, nominal_unfolded, stat_unfolded, math.sqrt(sumsq_total))
        print ""
        if nomatrix == False:
            print "%s covariance matrix:" % plot
            print binlist
            print covar_total
            print ""
            print "%s correlation matrix:" % plot
            print binlist
            print corr_total
            print ""
            #print "%s bin systematics:" % plot
            #print binlist
            #print binsyst_total
            #print ""
        print ""
        print "code fragment to paste in AfbFinalUnfold.h:"
        for row in range(nbins): print "syst_corr[%i] = %2.6f;" % (row, binsyst_total[row])

    #end loop over asymmetries
    

if __name__ == '__main__':
    main()
