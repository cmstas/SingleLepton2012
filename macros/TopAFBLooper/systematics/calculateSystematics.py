#!/usr/bin/env python

import sys,getopt,urllib2,json,math,fnmatch
from optparse import OptionParser
from numpy import *


def calculateVariation(refNominal, systematic, binname):
    """
    Calculates the variation due to a particular systematic/subtype
    """

    vartype = systematic['vartype']

    if vartype == 'updown_newnom': nominal = systematic['nominal'][binname][0]
    else: nominal = refNominal

    #Calculate the contribution from this systematic.
    #Calculation method depends on the vartype flag, determined in parseSystematics.py
    if vartype == 'updown' or vartype == 'updown_newnom':
        upvar   = systematic['up'][binname][0]   - nominal
        downvar = systematic['down'][binname][0] - nominal
        varsign = sign(systematic['up'][binname][0] - systematic['down'][binname][0])
        if abs(upvar) >= abs(downvar): maxvar  = abs(upvar)*varsign
        elif abs(upvar) < abs(downvar): maxvar = abs(downvar)*varsign
    elif vartype == 'half':
        maxvar = ( systematic['up'][binname][0] - systematic['down'][binname][0] ) /2
    elif vartype == 'diffnom':
        maxvar = systematic['diffnom'][binname][0] - nominal
    else:
        print ""
        print "I don't know what to do with this variation type:", vartype
        #print systematic, subtype
        sys.exit(1)

    return maxvar


def main():

    usage  = "Usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="verbose output")
    parser.add_option("-f", "--file", action="store", type="string", default=None, dest="systematicsjsonfilename", help="Name of systematics file in json format")

    (opts, args) = parser.parse_args()
    
    if ( opts.systematicsjsonfilename == None ) :
        print ""
        print "Please specify systematics file in json format!"
        print ""
        parser.print_help()
        sys.exit(2)

    verbose = opts.verbose
    systematicsjsonfilename = opts.systematicsjsonfilename

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
    for plot in systematics.keys():

        #if plot not in results.keys(): results[plot] = {}

        print plot

        if 'Nominal' not in systematics[plot].keys(): 
            print 'Could not find Nominal values'
            sys.exit(1)

        #Count how many bins we're dealing with
        typelist = list( systematics[plot]['Nominal']['default']['nominal'].keys() )
        binlist = fnmatch.filter(typelist, 'bin*')
        binlist.sort()
        nbins = len(binlist)

        sumsq_total = 0
        covar_total = zeros( [nbins,nbins] )
        corr_total  = zeros( [nbins,nbins] )
        bin_nominals = {}

        #Get the nominal values for each bin, and overall
        for i in binlist: bin_nominals[i] = systematics[plot]['Nominal']['default']['nominal'][i][0]
        (nominal_unfolded, stat_unfolded) = systematics[plot]['Nominal']['default']['nominal']['Unfolded'] 

        #Loop over the different systematics...
        for systematic in systematics[plot].keys():
            if systematic == 'Nominal': continue
            if systematic == 'name': continue

            sumsq_syst = 0
            covar_syst = zeros( [nbins,nbins] )

            #Loop over the subtypes within each systematic...
            for subtype in systematics[plot][systematic].keys():
                sumsq_subtype = 0
                covar_subtype = zeros( [nbins,nbins] )
                bin_variations = {}

                #Calculate the variation in the overall asymmetry
                var_overall = calculateVariation( nominal_unfolded, systematics[plot][systematic][subtype], 'Unfolded' )
                sumsq_subtype = var_overall*var_overall

                #Calculate the variation in individual bins
                for i in binlist: bin_variations[i] = calculateVariation( bin_nominals[i], systematics[plot][systematic][subtype], i )

                #Calculate the covariance matrix using the individual bin variations:
                #covar_ij = var(bin i) * var(bin j)
                for row in range(nbins):
                    for col in range(nbins):
                        covar_subtype[row, col] = bin_variations[binlist[row]] * bin_variations[binlist[col]]

                #Add to the running total of the variance and covariance
                sumsq_syst += sumsq_subtype
                covar_syst += covar_subtype

            #end loop over subtypes
            print "%15s systematic: %2.6f" % (systematic, math.sqrt(sumsq_syst))
            sumsq_total += sumsq_syst
            covar_total += covar_syst

        #end loop over systematics
        #Now calculate correlation matrix
        for row in range(nbins):
            for col in range(nbins):
                corr_total[row,col] = covar_total[row,col] / math.sqrt( covar_total[row,row] * covar_total[col,col] )

        print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)" % (plot, nominal_unfolded, stat_unfolded, math.sqrt(sumsq_total))
        print ""
        print "%s covariance matrix:" % plot
        print covar_total
        print ""
        print "%s correlation matrix:" % plot
        print corr_total
        print ""
        print ""

    #end loop over asymmetries
    

if __name__ == '__main__':
    main()
