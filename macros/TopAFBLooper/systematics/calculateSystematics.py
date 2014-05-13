#!/usr/bin/env python

import sys,getopt,urllib2,json,math
from optparse import OptionParser
            
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

    #Loop over the different asymmetry variables...
    for plot in systematics.keys():

        #if plot not in results.keys(): results[plot] = {}

        print plot

        if 'Nominal' not in systematics[plot].keys(): 
            print 'Could not find Nominal values'
            sys.exit(1)

        sumsq_total = 0
        (default_asym, default_stat) = systematics[plot]['Nominal']['default']['nominal']['Unfolded']

        #Loop over the different systematics...
        for systematic in systematics[plot].keys():
            if systematic == 'Nominal': continue
            if systematic == 'name': continue

            sumsq_syst = 0

            #Loop over the different subtypes...
            for subtype in systematics[plot][systematic].keys():
                sumsq_subtype = 0
                maxvar = 0
                vartype = systematics[plot][systematic][subtype]['vartype']

                if vartype == 'updown_newnom': nominal = systematics[plot][systematic][subtype]['nominal']['Unfolded'][0]
                else: nominal = default_asym

                #Calculate the contribution from this systematic.
                #Calculation method depends on the vartype flag, provided in parseSystematics.py
                if vartype == 'updown' or vartype == 'updown_newnom':
                    upvar =   abs( nominal - systematics[plot][systematic][subtype]['up']['Unfolded'][0] )
                    downvar = abs( nominal - systematics[plot][systematic][subtype]['down']['Unfolded'][0] )
                    maxvar =  max( upvar, downvar )
                elif vartype == 'half':
                    maxvar = abs(systematics[plot][systematic][subtype]['up']['Unfolded'][0] -
                                 systematics[plot][systematic][subtype]['down']['Unfolded'][0] ) /2
                elif vartype == 'diffnom':
                    maxvar = abs( nominal - systematics[plot][systematic][subtype]['diffnom']['Unfolded'][0] )
                else:
                    print ""
                    print "I don't know what to do with this variation type:", vartype
                    print systematic, subtype
                    sys.exit(1)

                #Add this variation to the running total(s), in quadrature
                sumsq_subtype += maxvar*maxvar
                sumsq_syst    += maxvar*maxvar
                sumsq_total   += maxvar*maxvar

            print "%15s systematic: %2.6f" % (systematic, math.sqrt(sumsq_syst))

        print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)" % (plot, default_asym, default_stat, math.sqrt(sumsq_total))
        print ""
    

if __name__ == '__main__':
    main()
