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

    for plot in systematics.keys():

        #if plot not in results.keys(): results[plot] = {}

        print plot

        if 'Default' not in systematics[plot].keys(): 
            print 'Could not find Default values'
            sys.exit(1)

        sumsq_total = 0
        (default_asym, default_stat) = systematics[plot]['Default']['default']['default']['Unfolded']

        for systematic in systematics[plot].keys():
            if systematic == 'Default': continue
            if systematic == 'name': continue

            sumsq_syst = 0

            for subtype in systematics[plot][systematic].keys():
                sumsq_subtype = 0
                maxdiff = 0

                #This block deals with max(varup,vardown)-type systematics. It's only one case out of many more to come.
                #Once more cases are added, the following should be enclosed in an "if flag == updown" statement, or similar.
                if 'up' not in systematics[plot][systematic][subtype].keys(): 
                    print 'Could not find upward variation'
                    sys.exit(1)
                if 'down' not in systematics[plot][systematic][subtype].keys(): 
                    print 'Could not find downward variation'
                    sys.exit(1)

                updiff =   abs( default_asym - systematics[plot][systematic][subtype]['up']['Unfolded'][0] )
                downdiff = abs( default_asym - systematics[plot][systematic][subtype]['down']['Unfolded'][0] )
                maxdiff =  max( updiff, downdiff )

                sumsq_subtype += maxdiff*maxdiff
                sumsq_syst    += maxdiff*maxdiff
                sumsq_total   += maxdiff*maxdiff

            print "%15s systematic: %2.6f" % (systematic, math.sqrt(sumsq_syst))

        print "%s = %2.6f +/- %2.6f (stat) +/- %2.6f (syst)" % (plot, default_asym, default_stat, math.sqrt(sumsq_total))
        print ""
    

if __name__ == '__main__':
    main()
