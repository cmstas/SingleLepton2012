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
        
    results = {}

    # example: use JES up and down
    for plot in systematics.keys():
        if plot not in results.keys(): results[plot] = {}
        if 'Default' not in systematics[plot].keys(): 
            print 'Could not find Default values'
            sys.exit(1)
        if 'JESup' not in systematics[plot].keys(): 
            print 'Could not find JESup values'
            sys.exit(1)
        if 'JESdown' not in systematics[plot].keys(): 
            print 'Could not find JESdown values'
            sys.exit(1)
        print "%10s: % 2.6f +/- %2.6f" % ('Default',systematics[plot]['Default']['Unfolded'][0],systematics[plot]['Default']['Unfolded'][1])
        print "%10s: % 2.6f +/- %2.6f" % ('JESup',systematics[plot]['JESup']['Unfolded'][0],systematics[plot]['JESup']['Unfolded'][1])
        print "%10s: % 2.6f +/- %2.6f" % ('JESdown',systematics[plot]['JESdown']['Unfolded'][0],systematics[plot]['JESdown']['Unfolded'][1])
        print ""
        max_difference = max(abs(systematics[plot]['Default']['Unfolded'][0]-systematics[plot]['JESup']['Unfolded'][0]),abs(systematics[plot]['Default']['Unfolded'][0]-systematics[plot]['JESdown']['Unfolded'][0]))
        
        print "% 2.6f +/- %2.6f +/- %2.6f" % (systematics[plot]['Default']['Unfolded'][0],systematics[plot]['Default']['Unfolded'][1],max_difference)    
    

    

if __name__ == '__main__':
    main()
