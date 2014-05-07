#!/usr/bin/env python

import sys,getopt,urllib2,json
from optparse import OptionParser

def parseSystFiles(systfiles):
    """

    data structure for systematics

    systematics dictionary: key: identifier (plotname)
    |->
       syst dictionary: key: systname
       |->
          type dictionary: key: type/bin name/etc
          |->
             value list: [asymmetry,error]

    """


    systematics = {}

    for syst in systfiles.keys():
        for subtype in systfiles[syst].keys():
            for direction in systfiles[syst][subtype].keys():
                try:
                    systfile = open(systfiles[syst][subtype][direction])
                except:
                    print ""
                    print "Syst file",systfiles[syst],"for systematics",syst,"cannot be opened"
                    print ""
                    sys.exit(1)
                for line in systfile.readlines():
                    array = line.split()
                    if len(array) != 6:
                        print "malformated line in systfile",systfiles[syst]
                        print "line:",line
                        sys,exit(1)
                    identifier = array[0]
                    name = array[1]
                    type = array[2].split(':')[0] #remove ':'
                    asymmetry = float(array[3])
                    error = float(array[5])
                    if identifier not in systematics.keys(): systematics[identifier] = {}
                    systematics[identifier]['name'] = name
                    if syst not in systematics[identifier].keys(): systematics[identifier][syst] = {}
                    if subtype not in systematics[identifier][syst].keys(): systematics[identifier][syst][subtype] = {}
                    if direction not in systematics[identifier][syst][subtype].keys(): systematics[identifier][syst][subtype][direction] = {}
                    systematics[identifier][syst][subtype][direction][type] = [asymmetry,error]
            
    return systematics
            
def main():

    usage  = "Usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="verbose output")
    parser.add_option("-f", "--file", action="store", type="string", default=None, dest="steeringfilename", help="Name of steering file")
    parser.add_option("-o", "--output", action="store", type="string", default="systematics.json", dest="outputfilename", help="Name of output file (default: systematics.json)")

    (opts, args) = parser.parse_args()
    
    if ( opts.steeringfilename == None ) :
        print ""
        print "Please specify steering file!"
        print ""
        parser.print_help()
        sys.exit(2)

    verbose = opts.verbose
    steeringfilename = opts.steeringfilename
    outputfilename = opts.outputfilename

    try:
        steeringfile = open(steeringfilename)
    except:
        print ""
        print "Steering file",steeringfilename,"cannot be opened"
        print ""
        sys.exit(1)
        
    systfiles = {}
    for line in steeringfile.readlines():
        if line.startswith('#'): continue        
        array = line.split()
        if array[0] not in systfiles.keys(): systfiles[array[0]] = {}
        if array[1] not in systfiles[array[0]].keys(): systfiles[array[0]][array[1]] = {}
        if array[2] in systfiles[array[0]][array[1]].keys():
            print 'systematics:',array[0],array[1],array[2],'was entered twice in steering file:',steeringfilename,'Please choose only one entry per systematics, using last entry.'
        systfiles[array[0]][array[1]][array[2]] = array[3]
        
    systematics = parseSystFiles(systfiles)
    
    outputfile = open(outputfilename,'w')
    json.dump(systematics,outputfile)
    outputfile.close
    print "Systematics dictionary parsed from systematics dump files defined in steering file",steeringfilename,'have been written to',outputfilename

    

if __name__ == '__main__':
    main()
