import sys, os

#_____________________________________________________________________________
def main(argv):
    if not (len(argv) >= 1):
        print "usage: python writeDatacards_WH.py <input.txt> <output_dir>"
        exit()
    print "input file:",argv[0]

    fin = open(argv[0])
    output_dir = 'datacards'
    if len(argv) > 1:
        output_dir = argv[1]
    for line in fin:
        # input file format: mchargino mlsp best_sr nsig staterr
        tokens = line.split()
        filename = "%s/datacard_Wino_DOW_%s_%s.txt" % (output_dir,tokens[0],tokens[1])
        best_sr = int(tokens[2])
        nsig = float(tokens[3])
        staterr = float(tokens[4])
        if best_sr == 0: 
            write_sr1(filename, nsig, staterr)
        elif best_sr == 1:
            write_sr2(filename, nsig, staterr)
        elif best_sr == 2:
            write_sr3(filename, nsig, staterr)
        elif best_sr == 3:
            write_sr4(filename, nsig, staterr)

    print "outputs written to:",output_dir

#_________________________________________________
# MET100 signal region
def write_sr1(filename, nsig, errsig):
    f = open(filename,'w')
    f.write("""imax 1  number of channels
jmax 6  number of backgrounds
kmax 16  number of nuisance parameters (sources of systematical uncertainties)
------------
bin         1
observation 7
------------
bin             1      1       1     1     1     1      1
process       signal  dilep  top1l  wzbb   wbb  wlight  rare
process         0      1       2      3     4      5    6
rate           %.1f    2.8     1.8   0.6   1.5    0.5    0.4 
------------			               
LUMI          lnN   1.026   -    -    -     -    -    -    lumi uncertainty
JES           lnN   1.05    -    -    -     -    -    -    jet energy scale
DILEPTRIG     lnN   1.05    -    -    -     -    -    -    (single) lepton trigger
LEPID         lnN   1.05    -    -    -     -    -    -    lep id and iso
BTAG          lnN   1.05    -    -    -     -    -    -    btagging eff
ISR           lnN   1.01    -    -    -     -    -    -    isr
DOWJER        lnN   1.05    -    -    -     -    -    -    jet energy resolution
FASTSIM       lnN   1.01    -    -    -     -    -    -    lep fastsim eff
DOWSIGSTAT    lnN   %.2f    -    -    -     -    -    -    stat uncert on signal
DOWDILEP      lnN     -   1.43   -    -     -    -    -    tot uncert on dilep bg
DOWLEPBSYST   lnN     -     -   1.57 1.36  1.43  -    -    syst uncert on lepb bg
DOWTOPSTAT    lnN     -     -   1.18  -     -    -    -    stat uncert on top1l bg
DOWWZBBSTAT   lnN     -     -    -   1.25   -    -    -    stat uncert on wbb bg
DOWWBBSTAT    lnN     -     -    -    -    1.38  -    -    stat uncert on wbb bg
DOWWLIGHT     lnN     -     -    -    -     -   1.40  -    tot uncert on wlight bg
RARE          lnN     -     -    -    -     -    -   1.50  syst uncert on rare bg
DOWRARESTAT   lnN     -     -    -    -     -    -   1.20  stat uncert on rare bg
""" % (nsig, errsig) )
    f.close()

#_________________________________________________
# MET125 signal region
def write_sr2(filename, nsig, errsig):
    f = open(filename,'w')
    f.write("""imax 1  number of channels
jmax 6  number of backgrounds
kmax 16  number of nuisance parameters (sources of systematical uncertainties)
------------
bin         1
observation 6
------------
bin             1      1       1     1     1     1      1
process       signal  dilep  top1l  wzbb   wbb  wlight  rare
process         0      1       2      3     4      5    6
rate           %.1f    2.3     0.9   0.4   1.0    0.3    0.3 
------------			               
LUMI          lnN   1.026   -    -    -     -    -    -    lumi uncertainty
JES           lnN   1.05    -    -    -     -    -    -    jet energy scale
DILEPTRIG     lnN   1.05    -    -    -     -    -    -    (single) lepton trigger
LEPID         lnN   1.05    -    -    -     -    -    -    lep id and iso
BTAG          lnN   1.05    -    -    -     -    -    -    btagging eff
ISR           lnN   1.01    -    -    -     -    -    -    isr
DOWJER        lnN   1.05    -    -    -     -    -    -    jet energy resolution
FASTSIM       lnN   1.01    -    -    -     -    -    -    lep fastsim eff
DOWSIGSTAT    lnN   %.2f    -    -    -     -    -    -    stat uncert on signal
DOWDILEP      lnN     -   1.43   -    -     -    -    -    tot uncert on dilep bg
DOWLEPBSYST   lnN     -     -   1.57 1.36  1.43  -    -    syst uncert on lepb bg
DOWTOPSTAT    lnN     -     -   1.23  -     -    -    -    stat uncert on top1l bg
DOWWZBBSTAT   lnN     -     -    -   1.29   -    -    -    stat uncert on wbb bg
DOWWBBSTAT    lnN     -     -    -    -    1.46  -    -    stat uncert on wbb bg
DOWWLIGHT     lnN     -     -    -    -     -   1.40  -    tot uncert on wlight bg
RARE          lnN     -     -    -    -     -    -   1.50  syst uncert on rare bg
DOWRARESTAT   lnN     -     -    -    -     -    -   1.24  stat uncert on rare bg
""" % (nsig, errsig) )
    f.close()

#_________________________________________________
# MET150 signal region
def write_sr3(filename, nsig, errsig):
    f = open(filename,'w')
    f.write("""imax 1  number of channels
jmax 6  number of backgrounds
kmax 16  number of nuisance parameters (sources of systematical uncertainties)
------------
bin         1
observation 3
------------
bin             1      1       1     1     1     1      1
process       signal  dilep  top1l  wzbb   wbb  wlight  rare
process         0      1       2      3     4      5    6
rate           %.1f    1.7     0.5   0.3   0.9    0.2    0.3 
------------			               
LUMI          lnN   1.026   -    -    -     -    -    -    lumi uncertainty
JES           lnN   1.05    -    -    -     -    -    -    jet energy scale
DILEPTRIG     lnN   1.05    -    -    -     -    -    -    (single) lepton trigger
LEPID         lnN   1.05    -    -    -     -    -    -    lep id and iso
BTAG          lnN   1.05    -    -    -     -    -    -    btagging eff
ISR           lnN   1.01    -    -    -     -    -    -    isr
DOWJER        lnN   1.05    -    -    -     -    -    -    jet energy resolution
FASTSIM       lnN   1.01    -    -    -     -    -    -    lep fastsim eff
DOWSIGSTAT    lnN   %.2f    -    -    -     -    -    -    stat uncert on signal
DOWDILEP      lnN     -   1.41   -    -     -    -    -    tot uncert on dilep bg
DOWLEPBSYST   lnN     -     -   1.57 1.36  1.43  -    -    syst uncert on lepb bg
DOWTOPSTAT    lnN     -     -   1.32  -     -    -    -    stat uncert on top1l bg
DOWWZBBSTAT   lnN     -     -    -   1.35   -    -    -    stat uncert on wbb bg
DOWWBBSTAT    lnN     -     -    -    -    1.51  -    -    stat uncert on wbb bg
DOWWLIGHT     lnN     -     -    -    -     -   1.40  -    tot uncert on wlight bg
RARE          lnN     -     -    -    -     -    -   1.50  syst uncert on rare bg
DOWRARESTAT   lnN     -     -    -    -     -    -   1.31  stat uncert on rare bg
""" % (nsig, errsig) )
    f.close()

#_________________________________________________
# MET175 signal region
def write_sr4(filename, nsig, errsig):
    f = open(filename,'w')
    f.write("""imax 1  number of channels
jmax 6  number of backgrounds
kmax 16  number of nuisance parameters (sources of systematical uncertainties)
------------
bin         1
observation 3
------------
bin             1      1       1     1     1     1      1
process       signal  dilep  top1l  wzbb   wbb  wlight  rare
process         0      1       2      3     4      5    6
rate           %.1f    1.2     0.2   0.3   0.2    0.2    0.2 
------------			               
LUMI          lnN   1.026   -    -    -     -    -    -    lumi uncertainty
JES           lnN   1.05    -    -    -     -    -    -    jet energy scale
DILEPTRIG     lnN   1.05    -    -    -     -    -    -    (single) lepton trigger
LEPID         lnN   1.05    -    -    -     -    -    -    lep id and iso
BTAG          lnN   1.05    -    -    -     -    -    -    btagging eff
ISR           lnN   1.01    -    -    -     -    -    -    isr
DOWJER        lnN   1.05    -    -    -     -    -    -    jet energy resolution
FASTSIM       lnN   1.01    -    -    -     -    -    -    lep fastsim eff
DOWSIGSTAT    lnN   %.2f    -    -    -     -    -    -    stat uncert on signal
DOWDILEP      lnN     -   1.42   -    -     -    -    -    tot uncert on dilep bg
DOWLEPBSYST   lnN     -     -   1.57 1.36  1.43  -    -    syst uncert on lepb bg
DOWTOPSTAT    lnN     -     -   1.48  -     -    -    -    stat uncert on top1l bg
DOWWZBBSTAT   lnN     -     -    -   1.34   -    -    -    stat uncert on wbb bg
DOWWBBSTAT    lnN     -     -    -    -    2.00  -    -    stat uncert on wbb bg
DOWWLIGHT     lnN     -     -    -    -     -   1.40  -    tot uncert on wlight bg
RARE          lnN     -     -    -    -     -    -   1.50  syst uncert on rare bg
DOWRARESTAT   lnN     -     -    -    -     -    -   1.35  stat uncert on rare bg
""" % (nsig, errsig) )
    f.close()

#_____________________________________________________________________________
if __name__ == '__main__':
    main(sys.argv[1:])
