#!/usr/bin/python

from utilsPoirel import *
import sys


#def getRavasiTFs():

def main(args):

    allTFs = set()

    # only include TFs classified as 'a', 'b', or 'other
    vInfile = './human-tf-vaquerizas-nature-2009.txt'
    validClasses = set(['a', 'b', 'other'])
    vTFs = set([tf.strip() for tfType, tf in readColumns(vInfile, 1,6) if tf!='' and tfType in validClasses])
    allTFs.update(vTFs)

    rInfile = './human-tf-ravasi-cell-2010.txt'
    rTFs = set([tf.strip() for tf in readItemSet(rInfile, 2) if tf!=''])
    allTFs.update(rTFs)
    
    filesTFmap = { vInfile : vTFs, \
                   rInfile : rTFs, \
                 }
    
    overlaps = {}

    for f1 in sorted(filesTFmap.keys()):
        for f2 in sorted(filesTFmap.keys()):
            if f2<f1:
                continue
            if f1==f2:
                overlaps[(f1,f2)] = len(filesTFmap[f1])
            else: # f1<f2
                overlaps[(f1,f2)] = len(filesTFmap[f1].intersection(filesTFmap[f2]))
            
    print '\nTF overlap anlysis:'
    print '\t#\t' + '\t'.join(sorted(filesTFmap.keys()))
    for f1 in sorted(filesTFmap.keys()):
        print '\t%s' %(f1),
        for f2 in sorted(filesTFmap.keys()):
            k = tuple(sorted([f1,f2]))
            print '\t%d' %(overlaps[k]),
        print


    print '\nVaqueerizas: %d' %(len(vTFs))
    print 'Ravasi:      %d' %(len(rTFs))
    print 'Number of total TFs: %d' %(len(allTFs))

    outfile = 'human-tfs.txt'
    f = open(outfile, 'w')
    f.write('#human_transcription_factor\n')
    for tf in sorted(allTFs):
        f.write('%s\n' %(tf))
    f.close()
    print '\nAll TFs written to "%s"' %(outfile)


if __name__=='__main__':
    main(sys.argv)
