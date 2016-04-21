'''
Creates a new pulse with shot number = sys.argv[1] based on the model tree.

Created on 3/13/2015 by Jason Milhone
'''
import MDSplus as mds
import sys

shotnum = int(sys.argv[1])


try:                                                                    
    fpTree = mds.Tree("fp_raw",shotnum)                             
except mds.TreeException:                                               
    #Tree does not exist yet.  Try creating a new pulse.            
    try:                                                            
        fpModelTree = mds.Tree("fp_raw",-1)                     
        fpModelTree.createPulse(shotnum)                        
    except mds.TreeException as e:                                  
        print "Unable to create Tree"                             
        print e                                                 
        print " "                                               
        #Exit the function                                      
