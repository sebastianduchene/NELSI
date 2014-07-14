import os, sys, subprocess, re, copy

tr_names_old = subprocess.check_output("ls | awk '/sc_/'", shell = True)
tr_names_old = re.split('\n', tr_names_old)
tr_names_old.pop(-1)
tr_names_new = copy.copy(tr_names_old) 
tr_names_new = [re.sub('sc_', '', i) for i in tr_names_new]
tr_names_new = [re.sub('[a-z]tr', 'e', i) for i in tr_names_new]
for i in range(len(tr_names_old)):
    print "mv "+tr_names_old[i]+" "+tr_names_new[i]
    comm_out = subprocess.check_output("mv "+tr_names_old[i]+" "+tr_names_new[i], shell = True)
    print comm_out
