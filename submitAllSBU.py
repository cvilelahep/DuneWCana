#!/usr/bin/python

import glob, os, subprocess

duneWCanaDir = "/storage/shared/cvilela/DuneWC/DuneWCana"
headerDir = duneWCanaDir+"/Headers/"
setupFile = duneWCanaDir+"/setupSBU.sh"
logDir = duneWCanaDir+"/Logs/"
scriptDir = duneWCanaDir+"/Scripts/"

for dirs in [logDir, scriptDir] :
    if not os.path.isdir(dirs) :
        os.makedirs(dirs)

headerFiles = glob.glob(headerDir+"/*.py")

for job in headerFiles :
    jobName = os.path.basename(job)
    jobName = jobName.split(".")[0]
    jobName = jobName.split("_")[1]+"_"+jobName.split("_")[2]

    print jobName

    script = open(scriptDir+"/"+jobName+".sh", 'w')
    script.write("#!/bin/bash\n")
    script.write("echo "+jobName+" starts\n")
    script.write("date\n")
    script.write("source "+setupFile+"\n")
    script.write("python "+duneWCanaDir+"/DuneWCana.py "+job+"\n")
    script.write("echo "+jobName+" ends\n")
    script.write("date\n")
    script.close()

    subprocess.call(["qsub", scriptDir+"/"+jobName+".sh", "-e", logDir+"/"+jobName+".stderr", "-o", logDir+"/"+jobName+".stdout"])
    
