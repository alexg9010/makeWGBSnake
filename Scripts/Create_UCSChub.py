#!/usr/bin/env python

import sys,os,json

DIR_bigwig=sys.argv[1]
DIR_ucschub=sys.argv[2]
assembly=sys.argv[3]
email=sys.argv[4]
SAMPLES_TREATMENT=json.loads(sys.argv[5])

# create {assembly} directory
os.makedirs(DIR_ucschub + assembly, exist_ok=True)

# # copy BigWig files into the UCSC_hub/{assembly} directory
bw_files = [DIR_bigwig+sample+"/"+sample+".bw" for sample in SAMPLES_TREATMENT.keys()]
for bw in bw_files:
  os.system('cp '+bw+" "+DIR_ucschub+assembly+"/")  



# tutorial: https://genome.ucsc.edu/goldenPath/help/hubQuickStart.html
# # genomes.txt
f=open(DIR_ucschub+"genomes.txt","w+")
f.write("genome "+ assembly+"\n")
f.write("trackDb "+ assembly+"/trackDb.txt"+"\n")
f.close() 

# hub.txt
f=open(DIR_ucschub+"hub.txt","w+")
hubtxt_content = f"""hub WGBS_{assembly}
shortLabel WGBS ({assembly}) 
longLabel WGBS ({assembly})
genomesFile genomes.txt
email {email}
descriptionUrl hub.html
"""
f.write(hubtxt_content+"\n")
f.close() 


# {assembly}/trackDB.txt
f=open(DIR_ucschub+ assembly+"/trackDb.txt","w+")
# write a separate track for each treatment group and 
# mark it with a different color
def getKey(dict,value):
     return [key for key in dict.keys() if (dict[key] == value)]
treatments=SAMPLES_TREATMENT.values()
for treat in treatments:
  trackDBtxt_content = f"""track {treat}
  superTrack on show
  group {treat}
  shortLabel {treat}
  longLabel {treat}
  """
  f.write(trackDBtxt_content)
  samples_treat=getKey(SAMPLES_TREATMENT,treat)
  for samp in samples_treat:
    trackDBtxt_content = f"""track {treat} {samp}
    parent {treat}
    bigDataUrl {samp}.bw
    shortLabel {treat}_{samp}
    longLabel {treat}
    visibility full
    type bigWig 0 1
    color 77,175,74
    maxHeightPixels 40:20:8
    viewLimitsMax 0:1
    """
    f.write(trackDBtxt_content)
f.close()