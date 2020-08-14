# Author: Sudheer Ganisetti
# Date: Fr 14. Aug 12:50:34 CEST 2020

import subprocess
import sys
import os
from time import sleep
from datetime import datetime

if os.geteuid() == 0:
    print("\nWe're root!\n")
else:
    #print("We're not root.")
    subprocess.call(['sudo', 'python3', *sys.argv])
    sys.exit()

print("\n\n#######################################################################")
print("#                                                                     #")
print("#                      Author: Sudheer Ganisetti                      #")
print("#  You have submitted this job for checking your internet connection  #")
print("#  you can find the log file at                                       #")
print("#       /var/log/checking_internet_connection_by_sudheer              #")
print("#  If the internet is up, the log is entered with 1 hour frequency    #")
print("#  By default, it can check a maximum of 5 days                       #")
print("#  If it is down, the last few minutes data is entered in detail      #")
print("#  also, it will attempt to check 5 times (total of 1 hour)           #")
print("#                                                                     #")
print("#######################################################################\n\n")

reconnections=0
now=datetime.now()
timestamp=now.timestamp()
date_and_time1=now.ctime()
log_file_path="/var/log/checking_internet_connection_by_sudheer/internet_connection_stability_checking_with_ping_"+str(timestamp)+".log"
output1=open(log_file_path,'w')
output1.write("# hours_and_minutes \t time_for_receiving_the_last_64bits_data \t time_and_date")
while reconnections < 5:
  hours=0
  while hours < 24:
    minutes=0
    every_our_data={}
    every_our_date={}
    while minutes < 360: # each loop takes roughly 10 sec, i.e, minutes loop changes 360*10 sec = 3600/60 min = 60 minutes = roughly 1 hour, however it depends on the network speed
      try:
        ping              = subprocess.check_output(['ping', '-c', '1', '-D', '8.8.8.8'])
        sleep(9)
        internet_status   = "up"
      except subprocess.CalledProcessError:
        internet_status   = "down"
        now=datetime.now()
        timestamp=now.timestamp()
        date_and_time1=now.ctime()
      if internet_status == "up": 
        ping              = ping.decode("utf-8").split(' ')
        timestamp         = float(ping[6].split()[1].split("[")[1].split("]")[0])
        performance       = float(ping[13].split("time=")[1])

        date_and_time1    = datetime.fromtimestamp(timestamp).ctime()
        date_and_time2    = date_and_time1.split(" ")[3].split(":")
        hours_and_minutes = float(date_and_time2[0])+0.01*float(date_and_time2[1])
        temp1={hours_and_minutes:performance}
        every_our_data.update(temp1)
        temp1={hours_and_minutes:date_and_time1}
        every_our_date.update(temp1)
        minutes = minutes + 1
      else:
        break
    if internet_status == "up":
      output1.write("%.2lf\t%.2lf\t%s\n" %(hours_and_minutes,performance,date_and_time1))
      hours=hours+1
    else:
      for i in every_our_data.keys():
        output1.write("%.2lf\t%.2lf\t%s\n" %(i,every_our_data[i],every_our_date[i]))
      output1.write("\n# !WARNING! internet is down at %s\n\n" %(date_and_time1))
      sleep(600)
      break
  reconnections=reconnections+1
output1.write("\n\n# I am done with checking so I am exiting now\n")
output1.close()
  
