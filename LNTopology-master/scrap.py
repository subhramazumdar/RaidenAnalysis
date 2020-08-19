import requests, json
import pprint
import time
import sys
from datetime import datetime

n = 1 #First token network (DAI)

while True:
    url = "https://explorer.raiden.network/json"
    response = requests.get(url)
    data = json.loads(response.text)
    l1 = []
    s = "file"
    now = datetime.now()
    dt_string = now.strftime("%d_%m_%Y_%H:%M:%S")
    l1.append(sys.argv[1]+s+dt_string+".json")
    f = open(l1[0],"w+")
    #str1=pprint.pformat(data,indent=4)
    channelsnum = len(data['networks'][n]['channels'])
    i = 0
    while True:
        i += 1
        if data['networks'][n]['channels'][i-1]['status'] != 'opened':
            data['networks'][n]['channels'].remove(data['networks'][n]['channels'][i-1])
            i -= 1
            channelsnum -= 1
        if i >= channelsnum:
            break
    print(data['networks'][n]['channels'])
    f.write(json.dumps(data['networks'][n], indent=4, sort_keys=True))
    f.close()
    print(l1)
    time.sleep(12.0 * 60.0 * 60.0)

