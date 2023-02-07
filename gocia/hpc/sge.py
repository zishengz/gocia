import os, sys

class SGE:
    def __init__(
        self,
        username=None,
        qinfoFile = 'q.info'
        ):
        if username is not None:
            self.username = username
        else:
            self.username = os.popen('whoami').read()
        self.mainDir = os.getcwd() 
        self.jidList = []
        self.qinfoFile = qinfoFile

    def get_allJobID(self):
        rawInfo = os.popen('qstat -u %s'%(self.username)).readlines()[2:]
        if len(rawInfo) > 1:
            rawInfo = [l.split()[0] for l in rawInfo]
        return rawInfo

    def get_run(self):
        myRun = []
        for j in self.jidList:
            myInfo = self.jobInfo(j)
            if myInfo['state'] == 'r':
                myRun.append(j)
        return myRun

    def get_wait(self):
        myWait = []
        for j in self.jidList:
            myInfo = self.jobInfo(j)
            if myInfo['state'] == 'q':
                myWait.append(j)
        return myWait

    def update(self):
        if len(self.jidList) > 0:
            self.jidList = [
                j for j in self.jidList\
                if j in self.get_allJobID()
            ]

    def write(self):
        with open(self.qinfoFile, 'w') as f:
            f.write(str(self.jidList))

    def read(self):
        if self.qinfoFile in os.listdir('.'):
            self.jidList = eval(open(self.qinfoFile).read())

    def __len__(self):
        return len(self.jidList)

    def submit(self, jobCommand=''):
        subOut = os.popen(jobCommand).read()
        if 'has been submitted' in subOut:
            myID = subOut.split()[2]
            self.add(myID)

    def add(self, addID):
        if type(addID) is int:
            addID = str(addID)
        if addID not in self.jidList:
            self.jidList.append(addID)

    def kill(self, killID):
        if killID in self.jidList:
            os.system('qdel %s'%killID)
            self.jidList.remove(killID)

    def jobInfo(self, jobid):
        rawInfo = os.popen('qstat -j %s'%jobid).readlines()[1:-1]
        jobInfo = {
            'state': 'd',
            'name': '',
            'wkdir': '',
        }
        if len(rawInfo) > 0:
            # update job state
            tmp = [l for l in rawInfo if 'job_state' in l]
            if len(tmp) > 0:
                jobInfo.update({'state': 'r'})
            else:
                jobInfo.update({'state': 'q'})
            # update job rundir
            tmp = [l for l in rawInfo if 'cwd' in l]
            jobInfo.update({'wkdir': tmp[0].split()[1]})
            # update job name
            tmp = [l for l in rawInfo if 'job_name' in l]
            jobInfo.update({'name': tmp[0].split()[1]})
        return jobInfo

    



