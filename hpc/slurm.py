import os, sys

class SLURM:
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
        rawInfo = os.popen('squeue -u %s'%(self.username)).readlines()[1:]
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
        subOut = os.popen(jobCommand).readlines()[-1]
        if 'Submitted batch job' in subOut:
            myID = subOut.split()[-1]
            self.add(myID)

    def add(self, addID):
        if type(addID) is int:
            addID = str(addID)
        if addID not in self.jidList:
            self.jidList.append(addID)

    def kill(self, killID):
        if killID in self.jidList:
            os.system('scancel %s'%killID)
            self.jidList.remove(killID)

# TODO: rewrite for SLURM
# this is for perlmutter
    def jobInfo(self, jobid):
        rawInfo = os.popen('squeue -j %s'%jobid).readlines()[1:]
        jobInfo = {
            'state': 'd',
            'name': '',
            'wkdir': '',
        }
        if len(rawInfo) > 0:
            if rawInfo.split()[4] == 'R':
                jobInfo.update({'state': 'r'})
            elif rawInfo.split()[4] =='PD':
                jobInfo.update({'state': 'q'})
            jobInfo.update({'name': rawInfo.split()[2]})
        return jobInfo

    



