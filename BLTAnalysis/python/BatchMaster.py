import sys, os, subprocess, fileinput, math, datetime

def get_current_time():
    now = datetime.datetime.now()
    currentTime = '{0:02d}{1:02d}{2:02d}_{3:02d}{4:02d}{5:02d}'.format(now.year, now.month, now.day, now.hour, now.minute, now.second)
    return currentTime

def make_directory(filePath, clear = True):
    if not os.path.exists(filePath):
        os.system('mkdir -p '+filePath)
    if clear and len(os.listdir(filePath)) != 0:
        os.system('rm '+filePath+'/*')

class JobConfig():
    '''Class for storing configuration for each dataset'''
    def __init__(self, data_name, path, suffix, nJobs=1):
        self._data_name  = data_name
        self._path      = path
        self._suffix    = suffix
        self._nJobs     = nJobs

class BatchMaster():
    '''A tool for submitting batch jobs'''
    def __init__(self, config_list, stage_dir, selection, period, executable='execBatch.sh', location='lpc'):
        self._current     = os.path.abspath('.')
        self._stage_dir   = stage_dir
        self._config_list = config_list
        self._selection   = selection
        self._period      = period
        self._executable  = executable
        self._location    = location
    
    def split_jobs_by_dataset(self, directory, nJobs):
        fileList = os.listdir(directory)
        nFiles = len(fileList)
        
        # Split files to requested number.  Cannot exceed
        # the number of files being run over.
        if nJobs > nFiles:
            nJobs = nFiles

        nFilesPerJob = int(math.ceil(float(nFiles)/float(nJobs)))
        fileSplit = [fileList[i:i+nFilesPerJob] for i in range(0, len(fileList), nFilesPerJob)]

        return fileSplit

    def make_batch_lpc(self, cfg, sources):
        '''
        Prepares for submission to lpc.  Does the following:

        1. Generates input_files.txt with files to run over
        2. Write batch configuration file
        '''

        ## Writing the batch config file
        batch_tmp = open('.batch_tmp_{0}'.format(cfg._data_name,), 'w')
        batch_tmp.write('Universe              = vanilla\n')
        batch_tmp.write('Should_Transfer_Files = YES\n')
        batch_tmp.write('WhenToTransferOutput  = ON_EXIT\n')
        batch_tmp.write('Notification          = Never\n')
        #batch_tmp.write('notify_user           = brian.lee.pollack@cern.ch\n')
        if self._location == 'nut3':
            batch_tmp.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
        elif self._location == 'lpc':
            batch_tmp.write('Requirements          = OpSys == "LINUX"&& (Arch != "DUMMY" )\n')
            batch_tmp.write('request_disk          = 2000000\n')
            batch_tmp.write('request_memory        = 2048\n')
        batch_tmp.write('\n')

        for i, source in enumerate(sources):

            ## make file with list of inputs ntuples for the analyzer
            input = open('input_{1}_{2}.txt'.format(self._stage_dir, cfg._data_name, str(i+1)), 'w')
            path = cfg._path
            if self._location == 'lpc':
                path = 'root://cmseos.fnal.gov/' + cfg._path[10:]

            for file in source:
                input.write(path+'/'+file+'\n')
            input.close()

            batch_tmp.write('Arguments             = {0} {1} {2} {3} {4}\n'.format(cfg._data_name, i+1, cfg._suffix, self._selection, self._period))
            batch_tmp.write('Executable            = {0}\n'.format(self._executable))
            batch_tmp.write('Transfer_Input_Files  = source.tar.gz, input_{0}_{1}.txt\n'.format(cfg._data_name, i+1))
            batch_tmp.write('Output                = reports/{0}_{1}_$(Cluster)_$(Process).stdout\n'.format(cfg._data_name, i+1))
            batch_tmp.write('Error                 = reports/{0}_{1}_$(Cluster)_$(Process).stderr\n'.format(cfg._data_name, i+1))
            batch_tmp.write('Log                   = reports/{0}_{1}_$(Cluster)_$(Process).log   \n'.format(cfg._data_name, i+1))
            batch_tmp.write('Queue\n\n')

        batch_tmp.close()
        

    def submit_to_batch(self):
        '''
        Submits batch jobs to batch.  Currently only works
        for lpc batch system, but should be updated for more 
        general use
        '''

        print 'Running on {0}'.format(self._location)
        print 'Setting up stage directory...'
        self._stage_dir  = '{0}/{1}_{2}_{3}'.format(self._stage_dir, self._selection, self._period, get_current_time())
        make_directory(self._stage_dir, clear=False)

        print 'Creating tarball of current workspace in {0}'.format(self._stage_dir)
        #os.system('tar czf {0}/source.tar.gz -C $CMSSW_BASE/src . 2> /dev/null'.format(self._stage_dir))
        if os.getenv('CMSSW_BASE') == '':
            print 'You must source the CMSSW environment you are working in...'
            exit()
        else:
            cmssw_version = os.getenv('CMSSW_BASE').split('/')[-1]
            os.system('tar czf {0}/source.tar.gz -C $CMSSW_BASE/.. {1}'.format(self._stage_dir, cmssw_version))

        subprocess.call('cp {0} {1}'.format(self._executable, self._stage_dir), shell=True)
        os.chdir(self._stage_dir)
        make_directory('reports', clear=False)
        
        print 'Ready to submit to batch system {0}!'.format(self._location)
        if self._location in ['lpc', 'nut3']:
            for cfg in self._config_list:
                print cfg._data_name
                sourceFiles = self.split_jobs_by_dataset(cfg._path, cfg._nJobs)
                self.make_batch_lpc(cfg, sourceFiles)
                subprocess.call('condor_submit .batch_tmp_{0}'.format(cfg._data_name), shell=True)
