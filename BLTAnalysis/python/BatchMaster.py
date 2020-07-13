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
    def __init__(self, dataname, path, suffix, n_jobs=1):
        self._dataname  = dataname
        self._path      = path
        self._suffix    = suffix
        self._n_jobs    = n_jobs

class BatchMaster():
    '''A tool for submitting batch jobs'''
    def __init__(self, config_list, stage_dir, output_dir, selection, period, executable='execBatch.sh', location='lpc'):
        self._current     = os.path.abspath('.')
        self._stage_dir   = stage_dir
        self._output_dir  = output_dir
        self._config_list = config_list
        self._selection   = selection
        self._period      = period
        self._executable  = executable
        self._location    = location
    
    def split_jobs_by_dataset(self, directory, n_jobs):

        file_list = []
        if directory[-1] == '*':
            directory = directory[:-1]
            subdirs = os.listdir(directory)
            for subdir in subdirs:
                new_path = '{0}/{1}'.format(directory, subdir)
                file_names = ['{0}/{1}'.format(subdir, fn) for fn in os.listdir(new_path)]
                file_list.extend(file_names)
        else:
            file_list = os.listdir(directory)
        n_files = len(file_list)
        
        # Split files to requested number.  Cannot exceed
        # the number of files being run over.
        if n_jobs > n_files:
            n_jobs = n_files

        n_files_per_job = int(math.ceil(float(n_files)/float(n_jobs)))
        file_split = [file_list[i:i+n_files_per_job] for i in range(0, len(file_list), n_files_per_job)]

        return file_split

    def make_batch_lpc(self, cfg, sources):
        '''
        Prepares for submission to lpc.  Does the following:

        1. Generates input_files.txt with files to run over
        2. Write batch configuration file
        '''

        if self._location == 'lpc':
            output_dir = 'root://cmseos.fnal.gov/' + self._output_dir
        else:
            output_dir = self._output_dir

        ## Writing the batch config file
        batch_tmp = open('.batch_tmp_{0}'.format(cfg._dataname,), 'w')
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
            batch_tmp.write('request_memory        = 4096\n')
        batch_tmp.write('\n')

        for i, source in enumerate(sources):

            ## make file with list of inputs ntuples for the analyzer
            input_file = open('input_{1}_{2}.txt'.format(self._stage_dir, cfg._dataname, str(i+1)), 'w')
            if self._location == 'lpc':
                path = 'root://cmseos.fnal.gov/' + cfg._path[10:]
            else:
                path = cfg._path

            if path[-1] == '*':
                path = path[:-2]

            for filename in source:
                input_file.write('{0}/{1}\n'.format(path, filename))
            input_file.close()

            ### set output directory

            batch_tmp.write('Arguments             = {0} {1} {2} {3} {4} {5}\n'.format(cfg._dataname, i+1, cfg._suffix, self._selection, self._period, output_dir))
            batch_tmp.write('Executable            = {0}\n'.format(self._executable))
            batch_tmp.write('Transfer_Input_Files  = source.tar.gz, input_{0}_{1}.txt\n'.format(cfg._dataname, i+1))
            batch_tmp.write('Output                = reports/{0}_{1}_$(Cluster)_$(Process).stdout\n'.format(cfg._dataname, i+1))
            batch_tmp.write('Error                 = reports/{0}_{1}_$(Cluster)_$(Process).stderr\n'.format(cfg._dataname, i+1))
            batch_tmp.write('Log                   = reports/{0}_{1}_$(Cluster)_$(Process).log   \n'.format(cfg._dataname, i+1))
            batch_tmp.write('Queue\n\n')

        batch_tmp.close()
        

    def submit_to_batch(self):
        '''
        Submits batch jobs to scheduler.  Currently only works
        for condor-based batch systems.
        '''

        print 'Running on {0}'.format(self._location)
        print 'Setting up stage directory...'
        self._stage_dir  = '{0}/{1}_{2}_{3}'.format(self._stage_dir, self._selection, self._period, get_current_time())
        make_directory(self._stage_dir, clear=False)

        print 'Setting up output directory...'
        self._output_dir  = '{0}/{1}_{2}_{3}'.format(self._output_dir, self._selection, self._period, get_current_time())
        make_directory('/eos/uscms/' + self._output_dir, clear=False)

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
                print cfg._dataname
                sourceFiles = self.split_jobs_by_dataset(cfg._path, cfg._n_jobs)
                self.make_batch_lpc(cfg, sourceFiles)
                subprocess.call('condor_submit .batch_tmp_{0}'.format(cfg._dataname), shell=True)

