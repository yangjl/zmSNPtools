import uuid
from IPython.parallel.apps import launcher
from IPython.utils import traitlets
from IPython.utils.traitlets import (List, Unicode, CRegExp)

DEFAULT_MEM_PER_CPU = 4000 # Mb

class SLURMLauncher(launcher.BatchSystemLauncher):
    """A BatchSystemLauncher subclass for SLURM
"""
submit_command = List(['sbatch'], config=True,
                      help="The SLURM submit command ['sbatch']")
# Send SIGKILL instead of term, otherwise the job is "CANCELLED", not
# "FINISHED"
delete_command = List(['scancel', '--signal=KILL'], config=True,
                      help="The SLURM delete command ['scancel']")
job_id_regexp = CRegExp(r'\d+', config=True,
                        help="A regular expression used to get the job id from the output of 'sbatch'")

batch_file = Unicode(u'', config=True,
                     help="The string that is the batch script template itself.")

queue_regexp = CRegExp('#SBATCH\W+-p\W+\w')
queue_template = Unicode('#SBATCH -p {queue}')

class SLURMEngineSetLauncher(SLURMLauncher, launcher.BatchClusterAppMixin):
    """Custom launcher handling heterogeneous clusters on SLURM
"""
batch_file_name = Unicode("SLURM_engine" + str(uuid.uuid4()))
mem = traitlets.Unicode("", config=True)
timelimit = traitlets.Unicode("", config=True)
resources = traitlets.Unicode("", config=True)
default_template = traitlets.Unicode("""#!/bin/sh
                                     #SBATCH --partition=defq
                                     #SBATCH --output=/dev/null
                                     #SBATCH --ntasks={n}
                                     #SBATCH --cpus-per-task={1}
                                     #SBATCH --time={timelimit}
                                     #SBATCH --export=ALL
                                     {mem}
                                     {resources}
                                     cd {profile_dir}
                                     mpirun -n {n} ipengine --profile-dir={profile_dir}
                                     """)

def start(self, n):
    if self.mem:
        self.context["mem"] = "#SBATCH --mem-per-cpu=%s" % int(float(self.mem))
    else:
        self.context["mem"] = "#SBATCH --mem-per-cpu=%d" % int(DEFAULT_MEM_PER_CPU)
    self.context["timelimit"] = self.timelimit
    self.context["resources"] = "\n".join(["#SBATCH --%s" % r.strip()
                                       for r in str(self.resources).split(";")
                                       if r.strip()])
    return super(SLURMEngineSetLauncher, self).start(n)
