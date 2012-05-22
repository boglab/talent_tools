cdef extern from "talesf.h":
        int run_talesf_task(char *seqfilename, char* rvdstring, char *output_filepath, char *log_filepath, double weight, double cutoff, int forwardonly, int c_upstream, int numprocs, char *organism_name)
        
def ScoreTalesfTask(char *seqfilename, char* rvdstring, char *output_filepath, char *log_filepath, bint forwardonly, int c_upstream, double cutoff, int numprocs, char *organism_name):
        return run_talesf_task(seqfilename, rvdstring, output_filepath, log_filepath, 0.9, cutoff, forwardonly, c_upstream, numprocs, organism_name)
