ctypedef void (*valuefreefunc)(void *)

cdef extern from "Hashmap.h":

        ctypedef struct Hashmap:
                pass

        Hashmap *hashmap_new(int size)
        int hashmap_add(Hashmap *hash, char *key, void *value)
        void hashmap_delete(Hashmap *hash, valuefreefunc)

cdef extern from "pairedtalesf.h":
        int run_paired_talesf_task(Hashmap *kwargs)
        
def ScorePairedTalesfTask(char *seqfilename, char *rvd_string, char *rvd_string2, char *output_filepath, char *log_filepath, int c_upstream, double cutoff, int spacer_min, int spacer_max, int numprocs, char *organism_name):

        cdef Hashmap *paired_talesf_kwargs = hashmap_new(32)
        
        cdef double weight = 0.9
        
        hashmap_add(paired_talesf_kwargs, "seq_filename", seqfilename)
        hashmap_add(paired_talesf_kwargs, "rvd_string", rvd_string)
        hashmap_add(paired_talesf_kwargs, "rvd_string2", rvd_string2)
        hashmap_add(paired_talesf_kwargs, "output_filepath", output_filepath)
        hashmap_add(paired_talesf_kwargs, "log_filepath", log_filepath)
        hashmap_add(paired_talesf_kwargs, "weight", &weight)
        hashmap_add(paired_talesf_kwargs, "cutoff", &cutoff)
        hashmap_add(paired_talesf_kwargs, "c_upstream", &c_upstream)
        hashmap_add(paired_talesf_kwargs, "spacer_min", &spacer_min)
        hashmap_add(paired_talesf_kwargs, "spacer_max", &spacer_max)
        hashmap_add(paired_talesf_kwargs, "num_procs", &numprocs)
        hashmap_add(paired_talesf_kwargs, "organism_name", organism_name)
        
        # add a way to pass in an error string
        
        cdef int task_result = run_paired_talesf_task(paired_talesf_kwargs)
        
        hashmap_delete(paired_talesf_kwargs, NULL)
        
        return task_result
