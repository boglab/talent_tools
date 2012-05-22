ctypedef void (*valuefreefunc)(void *)

cdef extern from "Hashmap.h":

        ctypedef struct Hashmap:
                pass

        Hashmap *hashmap_new(int size)
        int hashmap_add(Hashmap *hash, char *key, void *value)
        void hashmap_delete(Hashmap *hash, valuefreefunc)

cdef extern from "talesf.h":
        int run_talesf_task(Hashmap *kwargs)
        
def ScoreTalesfTask(char *seqfilename, char* rvdstring, char *output_filepath, char *log_filepath, bint forwardonly, int c_upstream, double cutoff, int numprocs, char *organism_name):

        cdef Hashmap *talesf_kwargs = hashmap_new(32)
        
        cdef double weight = 0.9
        
        hashmap_add(talesf_kwargs, "seq_filename", seqfilename)
        hashmap_add(talesf_kwargs, "rvd_string", rvdstring)
        hashmap_add(talesf_kwargs, "output_filepath", output_filepath)
        hashmap_add(talesf_kwargs, "log_filepath", log_filepath)
        hashmap_add(talesf_kwargs, "weight", &weight)
        hashmap_add(talesf_kwargs, "cutoff", &cutoff)
        hashmap_add(talesf_kwargs, "forward_only", &forwardonly)
        hashmap_add(talesf_kwargs, "c_upstream", &c_upstream)
        hashmap_add(talesf_kwargs, "num_procs", &numprocs)
        hashmap_add(talesf_kwargs, "organism_name", organism_name)
  
        cdef int task_result = run_talesf_task(talesf_kwargs)
        
        hashmap_delete(talesf_kwargs, NULL)
        
        return task_result
