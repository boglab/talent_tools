from libc.stdlib cimport malloc, free, calloc

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
        cdef int count_only = 0
        
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
        
        hashmap_add(paired_talesf_kwargs, "count_only", &count_only)
        
        # add a way to pass in an error string
        
        cdef int task_result = run_paired_talesf_task(paired_talesf_kwargs)
        
        hashmap_delete(paired_talesf_kwargs, NULL)
        
        return task_result

def CountOffTargetsTask(char *seqfilename, char *rvd_string, char *rvd_string2, int c_upstream, double cutoff, int spacer_min, int spacer_max, int numprocs):

        cdef Hashmap *paired_talesf_kwargs = hashmap_new(32)
        
        cdef double weight = 0.9
        cdef int count_only = 1
        
        cdef int *count_results_array_f = <int *> calloc(2, sizeof(int))
        cdef int *count_results_array_r = <int *> calloc(2, sizeof(int))
        
        cdef int **count_results_array = <int **> calloc(2, sizeof(int*))
        
        cdef int i
        cdef int j
        
        hashmap_add(paired_talesf_kwargs, "seq_filename", seqfilename)
        hashmap_add(paired_talesf_kwargs, "rvd_string", rvd_string)
        hashmap_add(paired_talesf_kwargs, "rvd_string2", rvd_string2)
        hashmap_add(paired_talesf_kwargs, "output_filepath", "/dev/null")
        hashmap_add(paired_talesf_kwargs, "log_filepath", "/dev/null")
        hashmap_add(paired_talesf_kwargs, "weight", &weight)
        hashmap_add(paired_talesf_kwargs, "cutoff", &cutoff)
        hashmap_add(paired_talesf_kwargs, "c_upstream", &c_upstream)
        hashmap_add(paired_talesf_kwargs, "spacer_min", &spacer_min)
        hashmap_add(paired_talesf_kwargs, "spacer_max", &spacer_max)
        hashmap_add(paired_talesf_kwargs, "num_procs", &numprocs)
        
        hashmap_add(paired_talesf_kwargs, "count_only", &count_only)
        hashmap_add(paired_talesf_kwargs, "count_results_array", count_results_array)
        
        # add a way to pass in an error string
        
        cdef int task_result = run_paired_talesf_task(paired_talesf_kwargs)
        
        hashmap_delete(paired_talesf_kwargs, NULL)
        
        count_results_list = []
        cdef int count_results_total = 0
        
        for i in range(2):
                for j in range(2):
                        count_results_total += count_results_array[i][j]
                        count_results_list.append(count_results_array[i][j])
        
        count_results_list.append(count_results_total)
        
        free(count_results_array_f)
        free(count_results_array_r)
        free(count_results_array)
        
        return count_results_list

