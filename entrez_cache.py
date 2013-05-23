import pickle
import uuid
import time
import os
import tempfile
import urllib2
import math
from Bio import Entrez, SeqIO
from talconfig import BASE_DIR, REDIS_SERVER_HOSTNAME, REDIS_SERVER_PORT
from talutil import TaskError

redis_found = True
try:
    import redis
except ImportError:
    redis_found = False

# modified version of Entrez.efetch
def efetch_post(db, **keywds):
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    variables = {'db' : db}
    keywords = keywds
    if "id" in keywds and isinstance(keywds["id"], list):
        #Fix for NCBI change (probably part of EFetch 2,0, Feb 2012) where
        #a list of ID strings now gives HTTP Error 500: Internal server error
        #This was turned into ...&id=22307645&id=22303114&... which used to work
        #while now the NCBI appear to insist on ...&id=22301129,22299544,...
        keywords = keywds.copy() #Don't alter input dict!
        keywords["id"] = ",".join(keywds["id"])
    variables.update(keywords)
    return Entrez._open(cgi, variables, post=True)


# redis locking
# http://www.dr-josiah.com/2012/01/creating-lock-with-redis.html

if redis_found:
    
    class CachedEntrezFile(object):
        
        WRITE_LOCK_DURATION = 3600
        READ_LOCK_DURATION = 36000
        
        CACHED_FILE_DIR = BASE_DIR + "/cached_gb_files"
        CACHED_FILE_PATTERN = CACHED_FILE_DIR + "/%s.fasta"
        
        def __init__(self, logger, assembly_id):
            
            self.conn = redis.StrictRedis(host=REDIS_SERVER_HOSTNAME, port=REDIS_SERVER_PORT, db=1)
            
            try:
                
                Entrez.email = "6e6a62393840636f726e656c6c2e656475".decode("hex")
                Entrez.tool = "https://tale-nt.cac.cornell.edu"
                
                ehandle = Entrez.elink(db="nuccore", dbfrom="assembly", id=assembly_id, linkname="assembly_nuccore_refseq")
                erecord = Entrez.read(ehandle)
                ehandle.close()
                
                canonical_assembly_ids = [x for x in erecord[0]["IdList"] if x != "-1"]
                
                if len(canonical_assembly_ids) == 0:
                    logger("Warning: NCBI search for assembly ID returned more than 1 record, choosing the first")
                
                self.assembly_id = canonical_assembly_ids[0]
                
                nuc_seq_ids = []
                
                if len(erecord[0]["LinkSetDb"]) != 0:
                    nuc_seq_ids = [link["Id"] for link in erecord[0]["LinkSetDb"][0]["Link"]]
                
                if len(nuc_seq_ids) == 0:
                    
                    logger("Provided NCBI ID did not match an assembly, or matched an assembly with no associated sequences")
                    logger("Checking if NCBI ID matches nucleotide sequence")
                    
                    ehandle = Entrez.esearch(db="nucleotide", term=assembly_id)
                    erecord = Entrez.read(ehandle)
                    ehandle.close()
                    
                    if int(erecord["Count"]) != 0:
                        
                        if int(erecord["Count"]) > 1:
                            logger("Warning: NCBI search for nucleotide sequence ID returned more than 1 record, choosing the first")
                        
                        self.assembly_id = "nucleotide:" + erecord["IdList"][0]
                        nuc_seq_ids.append(erecord["IdList"][0])
                        
                    else:
                        raise TaskError("No assemblies or nucleotide sequences found for given ID")
                
                self.nuc_seq_ids = nuc_seq_ids
                
                #search_handle = Entrez.esearch(db="nucleotide", term=assembly_id)
                #search_record = Entrez.read(search_handle)
                #search_handle.close()
                #if int(search_record["Count"]) < 1:
                #    raise TaskError("Invalid assembly ID provided")
                
            except (IOError, urllib2.HTTPError):
                raise TaskError("Invalid NCBI ID provided")
            
            self.lock_name = "gb_cache:%s" % self.assembly_id
            self.uuid = str(uuid.uuid4())
            self.filepath = CachedEntrezFile.CACHED_FILE_PATTERN % self.assembly_id
            self.file = None
        
        def __enter__(self):
            
            pipe = self.conn.pipeline(True)
            
            while True:
                
                try:
                    
                    pipe.watch(self.lock_name)
                    
                    if not pipe.exists(self.lock_name):
                        
                        print("No entry for key %s, creating" % self.lock_name)
                        
                        timestamp = time.time()
                        
                        pipe.multi()
                        
                        pipe.hmset(self.lock_name, {
                            "last-accessed": timestamp,
                            "readers": pickle.dumps({}),
                            "writers": pickle.dumps({
                                self.uuid: timestamp + CachedEntrezFile.WRITE_LOCK_DURATION
                            })
                        })
                        
                        pipe.execute()
                        
                        print("Key %s created" % self.lock_name)
                        print("Downloading file for newly created key %s" % self.lock_name)
                        
                        try:
                            os.remove(self.filepath)
                        except OSError:
                            pass
                        
                        try:
                            
                            with open(self.filepath, "w") as cached_file:
                                
                                for i in range(int(math.ceil(float(len(self.nuc_seq_ids)) / 5000))):
                                    
                                    
                                    start = i * 5000
                                    remaining = len(self.nuc_seq_ids) - start
                                    end = start + remaining if remaining < 5000 else start + 5000
                                    
                                    print("start %d remaining %d end %d" % (start, remaining, end))
                                    
                                    fetch_handle = efetch_post(db="nucleotide", id=self.nuc_seq_ids[start:end], rettype="fasta", retmode="text")
                                    
                                    for line in fetch_handle:
                                        cached_file.write(line)
                                    
                                    fetch_handle.close()
                            
                        except IOError:
                            # I dunno
                            raise
                            
                        finally:
                            
                            print("File for key %s finished downloading, removing writer")
                            
                            pipe.watch(self.lock_name)
                            
                            lock_writers = pickle.loads(pipe.hget(self.lock_name, "writers"))
                            
                            if self.uuid in lock_writers:
                                del lock_writers[self.uuid]
                            
                            pipe.multi()
                            pipe.hset(self.lock_name, "writers", pickle.dumps(lock_writers))
                            pipe.execute()
                            
                            print("Writer for key %s removed" % self.lock_name)
                        
                    else:
                        
                        print("Key %s exists, checking for writer" % self.lock_name)
                        
                        timestamp = time.time()
                        
                        lock_writers = pickle.loads(pipe.hget(self.lock_name, "writers"))
                        
                        if len(lock_writers.keys()) > 0:
                            
                            print("Writer for key %s exists, checking expiration" % self.lock_name)
                            
                            writers_changed = False
                            
                            for l_uuid, l_timestamp in lock_writers.copy().items():
                                
                                print("Writer for key %s with UUID %s, timestamp %s, current time %s" % (self.lock_name, l_uuid, str(l_timestamp), str(timestamp)))
                                
                                if l_uuid == self.uuid:
                                    
                                    print("Writer for key %s has our UUID, previous attempt to remove writer failed, retrying" % self.lock_name)
                                    del lock_writers[l_uuid]
                                    writers_changed = True
                                    
                                elif l_timestamp < timestamp:
                                    
                                    print("Writer for %s with UUID %s expired, removing" % (self.lock_name, l_uuid))
                                    del lock_writers[l_uuid]
                                    writers_changed = True
                            
                            if writers_changed:
                                
                                pipe.multi()
                                #pipe.hset(self.lock_name, "writers", pickle.dumps(lock_writers))
                                # since file will only be downloaded if key doesn't exist (because there's no way to check if local sequence
                                # file is corrupted), and there's only one writer at a time, removing the only writer means the local file
                                # is likely corrupted, so we remove the key so that redownload of file is attempted
                                pipe.delete(self.lock_name)
                                pipe.execute()
                                
                                print("Expired writers for key %s removed" % self.lock_name)
                            
                            time.sleep(5)
                        
                        else:
                            
                            print("No writers for key %s, adding reader" % self.lock_name)
                            
                            timestamp = time.time()
                            
                            lock_readers = pickle.loads(pipe.hget(self.lock_name, "readers"))
                            
                            lock_readers[self.uuid] = timestamp
                            
                            pipe.multi()
                            pipe.hset(self.lock_name, "last-accessed", timestamp)
                            pipe.hset(self.lock_name, "readers", pickle.dumps(lock_readers))
                            pipe.execute()
                            
                            print("Reader for key %s acquired, returning file" % self.lock_name)
                            
                            break
                    
                except redis.exceptions.WatchError:
                    pass
            
            
            self.file = open(self.filepath, "r")
            
            return self
        
        def __exit__(self, exc_type, exc_value, traceback):
            
            if self.file is not None:
                self.file.close()
            
            print("Releasing reader for key %s" % self.lock_name)
            
            pipe = self.conn.pipeline(True)
            
            while True:
                
                try:
                    
                    pipe.watch(self.lock_name)
                    
                    if pipe.exists(self.lock_name):
                        
                        lock_readers = pickle.loads(pipe.hget(self.lock_name, "readers"))
                        
                        if self.uuid in lock_readers:
                            del lock_readers[self.uuid]
                        
                        pipe.multi()
                        pipe.hset(self.lock_name, "readers", pickle.dumps(lock_readers))
                        pipe.execute()
                        
                    else:
                        
                        # read lock was broken by garbage collector
                        
                        pipe.unwatch()
                    
                    break
                    
                except redis.exceptions.WatchError:
                    pass
            
            print("Reader for key %s released" % self.lock_name)
            
else:
    
    class CachedEntrezFile(object):
        
        def __init__(self, logger, assembly_id):
            
            try:
                
                Entrez.email = "6e6a62393840636f726e656c6c2e656475".decode("hex")
                Entrez.tool = "https://tale-nt.cac.cornell.edu"
                
                ehandle = Entrez.elink(db="nuccore", dbfrom="assembly", id=assembly_id, linkname="assembly_nuccore_refseq")
                erecord = Entrez.read(ehandle)
                ehandle.close()
                
                canonical_assembly_ids = [x for x in erecord[0]["IdList"] if x != "-1"]
                
                if len(canonical_assembly_ids) == 0:
                    logger("Warning: NCBI search for assembly ID returned more than 1 record, choosing the first")
                
                self.assembly_id = canonical_assembly_ids[0]
                
                nuc_seq_ids = []
                
                if len(erecord[0]["LinkSetDb"]) != 0:
                    nuc_seq_ids = [link["Id"] for link in erecord[0]["LinkSetDb"][0]["Link"]]
                
                if len(nuc_seq_ids) == 0:
                    
                    logger("Provided NCBI ID did not match an assembly, or matched an assembly with no associated sequences")
                    logger("Checking if NCBI ID matches nucleotide sequence")
                    
                    ehandle = Entrez.esearch(db="nucleotide", term=assembly_id)
                    erecord = Entrez.read(ehandle)
                    ehandle.close()
                    
                    if int(erecord["Count"]) != 0:
                        
                        if int(erecord["Count"]) > 1:
                            logger("Warning: NCBI search for nucleotide sequence ID returned more than 1 record, choosing the first")
                        
                        self.assembly_id = "nucleotide:" + erecord["IdList"][0]
                        nuc_seq_ids.append(erecord["IdList"][0])
                        
                    else:
                        raise TaskError("No assemblies or nucleotide sequences found for given ID")
                
                self.nuc_seq_ids = nuc_seq_ids
                
            except (IOError, urllib2.HTTPError):
                raise TaskError("Invalid NCBI ID provided")
            
            self.file = None
            self.filepath = ""
        
        def __enter__(self):
            
            try:
                
                fetch_handle = efetch_post(db="nucleotide", id=self.nuc_seq_ids, rettype="fasta", retmode="text")
                self.file = tempfile.NamedTemporaryFile()
                
                for i in range(int(math.ceil(float(len(self.nuc_seq_ids)) / 5000))):
                    
                    start = i * 5000
                    remaining = len(self.nuc_seq_ids) - start
                    end = start + remaining if remaining < 5000 else start + 5000
                    
                    fetch_handle = efetch_post(db="nucleotide", id=self.nuc_seq_ids[start:end], rettype="fasta", retmode="text")
                    
                    for line in fetch_handle:
                        self.file.write(line)
                    
                    fetch_handle.close()
                
                self.file.seek(0)
                
            except IOError:
                raise TaskError("Could not download sequence from NCBI")
            
            self.filepath = self.file.name
            
            return self
        
        def __exit__(self, exc_type, exc_value, traceback):
            
            if self.file is not None:
                self.file.close()
    