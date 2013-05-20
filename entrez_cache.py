import redis
import pickle
import uuid
import time
import os
from Bio import Entrez, SeqIO
from talconfig import BASE_DIR
from talutil import TaskError


# redis locking
# http://www.dr-josiah.com/2012/01/creating-lock-with-redis.html

class CachedEntrezFile(object):
    
    WRITE_LOCK_DURATION = 3600
    READ_LOCK_DURATION = 36000
    
    CACHED_FILE_DIR = BASE_DIR + "/cached_gb_files"
    CACHED_FILE_PATTERN = CACHED_FILE_DIR + "/%s.fasta"
    
    def __init__(self, seq_id):
        
        self.conn = redis.StrictRedis(host='localhost', port=6379, db=1)
        
        try:
            Entrez.email = "6e6a62393840636f726e656c6c2e656475".decode("hex")
            Entrez.tool = "https://tale-nt.cac.cornell.edu"
            search_handle = Entrez.esearch(db="nucleotide", term=seq_id)
            search_record = Entrez.read(search_handle)
            search_handle.close()
            self.seq_id = search_record["IdList"][0]
        except IOError:
            raise TaskError("Invalid sequence ID provided")
        
        self.lock_name = "gb_cache:%s" % self.seq_id
        self.uuid = str(uuid.uuid4())
        self.filepath = CachedEntrezFile.CACHED_FILE_PATTERN % self.seq_id
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
                        
                        fetch_handle = Entrez.efetch(db="nucleotide", id=self.seq_id, rettype="fasta", retmode="text")
                        fetch_record = SeqIO.read(fetch_handle, "fasta")
                        fetch_handle.close()
                        
                        with open(self.filepath, "w") as cached_file:
                            SeqIO.write([fetch_record], cached_file, "fasta")
                        
                    except IOError:
                        # I dunno
                        pass
                    
                    #subprocess.check_call("rm -f %s/chromFa.tar.gz" % CachedEntrezFile.CACHED_FILE_DIR, shell=True)
                    #subprocess.check_call("wget -O %s/chromFa.tar.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/chromFa.tar.gz" % CachedEntrezFile.CACHED_FILE_DIR, shell = True)
                    
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
                            
                            if l_timestamp < timestamp:
                                print("Writer for %s with UUID %s expired, removing" % (self.lock_name, l_uuid))
                                del lock_writers[l_uuid]
                                writers_changed = True
                        
                        if writers_changed:
                            
                            pipe.multi()
                            #pipe.hset(self.lock_name, "writers", pickle.dumps(lock_writers))
                            # since file will only be downloaded if key doesn't exist (because there's no way to check if local sequence
                            # file is corrupted), and there's only one writer at a time, remove key so that redownload of file is attempted
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
        
        return self.file
    
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
