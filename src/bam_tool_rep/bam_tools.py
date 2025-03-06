import pysam
import tqdm
import numpy as np
import re
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def smooth(ser, sc):
    return np.array(pd.Series(ser).rolling(sc, min_periods=1, center=True).mean())

def find1(stro, ch):
  # 0.100 seconds for 1MB str
  npbuf = np.frombuffer(bytes(stro,'utf-8'), dtype=np.uint8) # Reinterpret str as a char buffer
  return np.where(npbuf == ord(ch))[0]


def convert_to_coordinate(seq, Ml, Mm, which=u"T"):
    # Ml base proba
    # Mm number of base to skip
    # print(sum(Mm),sum(np.array(list(seq))=="T"))
    # print(Mm)
    # print(seq)
    assert (len(Ml) == len(Mm))
    Mm = Mm + 1
    result = np.zeros((len(seq))) + np.nan
    n_which = 0
    # print(Mmc)
    cum = np.cumsum(Mm)
    cum -= 1
    #print(cum, cum.dtype)
    # which = r
    # array_s = np.fromiter(seq,dtype=np.char)
    # array_s=np.array(list(seq), dtype=np.unicode)
    #pos = np.array([m.start() for m in re.finditer(which, seq)])
    if which != "N":
        pos = find1(seq,which)
    else:
        pos = np.arange(len(result))
    #print(which,pos)
    #print()
    # shold not append
    if len(cum) != 0:
        if cum[-1] > len(pos) - 1:
            # truncate
            cum = cum[:np.argmax(cum > len(pos) - 1)]
            Ml = Ml[:len(cum)]
        # print(pos)

        result[pos[cum]] = np.array(Ml) / 255

    return result



def get_longest_low(v_mono):
    #print("Here")
    monov = smooth(v_mono, 1000)
    sum_b = np.nansum((smooth(v_mono, 100)) > 0.5)
    if np.isnan(sum_b) or (not sum_b > 30):
        return None,None
    selected = (monov > 0.015) & (monov < 0.35)
    found=False
    while np.sum(selected) > 1000:
        st = "".join([str(int(s)) for s in selected])
        sp = st.split("0")
        longest = np.argmax([len(ss) for ss in sp])

        def getl(ss):
            if ss == "":
                return 1
            else:
                return len(ss)

        start = int(sum([getl(ss) for ss in sp[:longest]]))
        end = int(start + len(sp[longest]))
        #print(start,end)
        if np.any(monov[start:end] > 0.1):
            found=True
            break
        else:
            selected[start:end] = 0
    if found:
        return start, end
    else:
        return None, None



class SerializableRead:
    """Lightweight, serializable version of pysam.AlignedSegment."""
    def __init__(self, read,sequence=""):
        # Basic attributes
        self.query_name = read.query_name
        self.reference_name = read.reference_name
        self.is_reverse = read.is_reverse
        self.seq = sequence  # Get the sequence directly

        # Positions (handle empty positions gracefully)
        pos = read.get_reference_positions()
        try:
            self.mapped_start = pos[0] if pos else None
            self.mapped_end = pos[-1] if pos else None
        except:
            pass
           

        # Tags (extract Ml/ML and Mm/MM)
        self.tags = {}
        for tag in ["Ml", "ML"]:
            if read.has_tag(tag):
                self.tags["Ml"] = read.get_tag(tag)
                break
        for tag in ["Mm", "MM"]:
            if read.has_tag(tag):
                self.tags["Mm"] = read.get_tag(tag)
                break

    def has_tag(self, tag):
        return tag in self.tags

    def get_tag(self, tag):
        return self.tags.get(tag)
    def get_forward_sequence(self):
        return self.seq
    
    def get_reference_positions(self):
        return [self.mapped_start,self.mapped_end]
    

def load_read_bam_multi(bam, threads=4,maxi=None, **kwargs):
    """Track progress correctly with tqdm on the read-fetching loop."""
    samfile = pysam.AlignmentFile(bam, "r", threads=threads, check_sq=False)
    Read = {}
    remove_shorter_than = kwargs.get("remove_shorter_that",None)

    with ProcessPoolExecutor(max_workers=threads) as executor:
        # Use tqdm on the loop that iterates over reads, not futures
        monitor = tqdm.tqdm(total=maxi, desc="Processing reads")
        futures = []

        for ir, read in enumerate(samfile):
            if maxi is not None and ir >= maxi:
                break

            seq = read.get_forward_sequence()
            
            if remove_shorter_than is not None and len(seq)<remove_shorter_than:
                continue
            read = SerializableRead(read,seq)
            futures.append(
                executor.submit(
                    process_single_read,
                    read=read,
                    verbose = kwargs.get("verbose",False),
                    chs=kwargs.get('chs',None),
                    res=kwargs.get('res', 1),
                    no_seq=kwargs.get('no_seq', False),
                    remove_less_than=kwargs.get('remove_less_than',None),
                    allready_mod=kwargs.get('allready_mod', True)
                )
            )
            monitor.update(1)  # Update progress bar here

        monitor.close()

        # Collect results (no tqdm here)
        for future in futures:
            result = future.result()
            if result:
                Read[result[0]] = result[1]

    samfile.close()
    return Read


def process_single_read(read,verbose=False,res=1,chs=None,
                        allready_mod=True,
                        remove_less_than=None,no_seq=False,remove_shorter_than=None):
    """
    filter_b and n_b are here to select signal with Brdu
    it select the signal if n_b points are higher that filter_b
    If you want to keep everything you can set n_b to 0
    res is the resolution of the final signal
    for example at 100, it does a smoothing average over 100 points and then select one point every 100 points

    chs can be a list of chromosome that you want to keep
    for example ["chr1","chr2"]
    the nomenclature has to be the same as the one of your reference file

    maxi is the maximum number of read that you want to process

    it returns an array for each read with [attr,b_val]
    The x coordinate is computed as x = np.arange(len(b_val)) * res + attr["mapped_start"]
    If the strand mapped to is "-" I already flipped the signal.
    (it means that the x coordinate are always increasing)
    remove_shorter_than is in bp
    """


    if verbose:
        print(read)

    seq = read.get_forward_sequence()   # return the actual sequence of the mapped read

    #print(Mm2)
    #print(len(Mm),len(Mm2))
    attr={}

    if no_seq:
        attr["mapped_strand"] = "+"
    else:
        if read.is_reverse:
            attr["mapped_strand"] = "-"
        else:
            attr["mapped_strand"] = "+"


        attr["mapped_chrom"] = read.reference_name


        pos = read.get_reference_positions()
        try:
            attr["mapped_start"] = pos[0]
            attr["mapped_end"] = pos[-1]
            attr["seq"]=seq
            

        except:
            pass	    
        #print(read.reference_name, read.reference_id)
        #print(read.header)
        #for attrn in dir(read):
    #        print(attrn,getattr(read,attrn))


        if chs is not None and attr["mapped_chrom"] not in chs:
            return None

    for tag in ["Ml","ML"]:
        if read.has_tag(tag):
            Ml = read.get_tag(tag)
            break

    for tag in ["Mm","MM"]:
        if read.has_tag(tag):
            Mmt = read.get_tag(tag).split(";")[:-1]
    
    Mm = {}
    base_ref={}
    for Smm in Mmt:
        base = Smm[2:3]
        #shift = [int(v) for v in Smm.split(",")[1:]]
        try:
            shift= np.fromstring(Smm[Smm.index(",")+1:], dtype=np.int, sep=',')
        except:
            print("Strange read skipping")
            return None

        Mm[base]=shift
        base_ref[base]=Smm[:1]
        
    
    Nn ={}
    start = 0
    for mod in Mm.keys():
        #if mod == "b":
        #    print(Ml[start:start+len(Mm[mod])])
        #    print(Mm[mod])
        #    print(seq)
        val = convert_to_coordinate(seq,Ml[start:start+len(Mm[mod])],Mm[mod],which=base_ref[mod])
        #if mod == "b":
        #    print(val[:15])
        #val[np.isnan(val) & (np.array(list(seq))=="T")]=0
        start += len(Mm[mod])
        #print(Ml)
        #print("Mm",Mm)
        #raise


        if attr["mapped_strand"] == "-":
            val = val[::-1]

        if not allready_mod:
            val[np.isnan(val) & (np.array(list(seq)) == base_ref[mod])] = 0

        if res != 1:
            #val = smooth(val,res)
            #val = np.array(val[::res],dtype=np.float16)
            n = res * (len(val)//res)
            val = val[:n].reshape(-1,res)
            npts =  np.nansum(~np.isnan(val),axis=1,dtype=np.float16)
            val= np.nanmean(val,axis=1,dtype=np.float16)
        else:
            npts =  np.array(~np.isnan(val),dtype=int)
        Nn[mod]=val
        Nn[mod+"_npoints"]=npts
    skip=False
    

    if remove_less_than is not None:
        if Nn == {}:
            return None
        else:
            for m in Mm.keys():
                if type(remove_less_than) == float:
                    th = remove_less_than
                else:
                    th = remove_less_than.get(m,0)
                if np.nanmean(Nn[m]) <= th:
                    return None



    return (read.query_name, (attr, Nn))

