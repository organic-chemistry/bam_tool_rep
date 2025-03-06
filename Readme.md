Install
==============

```bash
conda create -n bam_tool_rep -c bioconda pysam
conda activate bam_tool_rep
git clone git@github.com:organic-chemistry/bam_tool_rep.git
pip install -e .
```
#eventually install notebook : conda install jupyterlab



Usage
================


For DNAscent v4
```python
from bam_tool_rep.bam_tools import load_read_bam_multi
bam = "./test_data/first_1000_reads_DNAscent_v4.bam"
remove_read_shorter_than=5_000 # in kb
maximum_of_read_processed=None
lower_mean_brdu_threshold=0.05 # remove reads whom median brdu content are lower that this value
r = load_read_bam_multi(bam,res=100,remove_less_than={"b":lower_mean_brdu_threshold},
                        remove_shorter_than=remove_read_shorter_than,
                        maxi=maximum_of_read_processed,threads=3)

# r is a dictionnary where each item is a read_id / processed read
# r[k][1]["b"] contain the binned read, where k is a key.
# r[k][0] contain information about sequence chromosomes and so on
```


Notebook
=================
See the notebooks repertory for example of usage

