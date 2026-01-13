import jawm
import os

geneid=jawm.Process(
    name="geneid",
    when=lambda p: not os.path.isfile( os.path.join( p.var["featurecounts_output"] , p.var["pair_id"]+".gene.featureCounts.txt" ) ) ,
    script="""\
if [ -f {{strand}} ] ; then strand=$(cat {{strand}} ) ; else strand="{{strand}}" ; fi
echo ${strand}
featureCounts -a {{gtf}} -T {{cpus}} -g gene_id -o {{featurecounts_output}}/{{pair_id}}.gene.featureCounts.txt {{paired}} -s ${strand} {{bam}}
""",
    var={
        "paired":""
        },
    desc={
        "strand":"",
        "gtf":"",
        "cpus": "",
        "featurecounts_output":"",
        "pair_id":"",
        "pair_id":"",
        "bam":"",
        "paired":"",
        "mapping_folder":"Required to find the bams"
    },
    container="mpgagebioinformatics/subread:2.0.3",
    manager_slurm={ "-c": 8, "--mem": "80GB", "-t": "6:00:00" }
)

biotype=jawm.Process(
    name="biotype",
    when=lambda p: not os.path.isfile( os.path.join( p.var["featurecounts_output"] , p.var["pair_id"]+".biotype.featureCounts.txt" ) ) ,
    script="""\
cat << 'EOF' > {{featurecounts_output}}/biotypes_header.txt
# id: 'biotype-counts'
# section_name: 'Biotype Counts'
# description: "shows reads overlapping genomic features of different biotypes,
#     counted by <a href='http://bioinf.wehi.edu.au/featureCounts'>featureCounts</a>."
# plot_type: 'bargraph'
# anchor: 'featurecounts_biotype'
# pconfig:
#     id: "featureCounts_biotype_plot"
#     title: "featureCounts: Biotypes"
#     xlab: "# Reads"
#     cpswitch_counts_label: "Number of Reads
EOF
if [ -f {{strand}} ] ; then strand=$(cat {{strand}} ) ; else strand="{{strand}}" ; fi

echo "featureCounts -a {{gtf}} -T {{cpus}} -g gene_biotype -o {{featurecounts_output}}/{{pair_id}}.biotype.featureCounts.txt {{paired}} -s ${strand} {{bam}}"q
featureCounts -a {{gtf}} -T {{cpus}} -g gene_biotype -o {{featurecounts_output}}/{{pair_id}}.biotype.featureCounts.txt {{paired}} -s ${strand} {{bam}}
cp {{featurecounts_output}}/biotypes_header.txt {{featurecounts_output}}/{{pair_id}}.biotype_counts_mqc.txt
cut -f 1,7 {{featurecounts_output}}/{{pair_id}}.biotype.featureCounts.txt | tail -n +3 | grep -v '^\\s' >> {{featurecounts_output}}/{{pair_id}}.biotype_counts_mqc.txt
""",
    var={
        "paired":""
    },
    desc={
        "cpus": "",
        "bam":"",
        "featurecounts_output":"",
        "gtf":"",
        "pair_id":"",
        "paired":"",
        "strand":"",
        "mapping_folder":"Required to find the bams"
    },
    container="mpgagebioinformatics/subread:2.0.3",
    manager_slurm={ "-c": 8, "--mem": "80GB", "-t": "6:00:00" }
)

headers=jawm.Process(
    name="headers",
    when=lambda p: not os.path.isfile( os.path.join( p.var["featurecounts_output"] , "headers.touch" ) ) ,
    script="""\
#!/usr/local/bin/python
from pathlib import Path
import os
import csv
import pandas as pd
import numpy as np
files_path="{{featurecounts_output}}/"
files=os.listdir(files_path)
files=[s for s in files if "summary" not in s and "mqc" not in s]
gene_files=[ s for s in files if s [-len(".gene.featureCounts.txt"):] == ".gene.featureCounts.txt" ]
biotype_files=[ s for s in files if s[-len(".biotype.featureCounts.txt"):] == ".biotype.featureCounts.txt" ]
files=gene_files+biotype_files
for f in files:
    file=pd.read_csv(files_path+f, sep="\\t")
    header=file.columns[0]
    file_=pd.read_csv(files_path+f, sep="\\t", comment="#")

    if ".gene.featureCounts.txt" in f:
        sample_name=f.split("_gene",1)[0]
    elif ".biotype.featureCounts.txt" in f:
        sample_name=f.split("_biotype",1)[0]

    file_=file_.rename(columns={"pseudoalignments.bam":sample_name})

    text_file = open(files_path+f, 'w')
    text_file.write(header+"\\n")
    text_file.close()

    file_.to_csv(files_path+f, mode='a', header=True, sep="\\t", index=None)

    file_sum=pd.read_csv(files_path+f+".summary", sep="\\t")
    file_sum=file_sum.rename(columns={"pseudoalignments.bam":f})
    file_sum.to_csv(files_path+f+".summary", index=None, sep="\\t")
Path(f"{files_path}/headers.touch").touch()
""",
    desc={
        "bam":"",
        "featureCounts_output":""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-1",
    manager_slurm={ "-c": 1, "--mem": "1GB", "-t": "1:00:00" }
)


if __name__ == "__main__":
    import sys
    from jawm.utils import workflow, id_files
    from pathlib import Path

    workflows, var, args, unknown_args = jawm.utils.parse_arguments([ 'main', 'featurecounts', 'test' ])

    if workflow([ "main", "featurecounts", "test" ], workflows):

        id_bams=id_files( geneid.var["mapping_folder"] )

        feature_counts_jobs=[]

        for pair_id in id_bams.keys() : 

            geneid_=geneid.clone()
            geneid_.var["pair_id"]=pair_id
            geneid_.var["map.bam"]=id_bams[pair_id][0]
            geneid_.execute()

            biotype_=biotype.clone()
            biotype_.var["pair_id"]=pair_id
            biotype_.var["map.bam"]=id_bams[pair_id][0]
            biotype_.execute()
            
            feature_counts_jobs = feature_counts_jobs + [ geneid_.hash, biotype_.hash ]

        jawm.Process.wait( feature_counts_jobs ) 

        if os.path.basename( list(id_bams.values())[0][0] ) == "pseudoalignments.bam" :

            headers.execute( )

            jawm.Process.wait( [ headers.hash ]  )

    if workflow("test", workflows):

        # for the test workflow we might also do something more
        print("Test completed.")
