# jawm_featurecounts

This is a jawm featurecounts module.

Installing jawm:
```
pip install git+ssh://git@github.com/mpg-age-bioinformatics/jawm.git
```
For more information on jawm please visit jawm's repo on [GitHub.com](https://github.com/mpg-age-bioinformatics/jawm/tree/main).

Example usage:
```
# clone this module
git clone git@github.com:mpg-age-bioinformatics/jawm_featurecounts.git

cd jawm_featurecounts

# download test data
jawm-test -r download

# docker
jawm featurecounts.py featurecounts -p ./yaml/docker.yaml

# slurm & apptainer with multiple yaml files
jawm featurecounts.py featurecounts -p ./yaml/vars.yaml ./yaml/hpc.yaml
```

Additional jawm workflows are available [here (GitHub.com)](https://github.com/mpg-age-bioinformatics?q=jawm_&type=all&language=&sort=).
