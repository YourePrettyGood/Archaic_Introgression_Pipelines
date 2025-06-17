# Repository of pipelines for detecting archaic introgression

This repository contains pipelines for running a variety of different methods
to detect archaic introgression, especially in the context of introgression
from archaic hominins (e.g. Neanderthals and Denisovans) into modern humans.
The code is relatively generalized, so it may be possible to run on other
species, although this possibility is untested.

Furthermore, this repository contains scripts for simulations (using msprime)
and evaluation of the accuracy of introgressed haplotype detection by each
archaic introgression method. Such simulations and evaluations are being used
for a manuscript currently in progress.

## Archaic introgression pipelines available

Pipelines are currently available for detecting archaic introgression using
the following methods:

1. [Sprime](https://github.com/Browning-lab/sprime)
1. [hmmix](https://github.com/LauritsSkov/Introgression-detection)
1. [ArchaicSeeker2.0](https://github.com/Shuhua-Group/ArchaicSeeker2.0)

We currently have plans to add pipelines for:

- [ ] [S*](https://github.com/bvernot/freezing-archer)
- [ ] [IBDmix](https://github.com/PrincetonUniversity/IBDmix)
- [ ] [Relate](https://github.com/MyersGroup/relate)

## Usage

At the moment, pipelines are written in Nextflow DSL2, though they are not
modularized, so must be run separately for the same input dataset. All
configurable parameters including input paths or path globs, parsing regexes,
output file prefixes, introgression method parameters, and process submission
parameters must be supplied in a configuration file that is passed to the
Nextflow invocation.

In general, Nextflow runs of each pipeline are done with a command like:

```bash
NXF_OPTS="-Xms500M -Xmx3500M" /usr/bin/time -v nextflow -c [path to .config file] run -ansi-log false [absolute path to .nf pipeline script] -profile [comma-separated list of profiles] -w [scratch/working directory] 2> [STDERR log] > [STDOUT log]
```

Note that `-Xmx` here is setting the maximum heap size allowed for the Nextflow
driver/head process. `3500M` should work for most invocations, even those
generating thousands of individual jobs. Just note that an out-of-memory (OOM)
kill of the Nextflow head process will result in a dirty pipeline exit that
does not resume correctly and truncates the log and cache.

All Nextflow pipelines have been tested with Nextflow 22.10 and 24.10, though
may not be compatible with the most recent versions of Nextflow due to
breaking changes with the `shell:` process block.

## Limitations

Current pipeline limitations include:

1. Methods generally expect to be run on autosomes/diploid chromosomes
1. Pipelines are parallelized by chromosome, so chromosomes are expected to
be named `chr[0-9]+` or just `[0-9]+` (e.g. chr1 through chr22, or 1 through 22)
1. Most regexes are evaluated against the part of the filename before the
first period `.`, so file naming schemes that use `.` as a delimiter should
be avoided or renamed. Underscores `_` are generally preferred here
1. Pipelines that spit out hundreds of jobs in parallel include the `array`
directive in the respective `process` blocks to avoid hitting SLURM job
submission rate limits, though this is not a perfect solution

## Contributors

- [Chang Liu](https://github.com/Chg-Liu/)
- [Patrick F. Reilly](https://github.com/YourePrettyGood/)
- [Aaron Pfennig](https://github.com/arpfennig/)
