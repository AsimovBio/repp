---
layout: default
title: features
parent: make
grand_parent: repp
nav_order: 1
---
## repp make features

Find or build a plasmid from its constituent features

### Synopsis

Find or build a plasmid from its constituent features

```
repp make features "[feature],...[featureN]" [flags]
```

### Examples

```
repp make features "BBa_R0062,BBa_B0034,BBa_C0040,BBa_B0010,BBa_B0012" --backbone pSB1C3 --enzymes "EcoRI,PstI" --igem
```

### Options

```
  -a, --addgene           use the Addgene repository
  -b, --backbone string   backbone to insert the fragments into. Can either be an entry 
                          in one of the dbs or a file on the local filesystem.
  -d, --dbs string        comma separated list of local fragment databases
  -u, --dnasu             use the DNASU repository
  -e, --enzymes string    comma separated list of enzymes to linearize the backbone with.
                          The backbone must be specified. 'repp ls enzymes' prints a list of
                          recognized enzymes.
  -x, --exclude string    keywords for excluding fragments
  -h, --help              help for features
  -p, --identity int      %-identity threshold (see 'blastn -help') (default 98)
  -g, --igem              use the iGEM repository
  -o, --out string        output file name
```

### Options inherited from parent commands

```
  -s, --settings string   build settings (default "~/.repp/config.yaml")
  -v, --verbose           whether to log results to stdout
```

### SEE ALSO

* [repp make](repp_make)	 - Make a plasmid from its expected sequence, features or fragments

###### Auto generated by spf13/cobra on 29-Jan-2020
