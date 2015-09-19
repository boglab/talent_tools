# License

All source code is available under an ISC license.

Copyright (c) 2011-2015, Nick Booher <njbooher@gmail.com> and Erin Doyle <edoyle@iastate.edu>.

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

# Description

This directory contains the source code for the TALEN Targeter, Old TALEN Targeter, and TAL Effector Targeter tools from the [TALE-NT](https://tale-nt.cac.cornell.edu) website.

# Dependencies and Set Up

boglab_tools requires Python 2.6+ and [biopython](http://pypi.python.org/pypi/biopython).

By default the tools expect to be located at /opt/boglab/talent. Changing this requires editing 2 files:

1. In talconfig.py, change the value of BASE_DIR to the full path of the parent folder of boglab_tools
2. In findTAL.py find the line that says 'with open(BASE_DIR + "/talent/re_dict_dump", "rb") as re_dict_file:' and change 'talent' to 'boglab_tools'

Note that off-target counting with findTAL.py won't work unless [tfcount](https://github.com/boglab/tfcount) is installed, and findRvdTAL.py and findPairedRvdTALs.py won't work unless you install the C libraries and cython wrappers from [talesf](https://github.com/boglab/talesf) and [talesf/paired](https://github.com/boglab/talesf/tree/paired).

# Usage

Here's a mapping between the the tools on the site and the scripts in boglab_tools:

* TALEN Targeter - findTAL.py
* TALEN Targeter (old version with design guidelines) - findTAL_old.py
* TAL Effector Targeter - findSingleTALsite.py
* Target Finder - findRvdTAL.py
* Paired Target Finder - findPairedRvdTALs.py

You can run any of the scripts by doing:

```
python findTAL.py OPTIONS
```

You can find out what options exist for any of them by doing:

```
python findTAL.py --help
```

The command for running findTAL.py with the default parameters from the site is:

```
python findTAL.py --filter '0' --gspec --min '15' --max '30' --arraymin '15' --arraymax '20' --fasta 'filemx4P32.fasta' --outpath 'fileOgILYO.txt' --logpath '79637.log'
```

* --min and --max are the minimum and maximum spacer length
* --arraymin  and --arraymax are the minimum and maximum RVD sequence length
* --fasta is the path to the file containing your input sequences in [FASTA](http://en.wikipedia.org/wiki/FASTA_format) format
* --gspec tells the script to use NH instead of NN for G
* --outpath is where you want the output file to go
* --logpath is the path to write the process log to; if you leave this out it will print this to stdout on your terminal

To search for off-target sites in a custom genome, run findTAL.py as above, but add the options '--offtargets --offtargets-fasta PATH_TO_CUSTOM_GENOME'
