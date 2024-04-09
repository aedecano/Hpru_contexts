#!/usr/bin/env python

import sys
import subprocess

samples = sys.argv[1]

with open(samples) as sp:
    for line in sp:
        s = line.strip()
        subprocess.run(["esummary", "-db", "biosample", "-id", s, ">", f"{s}.html"])
        subprocess.run(["cp", f"{s}.html", "CTM15positive/"])

