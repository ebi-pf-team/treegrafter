#!/usr/bin/env python

import argparse
import json
import os
import re


def main():
    parser = argparse.ArgumentParser(description="""\
Pre-process TreeGrafter data files.
    """)
    parser.add_argument("datadir", help="TreeGrafter data directory")
    args = parser.parse_args()
    datadir = args.datadir

    paintdir = os.path.join(datadir, "PAINT_Annotations")
    paintfile = os.path.join(paintdir, "PAINT_Annotatations_TOTAL.txt")

    families = {}
    with open(paintfile, "rt") as fh:
        for line in fh:
            cols = line.split('\t')
            fam_id, an = cols[0].split(':')

            try:
                fam = families[fam_id]
            except KeyError:
                fam = families[fam_id] = {}

            # May contain: subfamily ID, GO terms, PANTHER Protein class
            # separated by spaces
            annotation = cols[1].strip()

            match = re.match(r"PTHR\d+:(SF\d+)", annotation)
            if match:
                subfam_id = match.group(1)
                annotation = annotation[match.end(1):].strip()
            else:
                subfam_id = ""

            # sub -> replace multiple spaces
            fam[an] = (subfam_id, re.sub(r"\s+", " ", annotation))

    os.unlink(paintfile)

    for fam_id, obj in families.items():
        with open(os.path.join(paintdir, fam_id + ".json"), "wt") as fh:
            json.dump(obj, fh)

    msfdir = os.path.join(datadir, "Tree_MSF")
    for name in os.listdir(msfdir):
        if name.endswith(".fasta"):
            src = os.path.join(msfdir, name)
            dst = src + ".tmp"
            with open(src, "rt") as fh1, open(dst, "wt") as fh2:
                for line in fh1:
                    if line[0] == ">":
                        fh2.write(line)
                    else:
                        # Replace Selenocysteine (U) and Pyrrolysine (P) AA
                        # by the undetermined AA character (X).
                        fh2.write(re.sub(r"[UO]", r"X", line.upper()))

            os.unlink(src)
            os.rename(dst, src)


if __name__ == '__main__':
    main()
