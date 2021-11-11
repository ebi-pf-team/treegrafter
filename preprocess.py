#!/usr/bin/env python

import argparse
import os
import re
import shutil


def main():
    parser = argparse.ArgumentParser(description="""\
Pre-process TreeGrafter data files.
    """)
    parser.add_argument("datadir", help="TreeGrafter data directory")
    args = parser.parse_args()
    datadir = args.datadir

    paintdir = os.path.join(datadir, "PAINT_Annotations")
    paintfile = os.path.join(paintdir, "PAINT_Annotatations_TOTAL.txt")
    dirs = set()
    with open(paintfile, "rt") as fh:
        for line in fh:
            cols = line.split('\t')
            fam_id, an = cols[0].split(':')

            famdir = os.path.join(paintdir, fam_id)
            if fam_id not in dirs:
                try:
                    shutil.rmtree(famdir)
                except FileNotFoundError:
                    pass
                finally:
                    os.makedirs(famdir)
                    dirs.add(fam_id)

            with open(os.path.join(famdir, an), "wt") as fh2:
                # fh2.write("{}\n".format(cols[1].strip()))
                fh2.write("{}\t{}\n".format(cols[1].strip(), cols[2].strip()))

    os.unlink(paintfile)

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
