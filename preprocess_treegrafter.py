#!/usr/bin/env python

import os
import shutil

import argparse
import logging


def get_args():
    """
    Command line arguments parser.
    """

    ap = argparse.ArgumentParser(
        prog='preprocess_treegrafter.py', description="TreeGrafter data preprocessor",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ap.add_argument(
        '-d', '--data', required=True,
        help='panther data directory')

    ap.add_argument(
        '-v', '--verbose', default='ERROR', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], type=str.upper,
        help='if set, logs debug info at the provided level')

    args = vars(ap.parse_args())

    return args


if __name__ == '__main__':

    args = get_args()

    data_folder = os.path.abspath(args['data'])

    log_formatter_str = '%(asctime)s | %(levelname)-8s | %(message)s'
    log_formatter = logging.Formatter(log_formatter_str)

    handlers = list()
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    handlers.append(console_handler)

    logging.basicConfig(level=args['verbose'],
                        format=log_formatter_str,
                        handlers=handlers)

    logger = logging.getLogger('treegrafter')

    logger.info('Preprocessing tregrafter data on ' + data_folder)

    annot_dir = os.path.join(data_folder, 'PAINT_Annotations')
    print(annot_dir)
    annot_file = os.path.join(annot_dir, 'PAINT_Annotatations_TOTAL.txt')
    print(annot_file)

    with open(annot_file) as f:
        for line in f:
            line = line.strip()
            # print(line)
            line_array = line.split("\t")
            # print(line_array)

            with open(annot_dir + '/' + line_array[0], 'w') as outfile:
                outfile.write(line_array[1] + "\t" + line_array[2])

    open(os.path.join(annot_dir, 'preprocessed'), 'w')

    logger.info('Preprocessing completed')

