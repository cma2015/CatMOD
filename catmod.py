# -*- coding: utf-8 -*-

"""CatMOD is a .

This is a rough draft of the author's analysis of ***.
We will continue to add applications until we felt ready to say OK, and
we will release an official version through pip and bioconda.

"""

from sys import exit, version_info

from CatMOD import fullhelp_argumentparser
from CatMOD.sys_output import Output

# version control
if version_info[0] == 3 and version_info[1] >= 9:
    pass
else:
    output = Output()
    output.error('This program requires at least python3.9')
    exit()


def main():
    """Create subcommands and execute."""
    parser = fullhelp_argumentparser.FullHelpArgumentParser()
    subparser = parser.add_subparsers()
    data_process = fullhelp_argumentparser.DPArgs(
        subparser,
        'data_process',
        """.""")
    extract_features = fullhelp_argumentparser.EFArgs(
        subparser,
        'extract_features',
        """.""")
    train = fullhelp_argumentparser.TrainArgs(
        subparser,
        'train',
        """.""")
    predict = fullhelp_argumentparser.PredictArgs(
        subparser,
        'predict',
        """.""")

    def bad_args(args):
        """Print help on bad arguments."""
        parser.print_help()
        exit()

    parser.set_defaults(func=bad_args)
    arguments = parser.parse_args()
    arguments.func(arguments)


if __name__ == '__main__':
    main()
