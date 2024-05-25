import argparse
import sys
import rnavirhost


def main():
    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"""RNAVirusHost: a machine learning-based method for predicting hosts of RNA viruses through viral genomes
https://github.com/GreyGuoweiChen/VirHost.git

usage: rnavirhost <program> [options]

programs:
    classify_order          classifier query viruses at order level
    predict       predict hosts of the query viruses""",
    )

    subparsers = parser.add_subparsers(help=argparse.SUPPRESS)

    classify_order_parser = subparsers.add_parser(
        "classify_order",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Classifier query viruses at order level
\nusage: rnavirhost classify_order -i <input> -o <taxa.csv> [options]""",
    )
    rnavirhost.classify_order.fetch_arguments(classify_order_parser)

    predict_parser = subparsers.add_parser(
        "predict",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Predict hosts of the query viruses
\nusage: rnavirhost predict -i <input> --taxa <taxa.csv> [options] -o <output> [options]""",
    )
    rnavirhost.predict.fetch_arguments(predict_parser)


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "classify_order":
            classify_order_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "predict":
            predict_parser.print_help()
            sys.exit(0)
        else:
            parser.print_help()
            sys.exit(0)

    args = vars(parser.parse_args())
    args["func"](args)

