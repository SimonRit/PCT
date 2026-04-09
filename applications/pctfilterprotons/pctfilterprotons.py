#!/usr/bin/env python
import argparse
import itk
import numpy as np


def build_parser():

    parser = argparse.ArgumentParser(
        description="Filter protons according to various criteria"
    )
    parser.add_argument("-i", "--input", help="Input list-mode file", required=True)
    parser.add_argument("-o", "--output", help="Output list-mode file", required=True)
    parser.add_argument("--roi-radius", help="Radius of the ROI in mm", type=float)
    parser.add_argument(
        "--fluence", help="Fluence level as fraction of full fluence", type=float
    )
    parser.add_argument("--min-wepl", help="Minimum WEPL", type=float)
    parser.add_argument("--max-wepl", help="Maximum WEPL", type=float)
    parser.add_argument(
        "--verbose", "-v", help="Verbose execution", default=False, action="store_true"
    )

    return parser


def process(args_info: argparse.Namespace):

    if args_info.verbose:

        def verbose(message):
            print(message)

    else:

        def verbose(message):
            pass

    pairs_itk = itk.imread(args_info.input)
    pairs = itk.GetArrayFromImage(pairs_itk)

    u_in = np.s_[:, 0, 0]
    w_in = np.s_[:, 0, 2]
    u_out = np.s_[:, 1, 0]
    w_out = np.s_[:, 1, 2]
    e_in = np.s_[:, 4, 0]
    wepl = np.s_[:, 4, 1]

    if np.any(pairs[e_in] != 0.0) and (
        args_info.min_wepl is not None or args_info.max_wepl is not None
    ):
        raise ValueError(
            "Cannot filter WEPLs because input file does not contain WEPL (e_in != 0)! Aborting."
        )

    interception = np.full_like(pairs[wepl], True)

    if args_info.roi_radius is not None:
        # Circular ROI intersection
        verbose("Filtering ions outside of the ROI…")
        radius = args_info.roi_radius
        du = pairs[u_out] - pairs[u_in]
        dt = pairs[w_out] - pairs[w_in]
        dr_sq = (du * du) + (dt * dt)
        d = (pairs[u_in] * pairs[w_out]) - (pairs[u_out] * pairs[w_in])
        delta = (radius * radius * dr_sq) - (d * d)
        interception = np.logical_and(interception, delta >= 0.0)

    if args_info.fluence is not None:
        # Fluence filtering
        verbose("Applying fluence level…")
        rng = np.random.default_rng()
        fluence_filter = rng.random(interception.shape)
        interception = np.logical_and(interception, fluence_filter < args_info.fluence)

    if args_info.min_wepl is not None:
        verbose("Filtering low WEPLs…")
        interception = np.logical_and(interception, pairs[wepl] >= args_info.min_wepl)

    if args_info.max_wepl is not None:
        verbose("Filtering high WEPLs…")
        interception = np.logical_and(interception, pairs[wepl] <= args_info.max_wepl)

    verbose("Applying interception filter…")

    ComponentType = itk.ctype("float")
    PixelType = itk.Vector[ComponentType, 3]
    ImageType = itk.Image[PixelType, 2]

    output_itk = itk.GetImageFromArray(pairs[interception], ttype=ImageType)
    itk.imwrite(output_itk, args_info.output)
    verbose(f"Wrote file {args_info.output}.")


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
