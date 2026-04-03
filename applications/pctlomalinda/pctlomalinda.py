#!/usr/bin/env python
import argparse
import uproot
import itk

import numpy as np


def build_parser():

    parser = argparse.ArgumentParser(description="Convert Loma Linda data to PCT data")
    parser.add_argument(
        "-i", "--input", help="Root phase space file of particles", required=True
    )
    parser.add_argument("-o", "--output", help="Output file name", required=True)
    parser.add_argument(
        "--plane-in",
        help="Plane position of incoming protons",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--plane-out",
        help="Plane position of outgoing protons",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--min-run", help="Minimum run (inclusive)", default=0, type=int
    )
    parser.add_argument(
        "--max-run", help="Maximum run (exclusive)", default=1e6, type=int
    )
    parser.add_argument(
        "--verbose", "-v", help="Verbose execution", default=False, action="store_true"
    )
    parser.add_argument(
        "--ps", help="Name of tree in input phase space", default="PhaseSpace"
    )

    return parser


def process(args_info: argparse.Namespace):

    if args_info.verbose:

        def verbose(message):
            print(message)

    else:

        def verbose(message):
            pass

    tree = uproot.open(args_info.input)[args_info.ps]
    pairs = tree.arrays(library="np")
    verbose("Read input phase space:\n" + str(pairs))

    # should match what it in the file, no checks is performed
    pairs["u_hit1"] = np.full_like(pairs["u_hit1"], args_info.plane_in)
    pairs["u_hit2"] = np.full_like(pairs["u_hit2"], args_info.plane_out)

    verbose("Filter RunIDs…")
    # "projection_angle" acts as the run ID here
    interception = np.logical_and(
        pairs["projection_angle"] < args_info.max_run,
        pairs["projection_angle"] >= args_info.min_run,
    )
    for k, v in pairs.items():
        pairs[k] = v[interception]

    verbose("Calculating directions…")
    dir_in_t = pairs["t_hit1"] - pairs["t_hit0"]
    dir_in_v = pairs["v_hit1"] - pairs["v_hit0"]
    dir_in_u = pairs["u_hit1"] - pairs["u_hit0"]
    norm_in = np.sqrt(dir_in_t**2 + dir_in_v**2 + dir_in_u**2)
    dir_out_t = pairs["t_hit3"] - pairs["t_hit2"]
    dir_out_v = pairs["v_hit3"] - pairs["v_hit2"]
    dir_out_u = pairs["u_hit3"] - pairs["u_hit2"]
    norm_out = np.sqrt(dir_out_t**2 + dir_out_v**2 + dir_out_u**2)

    ComponentType = itk.ctype("float")
    PixelType = itk.Vector[ComponentType, 3]
    ImageType = itk.Image[PixelType, 2]

    runs = np.unique(pairs["projection_angle"])
    number_of_runs = len(runs)
    verbose("Identified number of runs: " + str(number_of_runs))
    run_range = range(args_info.min_run, min(number_of_runs, args_info.max_run))
    for r in run_range:
        ps_run = pairs["projection_angle"] == runs[r]
        len_ps_run = ps_run.sum()
        if len_ps_run == 0:
            continue

        ps_np = np.empty(shape=(len_ps_run, 5, 3), dtype=np.float32)
        ps_np[:, 0, 0] = pairs["t_hit1"][ps_run]
        ps_np[:, 0, 1] = pairs["v_hit1"][ps_run]
        ps_np[:, 0, 2] = pairs["u_hit1"][ps_run]
        ps_np[:, 1, 0] = pairs["t_hit2"][ps_run]
        ps_np[:, 1, 1] = pairs["v_hit2"][ps_run]
        ps_np[:, 1, 2] = pairs["u_hit2"][ps_run]
        ps_np[:, 2, 0] = dir_in_t[ps_run] / norm_in[ps_run]
        ps_np[:, 2, 1] = dir_in_v[ps_run] / norm_in[ps_run]
        ps_np[:, 2, 2] = dir_in_u[ps_run] / norm_in[ps_run]
        ps_np[:, 3, 0] = dir_in_t[ps_run] / norm_out[ps_run]
        ps_np[:, 3, 1] = dir_in_v[ps_run] / norm_out[ps_run]
        ps_np[:, 3, 2] = dir_in_u[ps_run] / norm_out[ps_run]
        ps_np[:, 4, 0] = 0.0  # WEPL is already present in the ROOT file
        ps_np[:, 4, 1] = pairs["calculated_WEPL"][ps_run]
        ps_np[:, 4, 2] = pairs["timestamp"][ps_run]

        df_itk = itk.GetImageFromArray(ps_np, ttype=ImageType)

        output_file = args_info.output.replace(".", f"{r:04d}.")
        itk.imwrite(df_itk, output_file)
        verbose(f"Wrote file {output_file}.")


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
