#!/usr/bin/env python
import argparse
import uproot
import itk

import numpy as np
import numpy.lib.recfunctions as rfn

def pctpairprotons(
    input_in,
    input_out,
    output,
    plane_in,
    plane_out,
    min_run=0,
    max_run=1e6,
    no_nuclear=False,
    verbose=False,
    psin='PhaseSpace',
    psout='PhaseSpace'
):
    if verbose:
        def verbose(message):
            print(message)
    else:
        def verbose(message):
            pass

    def load_tree_as_df(root_file, tree_name):

        branch_names = [
            'RunID',
            'EventID',
            'TrackID',
            'KineticEnergy',
            'Position_X',
            'Position_Y',
            'Position_Z',
            'Direction_X',
            'Direction_Y',
            'Direction_Z'
        ]

        tree = uproot.open(root_file)[tree_name]
        branches = tree.arrays(branch_names, library='np')

        dtype = [(name, branch.dtype) for name, branch in branches.items()]
        ps = np.rec.recarray((len(branches['RunID']), ), dtype=dtype)
        for branch_name in branch_names:
            ps[branch_name] = branches[branch_name]

        ps = rfn.rename_fields(ps, {
            'Position_X': 'u',
            'Position_Y': 'v',
            'Position_Z': 'w',
        })
        ps = rfn.rename_fields(ps, {
            'Direction_X': 'du',
            'Direction_Y': 'dv',
            'Direction_Z': 'dw',
        })

        ps = ps[(ps['RunID'] >= min_run) & (ps['RunID'] < max_run)]

        return ps

    ps_in = load_tree_as_df(input_in, psin)
    ps_in['w'] = plane_in
    verbose("Read input phase space:\n" + str(ps_in))

    ps_out = load_tree_as_df(input_out, psout)
    ps_out['w'] = plane_out
    verbose("Read output phase space:\n" + str(ps_out))

    merge_columns = ['RunID', 'EventID']
    if no_nuclear:
        merge_columns.append('TrackID')

    # Remove duplicates
    _, unique_index = np.unique(ps_in[merge_columns], return_index=True)
    ps_in = ps_in[unique_index]

    ps_in.dtype.names = [n if n in merge_columns else n + '_in' for n in ps_in.dtype.names]
    ps_out.dtype.names = [n if n in merge_columns else n + '_out' for n in ps_out.dtype.names]
    ps_in_uniques = [n for n in ps_in.dtype.names if n not in merge_columns]
    ps_out_uniques = [n for n in ps_out.dtype.names if n not in merge_columns]

    if no_nuclear:
        # Easy case, there should be at most one row in ps_in and ps_out for keys ['RunID', 'EventID', 'TrackID']
        intersect, intersect_in, intersect_out = np.intersect1d(ps_in[merge_columns], ps_out[merge_columns], return_indices=True)
        pairs = rfn.merge_arrays((intersect, ps_in[ps_in_uniques][intersect_in], ps_out[ps_out_uniques][intersect_out]), asrecarray=True, flatten=True)
    else:
        # More complicated, there can be more than one row in ps_out for keys ['RunID', 'EventID'] (because of 'TrackID')
        # The solution is to repeat the computation for many TrackIDs, then merge it all together
        track_max = ps_out['TrackID_out'].max()
        verbose("Identified maximum number of tracks: " + str(track_max))
        pairs_list = []
        for t in range(track_max + 1):
            ps_out_t = ps_out[ps_out['TrackID_out'] == t]
            intersect, intersect_in, intersect_out = np.intersect1d(ps_in[merge_columns], ps_out_t[merge_columns], return_indices=True)
            pairs = rfn.merge_arrays((intersect, ps_in[ps_in_uniques][intersect_in], ps_out_t[ps_out_uniques][intersect_out]), asrecarray=True, flatten=True)
            if len(pairs) > 0:
                pairs_list.append(pairs)
        pairs = rfn.stack_arrays(pairs_list, asrecarray=True)
        np.recarray.sort(pairs, order=['RunID', 'EventID', 'TrackID_in', 'TrackID_out'])
    verbose("Merged input and output phase spaces.")

    number_of_runs = pairs['RunID'].max() + 1
    verbose("Identified number of runs: " + str(number_of_runs))

    ComponentType = itk.ctype('float')
    PixelType = itk.Vector[ComponentType, 3]
    ImageType = itk.Image[PixelType, 2]

    run_range = range(min_run, min(number_of_runs, max_run))
    for r in run_range:
        ps_run = pairs[pairs['RunID'] == r]
        if len(ps_run) == 0:
            continue

        ps_np = np.empty(shape=(len(ps_run), 5, 3), dtype=np.float32)
        ps_np[:,0,0] = ps_run['u_in']
        ps_np[:,0,1] = ps_run['v_in']
        ps_np[:,0,2] = ps_run['w_in']
        ps_np[:,1,0] = ps_run['u_out']
        ps_np[:,1,1] = ps_run['v_out']
        ps_np[:,1,2] = ps_run['w_out']
        ps_np[:,2,0] = ps_run['du_in']
        ps_np[:,2,1] = ps_run['dv_in']
        ps_np[:,2,2] = ps_run['dw_in']
        ps_np[:,3,0] = ps_run['du_out']
        ps_np[:,3,1] = ps_run['dv_out']
        ps_np[:,3,2] = ps_run['dw_out']
        ps_np[:,4,0] = ps_run['KineticEnergy_in']
        ps_np[:,4,1] = ps_run['KineticEnergy_out']
        ps_np[:,4,2] = ps_run['TrackID'] if no_nuclear else ps_run['TrackID_out']

        df_itk = itk.GetImageFromArray(ps_np, ttype=ImageType)

        output_file = output.replace('.', f'{r:04d}.')
        itk.imwrite(df_itk, output_file)
        verbose(f"Wrote file {output_file}.")

def main():

    parser = argparse.ArgumentParser(description="Pair corresponding protons from GATE ROOT files")
    parser.add_argument('-i', '--input-in', help="Root phase space file of particles before object", required=True)
    parser.add_argument('-j', '--input-out', help="Root phase space file of particles after object", required=True)
    parser.add_argument('-o', '--output', help="Output file name", required=True)
    parser.add_argument('--plane-in', help="Plane position of incoming protons", required=True, type=float)
    parser.add_argument('--plane-out', help="Plane position of outgoing protons", required=True, type=float)
    parser.add_argument('--min-run', help="Minimum run (inclusive)", default=0, type=int)
    parser.add_argument('--max-run', help="Maximum run (exclusive)", default=1e6, type=int)
    parser.add_argument('--no-nuclear', help="Remove inelastic nuclear collisions", default=False, action='store_true')
    parser.add_argument('--verbose', '-v', help="Verbose execution", default=False, action='store_true')
    parser.add_argument('--psin', help="Name of tree in input phase space", default='PhaseSpace')
    parser.add_argument('--psout', help="Name of tree in output phase space", default='PhaseSpace')
    args_info = parser.parse_args()

    pctpairprotons(**vars(args_info))

if __name__ == '__main__':
  main()