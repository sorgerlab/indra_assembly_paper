import sys
import csv
from os.path import join, dirname, abspath
import numpy as np
from matplotlib import pyplot as plt
from indra.tools import assemble_corpus as ac
from bioexp.curation.process_curations import get_curations_for_reader, \
                                    reader_input, load_curated_pkl_files, \
                                    get_correctness_data
from bioexp.util import set_fig_params, format_axis, fontsize, reader_name_map


def plot_correctness_curve(reader, all_stmts, show_ylabel=True, allow_incomplete=False):
    ev_correct_by_num_ev = get_curations_for_reader(reader, all_stmts,
                                                    'evidence',
                                                    allow_incomplete=allow_incomplete)
    data = ev_correct_by_num_ev
    data_stmt = {}
    # Convert at the evidence level to data at the stmt level and store
    for num_ev in data.keys():
        stmt_corrects = list(np.array(np.array(data[num_ev]) >= 1,
                                      dtype=int))
        data_stmt[num_ev] = stmt_corrects

    num_evs = sorted(data.keys())
    means = [np.mean(data_stmt[n]) for n in num_evs]
    # Stderr of proportion is sqrt(pq/n)
    std = [2*np.sqrt((np.mean(data_stmt[n]) *
                      (1 - np.mean(data_stmt[n]))) /
                       len(data_stmt[n]))
           for n in num_evs]

    # -- Plot statement correctness curve --
    #set_fig_params()
    plt.figure(figsize=(1.15, 1.25), dpi=150)
    plt.errorbar(num_evs, means, yerr=std, fmt='bo-', ls='none',
                 label='Empirical mean correctness', linewidth=1, markersize=3)
    # Legends, labels, etc.
    plt.ylim(0, 1.1)
    #plt.grid(True)
    plt.xticks(num_evs)
    plt.xlabel('Mentions')
    if show_ylabel:
        plt.yticks([0, 0.25, 0.5, 0.75, 1.0])
        plt.ylabel('Empirical precision')
    else:
        plt.yticks([])
    #plt.legend(loc='lower right', fontsize=fontsize)
    reader_name = reader_name_map[reader]

    print("Setting title", reader_name)
    plt.title(reader_name, fontsize=fontsize)
    plt.subplots_adjust(left=0.3, right=0.94, bottom=0.23, top=0.826)
    ax = plt.gca()
    format_axis(ax)
    plt.savefig(join(output_dir, f'fig4_{reader}_curve.pdf'))


def dataset_table(
    suffix,
    allow_incomplete=False,
):
    # Load curations for all readers
    curations = {}
    stmts_by_reader = {}
    header = [''] + list(range(1, 11))
    table_rows = [header]
    for reader, rd_dict in reader_input.items():
        # Note: if the filter below is not needed, the calls
        # to load_curated_pkl_files and get_correctness_data can
        # be replaced with a call to get_curations_for_reader
        #curations[reader] = get_curations_for_reader(reader, all_stmts,
        #                                             'evidence',
        #                                             allow_incomplete=allow_incomplete)

        reader_stmts = load_curated_pkl_files(rd_dict['pkl_list'], all_stmts, reader)

        # Note: get_combined_curations uses a filter like the one below, if
        # counts from Table 2
        #filtered_stmts = []
        #for stmt in reader_stmts:
        #    agent_names = [ag.name for ag in stmt.agent_list() if
        #                   ag is not None]
        #    # Complex > 3, translocations, autophosphorylations will be skipped
        #    if len(agent_names) != 2:
        #        continue
        #    filtered_stmts.append(stmt)
        #reader_stmts = filtered_stmts

        #reader_stmts_dict = {stmt.get_hash(): stmt for stmt in reader_stmts}
        curations[reader] = \
            get_correctness_data(rd_dict['source_list'],
                             reader_stmts, aggregation='evidence',
                             allow_incomplete=allow_incomplete)
        row = [reader]
        for i in range(1, 11):
            counts = curations[reader].get(i, {})
            if not counts:
                row.append('0 (0)')
            else:
                total_curated = len(counts)
                num_correct = len([n for n in counts if n > 0])
                row.append(f'{num_correct} ({total_curated})')
        table_rows.append(row)
    with open(join(output_dir,
                    f'table2_curation_dataset_{suffix}.csv'), 'w') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerows(table_rows)


if __name__ == '__main__':
    # Load the pickle file with all assembled statements
    asmb_pkl = join(dirname(abspath(__file__)), '..', '..', '..', 'data',
                    'bioexp_asmb_preassembled.pkl')
    all_stmts = ac.load_statements(asmb_pkl)
    # Get output directory
    output_dir = sys.argv[1]
    plot_correctness_curve('reach', all_stmts, show_ylabel=True, allow_incomplete=True)
    plot_correctness_curve('sparser', all_stmts, show_ylabel=False, allow_incomplete=True)
    plot_correctness_curve('medscan', all_stmts, show_ylabel=False, allow_incomplete=True)
    # Make the dataset table 
    dataset_table('complete', allow_incomplete=False)
    dataset_table('incomplete', allow_incomplete=True)
