import sys
import csv
from os.path import join
import texttable
import numpy as np
from matplotlib import pyplot as plt
from bioexp.curation.group_curations import reader_input
from bioexp.curation.process_curations import get_curations_for_reader
from bioexp.util import set_fig_params, format_axis, fontsize


def plot_correctness_curve(reader):
    ev_correct_by_num_ev = get_curations_for_reader(reader, 'evidence')
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
    plt.figure(figsize=(2, 2), dpi=150)
    plt.errorbar(num_evs, means, yerr=std, fmt='bo-', ls='none',
                 label='Empirical mean correctness', linewidth=1, markersize=3)
    # Legends, labels, etc.
    plt.ylim(0, 1)
    plt.grid(True)
    plt.xticks(num_evs)
    plt.xlabel('Number of evidence per INDRA Statement')
    plt.ylabel('Average statement correctness')
    #plt.legend(loc='lower right', fontsize=fontsize)
    reader_name = reader[0].upper() + reader[1:]
    plt.title('Statement Correctness', fontsize=fontsize)
    plt.subplots_adjust(left=0.15, right=0.94)
    ax = plt.gca()
    format_axis(ax)
    plt.savefig(join(output_dir, f'fig4_{reader}_curve.pdf'))


def dataset_table():
    # Load curations for all readers
    curations = {}
    stmts_by_reader = {}
    header = ['', list(range(1, 11))]
    table_rows = [header]
    for reader, rd_dict in reader_input.items():
        reader_stmts = load_curated_pkl_files(rd_dict['pkl_list'])
        #reader_stmts_dict = {stmt.get_hash(): stmt for stmt in reader_stmts}
        curations[reader] = get_correctness_data(rd_dict['source_list'],
                                   reader_stmts, aggregation='evidence')
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
    with open(join(output_dir, 'table1_curation_dataset.csv'), 'w') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerows(table_rows)


if __name__ == '__main__':
    # Get output directory
    output_dir = sys.argv[1]
    plot_correctness_curve('reach')
    plot_correctness_curve('sparser')
    # Make the dataset table 
    dataset_table()
