from os.path import join, dirname
from collections import Counter, defaultdict
from matplotlib import pyplot as plt
from bioexp.util import set_fig_params, fontsize, format_axis, red, blue, \
                        green, purple, black, yellow, orange, gray
from indra.util import write_unicode_csv


def plot_steps_before_preassembly(data):
    # The sequential order of the assembly steps, with plot labels
    step_order = (
        ('all_raw', 'Raw Statements'),
        ('filter_no_hypothesis', 'Filter out hypotheses'),
        ('filter_grounded_only', 'Filter ungrounded'),
        ('filter_genes_only', 'Filter to genes only'),
        ('filter_human_only', 'Filter to human only'),
        ('asmb_map_sequence', 'Map sites of PTMs'),
        )
    # Put source of greatest numbers of statements at the bottom of the stack
    reader_sources = ['reach', 'sparser', 'medscan', 'other_reader']
    other_readers = ['isi', 'trips', 'rlimsp']
    db_sources = ['biopax', 'bel', 'signor', 'other_db']
    other_dbs = ['hprd', 'phosphosite', 'cbn', 'trrust']
    source_order = reader_sources + db_sources
    # Colors for each source type
    source_colors = {'reach': red,
                     'sparser': blue,
                     'medscan': green,
                     'other_reader': black,
                     'biopax': orange,
                     'bel': purple,
                     'signor': yellow,
                     'other_db': gray}

    # Make the plot
    fig = plt.figure(figsize=(1.5, 2), dpi=150)
    for ix, (step, step_label) in enumerate(step_order):
        # Plot distribution of statement sources for steps prior to preassembly
        last_val = 0
        for source in source_order:
            col = source_colors[source]
            if source == 'other_db':
                count = sum(int(data[step][(db,)]) for db in other_dbs)
            elif source == 'other_reader':
                count = sum(int(data[step][(rd,)]) for rd in other_readers)
            else:
                count = int(data[step][(source,)])
            plt.bar(ix, count, bottom=last_val, color=col)
            last_val += count
    _, step_labels = zip(*step_order)
    plt.xticks(range(len(step_labels)), step_labels, rotation='vertical')
    plt.ylabel(r'Statements ($\times10^{6}$)')
    ax = fig.gca()
    format_axis(ax, tick_padding=2)
    plt.subplots_adjust(left=0.18, bottom=0.53, right=0.96, top=0.97)
    plt.yticks([0, 2000000, 4000000, 6000000, 8000000],
               [0, 2, 4, 6, 8])
    plt.show()
    filename = join(build_dir, 'fig2_stmt_counts_before_pa.pdf')
    plt.savefig(filename)


def plot_steps_after_preassembly(data):
    # The sequential order of the assembly steps, with plot labels
    step_order = (
        ('preassembled', 'Combine duplicates'),
        ('filter_belief', r'Filter belief $>$ 0.95'),
        ('filter_top_level', 'Filter to top-level'),
        ('filter_mod_nokinase', 'Filter phos. to kinases'),
        ('reduce_activities', 'Reduce activities'),
        ('reduce_mods', 'Reduce modifications'),
        )
    db_sources = set(['biopax', 'bel', 'signor'])
    reading_sources = set(['reach', 'sparser'])
    # Put source of greatest numbers of statements at the bottom of the stack
    # Make the plot
    fig = plt.figure(figsize=(1, 2), dpi=150)
    rows = []
    for ix, (step, step_label) in enumerate(step_order):
        last_val = 0
        # Source list is a list of tuples, with each tuple representing a
        # combination of individual sources
        source_list = data[step]
        reading_only = 0
        db_only = 0
        both = 0
        for sources, count in source_list.items():
            has_db = False
            has_reading = False
            count = int(count)
            # See if we have db, reading, or both
            if sources == ('total',):
                continue
            if set(sources).intersection(db_sources):
                has_db = True
            if set(sources).intersection(reading_sources):
                has_reading = True
            # Add to the appropriate counters
            if has_db and has_reading:
                both += count
            elif has_db and not has_reading:
                db_only += count
            elif not has_db and has_reading:
                reading_only += count
        plt.bar(ix, reading_only, bottom=0, color=red)
        plt.bar(ix, db_only, bottom=reading_only, color=blue)
        plt.bar(ix, both, bottom=db_only+reading_only, color=purple)
        print(step_label)
        print("Reading", reading_only)
        print("DB", db_only)
        print("Both", both)
        rows.append((step, reading_only, db_only, both))

    _, step_labels = zip(*step_order)
    plt.xticks(range(len(step_labels)), step_labels, rotation='vertical')
    plt.ylabel(r'Statements ($\times10^{-5}$)')
    ax = fig.gca()
    plt.yticks([0, 200000, 400000, 600000],
               [0, 2, 4, 6])
    format_axis(ax, tick_padding=2)
    plt.subplots_adjust(left=0.18, bottom=0.53, right=0.76, top=0.97)
    plt.show()
    fig_filename = join(build_dir, 'fig2_stmt_counts_after_pa.pdf')
    plt.savefig(fig_filename)
    csv_filename = join(build_dir, 'fig2_stmt_sources_after_pa.csv')
    write_unicode_csv(csv_filename, rows)


if __name__ == '__main__':
    data_dir = join(dirname(__file__), '..', '..', '..', 'data')
    build_dir = join(dirname(__file__), '..', '..', '..', 'output')

    # Read the statement source data from the TSV file
    data = defaultdict(dict)
    tsv_file = join(build_dir, 'fig2_stmt_counts.txt')
    with open(tsv_file, 'rt') as f:
        for line in f.readlines():
            entries = line.strip().split('\t')
            assembly_step = entries[0]
            sources = tuple(entries[1].split(','))
            count = entries[2]
            data[assembly_step][sources] = count

    plt.ion()
    set_fig_params()
    # Make plot for steps prior to preassembly, where bars are stratified
    # by individual sources
    plot_steps_before_preassembly(data)

    # Make plot for steps after preassembly, where bars are stratified
    # by combinations of sources
    # plot_steps_after_preassembly(data)

