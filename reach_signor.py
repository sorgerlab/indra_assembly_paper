"""
Compare statements extracted by REACH to Activation/Inhibition statements
in SIGNOR. Collect REACH extractions that match SIGNOR except for incorrect
polarity.
"""
import csv
import pickle
import random
from indra.sources import signor
from indra.tools import assemble_corpus as ac
from indra.statements import *
from indra.preassembler import *
from indra.preassembler.hierarchy_manager import hierarchies
from indra_db import get_primary_db, client
from indra.util import batch_iter

def _opposite_polarity(s1, s2):
    if type(s1) == type(s2):
        return False
    elif type(s1) != type(s2) and \
        s1.entities_match(s2):
        return True
    else:
        return False


def signor_ref_set(filename):
    """Get opposing and non-opposing statements from SIGNOR."""
    # Get SIGNOR statements
    sp = signor.process_from_web()

    # Filter to specific human genes only
    stmts = ac.filter_human_only(sp.statements)
    stmts = ac.filter_genes_only(stmts)
    act = ac.filter_by_type(stmts, Activation)
    inh = ac.filter_by_type(stmts, Inhibition)
    phos = ac.filter_by_type(stmts, Phosphorylation)
    dephos = ac.filter_by_type(stmts, Dephosphorylation)

    # Combine any duplicates in SIGNOR
    reg_stmts = act + inh
    phos_stmts = phos + dephos
    non_opposing = []
    for stmt_grp in (reg_stmts, phos_stmts):
        pa = Preassembler(hierarchies, stmt_grp)
        uniq = pa.combine_duplicates()

        # Compare all pairs of statements in the combined pool
        opposing = []
        combos = list(itertools.combinations(uniq, 2))
        for ix, (s1, s2) in enumerate(combos):
            # Progress msg
            if (ix + 1) % 1000000 == 0:
                print("%d of %d" % (ix+1, len(combos)))
            # If the statements are of opposing polarity, save them as a pair
            if _opposite_polarity(s1, s2):
                opposing.append((s1, s2))
                print(str((s1, s2)))
        # Flatten pairs of opposing statements into a single list of statements
        opposing_flat = list(set([s for op_tuple in opposing for s in op_tuple]))
        # Now get the statements that are *NOT* opposing--these are the reference
        # set, because they suggest that there is no controversy about the correct
        # polarity
        non_opposing += [s for s in uniq if s not in opposing_flat]
        # Now dump the non-opposing statements to a file
    ac.dump_statements(non_opposing, filename)


def get_matching_db_stmts(source_stmts, filename):
    # Make synthetic statements with inverted polarity
    opp_stmts = []
    for s in source_stmts:
        if isinstance(s, Activation):
            opp_stmts.append(Inhibition(s.subj, s.obj))
        elif isinstance(s, Inhibition):
            opp_stmts.append(Activation(s.subj, s.obj))
        elif isinstance(s, Modification):
            opposing_class = modclass_to_inverse[type(s)]
            opp_stmts.append(opposing_class(s.enz, s.sub, s.residue,
                             s.position))
    # Make hash lists
    match_hash_list = [s.get_hash(shallow=True) for s in source_stmts]
    opp_hash_list = [s.get_hash(shallow=True) for s in opp_stmts]
    # Get matching and opposing statements from the DB
    db = get_primary_db()
    # Function to run for both matching and opposing stmts
    def get_corresponding_stmts(hash_list, batch_size=50):
        stmts = []
        for ix, batch in enumerate(batch_iter(hash_list, batch_size)):
            print(f"Loading batch {ix}")
            result = client.get_statement_jsons_from_hashes(batch,
                                                            ev_limit=100000)
            jsons = list(result['statements'].values())
            stmts += stmts_from_json(jsons)
        return stmts
    # Run the queries
    print(f"Querying for matching JSONs for {len(match_hash_list)} statements")
    match_stmts = get_corresponding_stmts(match_hash_list)
    print(f"Querying for opposing JSONs for {len(opp_hash_list)} statements")
    opp_stmts = get_corresponding_stmts(opp_hash_list)
    # Save results in a single dict
    result = {'matching': match_stmts, 'opposing': opp_stmts}
    with open(filename, 'wb') as f:
        pickle.dump(result, f)
    return result


def dump_sentences(db_stmts, filename, reader=None):
    header = ['Statement', 'Source', 'SubjText', 'ObjText', 'PMID', 'Rule',
              'StmtHash', 'Sentence']
    rows = [header]
    for s in db_stmts:
        for ev in s.evidence:
            if reader is not None and ev.source_api != reader:
                continue
            subj_text = ev.annotations['agents']['raw_text'][0]
            obj_text = ev.annotations['agents']['raw_text'][1]
            row = [str(s), ev.source_api, subj_text, obj_text,
                   ev.pmid, ev.annotations.get('found_by'),
                   s.get_hash(shallow=True), ev.text]
            rows.append(row)
    print(f"Dumping {len(rows)} sentences")
    with open(filename, 'wt') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerows(rows)


def train_naive_bayes(db_stmts_dict):
    from sklearn.naive_bayes import MultinomialNB
    random.seed(1)
    random.shuffle(non_opposing)
    stmts_match = []
    stmts_not_match = []
    matching_ct = 0
    nonmatching_ct = 0
    m_ev = []
    nm_ev = []

    # is activation
    # is inhibition
    # has_mod
    # has_complex
    # has_inc_amount
    # has_dec_amount
    # opposite polarity greater
    db_stmts_act = ac.filter_by_type(db_stmts, Activation)
    db_stmts_inh = ac.filter_by_type(db_stmts, Inhibition)
    act_reach = [e for s in db_stmts_act for e in s.evidence
                 if e.source_api == 'reach']
    inh_reach = [e for s in db_stmts_inh for e in s.evidence
                 if e.source_api == 'reach']
    if isinstance(s, Activation):
        if db_stmts_act:
            stmts_match.append((s, db_stmts_act))
        if db_stmts_inh:
            stmts_not_match.append((s, db_stmts_inh))
        if act_reach:
            m_ev.extend(act_reach)
        if inh_reach:
            nm_ev.extend(inh_reach)
        matching_ct += len(act_reach)
        nonmatching_ct += len(inh_reach)
    elif isinstance(s, Inhibition):
        if db_stmts_act:
            stmts_not_match.append((s, db_stmts_act))
        if db_stmts_inh:
            stmts_match.append((s, db_stmts_inh))
        if act_reach:
            nm_ev.extend(act_reach)
        if inh_reach:
            m_ev.extend(inh_reach)
        matching_ct += len(inh_reach)
        nonmatching_ct += len(act_reach)

    # Preprocess the data into a vector of counts for each evidence type,
    # plus a matching/opposing (correct/incorrect) class identifier
    data = list(zip([e.text for e in m_ev], [0] * len(m_ev))) + \
           list(zip([e.text for e in nm_ev], [1] * len(nm_ev)))
    random.shuffle(data)
    corpus, classes = zip(*data)
    y = np.array(classes)
    #cv = CountVectorizer()
    cv = HashingVectorizer(n_features=1000)
    #cv = TfidfVectorizer()
    X = cv.fit_transform(corpus)
    X_train = X[0:4000]
    y_train = y[0:4000]
    X_test = X[4000:]
    y_test = y[4000:]
    clf = MultinomialNB()
    clf.fit(X_train, y_train)
    train_score = clf.score(X_train, y_train)
    test_score = clf.score(X_test, y_test)
    (train_score, test_score)


if __name__ == '__main__':
    # RELOAD?
    reload_signor = False
    signor_nonopp_file = 'signor_non_opposing.pkl'
    if reload_signor:
        signor_stmts = signor_ref_set(signor_nonopp_file)
    else:
        signor_stmts = ac.load_statements(signor_nonopp_file)

    # RELOAD?
    reload_db = False
    db_stmts_file = 'db_stmts_dict.pkl'
    if reload_db == True:
        db_stmts = get_matching_db_stmts(signor_stmts, db_stmts_file)
    else:
        with open(db_stmts_file, 'rb') as f:
            db_stmts = pickle.load(f)
    #print("Dumping opposing sentences")
    #dump_sentences(db_stmts['opposing'], 'opposing_reach_stmts.csv',
                    #reader='reach')
    #print("Dumping concurring sentences")
    #dump_sentences(db_stmts['matching'], 'concurring_reach_stmts.csv',
                    #reader='reach')

    belief_mek_erk_example(db_stmts)





