"""
Compare statements extracted by REACH to Activation/Inhibition statements
in SIGNOR. Collect REACH extractions that match SIGNOR except for incorrect
polarity.
"""
import random
from indra.sources import signor
from indra.tools import assemble_corpus as ac
from indra.statements import *
from indra.preassembler import *
from indra.preassembler.hierarchy_manager import hierarchies
#from indra.sources.indra_db_rest import get_statements
from indra.db import get_primary_db, client

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
    act = ac.filter_by_type(stmts, Activation)
    inh = ac.filter_by_type(stmts, Inhibition)

    # Combine any duplicates in SIGNOR
    pa = Preassembler(hierarchies, act + inh)
    uniq = pa.combine_duplicates()

    # Compare all pairs of statements in the combined pool of Act/Inh
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
    non_opposing = [s for s in uniq if s not in opposing_flat]
    # Now dump the non-opposing statements to a file
    ac.dump_statements(non_opposing, filename)


def get_matching_db_stmts(source_stmts, filename):
    all_db_stmts = []
    db = get_primary_db()
    hash_list = [s.get_hash(shallow=True) for s in source_stmts]
    stmts = client.get_statements([db.PAStatements.mk_hash.in_(hash_list)])
    ac.dump_statements(stmts, filename)
    return stmts

    """
    client.get_statements([
    for ix, s in enumerate(signor_stmts):
        if ix % 10 == 0:
            print(ix)
        if ix % 100 ==0:
            print("Saving...")
            ac.dump_statements(all_db_stmts, '%s_%d.pkl' % (filename, ix))
        subj = s.agent_list()[0]
        obj = s.agent_list()[1]
        subj_hgnc = subj.db_refs.get('HGNC')
        obj_hgnc = obj.db_refs.get('HGNC')
        if subj_hgnc is None or obj_hgnc is None:
            continue
        db_stmts = get_statements(subject=subj.name, object=obj.name)
        all_db_stmts.append((s, db_stmts))
    ac.dump_statements(all_db_stmts, '%s_all.pkl' % filename)
    return all_db_stmts
    """

if __name__ == '__main__':
    # RELOAD?
    reload_signor = False
    signor_nonopp_file = 'signor_non_opposing_regulate.pkl'
    if reload_signor:
        signor_stmts = signor_ref_set(signor_nonopp_file)
    else:
        signor_stmts = ac.load_statements(signor_nonopp_file)

    # RELOAD?
    reload_db = True
    db_stmts_file = 'non_opp_db_stmts.pkl'
    if reload_db == True:
        db_stmts = get_matching_db_stmts(signor_stmts, db_stmts_file)
    else:
        db_stmts = ac.load_statements(db_stmts_file)



"""
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
"""









"""
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer, \
                                            HashingVectorizer
from sklearn.naive_bayes import MultinomialNB

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
"""

# For each statement where there is a DB stmt with a polarity that doesn't
# match the polarity of the SIGNOR statement, create an entry marked "incorrect"


#contra = pa.find_contradicts()

# We get statement from REACH; we run it through our classifier and get a
# probability estimate of the statement being correct
# Features:
# - length of sentence,
# - num of words between entities
# - do we know if there is a phosphorylation or binding event between the two
#   proteins?
# - number of evidences
# - number of opposing evidences from reading
# - ratio between supporting and opposing evidences
# Stemmed word vector?

# Of course, could use the sentences to create a training corpus 

# What fraction of the time do we get statements from the database 
# Purpose: predict when the statement is likely to be incorrect, at least
# due to polarity

