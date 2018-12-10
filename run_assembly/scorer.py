import numpy
from copy import deepcopy
from indra.belief import BeliefScorer, load_default_probs
from indra.statements import *

# Probabilities of correct grounding given a correct statement
REACH_GND = 0.90
SPARSER_GND_TWO = 0.71 # Stmts with two participants
SPARSER_GND_ONE = 0.86 # Stmts with one participant
MEDSCAN_GND = 0.96

precision_scores = {
                'reach':
                   {'regulation': 0.41 * REACH_GND,
                    'complex': 0.50 * REACH_GND,
                    'mod_two': 0.66 * REACH_GND,
                    'mod_one': 0.66 * REACH_GND, # Assumed from mod_two
                    'other': 0.46 * REACH_GND}, # Weighted avg precision
                'sparser':
                   {'regulation': 0.45 * SPARSER_GND_TWO,
                    'complex': 0.64 * SPARSER_GND_TWO, # Weighted avg precision
                    'mod_two': 0.60 * SPARSER_GND_TWO,
                    'mod_one': 0.68 * SPARSER_GND_ONE,
                    'other': 0.64 * SPARSER_GND_TWO}, # Weighted avg precision
                'medscan':
                   {'regulation': 0.48 * MEDSCAN_GND,
                    'complex': 0.58 * MEDSCAN_GND,
                    'mod_two': 0.52 * MEDSCAN_GND, # Weighted avg precision
                    'mod_one': 0.52 * MEDSCAN_GND, # Weighted avg precision
                    'other': 0.52 * MEDSCAN_GND} # Weighted avg precision
}


default_probs = load_default_probs()


# Set systematic errors for readers to 0 for now
local_probs = deepcopy(default_probs)
local_probs['syst']['reach'] = 0
local_probs['syst']['sparser'] = 0
local_probs['syst']['medscan'] = 0


class CuratedScorer(BeliefScorer):
    def score_statement(self, st, extra_evidence=None):
        if extra_evidence is None:
            extra_evidence = []
        # Collect random errors values for each
        all_evidence = st.evidence + extra_evidence
        # This part is the same as in the SimpleScorer of the BeliefEngine
        sources = [ev.source_api for ev in all_evidence]
        uniq_sources = numpy.unique(sources)
        syst_factors = {s: local_probs['syst'][s] for s in uniq_sources}
        rand_factors = {k: [] for k in uniq_sources}
        # Get the relevant precision category for this statement
        if isinstance(st, RegulateActivity) or isinstance(st, RegulateAmount):
            stmt_type = 'regulation'
        elif isinstance(st, Modification):
            if st.enz is None:
                stmt_type = 'mod_one'
            else:
                stmt_type = 'mod_two'
        elif isinstance(st, Translocation):
            stmt_type = 'mod_two'
        elif isinstance(st, Complex):
            stmt_type = 'complex'
        else:
            stmt_type = 'other'
        for ev in all_evidence:
            if ev.source_api in ('reach', 'sparser', 'medscan'):
                prob_incorrect = 1 - precision_scores[ev.source_api][stmt_type]
            else:
                prob_incorrect = local_probs['rand'][ev.source_api]
            rand_factors[ev.source_api].append(prob_incorrect)
        # Calculate the probability of the statement being incorrect
        neg_prob_prior = 1
        for s in uniq_sources:
            neg_prob_prior *= (syst_factors[s] + numpy.prod(rand_factors[s]))
        prob_prior = 1 - neg_prob_prior
        return prob_prior

    def check_prior_probs(self, statements):
        """Make sure the scorer has all the information needed to compute
        belief scores of each statement in the provided list, and raises an
        exception otherwise.

        Parameters
        ----------
        statements : list<indra.statements.Statement>
            List of statements to check
        """
        sources = set()
        for stmt in statements:
            sources |= set([ev.source_api for ev in stmt.evidence])
        for source in sources:
            # Check random error entries
            if source in ('reach', 'sparser', 'medscan') and \
               source not in precision_scores:
                msg = 'CurationScorer missing precision value' + \
                    ' for reader: %s' % source
                raise Exception(msg)
            elif source not in ('reach', 'sparser', 'medscan') and \
                 source not in local_probs['rand']:
                msg = 'CurationScorer missing precision value' + \
                    ' for source: %s' % source
                raise Exception(msg)
            # Check systematic error entries
            if source not in local_probs['syst']:
                msg = 'CurationScorer missing systematic error estimate ' + \
                    ' for source: %s' % source
                raise Exception(msg)
