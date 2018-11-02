"""This script does additional Statement assembly specific to PySB models,
and then assembles and contextualizes a PySB model from the Statements."""

from indra.util import _require_python3
import re
from pysb import Observable, ReactionPattern, ComplexPattern
from indra.statements import *
from indra.mechlinker import MechLinker
import indra.tools.assemble_corpus as ac
from indra.assemblers.pysb import PysbAssembler
import sys
sys.path.append('..')
sys.path.append('../run_assembly')
import process_data
from util import *
from run_assembly.util import pklload, pkldump, prefixed_file


def assemble_pysb(stmts, data_genes):
    """Return an assembled PySB model."""
    # Save a version of statements with no evidence for faster loading
    #for s in stmts:
    #    s.evidence = []
    #    for ss in s.supports + s.supported_by:
    #        ss.evidence = []
    #pkldump(stmts, 'no_evidence')

    # Assemble model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model(reverse_effects=False)
    # Set context
    set_context(pa)
    # Add observables
    add_observables(pa.model)
    pa.save_model(prefixed_file('pysb', 'py'))
    pkldump(pa.model, 'pysb')
    return pa.model


def preprocess_stmts(stmts, data_genes):
    """Preprocess Statements specifically for PySB assembly."""
    # Filter the INDRA Statements to be put into the model
    stmts = ac.filter_mutation_status(stmts,
                                      {'BRAF': [('V', '600', 'E')]}, ['PTEN'])
    stmts = ac.filter_by_type(stmts, Complex, invert=True)
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_gene_list(stmts, data_genes, 'all')
    # Simplify activity types
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.reduce_activities()
    ml.gather_modifications()
    ml.reduce_modifications()
    af_stmts = ac.filter_by_type(ml.statements, ActiveForm)
    non_af_stmts = ac.filter_by_type(ml.statements, ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(af_stmts)
    stmts = af_stmts + non_af_stmts
    # Replace activations when possible
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.replace_activations()
    # Require active forms
    ml.require_active_forms()
    num_stmts = len(ml.statements)
    while True:
        # Remove inconsequential PTMs
        ml.statements = ac.filter_inconsequential_mods(ml.statements,
                                                       get_mod_whitelist())
        ml.statements = ac.filter_inconsequential_acts(ml.statements,
                                                       get_mod_whitelist())
        if num_stmts <= len(ml.statements):
            break
        num_stmts = len(ml.statements)
    stmts = ml.statements
    return stmts


def set_context(pa):
    """Set expression amounts for the cell line and make sure BRAF is
    in a mutated form."""
    pa.set_context('SKMEL28_SKIN')
    # Set BRAF V600E
    for ic in pa.model.initial_conditions:
        if str(ic[0]).startswith('BRAF'):
            ic[0].monomer_patterns[0].site_conditions['V600'] = 'E'


def add_observables(model):
    """Get the antibody targets and add observables for them to the model."""
    data = process_data.read_data()
    ab_map = process_data.get_antibody_map(data)
    for ab_name, agents in ab_map.items():
        patterns = []
        for agent in agents:
            try:
                monomer = model.monomers[agent.name]
            except KeyError:
                continue
            if agent.mods:
                mc = agent.mods[0]
                site_names = ['phospho', mc.residue]
                if mc.position is not None:
                    site_names.append(mc.residue + mc.position)
                for site_name in site_names:
                    try:
                        pattern = monomer(**{site_name: 'p'})
                        patterns.append(ComplexPattern([pattern], None))
                    except Exception:
                        pass
            else:
                patterns.append(ComplexPattern([monomer()], None))
        if patterns:
            if model.monomers.get(ab_name) is not None:
                obs_name = ab_name + '_obs'
            else:
                obs_name = ab_name
            if not re.match(r'[_a-z][_a-z0-9]*\Z', obs_name, re.IGNORECASE):
                obs_name = obs_name.replace('-', '_')
            if not re.match(r'[_a-z][_a-z0-9]*\Z', obs_name, re.IGNORECASE):
                obs_name = 'p' + obs_name
            o = Observable(obs_name, ReactionPattern(patterns))
            model.add_component(o)
    '''
    o = Observable(b'MAPK1p', model.monomers['MAPK1'](T185='p', Y187='p'))
    model.add_component(o)
    o = Observable(b'MAPK3p', model.monomers['MAPK3'](T202='p', Y204='p'))
    model.add_component(o)
    o = Observable(b'MAPK14p', model.monomers['MAPK14'](T180='p'))
    model.add_component(o)
    o = Observable(b'GSK3Ap', model.monomers['GSK3A'](S21='p'))
    model.add_component(o)
    o = Observable(b'GSK3Bp', model.monomers['GSK3B'](S9='p'))
    model.add_component(o)
    o = Observable(b'RPS6pS235', model.monomers['RPS6'](S235='p'))
    model.add_component(o)
    o = Observable(b'RPS6pS240', model.monomers['RPS6'](S240='p'))
    model.add_component(o)
    o = Observable(b'EIF4EBP1p', model.monomers['EIF4EBP1'](S65='p'))
    model.add_component(o)
    o = Observable(b'JUNp', model.monomers['JUN'](S73='p'))
    model.add_component(o)
    o = Observable(b'FOXO3p', model.monomers['FOXO3'](S315='p'))
    model.add_component(o)
    o = Observable(b'AKT1p', model.monomers['AKT1'](S473='p'))
    model.add_component(o)
    o = Observable(b'AKT2p', model.monomers['AKT2'](S474='p'))
    model.add_component(o)
    o = Observable(b'AKT3p', model.monomers['AKT3'](phospho='p'))
    model.add_component(o)
    o = Observable(b'ELK1p', model.monomers['ELK1'](S383='p'))
    model.add_component(o)
    o = Observable(b'RB1p', model.monomers['RB1'](S807='p'))
    model.add_component(o)
    o = Observable(b'RPS6KA1p', model.monomers['RPS6KA1'](T359='p'))
    model.add_component(o)
    o = Observable(b'RPS6KB1p', model.monomers['RPS6KB1'](phospho='p'))
    model.add_component(o)
    o = Observable(b'PDPK1p', model.monomers['PDPK1'](S241='p'))
    model.add_component(o)
    o = Observable(b'PTK2p', model.monomers['PTK2'](Y397='p'))
    model.add_component(o)
    o = Observable(b'STAT3p', model.monomers['STAT3'](S727='p'))
    model.add_component(o)
    o = Observable(b'IRS1p', model.monomers['IRS1'](S307='p'))
    model.add_component(o)
    o = Observable(b'ESR1p', model.monomers['ESR1'](S118='p'))
    model.add_component(o)
    '''


def get_mod_whitelist():
    """Return a list of modifications that should not be filtered out

    These sites should be preserved in the model so that they can be
    observed, even if they are otherwise considered inconsequential.
    """
    mod_whitelist = {}
    # Here we take the targets of each antybody to make sure we
    # don't remove these observables from the model
    ab_map = process_data.get_phospho_antibody_map()
    for k, v in ab_map.items():
        for agent in v:
            mod = ('phosphorylation', agent.mods[0].residue,
                   agent.mods[0].position)
            try:
                mod_whitelist[agent.name].append(mod)
            except KeyError:
                mod_whitelist[agent.name] = [mod]
            # Add generic mods to make sure we keep them
            mod = ('phosphorylation', agent.mods[0].residue, None)
            mod_whitelist[agent.name].append(mod)
            mod = ('phosphorylation', None, None)
            mod_whitelist[agent.name].append(mod)
    return mod_whitelist


if __name__ == '__main__':
    data = process_data.read_data()
    data_genes = process_data.get_all_gene_names(data)

    cmd = sys.argv[1]
    cmds = ('preprocess_stmts', 'assemble_pysb')
    if cmd in cmds:
        input_file = sys.argv[2]
        output_file = sys.argv[3]
    # Break each assembly cmd down into intermediate steps
    if cmd == 'preprocess_stmts':
        stmts_in = pklload(input_file)
        stmts_out = preprocess_stmts(stmts_in, data_genes)
        pkldump(stmts_out, output_file)
    elif cmd == 'assemble_pysb':
        stmts_in = pklload(input_file)
        model = assemble_pysb(stmts_in, data_genes)
        pkldump(model, output_file)
