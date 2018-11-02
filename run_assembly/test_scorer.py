from indra.statements import *
from indra.tools import assemble_corpus as ac
from scorer import CuratedScorer

sparser_ev = Evidence(source_api='sparser', text='foo1')
reach_ev1 = Evidence(source_api='reach', text='reach1')
reach_ev2 = Evidence(source_api='reach', text='reach2')
medscan_ev = Evidence(source_api='medscan', text='foo3')
cs = CuratedScorer()

def test_reach_activation():
    s1 = Activation(Agent('MAP2K1'), Agent('MAPK1'), evidence=[reach_ev1])
    s2 = Activation(Agent('MAP2K1'), Agent('MAPK1'), evidence=[reach_ev2])
    # 1 stmt
    stmts = ac.run_preassembly([s1], belief_scorer=cs)
    assert stmts[0].belief == 1 - (1 - (0.41 * 0.9))
    # 2 stmts
    stmts = ac.run_preassembly([s1, s2], belief_scorer=cs)
    assert len(stmts) == 1
    assert stmts[0].belief == 1 - (1 - (0.41 * 0.9))**2


def test_reach_sparser_activation():
    s1 = Activation(Agent('MAP2K1'), Agent('MAPK1'), evidence=[reach_ev1])
    s2 = Activation(Agent('MAP2K1'), Agent('MAPK1'), evidence=[sparser_ev])
    # 2 stmts
    stmts = ac.run_preassembly([s1, s2], belief_scorer=cs)
    assert len(stmts) == 1
    assert stmts[0].belief == 1 - ((1 - (0.41 * 0.9)) * (1 - (0.45 * 0.71)))


if __name__ == '__main__':
    test_reach_activation()
    test_reach_sparser_activation()
