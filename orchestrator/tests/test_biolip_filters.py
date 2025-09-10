from orchestrator.curation.ccd import is_common_additive
from orchestrator.curation.biolip import BiolipLigand, select_relevant_ligand


def test_is_common_additive_basic():
    assert is_common_additive("HOH") is True
    assert is_common_additive("EDO") is True
    assert is_common_additive("ATP") is False


def test_select_relevant_ligand():
    ligs = [
        BiolipLigand(ligand_id="ATP", chain="A", contacts=["A:10", "A:11"], notes={}),
        BiolipLigand(ligand_id="GOL", chain="A", contacts=[], notes={}),
    ]
    sel = select_relevant_ligand(ligs)
    assert sel is not None
    assert sel.ligand_id == "ATP"

