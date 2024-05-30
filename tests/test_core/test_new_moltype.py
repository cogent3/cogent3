import pytest

from cogent3.core import moltype, new_moltype


def test_make_pairs():
    orig = moltype.make_pairs(
        pairs=moltype.RnaStandardPairs,
        monomers=moltype.IUPAC_RNA_chars,
        gaps=moltype.IUPAC_gap,
        degenerates=moltype.IUPAC_RNA_ambiguities,
    )
    new = new_moltype.make_pairs(
        pairs=new_moltype.RNA_STANDARD_PAIRS,
        monomers=new_moltype.IUPAC_RNA_chars,
        gaps=new_moltype.IUPAC_gap,
        degenerates=new_moltype.IUPAC_RNA_ambiguities,
    )
    # convert the old object to a frozen set keyed one for comparison
    orig = {frozenset(k): v for k, v in orig.items()}
    assert new == orig


def test_is_compatible_alphabet():
    from cogent3.core.new_alphabet import CharAlphabet

    dna = new_moltype.get_moltype("dna")
    alpha = CharAlphabet("TCAG")
    assert dna.is_compatible_alphabet(alpha)
    rna = new_moltype.get_moltype("rna")
    assert not rna.is_compatible_alphabet(alpha)
    alpha = CharAlphabet("".join(dna.ambiguities))
    prot = new_moltype.get_moltype("protein")
    assert not prot.is_compatible_alphabet(alpha)


def test_is_compatible_alphabet_strict():
    from cogent3.core.alphabet import CharAlphabet

    dna = new_moltype.get_moltype("dna")
    alpha1 = CharAlphabet("TCAG")
    assert dna.is_compatible_alphabet(alpha1, strict=True)
    # returns False if the order is not exactly the same
    alpha1 = CharAlphabet("CTAG")
    assert not dna.is_compatible_alphabet(alpha1, strict=True)


@pytest.mark.parametrize(
    "name", ("dna", "rna", "protein", "protein_with_stop", "bytes", "text")
)
def test_get_moltype(name):
    """correctly return a moltype by name"""
    mt = new_moltype.get_moltype(name)
    assert mt.name == name
    mt = new_moltype.get_moltype(name.upper())
    assert mt.name == name
    got = new_moltype.get_moltype(mt)
    assert got is mt


def test_available_moltypes():
    t = new_moltype.available_moltypes()
    assert t.shape[0] == 6


def test_str_moltype():
    dna = new_moltype.get_moltype("dna")
    text = str(dna)
    assert isinstance(text, str)
    assert text == f"MolType({tuple('TCAG')})"


@pytest.mark.parametrize("seq", ("ACCCG", b"ACCCG"))
@pytest.mark.parametrize("name", ("dna", "rna"))
def test_complement(name, seq):
    dna = new_moltype.get_moltype(name)
    expect = "TGGGC" if name == "dna" else "UGGGC"
    got = dna.complement(seq)
    assert got == expect
