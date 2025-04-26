import pytest

from cogent3 import available_datasets, get_dataset


def test_available_data_sets():
    datasets = available_datasets()
    assert datasets.shape[0] > 1


def test_get_dataset_matched():
    aln = get_dataset("brca1")
    assert len(aln) > 3000
    tree = get_dataset("mammal-tree")
    assert aln.num_seqs == len(tree.get_tip_names())
    assert set(aln.names) == set(tree.get_tip_names())


def name_type():
    datasets = available_datasets()
    return datasets.to_list(columns=["name", "type"])


@pytest.mark.parametrize(("name", "type_"), name_type())
def test_dataset_name_type(name, type_):
    dset = get_dataset(name)
    assert dset.__class__.__name__ == type_


def test_invalid_name():
    with pytest.raises(ValueError):
        get_dataset("invalid-dataset-name")
