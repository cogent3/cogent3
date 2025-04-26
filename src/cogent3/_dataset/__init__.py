import dataclasses
import pathlib
import typing

if typing.TYPE_CHECKING:
    from cogent3.core.new_alignment import Alignment
    from cogent3.core.table import Table
    from cogent3.core.tree import PhyloNode

DATASET_DIR = pathlib.Path(__file__).parent

_brca1_primates = (
    "FlyingLem",
    "TreeShrew",
    "Galago",
    "HowlerMon",
    "Rhesus",
    "Orangutan",
    "Gorilla",
    "Chimpanzee",
    "Human",
)


def _mammal_tree(tree: "PhyloNode") -> "PhyloNode":
    tip_names = set(tree.get_tip_names()) - {"Chook"}
    return tree.get_sub_tree(tip_names)


def _select_primate_tree(tree: "PhyloNode") -> "PhyloNode":
    return tree.get_sub_tree(_brca1_primates)


def _select_primate_brca1(aln: "Alignment") -> "Alignment":
    return aln.take_seqs(_brca1_primates).omit_gap_pos()


@dataclasses.dataclass
class Dataset:
    name: str
    filename: str
    cogent3_type: str
    loader_func_name: str
    load_args: dict | None
    description: str
    post_process: typing.Callable | None = None

    def load(self) -> typing.Any:  # noqa: ANN401
        """Load the dataset using the provided loader."""
        import cogent3

        func = getattr(cogent3, self.loader_func_name)
        kwargs = self.load_args or {}
        data = func(DATASET_DIR / self.filename, **kwargs)

        if self.post_process:
            data = self.post_process(data)

        return data


_aln_kwargs = {"moltype": "dna", "new_type": True}

_datasets = [
    Dataset(
        name="brca1",
        cogent3_type="Alignment",
        filename="brca1.fasta",
        loader_func_name="load_aligned_seqs",
        load_args=_aln_kwargs,
        description="Multiple sequence alignment of mammal BRCA1 sequences.",
    ),
    Dataset(
        name="mammal-tree",
        cogent3_type="PhyloNode",
        filename="murphy.tree",
        loader_func_name="load_tree",
        load_args=None,
        description="Unrooted phylogeny for mammal data in brca1 dataset.",
        post_process=_mammal_tree,
    ),
    Dataset(
        name="murphy-tree",
        cogent3_type="PhyloNode",
        filename="murphy.tree",
        loader_func_name="load_tree",
        load_args=None,
        description="Unrooted phylogeny for mammal data with chicken from Murphy et al 2001.",
    ),
    Dataset(
        name="primate-brca1",
        filename="brca1.fasta",
        cogent3_type="Alignment",
        loader_func_name="load_aligned_seqs",
        load_args=_aln_kwargs,
        description="Alignment of BRCA1 for the primates.",
        post_process=_select_primate_brca1,
    ),
    Dataset(
        name="primate-tree",
        cogent3_type="PhyloNode",
        filename="murphy.tree",
        loader_func_name="load_tree",
        load_args=None,
        description="Unrooted phylogeny for primate data.",
        post_process=_select_primate_tree,
    ),
]

_dataset_mapping = {ds.name: ds for ds in _datasets}


def get_dataset(name: str) -> typing.Any:  # noqa: ANN401
    """returns the named dataset as cogent3 type"""
    ds = _dataset_mapping.get(name)
    if ds is None:
        msg = f"Dataset '{name!r}' not found."
        raise ValueError(msg)
    return ds.load()


def available_datasets() -> "Table":
    """Return table of the available datasets."""
    import cogent3

    data = [[d.name, d.cogent3_type, d.description] for d in _datasets]
    table = cogent3.make_table(
        header=["name", "type", "description"],
        data=data,
        title="Available sample data sets.",
    )
    table.set_repr_policy(show_shape=False)
    return table
