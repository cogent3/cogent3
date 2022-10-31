#!/usr/bin/env python
"""
This takes doctest files and turns them into standalone scripts.
"""
import pathlib
import re

import click


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__contributors__ = ["Gavin Huttley", "Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

not_wsp = re.compile(r"^\S+")
code_eval = re.compile(r"(?<=\s{4})([>]{3}|[.]{3})\s")
code_block = re.compile(r".. (jupyter-execute|doctest)::")
block_option = re.compile(r"\s+:[a-z\-]+:")
raise_option = re.compile(r"\s+:raises:")


def get_error_type(line):
    if raise_option.search(line) is None:
        return None
    return raise_option.sub("", line).strip()


def get_path_update():
    """returns code block to allow import set_working_directory"""
    swd_script = pathlib.Path("set_working_directory.py").absolute()
    assert swd_script.exists()
    rootdir = str(swd_script.parent)
    block = ["import sys", f"sys.path.append({rootdir!r})", "import os", f"os.chdir({rootdir!r})"]
    return "\n".join(block)


def deindent(line):
    if line.startswith(" " * 4):
        line = line[4:]
    return line


def get_code_block_line_numbers(doc):
    """returns the (start, end) of codeblock sections"""
    lines = []
    in_code_block = False
    start = None
    for i, line in enumerate(doc):
        if code_block.search(line):
            if in_code_block:
                lines.append((start, i - 1))
            start = i
            in_code_block = True
            continue

        if in_code_block and not_wsp.search(line):
            lines.append((start, i))
            in_code_block = False
            continue

    if in_code_block:
        if i == len(doc) - 1:
            i += 1
        lines.append((start, i))

    return lines


def format_block(block):
    """handles exceptions, de-indent, etc..."""
    error_type = get_error_type(block[0])
    format_line = (lambda x: x) if error_type else deindent
    code = [format_line(l) for l in block if not block_option.search(l)]
    if error_type:
        code.insert(0, "try:")
        code.extend([f"except {error_type}:", "    pass"])
    return code


def get_code_blocks(doc: list[str]) -> str:
    coords = get_code_block_line_numbers(doc)
    refactored = [get_path_update()]
    for start, end in coords:
        code = format_block(doc[start + 1 : end])
        refactored.extend([""] + code)

    return "\n".join(refactored)


def _rst_path(*args):
    path = pathlib.Path(args[-1])
    assert path.suffix == ".rst"
    return path


@click.command(no_args_is_help=True, context_settings={"show_default": True})
@click.argument("rst_path", callback=_rst_path)
@click.option("-t", "--test", is_flag=True, help="don't write script, print it")
def main(rst_path, test):
    """extracts code under jupyter_execute or doctest blocks to pythomn script"""
    outpath = rst_path.parent / f"{rst_path.stem}.py"

    doc = rst_path.read_text()
    doc = doc.splitlines()

    just_code = get_code_blocks(doc)
    if test:
        print(just_code)
        exit()

    outpath.write_text(just_code)


if __name__ == "__main__":
    main()
