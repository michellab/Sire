"""
.. currentmodule:: sire.search

"""

from ..legacy import Search as _Search

from .. import use_new_api as _use_new_api
_use_new_api()

from ..legacy.Search import approx_equal, \
    set_approx_epsilon, get_approx_epsilon, \
    get_min_protein_residues, set_min_protein_residues, \
    get_protein_residue_names, set_protein_residue_names, \
    set_token, get_token, delete_token, delete_all_tokens
