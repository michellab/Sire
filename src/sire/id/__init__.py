
__all__ = ["CaseInsensitive",
           "CaseSensitive"]

from ..legacy import ID as _ID

from .. import use_new_api as _use_new_api
_use_new_api()

CaseSensitive = _ID.CaseSensitive
CaseInsensitive = _ID.CaseInsensitive
