diff --git a/src/cogent3/core/sequence.py b/src/cogent3/core/sequence.py
index 48434579b..831a21bec 100644
--- a/src/cogent3/core/sequence.py
+++ b/src/cogent3/core/sequence.py
@@ -20,11 +20,12 @@ import os
 import re
 import warnings
 
+from abc import ABC, abstractmethod
 from collections import defaultdict
 from functools import singledispatch, total_ordering
 from operator import eq, ne
 from random import shuffle
-from typing import Any, Generator, Iterable, List, Optional, Tuple
+from typing import Any, Generator, Iterable, List, Optional, Tuple, Union
 
 from numpy import (
     arange,
@@ -43,8 +44,6 @@ from numpy import (
 )
 from numpy.random import permutation
 
-import cogent3.util.warning as c3warn
-
 from cogent3._version import __version__
 from cogent3.core.alphabet import AlphabetError
 from cogent3.core.annotation import Feature
@@ -61,6 +60,7 @@ from cogent3.core.location import FeatureMap, IndelMap, LostSpan
 from cogent3.format.fasta import alignment_to_fasta
 from cogent3.maths.stats.contingency import CategoryCounts
 from cogent3.maths.stats.number import CategoryCounter
+from cogent3.util import warning as c3warn
 from cogent3.util.dict_array import DictArrayTemplate
 from cogent3.util.misc import (
     DistanceFromMatrix,
@@ -1037,7 +1037,7 @@ class Sequence(SequenceI):
         We assume that spans represent the coordinates for this instance!
         """
         feature = dict(feature)
-        seq_rced = self._seq.reverse
+        seq_rced = self._seq.is_reversed
         spans = feature.pop("spans", None)
         revd = feature.pop("strand", None) == "-"
         feature["strand"] = "+" if revd == seq_rced else "-"
@@ -1160,7 +1160,7 @@ class Sequence(SequenceI):
 
         s = moltype.coerce_str(self._seq.value)
         moltype.verify_sequence(s, gaps_allowed=True, wildcards_allowed=True)
-        sv = SeqView(s)
+        sv = SeqView(seq=s)
         new = moltype.make_seq(sv, name=self.name, info=self.info)
         new.annotation_db = self.annotation_db
         return new
@@ -1363,7 +1363,7 @@ class Sequence(SequenceI):
 
     def __str__(self):
         result = str(self._seq)
-        if self._seq.reverse:
+        if self._seq.is_reversed:
             with contextlib.suppress(TypeError):
                 result = self.moltype.complement(result)
         return result
@@ -1640,7 +1640,7 @@ class Sequence(SequenceI):
         seqid, start, end, strand of this sequence on the parent. strand is either
         -1 or 1.
         """
-        strand = -1 if self._seq.reverse else 1
+        strand = -1 if self._seq.is_reversed else 1
         return self._seq.seqid, self._seq.parent_start, self._seq.parent_stop, strand
 
 
@@ -1906,39 +1906,39 @@ def _input_vals_neg_step(seqlen, start, stop, step):
     return (0, 0, 1) if start < stop else (start, stop, step)
 
 
-class SeqView:
-    __slots__ = ("seq", "start", "stop", "step", "_offset", "_seqid", "_seq_len")
+class SliceRecordABC(ABC):
+    """Abstract base class for recording the history of operations to be applied
+    to some underlying data
 
-    def __init__(
-        self,
-        seq,
-        *,
-        start: Optional[int] = None,
-        stop: Optional[int] = None,
-        step: Optional[int] = None,
-        offset: int = 0,
-        seqid: Optional[str] = None,
-        seq_len: Optional[int] = None,
-    ):
-        if step == 0:
-            raise ValueError("step cannot be 0")
-        step = 1 if step is None else step
+    Notes
+    -----
+    seq_len refers to a Python typing.Sequence object, e.g. array, str, list.
+    """
 
-        self._seq_len = self._checked_seq_len(seq, seq_len)
-        func = _input_vals_pos_step if step > 0 else _input_vals_neg_step
-        start, stop, step = func(self._seq_len, start, stop, step)
-        self.seq = seq
-        self.start = start
-        self.stop = stop
-        self.step = step
-        self._offset = offset
-        self._seqid = seqid
+    start: int
+    stop: int
+    step: int
+    seq_len: int
+    _offset: int
 
-    def _checked_seq_len(self, seq, seq_len) -> int:
-        if seq_len is not None and seq_len != len(seq):
-            raise AssertionError(f"{seq_len} != {len(seq)})")
-        seq_len = seq_len or len(seq)
-        return seq_len
+    @abstractmethod
+    def _get_init_kwargs(self) -> dict:
+        """return required arguments for construction that are unique to the
+        subclass"""
+        ...
+
+    @abstractmethod
+    def copy(self): ...
+
+    @property
+    @abstractmethod
+    def _zero_slice(self): ...
+
+    @abstractmethod
+    def to_rich_dict(self): ...
+
+    @abstractmethod
+    def from_rich_dict(self): ...
 
     @property
     def offset(self) -> int:
@@ -1949,14 +1949,6 @@ class SeqView:
         value = value or 0
         self._offset = int(value)
 
-    @property
-    def seqid(self) -> str:
-        return self._seqid
-
-    @property
-    def seq_len(self) -> int:
-        return self._seq_len
-
     @property
     def parent_start(self) -> int:
         """returns the start on the parent plus strand
@@ -1966,7 +1958,7 @@ class SeqView:
         offset + start, taking into account whether reversed. Result
         is positive.
         """
-        if self.reverse:
+        if self.is_reversed:
             # self.stop becomes the start, self.stop will be negative
             assert self.stop < 0, "expected stop on reverse strand SeqView < 0"
             start = self.stop + self.seq_len + 1
@@ -1975,6 +1967,10 @@ class SeqView:
 
         return self.offset + start
 
+    @property
+    def is_reversed(self):
+        return self.step < 0
+
     @property
     def parent_stop(self) -> int:
         """returns the stop on the parent plus strand
@@ -1984,7 +1980,7 @@ class SeqView:
         offset + stop, taking into account whether reversed. Result
         is positive.
         """
-        if self.reverse:
+        if self.is_reversed:
             # self.start becomes the stop, self.start will be negative
             assert self.start < 0, "expected start on reverse strand SeqView < 0"
             stop = self.start + self.seq_len + 1
@@ -1992,13 +1988,9 @@ class SeqView:
             stop = self.stop
         return self.offset + stop
 
-    @property
-    def reverse(self):
-        return self.step < 0
-
-    def absolute_position(self, rel_index: int, include_boundary=False):
+    def absolute_position(self, rel_index: int, include_boundary: bool = False):
         """Converts an index relative to the current view to be with respect
-        to the coordinates of the sequence's annotations
+        to the coordinates of the original "Python sequence".
 
         Parameters
         ----------
@@ -2007,8 +1999,8 @@ class SeqView:
 
         Returns
         -------
-        the absolute index with respect to the coordinates of the self
-        including offset
+        the absolute index with respect to the coordinates of the original
+        sequence (including offset if present).
         """
         if not self:
             return 0
@@ -2019,20 +2011,22 @@ class SeqView:
         # _get_index return the absolute position relative to the underlying sequence
         seq_index, _, _ = self._get_index(rel_index, include_boundary=include_boundary)
 
-        # add offset and handle reversed views, now absolute relative to annotation coordinates
         offset = self.offset
-        if self.reverse:
+        if self.is_reversed:
             abs_index = offset + self.seq_len + seq_index + 1
         else:
             abs_index = offset + seq_index
 
         return abs_index
 
-    def relative_position(self, abs_index, stop=False):
-        """converts an index relative to annotation coordinates to be with respect to the current sequence view
+    def relative_position(self, abs_index: int, stop: bool = False):
+        """converts an index on the original "Python sequence" into an index
+        on this "view"
 
-        NOTE: the returned value DOES NOT reflect python indexing. Importantly, negative values represent positions that
-        precede the current view.
+        Notes
+        -----
+        The returned value DOES NOT reflect python indexing. Importantly,
+        negative values represent positions that precede the current view.
         """
         if not self:
             return 0
@@ -2040,7 +2034,7 @@ class SeqView:
         if abs_index < 0:
             raise IndexError("Index must be +ve and relative to the + strand")
 
-        if self.reverse:
+        if self.is_reversed:
             offset = self.offset
 
             if (
@@ -2060,7 +2054,41 @@ class SeqView:
 
         return rel_pos
 
-    def _get_index(self, val, include_boundary=False):
+    def __getitem__(self, segment: Union[int, slice]):
+        kwargs = self._get_init_kwargs()
+
+        if _is_int(segment):
+            start, stop, step = self._get_index(segment)
+            return self.__class__(
+                start=start,
+                stop=stop,
+                step=step,
+                offset=self.offset,
+                seq_len=self.seq_len,
+                **kwargs,
+            )
+
+        if segment.start is segment.stop is segment.step is None:
+            return self.copy()
+
+        if len(self) == 0:
+            return self
+
+        if segment.start is not None and segment.start == segment.stop:
+            return self._zero_slice
+
+        slice_step = 1 if segment.step is None else segment.step
+
+        if slice_step > 0:
+            return self._get_slice(segment, slice_step, **kwargs)
+        elif slice_step < 0:
+            return self._get_reverse_slice(segment, slice_step, **kwargs)
+        else:
+            raise ValueError(
+                f"{self.__class__.__name__} cannot be sliced with a step of 0"
+            )
+
+    def _get_index(self, val: int, include_boundary: bool = False):
         if len(self) == 0:
             raise IndexError(val)
 
@@ -2089,21 +2117,23 @@ class SeqView:
 
             return val, val - 1, -1
 
-    def _get_slice(self, segment, step):
+    def _get_slice(self, segment: slice, step: int, **kwargs):
         slice_start = segment.start if segment.start is not None else 0
         slice_stop = segment.stop if segment.stop is not None else len(self)
 
         if self.step > 0:
             return self._get_forward_slice_from_forward_seqview_(
-                slice_start, slice_stop, step
+                slice_start, slice_stop, step, **kwargs
             )
 
         elif self.step < 0:
             return self._get_forward_slice_from_reverse_seqview_(
-                slice_start, slice_stop, step
+                slice_start, slice_stop, step, **kwargs
             )
 
-    def _get_forward_slice_from_forward_seqview_(self, slice_start, slice_stop, step):
+    def _get_forward_slice_from_forward_seqview_(
+        self, slice_start: int, slice_stop: int, step: int, **kwargs
+    ):
         start = (
             self.start + slice_start * self.step
             if slice_start >= 0
@@ -2123,25 +2153,26 @@ class SeqView:
 
         # if -ve, it's an invalid slice
         if start < 0 or stop < 0:
-            return _zero_slice
+            return self._zero_slice
 
         # checking for zero-length slice
         if stop < start:
-            return _zero_slice
+            return self._zero_slice
         if start > self.seq_len:
-            return _zero_slice
+            return self._zero_slice
 
         return self.__class__(
-            self.seq,
             start=start,
             stop=min(self.stop, stop),
             step=self.step * step,
             offset=self.offset,
-            seqid=self.seqid,
             seq_len=self.seq_len,
+            **kwargs,
         )
 
-    def _get_forward_slice_from_reverse_seqview_(self, slice_start, slice_stop, step):
+    def _get_forward_slice_from_reverse_seqview_(
+        self, slice_start: int, slice_stop: int, step: int, **kwargs
+    ):
         if slice_start >= 0:
             start = self.start + slice_start * self.step
         elif abs(slice_start) > len(self):
@@ -2156,32 +2187,33 @@ class SeqView:
 
         # if +ve, it's an invalid slice
         if start >= 0 or stop >= 0:
-            return _zero_slice
+            return self._zero_slice
 
         return self.__class__(
-            self.seq,
             start=start,
             stop=max(self.stop, stop),
             step=self.step * step,
             offset=self.offset,
-            seqid=self.seqid,
             seq_len=self.seq_len,
+            **kwargs,
         )
 
-    def _get_reverse_slice(self, segment, step):
+    def _get_reverse_slice(self, segment: slice, step: int, **kwargs):
         slice_start = segment.start if segment.start is not None else -1
         slice_stop = segment.stop if segment.stop is not None else -len(self) - 1
 
         if self.step < 0:
             return self._get_reverse_slice_from_reverse_seqview_(
-                slice_start, slice_stop, step
+                slice_start, slice_stop, step, **kwargs
             )
         elif self.step > 0:
             return self._get_reverse_slice_from_forward_seqview_(
-                slice_start, slice_stop, step
+                slice_start, slice_stop, step, **kwargs
             )
 
-    def _get_reverse_slice_from_forward_seqview_(self, slice_start, slice_stop, step):
+    def _get_reverse_slice_from_forward_seqview_(
+        self, slice_start: int, slice_stop: int, step: int, **kwargs
+    ):
         # "true stop" adjust for if abs(stop-start) % step != 0
         # max possible start is "true stop" - step, because stop is not inclusive
         # "true stop" - step is converted to -ve index via subtracting len(self)
@@ -2198,7 +2230,7 @@ class SeqView:
             )
 
         if slice_stop >= self.seq_len:
-            return _zero_slice
+            return self._zero_slice
 
         if slice_stop >= 0:
             stop = self.start + (slice_stop * self.step) - self.seq_len
@@ -2211,19 +2243,20 @@ class SeqView:
             )
 
         if start >= 0 or stop >= 0:
-            return _zero_slice
+            return self._zero_slice
 
         return self.__class__(
-            self.seq,
             start=start,
             stop=max(stop, self.start - self.seq_len - 1),
             step=self.step * step,
             offset=self.offset,
-            seqid=self.seqid,
             seq_len=self.seq_len,
+            **kwargs,
         )
 
-    def _get_reverse_slice_from_reverse_seqview_(self, slice_start, slice_stop, step):
+    def _get_reverse_slice_from_reverse_seqview_(
+        self, slice_start: int, slice_stop: int, step: int, **kwargs
+    ):
         # max start is "true stop" + abs(step), because stop is not inclusive
         # "true stop" adjust for if abs(stop-start) % step != 0
         if slice_start >= len(self):
@@ -2238,7 +2271,7 @@ class SeqView:
         if slice_stop >= 0:
             stop = self.seq_len + (self.start + slice_stop * self.step)
             if stop <= self.seq_len + self.stop:
-                return _zero_slice
+                return self._zero_slice
         else:
             stop = self.seq_len + (
                 self.start + len(self) * self.step + slice_stop * self.step
@@ -2249,60 +2282,72 @@ class SeqView:
         # if -ve, it's an invalid slice becomes zero
         # checking for zero-length slice
         if stop < start or start > self.seq_len or min(start, stop) < 0:
-            return _zero_slice
+            return self._zero_slice
 
         return self.__class__(
-            self.seq,
             start=start,
             stop=stop,
             step=self.step * step,
             offset=self.offset,
-            seqid=self.seqid,
             seq_len=self.seq_len,
+            **kwargs,
         )
 
-    def __getitem__(self, segment):
-        if _is_int(segment):
-            start, stop, step = self._get_index(segment)
-            return self.__class__(
-                self.seq,
-                start=start,
-                stop=stop,
-                step=step,
-                offset=self.offset,
-                seqid=self.seqid,
-                seq_len=self.seq_len,
-            )
 
-        if segment.start is segment.stop is segment.step is None:
-            return self.copy(sliced=False)
+class SeqView(SliceRecordABC):
+    __slots__ = ("seq", "start", "stop", "step", "_offset", "_seqid", "_seq_len")
 
-        if len(self) == 0:
-            return self
+    def __init__(
+        self,
+        *,
+        seq: str,
+        start: Optional[int] = None,
+        stop: Optional[int] = None,
+        step: Optional[int] = None,
+        offset: int = 0,
+        seqid: Optional[str] = None,
+        seq_len: Optional[int] = None,
+    ):
+        if step == 0:
+            raise ValueError("step cannot be 0")
+        step = 1 if step is None else step
 
-        if segment.start is not None and segment.start == segment.stop:
-            return _zero_slice
+        func = _input_vals_pos_step if step > 0 else _input_vals_neg_step
+        start, stop, step = func(len(seq), start, stop, step)
+        self.seq = seq
+        self.start = start
+        self.stop = stop
+        self.step = step
+        self._offset = offset
+        self._seqid = seqid
+        if seq_len is not None and seq_len != len(seq):
+            raise AssertionError(f"{seq_len} != {len(self.seq)})")
+        self._seq_len = seq_len or len(self.seq)
 
-        slice_step = 1 if segment.step is None else segment.step
+    @property
+    def _zero_slice(self):
+        return self.__class__(seq="")
 
-        if slice_step > 0:
-            return self._get_slice(segment, slice_step)
-        elif slice_step < 0:
-            return self._get_reverse_slice(segment, slice_step)
-        else:
-            raise ValueError(
-                f"{self.__class__.__name__} cannot be sliced with a step of 0"
-            )
+    @property
+    def seqid(self) -> str:
+        return self._seqid
+
+    @property
+    def seq_len(self) -> int:
+        return self._seq_len
+
+    def _get_init_kwargs(self):
+        return {"seq": self.seq, "seqid": self.seqid}
 
     @property
     def value(self):
         return self.seq[self.start : self.stop : self.step]
 
-    def replace(self, old, new):
+    def replace(self, old: str, new: str):
         new_seq = self.seq.replace(old, new)
         if len(old) == len(new):
             return self.__class__(
-                new_seq,
+                seq=new_seq,
                 start=self.start,
                 stop=self.stop,
                 step=self.step,
@@ -2311,7 +2356,7 @@ class SeqView:
                 seq_len=self.seq_len,
             )
 
-        return self.__class__(new_seq)
+        return self.__class__(seq=new_seq)
 
     def __len__(self):
         return abs((self.start - self.stop) // self.step)
@@ -2333,10 +2378,11 @@ class SeqView:
     def to_rich_dict(self):
         # get the current state
         data = {"type": get_object_provenance(self), "version": __version__}
+        data["init_args"] = self._get_init_kwargs()
         # since we will truncate the seq, we don't need start, stop,
         # step is sufficient
-        data["init_args"] = {"step": self.step}
-        if self.reverse:
+        data["init_args"]["step"] = self.step
+        if self.is_reversed:
             adj = self.seq_len + 1
             start, stop = self.stop + adj, self.start + adj
         else:
@@ -2344,7 +2390,6 @@ class SeqView:
 
         data["init_args"]["seq"] = self.seq[start:stop]
         data["init_args"]["offset"] = int(self.parent_start)
-        data["init_args"]["seqid"] = self.seqid
         return data
 
     @classmethod
@@ -2355,7 +2400,7 @@ class SeqView:
         sv = cls(**init_args)
         return sv
 
-    def copy(self, sliced=False):
+    def copy(self, sliced: bool = False):
         """returns copy
 
         Parameters
@@ -2366,7 +2411,7 @@ class SeqView:
         """
         if not sliced:
             return self.__class__(
-                self.seq,
+                seq=self.seq,
                 start=self.start,
                 stop=self.stop,
                 step=self.step,
@@ -2377,9 +2422,6 @@ class SeqView:
         return self.from_rich_dict(self.to_rich_dict())
 
 
-_zero_slice = SeqView(seq="")
-
-
 class DnaSequence(NucleicAcidSequence):
     """Holds the standard DNA sequence."""
 
@@ -2993,7 +3035,7 @@ def _(data: str, seqid, preserve_case, checker):
     if not preserve_case:
         data = data.upper()
     checker(data)
-    return SeqView(data, seqid=seqid)
+    return SeqView(seq=data, seqid=seqid)
 
 
 @_coerce_to_seqview.register
@@ -3002,7 +3044,7 @@ def _(data: bytes, seqid, preserve_case, checker):
         data = data.upper()
     data = data.decode("utf8")
     checker(data)
-    return SeqView(data, seqid=seqid)
+    return SeqView(seq=data, seqid=seqid)
 
 
 @_coerce_to_seqview.register
