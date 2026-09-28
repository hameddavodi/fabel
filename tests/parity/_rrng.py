"""A stand-in for R's default random number generator, for seeded parity cases.

R's permutation tests (``tperm.fd``, ``Fperm.fd``) draw their permutations with
``sample(n)`` after ``set.seed(seed)``.  fdatools draws them from a NumPy
``Generator``, so the two null distributions are different samples from the
same law.  To compare the *seeded* golden outputs value for value, the parity
tests hand fdatools an object that exposes the one method fdatools calls,
``permutation(n)``, and answers with exactly the permutation R would draw.

R's default generator is the standard Mersenne Twister (MT19937, Matsumoto and
Nishimura, 1998), which NumPy ships as :class:`numpy.random.MT19937`.  What is
specific to R is documented in its ``?RNG`` and ``?Random`` help pages:

* ``set.seed(s)`` scrambles the integer seed with the linear congruential map
  ``s -> 69069 s + 1 (mod 2^32)`` -- fifty times, then once more for every word
  of generator state -- and marks the state as exhausted, so that the first
  draw regenerates it.
* ``unif_rand()`` is the 32-bit output scaled by ``2^-32`` and nudged into the
  open interval ``(0, 1)``.
* Since R 3.6.0 ``sample()`` uses rejection sampling (``sample.kind =
  "Rejection"``): an index below ``m`` is drawn from ``ceil(log2 m)`` random
  bits, assembled 16 bits per uniform, and redrawn until it falls below ``m``.
  A permutation is built by the textbook "draw a position, move the last
  element into its slot" shuffle.

This module is test-support only; nothing under ``src/`` imports it.  Its
correctness is pinned by :func:`tests.parity.test_stats.test_r_generator_matches_r`
against values printed by R 4.6.1 itself.
"""

from __future__ import annotations

import math

import numpy as np

_LCG_MULTIPLIER = 69069
_MASK32 = 0xFFFFFFFF
_STATE_WORDS = 624
_TWO_POW_M32 = 2.3283064365386963e-10
_FIXUP = 2.328306437080797e-10


class RRandom:
    """Reproduce ``set.seed(seed)`` followed by ``runif`` / ``sample`` in R.

    Parameters
    ----------
    seed : int
        The argument given to R's ``set.seed``.
    """

    def __init__(self, seed: int) -> None:
        scrambled = seed & _MASK32
        for _ in range(50):
            scrambled = (_LCG_MULTIPLIER * scrambled + 1) & _MASK32
        words = []
        for _ in range(_STATE_WORDS + 1):
            scrambled = (_LCG_MULTIPLIER * scrambled + 1) & _MASK32
            words.append(scrambled)
        # The first word is R's position counter; set.seed overwrites it with
        # 624, i.e. "state exhausted, regenerate before the next draw".
        self._bits = np.random.MT19937()
        self._bits.state = {
            "bit_generator": "MT19937",
            "state": {"key": np.asarray(words[1:], dtype=np.uint32), "pos": _STATE_WORDS},
        }

    def unif(self) -> float:
        """Return the next value of R's ``unif_rand()``."""
        value = float(self._bits.random_raw()) * _TWO_POW_M32
        if value <= 0.0:
            return 0.5 * _FIXUP
        if 1.0 - value <= 0.0:
            return 1.0 - 0.5 * _FIXUP
        return value

    def _bits_value(self, bits: int) -> int:
        value = 0
        for _ in range(0, bits + 1, 16):
            value = 65536 * value + math.floor(self.unif() * 65536)
        return value & ((1 << bits) - 1)

    def index(self, size: int) -> int:
        """Return a uniform index in ``[0, size)`` as R's rejection sampler does."""
        if size <= 0:
            return 0
        bits = math.ceil(math.log2(size))
        while True:
            value = self._bits_value(bits)
            if value < size:
                return value

    def sample(self, n: int, size: int) -> np.ndarray:
        """Return ``sample(n, size) - 1``: ``size`` draws without replacement."""
        pool = list(range(n))
        out = []
        remaining = n
        for _ in range(size):
            j = self.index(remaining)
            out.append(pool[j])
            remaining -= 1
            pool[j] = pool[remaining]
        return np.asarray(out, dtype=np.int64)

    def permutation(self, n: int) -> np.ndarray:
        """Return ``sample(n) - 1``: R's random permutation, zero-based."""
        return self.sample(n, n)
