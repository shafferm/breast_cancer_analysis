"""
Tests for etl/liftover_mutations.py helper functions.

Uses an inline MockLO instead of the real pyliftover — no network access required.

Run with: python -m pytest tests/test_liftover.py -v
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import etl._bootstrap  # noqa: F401
from etl.liftover_mutations import _lift_position


# ---------------------------------------------------------------------------
# Inline LiftOver mock
# ---------------------------------------------------------------------------

class MockLO:
    """Minimal LiftOver stand-in that returns pre-configured results."""

    def __init__(self, results):
        # Each element is a return value for one convert_coordinate call
        self.results = list(results)
        self.calls: list[tuple] = []

    def convert_coordinate(self, chrom: str, pos: int):
        self.calls.append((chrom, pos))
        return self.results.pop(0)


# ---------------------------------------------------------------------------
# TestLiftPosition
# ---------------------------------------------------------------------------

class TestLiftPosition:
    def test_none_chrom(self):
        lo = MockLO([])
        result = _lift_position(lo, None, 100)
        assert result == (None, None)
        assert len(lo.calls) == 0

    def test_none_pos(self):
        lo = MockLO([])
        result = _lift_position(lo, "1", None)
        assert result == (None, None)
        assert len(lo.calls) == 0

    def test_adds_chr_prefix(self):
        lo = MockLO([[("chr7", 99, "+", 1.0)]])
        _lift_position(lo, "7", 100)
        assert lo.calls[0][0] == "chr7"

    def test_keeps_existing_chr(self):
        lo = MockLO([[("chr7", 99, "+", 1.0)]])
        _lift_position(lo, "chr7", 100)
        # Should not become "chrchr7"
        assert lo.calls[0][0] == "chr7"

    def test_empty_result(self):
        lo = MockLO([[]])
        assert _lift_position(lo, "1", 100) == (None, None)

    def test_none_result(self):
        lo = MockLO([None])
        assert _lift_position(lo, "1", 100) == (None, None)

    def test_valid_lift_strips_chr(self):
        lo = MockLO([[("chr7", 99, "+", 1.0)]])
        chrom, _ = _lift_position(lo, "7", 100)
        assert chrom == "7"

    def test_valid_lift_1based_output(self):
        # Converter returns 0-based position 99 → function should return 1-based 100
        lo = MockLO([[("chr7", 99, "+", 1.0)]])
        _, pos = _lift_position(lo, "7", 100)
        assert pos == 100

    def test_input_converted_to_0based(self):
        # 1-based input 100 → converter should be called with 99 (0-based)
        lo = MockLO([[("chr7", 99, "+", 1.0)]])
        _lift_position(lo, "7", 100)
        assert lo.calls[0][1] == 99

    def test_takes_first_mapping(self):
        # Two results returned; only first should be used
        lo = MockLO([[("chr7", 99, "+", 1.0), ("chr8", 200, "+", 0.5)]])
        chrom, pos = _lift_position(lo, "7", 100)
        assert chrom == "7"
        assert pos == 100
