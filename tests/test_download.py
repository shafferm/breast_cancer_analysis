"""
Tests for etl/download_data.py.

Network calls are monkeypatched — no real downloads occur.

Run with: python -m pytest tests/test_download.py -v
"""

import sys
import tarfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import etl._bootstrap  # noqa: F401
from etl.download_data import (
    CDR_FALLBACK,
    CDR_FILENAME,
    CDR_URL,
    _extract_tarball,
    _has_expected_files,
    _validate,
    download_cdr,
    download_metabric,
    download_tcga,
)


# ---------------------------------------------------------------------------
# _has_expected_files
# ---------------------------------------------------------------------------

class TestHasExpectedFiles:
    def test_nonexistent_dir(self, tmp_path):
        assert _has_expected_files(tmp_path / "no_such_dir", ["a.txt"]) is False

    def test_empty_dir(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        assert _has_expected_files(d, ["a.txt"]) is False

    def test_all_present(self, tmp_path):
        d = tmp_path / "full"
        d.mkdir()
        (d / "a.txt").write_text("x")
        (d / "b.txt").write_text("x")
        assert _has_expected_files(d, ["a.txt", "b.txt"]) is True

    def test_partial_returns_false(self, tmp_path):
        d = tmp_path / "partial"
        d.mkdir()
        (d / "a.txt").write_text("x")
        assert _has_expected_files(d, ["a.txt", "b.txt"]) is False

    def test_unrelated_file_not_enough(self, tmp_path):
        d = tmp_path / "wrong"
        d.mkdir()
        (d / "unrelated.xlsx").write_text("x")
        assert _has_expected_files(d, ["data_mutations.txt"]) is False


# ---------------------------------------------------------------------------
# _extract_tarball
# ---------------------------------------------------------------------------

class TestExtractTarball:
    def _make_tarball(self, tmp_path, arcname: str, content: bytes = b"hello") -> Path:
        """Create a .tar.gz with a single file at the given arcname path."""
        src = tmp_path / "src_file.txt"
        src.write_bytes(content)
        tarball = tmp_path / "test.tar.gz"
        with tarfile.open(tarball, "w:gz") as tf:
            tf.add(src, arcname=arcname)
        return tarball

    def test_strips_top_level_dir(self, tmp_path):
        tarball = self._make_tarball(tmp_path, "toplevel/file.txt")
        dest = tmp_path / "dest"
        _extract_tarball(tarball, dest, "test")
        assert (dest / "file.txt").exists()
        assert not (dest / "toplevel").exists()

    def test_file_content_preserved(self, tmp_path):
        content = b"test content 12345"
        tarball = self._make_tarball(tmp_path, "toplevel/data.txt", content)
        dest = tmp_path / "dest"
        _extract_tarball(tarball, dest, "test")
        assert (dest / "data.txt").read_bytes() == content

    def test_skips_bare_dir_entry(self, tmp_path):
        """Top-level directory entry (no slash) must not create a stray artifact."""
        src = tmp_path / "src_file.txt"
        src.write_text("content")
        tarball = tmp_path / "test.tar.gz"
        with tarfile.open(tarball, "w:gz") as tf:
            # Add top-level dir entry
            info = tarfile.TarInfo(name="toplevel")
            info.type = tarfile.DIRTYPE
            tf.addfile(info)
            # Add file inside it
            tf.add(src, arcname="toplevel/file.txt")
        dest = tmp_path / "dest"
        _extract_tarball(tarball, dest, "test")
        # file.txt should exist; 'toplevel' dir entry should not create artifact at dest root
        assert (dest / "file.txt").exists()


# ---------------------------------------------------------------------------
# _validate
# ---------------------------------------------------------------------------

class TestValidate:
    def test_all_present_logs_info(self, tmp_path, caplog):
        import logging
        (tmp_path / "a.txt").write_text("x")
        (tmp_path / "b.txt").write_text("x")
        with caplog.at_level(logging.INFO):
            _validate(tmp_path, ["a.txt", "b.txt"], "TEST")
        assert any("All expected files present" in r.message for r in caplog.records)

    def test_missing_warns(self, tmp_path, caplog):
        import logging
        (tmp_path / "present.txt").write_text("x")
        with caplog.at_level(logging.WARNING):
            _validate(tmp_path, ["present.txt", "missing.txt"], "TEST")
        warning_msgs = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("missing.txt" in m for m in warning_msgs)


# ---------------------------------------------------------------------------
# download_cdr
# ---------------------------------------------------------------------------

class TestDownloadCdr:
    def _patch_tcga_dir(self, monkeypatch, tmp_path) -> Path:
        tcga_dir = tmp_path / "tcga"
        tcga_dir.mkdir()
        monkeypatch.setattr("etl.download_data.TCGA_DIR", tcga_dir)
        return tcga_dir

    def test_skips_if_file_present(self, tmp_path, monkeypatch):
        tcga_dir = self._patch_tcga_dir(monkeypatch, tmp_path)
        (tcga_dir / CDR_FILENAME).write_bytes(b"x" * 60_000)

        calls = []
        monkeypatch.setattr("urllib.request.urlretrieve",
                            lambda url, dest, reporthook=None: calls.append(url))

        download_cdr(force=False)
        assert len(calls) == 0

    def test_force_redownloads(self, tmp_path, monkeypatch):
        tcga_dir = self._patch_tcga_dir(monkeypatch, tmp_path)
        dest = tcga_dir / CDR_FILENAME
        dest.write_bytes(b"x" * 60_000)

        def fake_urlretrieve(url, path, reporthook=None):
            Path(path).write_bytes(b"x" * 60_000)

        monkeypatch.setattr("urllib.request.urlretrieve", fake_urlretrieve)

        download_cdr(force=True)
        assert dest.exists()

    def test_primary_url_success(self, tmp_path, monkeypatch):
        tcga_dir = self._patch_tcga_dir(monkeypatch, tmp_path)
        urls_called = []

        def fake_urlretrieve(url, path, reporthook=None):
            urls_called.append(url)
            Path(path).write_bytes(b"x" * 60_000)

        monkeypatch.setattr("urllib.request.urlretrieve", fake_urlretrieve)

        download_cdr(force=False)
        assert (tcga_dir / CDR_FILENAME).exists()
        assert urls_called[0] == CDR_URL
        assert len(urls_called) == 1  # no fallback needed

    def test_falls_back_if_too_small(self, tmp_path, monkeypatch):
        tcga_dir = self._patch_tcga_dir(monkeypatch, tmp_path)
        urls_called = []

        def fake_urlretrieve(url, path, reporthook=None):
            urls_called.append(url)
            # Write <50 KB (simulating a JSON error response from GDC)
            Path(path).write_bytes(b"x" * 100)

        monkeypatch.setattr("urllib.request.urlretrieve", fake_urlretrieve)

        download_cdr(force=False)
        assert CDR_URL in urls_called
        assert CDR_FALLBACK in urls_called

    def test_falls_back_on_oserror(self, tmp_path, monkeypatch):
        tcga_dir = self._patch_tcga_dir(monkeypatch, tmp_path)
        urls_called = []

        def fake_urlretrieve(url, path, reporthook=None):
            urls_called.append(url)
            if url == CDR_URL:
                raise OSError("network failure")
            Path(path).write_bytes(b"x" * 60_000)

        monkeypatch.setattr("urllib.request.urlretrieve", fake_urlretrieve)

        download_cdr(force=False)
        assert CDR_URL in urls_called
        assert CDR_FALLBACK in urls_called
        assert (tcga_dir / CDR_FILENAME).exists()


# ---------------------------------------------------------------------------
# download_metabric / download_tcga — idempotency
# ---------------------------------------------------------------------------

class TestDownloadIdempotency:
    def _create_metabric_files(self, metabric_dir):
        from etl.download_data import METABRIC_EXPECTED
        metabric_dir.mkdir(exist_ok=True)
        for f in METABRIC_EXPECTED:
            (metabric_dir / f).write_text("existing")

    def _create_tcga_files(self, tcga_dir):
        from etl.download_data import TCGA_EXPECTED
        tcga_dir.mkdir(exist_ok=True)
        for f in TCGA_EXPECTED:
            (tcga_dir / f).write_text("existing")

    def test_metabric_skips_if_all_expected(self, tmp_path, monkeypatch):
        metabric_dir = tmp_path / "metabric"
        self._create_metabric_files(metabric_dir)
        monkeypatch.setattr("etl.download_data.METABRIC_DIR", metabric_dir)

        calls = []
        monkeypatch.setattr("etl.download_data._download",
                            lambda *a, **kw: calls.append(a))

        download_metabric(force=False)
        assert len(calls) == 0

    def test_metabric_downloads_if_partial(self, tmp_path, monkeypatch):
        """Only some expected files present — should re-download."""
        metabric_dir = tmp_path / "metabric"
        metabric_dir.mkdir()
        (metabric_dir / "data_clinical_patient.txt").write_text("existing")
        monkeypatch.setattr("etl.download_data.METABRIC_DIR", metabric_dir)

        calls = []
        monkeypatch.setattr("etl.download_data._download",
                            lambda *a, **kw: calls.append(a))
        monkeypatch.setattr("etl.download_data._extract_tarball",
                            lambda *a, **kw: None)

        download_metabric(force=False)
        assert len(calls) == 1

    def test_metabric_force_downloads(self, tmp_path, monkeypatch):
        metabric_dir = tmp_path / "metabric"
        self._create_metabric_files(metabric_dir)
        monkeypatch.setattr("etl.download_data.METABRIC_DIR", metabric_dir)

        calls = []
        monkeypatch.setattr("etl.download_data._download",
                            lambda *a, **kw: calls.append(a))
        monkeypatch.setattr("etl.download_data._extract_tarball",
                            lambda *a, **kw: None)

        download_metabric(force=True)
        assert len(calls) == 1

    def test_tcga_skips_if_all_expected(self, tmp_path, monkeypatch):
        tcga_dir = tmp_path / "tcga"
        self._create_tcga_files(tcga_dir)
        monkeypatch.setattr("etl.download_data.TCGA_DIR", tcga_dir)

        calls = []
        monkeypatch.setattr("etl.download_data._download",
                            lambda *a, **kw: calls.append(a))
        monkeypatch.setattr("etl.download_data.download_cdr",
                            lambda force=False: None)

        download_tcga(force=False)
        assert len(calls) == 0

    def test_tcga_downloads_if_only_cdr_present(self, tmp_path, monkeypatch):
        """CDR file present but tarball files missing — must not skip tarball."""
        tcga_dir = tmp_path / "tcga"
        tcga_dir.mkdir()
        (tcga_dir / CDR_FILENAME).write_bytes(b"x" * 60_000)
        monkeypatch.setattr("etl.download_data.TCGA_DIR", tcga_dir)

        calls = []
        monkeypatch.setattr("etl.download_data._download",
                            lambda *a, **kw: calls.append(a))
        monkeypatch.setattr("etl.download_data._extract_tarball",
                            lambda *a, **kw: None)
        monkeypatch.setattr("etl.download_data.download_cdr",
                            lambda force=False: None)

        download_tcga(force=False)
        assert len(calls) == 1  # tarball download should proceed
