"""No third-party infrastructure topology in a public repository.

The August sanitisation commit said the cluster topology had been removed. It
replaced the account name and the key path with placeholders and left both real
hostnames in the committed example file, where they stayed for three days. The
commit message was more thorough than the commit.

This asserts the property rather than trusting a commit message.
"""
from __future__ import annotations

import pathlib
import re
import subprocess

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent

#: Public API hosts and documentation placeholders that are meant to be here.
ALLOWED = re.compile(r"(ebi\.ac\.uk|\bexample\.ac\.uk|[a-z-]+\.example\.ac\.uk)")

HOSTNAME = re.compile(r"\b[a-zA-Z0-9][a-zA-Z0-9.-]*\.ac\.uk\b")


def _tracked() -> list[pathlib.Path]:
    out = subprocess.run(["git", "ls-files"], cwd=REPO,
                         capture_output=True, text=True).stdout.split("\n")
    return [REPO / f for f in out if f.strip()]


class TestNoRealHostnamesAreTracked:
    def test_working_tree_is_clean_of_topology(self):
        offenders = {}
        for p in _tracked():
            if not p.exists() or p.suffix in {".json", ".gz", ".bgz", ".tbi"}:
                continue
            try:
                text = p.read_text(errors="ignore")
            except (OSError, UnicodeDecodeError):
                continue
            hits = [h for h in HOSTNAME.findall(text) if not ALLOWED.search(h)]
            if hits:
                offenders[str(p.relative_to(REPO))] = sorted(set(hits))
        assert not offenders, (
            f"third-party hostnames are tracked in a public repository: "
            f"{ {k: len(v) for k, v in offenders.items()} }")

    def test_the_example_file_uses_placeholders(self):
        p = REPO / "scripts" / "gpu" / ".local.sh.example"
        if not p.exists():
            pytest.skip("example file absent")
        text = p.read_text()
        for var in ("GPU_USER", "JUMP_HOST", "GPU_NODE", "SSH_KEY"):
            m = re.search(rf"export {var}=(\S+)", text)
            assert m, f"{var} is missing from the example"
            value = m.group(1)
            assert ("example" in value or "your" in value or value.startswith("$")), (
                f"{var} in the committed example looks like a real value, not a "
                f"placeholder")

    def test_no_ssh_key_path_that_is_not_a_placeholder(self):
        p = REPO / "scripts" / "gpu" / ".local.sh.example"
        if not p.exists():
            pytest.skip("example file absent")
        m = re.search(r"export SSH_KEY=(\S+)", p.read_text())
        assert m and "your" in m.group(1), (
            "the example declares a specific private key path")


class TestHistoryIsCleanToo:
    """A public repository leaks its history, not just its tip."""

    def test_no_real_hostname_survives_in_any_reachable_commit(self):
        revs = subprocess.run(["git", "rev-list", "--all"], cwd=REPO,
                              capture_output=True, text=True).stdout.split()
        if not revs:
            pytest.skip("no history")
        for rev in revs:
            listing = subprocess.run(
                ["git", "ls-tree", "-r", "--name-only", rev], cwd=REPO,
                capture_output=True, text=True).stdout.split("\n")
            for path in listing:
                if not path.strip() or not path.startswith("scripts/gpu"):
                    continue
                blob = subprocess.run(["git", "show", f"{rev}:{path}"], cwd=REPO,
                                      capture_output=True, text=True).stdout
                hits = [h for h in HOSTNAME.findall(blob) if not ALLOWED.search(h)]
                assert not hits, (
                    f"commit {rev[:7]} still carries {len(hits)} third-party "
                    f"hostname(s) in {path}")
